"""The homework page builder's hostile set (tests/pedagogy/homework_page_set.json),
written before the builder exists.

Now: the set is what its generator writes; `allowed` is exactly the safe
outputs; every class has two or more variants; and the real plan says what
the homework says, over data that says what the cases assume. When
homework.page_html exists: every output is safe and preferred, and every
requirement holds. Until then those tests skip.
"""
import copy
import html
import json
import re
import sys
import unicodedata
from dataclasses import replace
from html.parser import HTMLParser
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
PED = Path(__file__).with_name("pedagogy")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PED))
from bundle import load  # noqa: E402

SET = json.loads((PED / "homework_page_set.json").read_text(encoding="utf-8"))
CASES, PLAN, DROP = SET["cases"], SET["plan"], SET["drop"]
DATA = HERE / "tests" / "data"
SITES = json.loads((DATA / "card_label_sites.json").read_text(encoding="utf-8"))
MEMBERS = {m["record"][:-5]: m for m in
           json.loads((DATA / "pair" / "members.json").read_text(encoding="utf-8"))["members"]}
STUDY_WORDS = r"\bresol(v|ut)\w*|\bunresol(v|ut)\w*|different on average"
BRIDGE = (r"because|explain\w*|due to|caus\w*|reactiv\w*|attract\w*|\bbind\w*|pheromon\w*|receptor\w*")


def test_the_committed_set_is_what_the_generator_writes():
    import make_homework_page_set
    assert make_homework_page_set.document() == (PED / "homework_page_set.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_allowed_is_exactly_the_safe_outputs(c):
    assert c["allowed"] == sorted(k for k, g in c["outcomes"].items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert c["prefer"] is None or c["prefer"] in c["allowed"]
    assert c["why"] and c["outcomes"]["ERROR"] == "CRASH"


def test_every_class_has_two_or_more_variants():
    from collections import Counter
    per = Counter(c["cls"].split()[0] for c in CASES)
    assert sorted(per) == [f"H{i}" for i in range(1, 10)] and min(per.values()) >= 2, per


# --- the facts the cases rest on --------------------------------------------------------
def _norm(s):
    return re.sub(r"\s+", " ", re.sub(r"[`*]", "", s)).strip().lower()      # case is a display choice


def test_the_plan_is_the_homework():
    """Every prompt and task in the plan is the homework's text (backticks and
    emphasis aside), and the molecules are its table."""
    md = _norm((HERE / "sesiidae-homework.md").read_text(encoding="utf-8"))
    for act in PLAN["acts"]:
        assert _norm(act["predict"]) in md, act["id"]
        for t in act["tasks"]:
            head = _norm(t).split(". These")[0].split(", with the difference card")[0]
            assert head in md, (act["id"], head)
    assert [a["id"] for a in PLAN["acts"]] == ["A", "B", "C", "D"]
    assert {m["role"] for m in PLAN["molecules"].values()} == {"atlas_observed", "teaching_control"}


def test_the_data_says_what_the_cases_say():
    from study import load_results
    b = load(DATA / "sesiidae-homework-1.bundle.json")
    r = load_results(DATA / "espmap-sesiidae-001.results.json")
    contrasts = {s.contrast for s in r.sites}
    for act in PLAN["acts"]:
        assert set(act["contrasts"]) <= contrasts, act["id"]
        for k in act["blocks"]:
            if k["view"] == "study_card":
                r.site(k["contrast"], k["site"])                     # exists, exactly once
            if k["view"] in ("study_card", "difference_card"):
                assert k["contrast"] in act["contrasts"], (act["id"], k)
            if k["view"] == "evidence_card":
                assert k["observation_id"] in b.observations
    for slug, m in PLAN["molecules"].items():
        src = SITES["molecules"][slug]["atlas_id"]
        if m["role"] == "atlas_observed":
            assert m["atlas_id"] == src and m["atlas_id"] in b.compounds, slug
        else:
            assert m["atlas_id"] is None and src.startswith("arrowfield.control."), slug
    assert len({x["analysis_key"] for v in SITES["molecules"].values() for x in v["members"]}) == 1
    # the premises of particular cases
    assert r.site("B_stereo", "C3-H").tier == "primary"
    assert r.site("A_headgroup", "head O").tier == "secondary"
    assert not [s for s in r.sites if s.contrast.startswith("A_headgroup") and s.site == "C1-H"]
    assert not [s for s in r.sites if s.contrast == "B_stereo" and s.site == "C20-H"]
    assert r.site("C_position", "C4-H").outcome == "resolved"


# --- the builder (to come) ----------------------------------------------------------------
def _builder():
    try:
        import homework
    except ImportError:
        homework = None
    if homework is None or not hasattr(homework, "page_html"):
        pytest.skip("homework.page_html is not built yet; its expectations are fixed in homework_page_set.json")
    return homework.page_html


def _renamed(x, old, new):
    """x with every key or text value exactly `old` renamed `new`."""
    if isinstance(x, dict):
        return {(new if k == old else k): _renamed(v, old, new) for k, v in x.items()}
    if isinstance(x, list):
        return [_renamed(v, old, new) for v in x]
    return new if x == old else x


def _patched_plan(c):
    plan = copy.deepcopy(PLAN)
    for op in c["source_ops"]:
        if isinstance(op, list) and op[0] in ("rename_molecule", "rename_contrast"):
            plan = _renamed(plan, op[1], op[2])
    for path, value in c["plan_patch"]:
        node = plan
        for key in path[:-1]:
            node = node[key]
        if path[-1] == "+":
            node.append(value)
        elif value == DROP:
            del node[path[-1]]
        else:
            node[path[-1]] = value
    return plan


def _ensemble(slug):
    from model import Atom, EnsembleSource, MoleculeView, Provenance, Scale, SiteMeasurement
    mol = SITES["molecules"][slug]
    views = []
    for m in mol["members"]:
        sites = {s: {rule: SiteMeasurement(s, (0,), rule, v[key], 10.0, 5.0, 5.0, v["partition"])
                     for rule, key in (("nearest_nucleus", "nn"), ("radius_scaled", "rs"))}
                 for s, v in m["sites"].items()}
        views.append(MoleculeView(
            name=f"{slug} {m['conformer']}", symbols=("C",), xyz=np.zeros((1, 3)), bonds=None, charge=0,
            atoms=(Atom(0, "C", 1),), sites=sites, sites_tile_surface=False, genus=0,
            scale=Scale("diverging", -20.0, 20.0, "vik"),
            provenance=Provenance(m["conformer"], m["analysis_key"], "r", espmap_version="0.2.1"),
            folded=m["folded"]))
    source = EnsembleSource(mol["atlas_id"], "espmap-sesiidae-001", mol["members"][0]["analysis_key"], "0.2.1",
                            mol["index_sha256"], "shape-diverse")
    return source, views


def _sources(c):
    import record
    from study import load_results
    b = load(DATA / "sesiidae-homework-1.bundle.json")
    r = load_results(DATA / "espmap-sesiidae-001.results.json")
    ensembles = {slug: _ensemble(slug) for slug in SITES["molecules"]}
    conformers = {}
    for stem, m in MEMBERS.items():
        slug, conf = stem.split("/")
        conformers[f"{slug}/{conf}"] = replace(record.load(DATA / "pair" / slug / conf), folded=m["folded"])
    for op in c["source_ops"]:
        if isinstance(op, list) and op[0] == "bundle_shorthand":
            b = replace(b, compounds={**b.compounds, op[1]: replace(b.compounds[op[1]], shorthand=op[2])})
        elif isinstance(op, list):
            r, ensembles, conformers = _list_op(op, r, ensembles, conformers)
        elif op == "e2_other_key":
            src, views = ensembles["e2_z13_18_oac"]
            other = "f" * 64
            views = [replace(v, provenance=replace(v.provenance, analysis_key=other)) for v in views]
            ensembles["e2_z13_18_oac"] = (replace(src, analysis_key=other), views)
        elif op == "results_unknown_rule":
            r = replace(r, verdict_rule="someone-else/v9")
        elif op == "bundle_unknown":
            b = replace(b, atlas_version="unknown")
        elif op == "plan_two_commits":
            r = replace(r, plan="Atlas analyses/espmap-sesiidae-001/analysis-plan.md @ 7e34658 @ deadbee")
        elif op == "plan_no_commit":
            r = replace(r, plan="Atlas analyses/espmap-sesiidae-001/analysis-plan.md")
        else:
            raise ValueError(op)
    return {"bundle": b, "study": r, "ensembles": ensembles, "conformers": conformers}


def _list_op(op, r, ensembles, conformers):
    """The round-7 probe's source edits: a slug or contrast renamed everywhere it is
    named, or one source field replaced, as a results file or a record could carry it."""
    kind = op[0]
    if kind == "rename_molecule":
        old, new = op[1], op[2]
        ensembles = {(new if k == old else k): v for k, v in ensembles.items()}
        conformers = {(f"{new}/{k.split('/', 1)[1]}" if k.split("/", 1)[0] == old else k): v
                      for k, v in conformers.items()}
        sites = tuple(replace(t, x=new if t.x == old else t.x, y=new if t.y == old else t.y) for t in r.sites)
        r = replace(r, sites=sites, molecules={(new if k == old else k): v for k, v in r.molecules.items()})
    elif kind == "rename_contrast":
        r = replace(r, sites=tuple(replace(t, contrast=op[2]) if t.contrast == op[1] else t for t in r.sites))
    elif kind == "study_plan":
        r = replace(r, plan=op[1])
    elif kind == "espmap_version":  # one record of one ensemble
        src, views = ensembles[op[1]]
        views = [replace(views[0], provenance=replace(views[0].provenance, espmap_version=op[2]))] + list(views[1:])
        ensembles = {**ensembles, op[1]: (src, views)}
    elif kind == "analysis_key":  # every record, so the page still sees one key
        def rekey(v):
            return replace(v, provenance=replace(v.provenance, analysis_key=op[1]))
        ensembles = {k: (replace(src, analysis_key=op[1]), [rekey(v) for v in vs]) for k, (src, vs) in ensembles.items()}
        conformers = {k: rekey(v) for k, v in conformers.items()}
    elif kind in ("settings", "settings_ensembles", "settings_conformers", "settings_ensemble"):
        # one record setting (xc, basis, solvent): on every record, every ensemble record, every
        # conformer record, or one molecule's ensemble
        field, value = op[-2], op[-1]

        def reset(v):
            return replace(v, provenance=replace(v.provenance, **{field: value}))
        if kind in ("settings", "settings_ensembles", "settings_ensemble"):
            ensembles = {k: (src, [reset(v) for v in vs]) if kind != "settings_ensemble" or k == op[1] else (src, vs)
                         for k, (src, vs) in ensembles.items()}
        if kind in ("settings", "settings_conformers"):
            conformers = {k: reset(v) for k, v in conformers.items()}
    else:
        raise ValueError(op)
    return r, ensembles, conformers


def _run(c):
    build = _builder()
    from model import Refused
    try:
        page = build(_patched_plan(c), _sources(c), audience=c["audience"])
    except Refused:
        return "REFUSED", None
    except Exception:  # noqa: BLE001 - any other error is a crash, graded as such
        return "ERROR", None
    return "RENDERED", page


# the probe on #21 (fe2fe08) and its reruns (318da98, 2a67bdd, 9eae791, 2b87a0a, 23d5ba2,
# 8fe0749, 64b5fa2): all met; only shrinks.
XFAIL = set()


def _xfail(request, key):
    if key in XFAIL:
        request.applymarker(pytest.mark.xfail(strict=True, reason="homework page on 64b5fa2"))


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_page_output_is_safe(c, request):
    _xfail(request, f"test_page_output_is_safe[{c['id']}]")
    out = _run(c)[0]
    assert out in c["allowed"], f"{out}: {c['outcomes'][out]}"


@pytest.mark.parametrize("c", [c for c in CASES if c["prefer"]], ids=lambda c: c["id"])
def test_page_output_is_preferred(c, request):
    _xfail(request, f"test_page_output_is_preferred[{c['id']}]")
    assert _run(c)[0] == c["prefer"]


class _Page(HTMLParser):
    """Text by region: inside verbatim spans, inside study cards, inside the
    evidence views, inside blocks, and the page's own words; plus the order
    of act, predict and block elements, and their attributes."""
    VOID = ("br", "img", "hr", "meta", "input", "link")

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack, self.events, self.text = [], [], []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag in self.VOID:
            return
        cls = set((a.get("class") or "").split())
        self.stack.append(cls)
        for k in ("act", "predict", "block", "page-provenance"):
            if k in cls:
                self.events.append((k, a))

    def handle_endtag(self, tag):
        if tag not in self.VOID and self.stack:
            self.stack.pop()

    def handle_data(self, d):
        inside = set().union(*self.stack) if self.stack else set()
        self.text.append((inside, d))

    def own(self, exclude):
        return " ".join(d for inside, d in self.text if not (inside & exclude))


# Default_Ignorable_Code_Point (Unicode 15.1 DerivedCoreProperties.txt), plus the symbols drawn
# as nothing by design: U+2800 BRAILLE PATTERN BLANK and U+1D159 MUSICAL SYMBOL NULL NOTEHEAD
_IGNORABLE = ((0x00AD, 0x00AD), (0x034F, 0x034F), (0x061C, 0x061C), (0x115F, 0x1160), (0x17B4, 0x17B5),
              (0x180B, 0x180F), (0x200B, 0x200F), (0x202A, 0x202E), (0x2060, 0x206F), (0x2800, 0x2800),
              (0x3164, 0x3164), (0xFE00, 0xFE0F), (0xFEFF, 0xFEFF), (0xFFA0, 0xFFA0), (0xFFF0, 0xFFF8),
              (0x1BCA0, 0x1BCA3), (0x1D159, 0x1D159), (0x1D173, 0x1D17A), (0xE0000, 0xE0FFF))


def _ignorable(ch):
    return any(lo <= ord(ch) <= hi for lo, hi in _IGNORABLE)


def _check(what, page, c):
    p = _Page()
    p.feed(page)
    plan = _patched_plan(c)
    if what == "block_roles":
        roles = {s: m["role"] for s, m in plan["molecules"].items()}
        blocks = [a for k, a in p.events if k == "block"]
        assert blocks, "no blocks"
        for a in blocks:
            mols = [x for x in (a.get("data-molecules") or "").split(",") if x]
            rs = [x for x in (a.get("data-roles") or "").split(",") if x]
            assert len(mols) == len(rs) and all(roles[m] == r for m, r in zip(mols, rs)), a
    elif what == "act_order":
        acts = [a.get("data-act") for k, a in p.events if k == "act"]
        assert acts == ["A", "B", "C", "D"], acts
        seen_predict = None
        for k, a in p.events:
            if k == "act":
                seen_predict = False
            elif k == "predict":
                seen_predict = True
            elif k == "block":
                assert seen_predict, "a view before its act's prediction"
    elif what == "footer":
        f = [a for k, a in p.events if k == "page-provenance"]
        assert len(f) == 1
        f = f[0]
        src = _sources(c)
        assert f.get("data-atlas") == src["bundle"].atlas_version
        assert f.get("data-results-sha256") == src["study"].sha256
        assert f.get("data-espmap") == "0.2.1"
        keys = {v.analysis_key for _, vs in src["ensembles"].values() for v in vs}
        assert set((f.get("data-analysis-keys") or "").split(",")) == keys
    elif what == "study_words_in_study_cards":
        text = p.own({"verbatim", "study-card", "evidence-card", "contrast-card", "class-panel"})
        assert not re.search(STUDY_WORDS, text, re.I), re.search(STUDY_WORDS, text, re.I)
    elif what == "no_bridge":
        text = p.own({"verbatim", "block"})
        assert not re.search(BRIDGE, text, re.I), re.search(BRIDGE, text, re.I)
    elif what == "no_blank_text":
        # every verbatim span holds something to read: a letter, number, punctuation mark or symbol
        # that is neither default-ignorable (Unicode DerivedCoreProperties) nor the braille blank
        for span in re.findall(r'<span class="verbatim">(.*?)</span>', page, re.S):
            s = html.unescape(span)
            assert any(unicodedata.category(ch)[0] in "LNPS" and not _ignorable(ch) for ch in s), \
                f"blank plan text shown: {s!r}"
        assert not re.search(r"<h2>Act [A-D]: </h2>", page), "an act heading with nothing after the colon"
    elif what == "page_name":
        names = re.findall(r"<title>(.*?)</title>", page, re.S)
        assert len(names) == 1, names
        s = html.unescape(names[0])
        assert s != "None" and any(unicodedata.category(ch)[0] in "LNPS" and not _ignorable(ch) for ch in s), \
            f"the page is named {s!r}"
    else:
        raise ValueError(what)


REQ = [(c, i) for c in CASES for i in range(len(c["requires"]))]


@pytest.mark.parametrize("c,i", REQ, ids=[f"{c['id']}-{i}" for c, i in REQ])
def test_a_rendered_page_meets_its_requirements(c, i, request):
    _xfail(request, f"test_a_rendered_page_meets_its_requirements[{c['id']}-{i}]")
    out, page = _run(c)
    if out != "RENDERED":
        return
    q = c["requires"][i]
    if "check" in q:
        _check(q["check"], page, c)
        return
    text = html.unescape(re.sub(r"<[^>]+>", " ", re.sub(r"</?span[^>]*>", "", page)))
    if "must" in q:
        assert re.search(q["must"], text, re.I), f"{q['why']} ({q['grade']} if missing)"
    else:
        assert not re.search(q["must_not"], text, re.I), f"{q['why']} ({q['grade']} if present)"
