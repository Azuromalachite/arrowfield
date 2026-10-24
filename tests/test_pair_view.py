"""The pair view's hostile set (tests/pedagogy/pair_view_set.json), written
before the view exists.

Now: the set is what its generator writes; every case's `allowed` is exactly
its safe outputs; all six classes have two or more variants; the fixtures in
tests/data/pair are the Atlas's records, byte for byte; and the facts each
case's reason cites are true of the records (a case that misdescribes its
input would test the wrong thing).
When model.relation and views.pair_view_html exist: every output is safe and
preferred, and every requirement holds. Until then each one's tests skip. The view's PR may add
strict-xfail lists here, never edit an expectation.
"""
import hashlib
import html
import json
import math
import re
import sys
from dataclasses import replace
from html.parser import HTMLParser
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
PED = Path(__file__).with_name("pedagogy")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PED))
import record  # noqa: E402

SET = json.loads((PED / "pair_view_set.json").read_text(encoding="utf-8"))
CASES = SET["cases"]
PAIR = HERE / "tests" / "data" / "pair"


def test_the_committed_set_is_what_the_generator_writes():
    import make_pair_view_set
    assert make_pair_view_set.document() == (PED / "pair_view_set.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_allowed_is_exactly_the_safe_outputs(c):
    assert c["allowed"] == sorted(k for k, g in c["outcomes"].items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert c["prefer"] is None or c["prefer"] in c["allowed"]
    assert c["why"] and c["outcomes"]["ERROR"] == "CRASH"


def test_every_class_has_two_or_more_variants():
    from collections import Counter
    per = Counter(c["cls"].split()[0] for c in CASES)
    assert sorted(per) == [f"P{i}" for i in range(1, 10)] + ["R1", "R2"] and min(per.values()) >= 2, per


def test_retired_cases_name_existing_replacements():
    ids = {c["id"] for c in CASES}
    old = {c["id"] for c in json.loads((PED / "hostile_set.json").read_text(encoding="utf-8"))["cases"]
           if c["component"] == "pair_view"}
    assert set(SET["retires"]) == old
    assert all(r in ids for rs in SET["retires"].values() for r in rs)


# --- the fixtures, and the facts the cases rest on -------------------------------------------
def _sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_pair_fixtures_are_the_atlas_records():
    members = json.loads((PAIR / "members.json").read_text(encoding="utf-8"))["members"]
    assert len(members) == 4
    for m in members:
        j = PAIR / m["record"]
        assert _sha(j) == m["record_sha256"] and _sha(j.with_suffix(".npz")) == m["npz_sha256"], m["record"]
        z = np.load(j.with_suffix(".npz"))
        r = json.loads(j.read_text(encoding="utf-8"))
        assert str(z["analysis_key"]) == r["analysis_key"] and str(z["render_key"]) == r["render_key"]


def _load(stem):
    return record.load(HERE / stem)


def _dist(v):
    return np.linalg.norm(v.xyz[:, None] - v.xyz[None], axis=-1)


def test_the_records_say_what_the_cases_say():
    import make_pair_view_set as g
    v = {s: _load(s) for s in (g.ETOH, g.ETOH2, g.ACAL, g.ACAL2, g.E2B, g.E2B2, g.Z2B, g.ACO,
                               g.E3F, g.E3X, g.Z3X, g.Z3F)}
    # P1: same method, different own scales
    assert v[g.E2B].analysis_key == v[g.Z2B].analysis_key and v[g.E2B].render_key != v[g.Z2B].render_key
    assert round(v[g.E2B].scale.vmax, 1) == 17.3 and round(v[g.Z2B].scale.vmax, 1) == 18.1
    assert round(v[g.ETOH].scale.vmax, 1) == 50.0 and round(v[g.ACAL].scale.vmax, 1) == 39.3
    assert v[g.ETOH].analysis_key == v[g.ACAL].analysis_key
    assert (round(v[g.E3F].scale.vmax, 2), round(v[g.E3X].scale.vmax, 2)) == (40.33, 40.77)
    assert v[g.E2B].render_key == v[g.E2B2].render_key
    # P2: the keys that differ, and the ones that do not
    assert v[g.ACO].analysis_key != v[g.ETOH].analysis_key and v[g.ACO].charge == -1
    assert v[g.E3X].analysis_key != v[g.ETOH].analysis_key
    assert v[g.E3X].provenance.method == v[g.ETOH].provenance.method == "b3lyp/def2-svp, gas"
    # P3: same geometry three ways, and true conformers
    assert np.array_equal(v[g.E2B].xyz, v[g.E2B2].xyz)
    assert np.abs(_dist(v[g.ACAL]) - _dist(v[g.ACAL2])).max() < 1e-5
    assert np.abs(v[g.ACAL].xyz - v[g.ACAL2].xyz).max() > 2.0
    assert round(float(np.abs(_dist(v[g.ETOH]) - _dist(v[g.ETOH2])).max()), 2) == 0.70
    assert np.abs(_dist(v[g.E3F]) - _dist(v[g.E3X])).max() > 1.0
    # P4: the producer's flags, and genus that does not follow them
    folded = {m["record"][:-5]: m["folded"] for m in
              json.loads((PAIR / "members.json").read_text(encoding="utf-8"))["members"]}
    assert folded == {"e3_z13_18_oac/conf_00": True, "e3_z13_18_oac/conf_02": False,
                      "z3_z13_18_oac/conf_00": False, "z3_z13_18_oac/conf_01": True}
    assert v[g.E3F].genus == 0 and v[g.Z3F].genus == 3
    assert all(v[s].folded is None for s in (g.E2B, g.Z2B))
    # the loader names both compounds' conformers alike
    assert v[g.E3X].name.startswith("conf_0") and v[g.Z3X].name.startswith("conf_0")


def _mol(v):
    import views
    return views.mol_from_view(v, v.mapped_smiles)


def _rms(a, b):
    """Symmetric RMSD after superposition, all atoms (2000 matches: an upper bound)."""
    from rdkit.Chem import rdMolAlign
    return rdMolAlign.GetBestRMS(_mol(a), _mol(b), maxMatches=2000)


def test_the_same_geometry_threshold_sits_where_it_was_measured():
    """The docstring's measurement, re-checked on the committed records: the
    equivalent seed pairs sit far below SAME_GEOMETRY_RMSD, the true conformers
    far above it, and the synthetic displacements on their stated sides."""
    import make_pair_view_set as g
    t = SET["same_geometry_rmsd"]
    assert t == g.SAME_GEOMETRY_RMSD == 0.01
    for a, b in ((g.E2B, g.E2B2), (g.ACAL, g.ACAL2), (g.ETOH, g.ETOH2)):
        assert _rms(_load(a), _load(b)) < 1e-5 < t
    assert _rms(_load(g.E3F), _load(g.E3X)) > 2.0 > 100 * t
    e3 = _load(g.E3X)
    assert _rms(e3, _side(g.side(g.E3X, geometry=["rmsd", 0.005, 1]))) <= 0.005 < t
    assert _rms(e3, _side(g.side(g.E3X, geometry=["rmsd", 0.04, 2]))) > 3 * t


def test_model_threshold_is_the_frozen_one():
    import model
    if not hasattr(model, "SAME_GEOMETRY_RMSD"):
        pytest.skip("model.SAME_GEOMETRY_RMSD is not defined yet")
    assert model.SAME_GEOMETRY_RMSD == SET["same_geometry_rmsd"]


# --- the view (to come) --------------------------------------------------------------------
def _view(component):
    import model
    import views
    mod, name = (views, "pair_view_html") if component == "pair_view" else (model, "relation")
    if not hasattr(mod, name):
        pytest.skip(f"{mod.__name__}.{name} is not built yet; its expectations are fixed in pair_view_set.json")
    return getattr(mod, name)


def _side(s):
    from model import Surface
    v = _load(s["record"])
    for field, value in s["patch"]:
        if "." in field:
            outer, inner = field.split(".")
            v = replace(v, **{outer: replace(getattr(v, outer), **{inner: value})})
        else:
            v = replace(v, **{field: value})
    geo = s.get("geometry")
    if geo == "zeros":
        v = replace(v, xyz=np.zeros_like(v.xyz))
    elif geo is not None:
        _, target, seed = geo
        d = np.random.default_rng(seed).standard_normal(v.xyz.shape)
        d -= d.mean(0)
        d *= target / math.sqrt(float((d ** 2).sum(1).mean()))
        v = replace(v, xyz=v.xyz + d)
    op = s["surface"]
    if op == "none":
        v = replace(v, surface_loader=None)
    elif op is not None:
        base = v.surface_loader

        def loader(rule, base=base, op=op):
            x = base(rule)
            V, faces = x.V.copy(), x.faces
            if op == "nan":
                V[0] = math.nan
            elif op == "inf":
                V[0] = math.inf
            elif op == "empty":
                faces = np.zeros((0, 3), dtype=x.faces.dtype)
            return Surface(x.verts, faces, V, x.owner[:len(faces)])
        v = replace(v, surface_loader=loader)
    return v


def _run(c):
    view = _view(c["component"])
    from model import Refused
    left = _side(c["left"])
    right = left if c["same_object"] else _side(c["right"])
    try:
        if c["component"] == "relation":
            return view(left, right).value.upper(), None, (left, right)
        page = view(left, right, dict(c["context"]))
    except Refused:
        return "REFUSED", None, (left, right)
    except Exception:  # noqa: BLE001 - any other error is a crash, graded as such
        return "ERROR", None, (left, right)
    return "RENDERED", page, (left, right)


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_view_output_is_safe(c):
    out = _run(c)[0]
    assert out in c["allowed"], f"{out}: {c['outcomes'][out]}"


@pytest.mark.parametrize("c", [c for c in CASES if c["prefer"]], ids=lambda c: c["id"])
def test_view_output_is_preferred(c):
    assert _run(c)[0] == c["prefer"]


class _Attrs(HTMLParser):
    def __init__(self):
        super().__init__()
        self.found = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if any(k.startswith("data-") for k in a):
            self.found.append(a)


def _check(what, page, sides):
    p = _Attrs()
    p.feed(page)
    classes = [(a.get("class") or "").split() for a in p.found]
    scale = [a for a, k in zip(p.found, classes) if "pair-scale" in k]
    panels = [a for a, k in zip(p.found, classes) if "pair-panel" in k]
    if what == "shared_scale":
        # the frozen rule (make_pair_view_set's docstring), exactly
        assert len(scale) == 1 and len(panels) == 2, "one pair-scale and two pair-panels"
        lo = min(v.scale.vmin for v in sides)
        hi = max(v.scale.vmax for v in sides)
        if any(v.scale.kind == "diverging" for v in sides) or lo < 0 < hi:
            m = max(abs(lo), abs(hi))
            kind, want = "diverging", (-m, m)
        else:
            kind, want = "sequential", (lo, hi)
        s = scale[0]
        assert s.get("data-unit") == "kcal/mol"
        assert s["data-kind"] == kind and s["data-cmap"] == SET["cmap"][kind], (s, kind)
        got = (float(s["data-vmin"]), float(s["data-vmax"]))
        assert all(math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-9) for x, y in zip(got, want)), (got, want)
        for q in panels:          # each picture's effective scale is the shared one
            assert q.get("data-scale") == "shared" and q.get("data-cmap") == s["data-cmap"]
            assert (float(q["data-vmin"]), float(q["data-vmax"])) == got
        return
    if what.startswith("relation="):
        rel = [a["data-relation"] for a in p.found if "data-relation" in a]
        assert rel == [what.split("=")[1]], rel
        return
    want = what.split("=")[1].split(",")
    assert [q.get("data-folded") for q in panels] == want


REQ = [(c, i) for c in CASES for i in range(len(c["requires"]))]


@pytest.mark.parametrize("c,i", REQ, ids=[f"{c['id']}-{i}" for c, i in REQ])
# (relation cases carry no requirements: their output is the judgement itself)
def test_a_rendered_view_meets_its_requirements(c, i):
    out, page, sides = _run(c)
    if out != "RENDERED":
        return
    q = c["requires"][i]
    if "check" in q:
        _check(q["check"], page, sides)
        return
    if "must" in q:          # inline spans add no space; tags go before entities are unescaped
        text = html.unescape(re.sub(r"<[^>]+>", " ", re.sub(r"</?span[^>]*>", "", page)))
        target = page if q["on"] == "html" else text
        assert re.search(q["must"], target, re.I), f"{q['why']} ({q['grade']} if missing)"
    else:                    # a ban reads the view's own words; verbatim text is exempt (the amendment)
        target = page if q["on"] == "html" else _split(page)[0]
        assert not re.search(q["must_not"], target, re.I), f"{q['why']} ({q['grade']} if present)"


@pytest.mark.parametrize("c", [c for c in CASES if c["component"] == "pair_view"], ids=lambda c: c["id"])
def test_verbatim_spans_are_names_sources_or_labels(c):
    """The exemption covers only text the view did not write: each span is
    a molecule's name or record source, or a caller's label, exactly."""
    out, page, sides = _run(c)
    if out != "RENDERED":
        return
    labels = (c["context"].get("labels") or {}) if isinstance(c["context"].get("labels"), dict) else {}
    allowed = {v.name for v in sides} | {v.provenance.source for v in sides} | set(labels.values())
    for span in _split(page)[1]:
        assert span in allowed, f"verbatim span is not a name, record source or label: {span!r}"


def test_bans_catch_word_families():
    import make_pair_view_set as g
    ban = g.NO_CAUSE["must_not"]
    for w in ("explains", "explained", "causes", "causal", "attractive", "attraction", "binds", "binding",
              "pheromonal", "receptors", "reactivity", "resolved", "resolution", "specificity", "because"):
        assert re.search(ban, w, re.I), w
    for w in ("unresolved", "specified", "rebinding"):
        assert not re.search(ban, w, re.I), w


class _Verbatim(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.depth, self.stack, self.own, self.spans = 0, [], [], []

    def handle_starttag(self, tag, attrs):
        if tag in ("br", "img", "hr", "meta", "input"):
            return
        v = "verbatim" in (dict(attrs).get("class") or "").split()
        self.stack.append(v)
        if v:
            self.depth += 1
            if self.depth == 1:
                self.spans.append("")

    def handle_endtag(self, tag):
        if tag in ("br", "img", "hr", "meta", "input") or not self.stack:
            return
        if self.stack.pop():
            self.depth -= 1

    def handle_data(self, d):
        if self.depth:
            self.spans[-1] += d
        else:
            self.own.append(d)


def _split(page):
    v = _Verbatim()
    v.feed(page)
    return " ".join(v.own), [x.strip() for x in v.spans]
