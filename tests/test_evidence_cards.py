"""The evidence views' hostile set (tests/pedagogy/evidence_card_set.json),
written before the views exist.

Now: the set is what its generator writes; every case's `allowed` is exactly
its safe outputs; all seven classes have two or more variants; and the
fixture facts each case's reason cites are true of the committed bundle (a
case that misdescribes its record would test the wrong thing).
When views.evidence_card_html, contrast_card_html and evidence_class_panel_html
exist: every output is safe and preferred, and every requirement holds.
Until then those tests skip. The views' PR may add strict-xfail lists here,
never edit an expectation. These cases supersede hostile_set.json's B1-B5.

The view is handed what the case builds: the bundle's own statement and
contrast, or, with `statement_patch`, ones built from a patched copy while
the view gets the unpatched bundle (E7, E8). The harness never refuses on
the view's behalf. XFAIL lists the (test, case) pairs the views did not yet
meet when the trust-boundary cases were added; it only shrinks, and a pair
that starts passing fails the run until it is removed.
"""
import copy
import html
import json
import re
import sys
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
PED = Path(__file__).with_name("pedagogy")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PED))
from bundle import InvalidBundle, load  # noqa: E402

SET = json.loads((PED / "evidence_card_set.json").read_text(encoding="utf-8"))
CASES = SET["cases"]
BUNDLE = HERE / SET["bundle"]
RAW = json.loads(BUNDLE.read_text(encoding="utf-8"))
VIEWS = {"evidence_card": "evidence_card_html", "contrast_card": "contrast_card_html",
         "class_panel": "evidence_class_panel_html"}


def test_the_committed_set_is_what_the_generator_writes():
    import make_evidence_card_set
    assert make_evidence_card_set.document() == (PED / "evidence_card_set.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_allowed_is_exactly_the_safe_outputs(c):
    assert c["allowed"] == sorted(k for k, g in c["outcomes"].items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert c["prefer"] is None or c["prefer"] in c["allowed"]
    assert c["why"] and c["outcomes"]["ERROR"] == "CRASH" and c["component"] in VIEWS


def test_every_class_has_two_or_more_variants():
    per = Counter(c["cls"].split()[0] for c in CASES)
    assert sorted(per) == sorted(f"E{i}" for i in range(1, 13)) and min(per.values()) >= 2, per


# --- the facts the cases rest on, checked against the committed bundle --------------------
def _obs(oid):
    return next(o for o in RAW["observations"] if o["observation_id"] == "spa.obs." + oid)


def _classes(comp, tax):
    c = next(c for c in RAW["compounds"] if c["atlas_id"] == comp)
    return {e["evidence_class"]: e for e in c["evidence_classes"] if e["taxon_id"] == tax}


def test_the_fixture_says_what_the_cases_say():
    oac, oh, e3 = "spa.cmpd.z3_z13_18_oac", "spa.cmpd.z3_z13_18_oh", "spa.cmpd.e3_z13_18_oac"
    tenuis, tipu = "spa.tax.synanthedon_tenuis", "spa.tax.synanthedon_tipuliformis"
    y = next(x for x in _obs("synanthedon_tenuis.yang_2012.gland")["components"] if x["atlas_id"] == oac)
    assert y["detected"] is False and y["detection_limit"] == {"unit": "picogram", "value": 100}
    n = next(x for x in _obs("synanthedon_tenuis.naka_2013.gland")["components"] if x["atlas_id"] == oh)
    assert n["detected"] is False and n["detection_limit"] is None
    t87 = _obs("synanthedon_tipuliformis.szocs_1990.trap_e3.tasmania_1987")
    assert t87["outcome"]["effect"] == "no_effect" and t87["outcome"]["significance"] == "not_significant"
    assert "Slateford" in t87["context"]["locality"]
    assert oh not in {b["atlas_id"] for b in t87["treatment"] + t87["comparator"]["blend"]}
    nz = _obs("synanthedon_tipuliformis.suckling_2005.trap_e3.new_zealand")["outcome"]
    assert (nz["effect"], nz["significance"]) == ("no_effect", "not_reported")
    r = _obs("synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1")["outcome"]
    assert (r["effect"], r["significance"]) == ("decreased", "significant")
    assert (r["response"]["treatment_value"], r["response"]["comparator_value"]) == (4.6, 17.0)
    assert "detection limit for the negatives" in _obs("synanthedon_tenuis.naka_2013.gland")["evidence"]["unresolved"]
    k = _classes(oac, tenuis)
    assert {"female_pheromone_component", "female_not_detected", "addition_decreased"} <= set(k)
    assert k["female_pheromone_component"]["population_match"] == "yes"
    assert k["female_not_detected"]["populations"] == ["Korean population"]
    obs = {o["observation_id"] for o in RAW["observations"]}
    fpc = k["female_pheromone_component"]
    assert fpc["n_observations"] == 8 and len(set(fpc["observation_ids"]) - obs) == 7
    k3 = _classes(e3, tipu)
    assert {"addition_increased", "addition_decreased"} <= set(k3)
    assert any("Tasmania" in p for p in k3["addition_decreased"]["populations"])
    assert any("Hungary" in p for p in k3["addition_increased"]["populations"])
    assert _classes(oh, tipu) == {}


# --- the views (to come) ---------------------------------------------------------------
def _view(component):
    import views
    name = VIEWS[component]
    if not hasattr(views, name):
        pytest.skip(f"views.{name} is not built yet; its expectations are fixed in evidence_card_set.json")
    return getattr(views, name)


def _patched(raw, patch, path_):
    raw = copy.deepcopy(raw)
    for path, value in patch:
        node = raw
        for key in path[:-1]:
            node = node[key]
        if path[-1] == "+":
            node.append(value)
        else:
            node[path[-1]] = value
    path_.write_text(json.dumps(raw), encoding="utf-8")
    return load(path_), raw


def _bundle(c, tmp_path):
    if not c["bundle_patch"]:
        return load(BUNDLE), RAW
    return _patched(RAW, c["bundle_patch"], tmp_path / "patched.json")


XFAIL = set()


def _xfail(request, c):
    if f"{request.function.__name__}[{c['id']}]" in XFAIL:
        request.applymarker(pytest.mark.xfail(strict=True, reason="trust-boundary case not yet met"))


class _Verbatim(HTMLParser):
    """Splits a page's text into the view's own words and the verbatim record
    text (elements with class "verbatim", and everything inside them)."""
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
    p = _Verbatim()
    p.feed(page)
    return " ".join(p.own), [x.strip() for x in p.spans]


def _strings(node, out):
    if isinstance(node, str):
        out.append(node)
    elif isinstance(node, dict):
        for v in node.values():
            _strings(v, out)
    elif isinstance(node, list):
        for v in node:
            _strings(v, out)
    return out


def _run(c, tmp_path):
    view = _view(c["component"])
    from model import Refused
    try:
        b, raw = _bundle(c, tmp_path)
    except InvalidBundle:            # the adapter refused the bundle: nothing is shown (E11)
        return "REFUSED", None
    src = _patched(raw, c["statement_patch"], tmp_path / "inputs.json")[0] if c["statement_patch"] else b
    c["_raw"] = raw
    ctx = dict(c["context"], bundle_object=b)
    if isinstance(ctx.get("bundle"), dict) and ctx["bundle"].get("sha256") == "from bundle.load":
        ctx["bundle"] = dict(ctx["bundle"], sha256=b.sha256)
    i = c["inputs"]
    # an input the case cannot build is the case's error, not a refusal
    if c["component"] == "evidence_card":
        args = (src.observations[i["observation_id"]], i["compound"], ctx)
    elif c["component"] == "contrast_card":
        k = next(k for k in src.contrasts if k.label == i["contrast_label"])
        args = (k, [src.observations[o] for o in k.observation_ids], ctx)
    else:
        args = (i["compound"], i["taxon_id"], ctx)
    try:
        page = view(*args)
    except Refused:
        return "REFUSED", None
    except Exception:  # noqa: BLE001 - any other error is a crash, graded as such
        return "ERROR", None
    return "RENDERED", page


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_view_output_is_safe(c, tmp_path, request):
    _xfail(request, c)
    out, _ = _run(c, tmp_path)
    assert out in c["allowed"], f"{out}: {c['outcomes'][out]}"


@pytest.mark.parametrize("c", [c for c in CASES if c["prefer"]], ids=lambda c: c["id"])
def test_view_output_is_preferred(c, tmp_path, request):
    _xfail(request, c)
    assert _run(c, tmp_path)[0] == c["prefer"]


REQ = [(c, i) for c in CASES for i in range(len(c["requires"]))]


@pytest.mark.parametrize("c,i", REQ, ids=[f"{c['id']}-{i}" for c, i in REQ])
def test_a_rendered_view_meets_its_requirements(c, i, tmp_path, request):
    if f"test_a_rendered_view_meets_its_requirements[{c['id']}-{i}]" in XFAIL:
        request.applymarker(pytest.mark.xfail(strict=True, reason="trust-boundary case not yet met"))
    out, page = _run(c, tmp_path)
    if out != "RENDERED":
        return
    q = c["requires"][i]
    own, _ = _split(page)
    if "must" in q:          # what must be shown may be shown in the record's words
        # inline spans add no space (as in a browser); block tags separate words
        text = html.unescape(re.sub(r"<[^>]+>", " ", re.sub(r"</?span[^>]*>", "", page)))
        target = page if q["on"] == "html" else text
        assert re.search(q["must"], target, re.I), f"{q['why']} ({q['grade']} if missing)"
    else:                    # a ban reads the view's own words; verbatim record text is exempt
        target = page if q["on"] == "html" else own
        assert not re.search(q["must_not"], target, re.I), f"{q['why']} ({q['grade']} if present)"


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_verbatim_spans_are_record_text(c, tmp_path, request):
    """The exemption cannot hide the view's own words: every verbatim span is
    text the record carries."""
    _xfail(request, c)
    out, page = _run(c, tmp_path)
    if out != "RENDERED":
        return
    record = _strings(c["_raw"], []) + _strings(c["statement_patch"], [])
    for span in _split(page)[1]:
        assert any(span in r for r in record), f"verbatim span not found in the record: {span!r}"


def test_bans_catch_word_families():
    """A ban is a rule too: its word families on one side, near-misses on the other."""
    import make_evidence_card_set as g
    over, esp = g.NO_OVERREACH["must_not"], g.NO_ESP["must_not"]
    for w in ("inhibit", "inhibited", "inhibitor", "Inhibition", "repel", "repellent", "antagonist",
              "antagonism", "inactive", "inactivity", "not produced"):
        assert re.search(over, w, re.I), w
    for w in ("resolved", "Resolves", "resolvable", "resolution", "ESP", "electrostatic", "because"):
        assert re.search(esp, w, re.I), w
    for w in ("unresolved", "Unresolved in the record", "uninhabited", "rebel"):
        assert not re.search(over, w, re.I) and not re.search(esp, w, re.I), w


def test_the_views_phrase_exactly_the_atlas_vocabularies():
    """Every controlled word bundle.load admits has a phrase in the views, and
    no phrase exists for a word it refuses: the two lists cannot drift."""
    import model
    import views
    assert set(views._EFFECT_WORDS) == set(model.EFFECTS)
    assert set(views._SIG_WORDS) == set(model.SIGNIFICANCES)
    assert set(views._ROLE_PHRASES) == set(model.ROLES)
    assert set(views._OBS_TYPE_WORDS) == set(model.OBSERVATION_TYPES)


HAND = [("detections", "one Detection, not a list"), ("roles", "a bare role string"),
        ("unresolved", "a bare string"), ("observation_type", "repellency_assay"),
        ("effect", "inhibited"), ("significance", "resolved"), ("limit", "0.5 as text"),
        ("limit_bool", "true"), ("response", "a bool value"), ("class_name", "antagonist_effect"),
        ("count_fraction", "3.5"), ("count_ids", "unequal to ids"), ("populations", "a bare string")]


@pytest.mark.parametrize("what,why", HAND, ids=[h[0] for h in HAND])
def test_a_hand_built_bundle_meets_the_same_validator(what, why):
    """The shared validator guards a bundle built without bundle.load: each
    broken field is a bundle_problem, and every view refuses it."""
    from dataclasses import replace
    import model
    import views
    b = load(BUNDLE)
    gy = "spa.obs.synanthedon_tenuis.yang_2012.gland"
    hu = "spa.obs.synanthedon_tipuliformis.szocs_1990.trap_e3.hungary"
    oac, e3, tipu = "spa.cmpd.z3_z13_18_oac", "spa.cmpd.e3_z13_18_oac", "spa.tax.synanthedon_tipuliformis"
    obs, comps = dict(b.observations), dict(b.compounds)
    st, oid, comp = obs[hu], hu, e3
    if what == "detections":
        st, oid, comp = replace(obs[gy], detections=obs[gy].detections[0]), gy, oac
    elif what == "roles":
        st = replace(st, roles={e3: "treatment"})
    elif what == "unresolved":
        st, oid, comp = replace(obs[gy], unresolved="detection limit"), gy, oac
    elif what in ("observation_type", "effect", "significance"):
        st = replace(st, **{what: why})
    elif what in ("limit", "limit_bool"):
        d = obs[gy].detections
        bad = replace(d[1], detection_limit={"value": "0.5" if what == "limit" else True, "unit": "picogram"})
        st, oid, comp = replace(obs[gy], detections=(d[0], bad) + d[2:]), gy, oac
    elif what == "response":
        st = replace(st, response=replace(st.response, treatment_value=True))
    else:
        c = comps[e3]
        k = next(i for i, e in enumerate(c.evidence_classes) if e.taxon_id == tipu)
        e = c.evidence_classes[k]
        e = {"class_name": replace(e, evidence_class="antagonist_effect"),
             "count_fraction": replace(e, n_observations=3.5),
             "count_ids": replace(e, n_observations=e.n_observations + 2),
             "populations": replace(e, populations="Tasmania")}[what]
        comps[e3] = replace(c, evidence_classes=c.evidence_classes[:k] + (e,) + c.evidence_classes[k + 1:])
    obs[oid] = st
    hb = replace(b, observations=obs, compounds=comps)
    assert model.bundle_problems(hb), why
    ctx = {"bundle_object": hb, "roles": {}}
    with pytest.raises(model.Refused):
        views.evidence_card_html(hb.observations[oid], comp, ctx)
    with pytest.raises(model.Refused):
        views.evidence_class_panel_html(e3, tipu, ctx)
