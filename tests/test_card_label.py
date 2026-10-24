"""The difference card's label (tests/pedagogy/card_label_set.json), written
before the fix.

The set is what its generator writes; `allowed` is exactly the safe
outputs; every class has two or more variants; the extracted site values
reproduce the verdict each case measured; and the study says what the cases
say. The card is measured as it is: until card_html takes a context it is
called without one, and what it gets wrong is listed below as strict
xfails, which only shrink.
"""
import inspect
import json
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
PED = Path(__file__).with_name("pedagogy")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PED))
from model import (Atom, MoleculeView, Provenance, Refused, Scale, SiteMeasurement,  # noqa: E402
                   compare)
from record import load_all  # noqa: E402

SET = json.loads((PED / "card_label_set.json").read_text(encoding="utf-8"))
CASES = SET["cases"]
SITES = json.loads((HERE / "tests" / "data" / "card_label_sites.json").read_text(encoding="utf-8"))
RESULTS = HERE / "tests" / "data" / "espmap-sesiidae-001.results.json"

# --- what the card gets wrong today (measured on the committed card) --------------------
# Empty since the label fix (the 53 pairs #15 measured on the card before it). A case added
# later that the card fails goes here, as "test[case]" or "test[case-i]", with why.
XFAIL = set()


def _xfail(request, key):
    if key in XFAIL:
        request.applymarker(pytest.mark.xfail(strict=True, reason="the card before its label fix"))


def test_the_committed_set_is_what_the_generator_writes():
    import make_card_label_set
    assert make_card_label_set.document() == (PED / "card_label_set.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_allowed_is_exactly_the_safe_outputs(c):
    assert c["allowed"] == sorted(k for k, g in c["outcomes"].items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert c["prefer"] is None or c["prefer"] in c["allowed"]
    assert c["why"] and c["outcomes"]["ERROR"] == "CRASH"


def test_every_class_has_two_or_more_variants():
    per = Counter(c["cls"].split()[0] for c in CASES)
    assert sorted(per) == ["L1", "L2", "L3", "L4"] and min(per.values()) >= 2, per


# --- inputs ------------------------------------------------------------------------------
def _views(slug):
    if slug.startswith("pilot:"):
        return load_all(HERE / "records", slug.split(":")[1])
    out = []
    for k, m in enumerate(SITES["molecules"][slug]["members"]):
        sites = {s: {rule: SiteMeasurement(s, (0,), rule, v[key], 10.0, 5.0, 5.0, v["partition"])
                     for rule, key in (("nearest_nucleus", "nn"), ("radius_scaled", "rs"))}
                 for s, v in m["sites"].items()}
        out.append(MoleculeView(
            name=f"{slug} {m['conformer']}", symbols=("C",), xyz=np.zeros((1, 3)), bonds=None, charge=0,
            atoms=(Atom(0, "C", 1),), sites=sites, sites_tile_surface=False, genus=0,
            scale=Scale("diverging", -20.0, 20.0, "vik"),
            provenance=Provenance(m["conformer"], m["analysis_key"], "r"), folded=m["folded"]))
    return out


def _context(c):
    spec = c["context"]
    if not isinstance(spec, dict):
        return spec
    ctx = dict(spec)
    if spec["study"] == "results":
        from study import load_results
        ctx["study"] = load_results(RESULTS)
    elif spec["study"] == "raw":
        ctx["study"] = json.loads(RESULTS.read_text(encoding="utf-8"))
    return ctx


def _cmp(c):
    return compare(_views(c["x"]), _views(c["y"]), c["site"], c["x"], c["y"])


# --- the facts the cases rest on ---------------------------------------------------------
@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_the_measured_verdict_reproduces(c):
    assert _cmp(c).verdict.name == c["measured"]


def test_the_extract_is_whole():
    for slug, m in SITES["molecules"].items():
        assert len(m["members"]) == 20, slug
        assert len({x["analysis_key"] for x in m["members"]}) == 1, slug
        assert all(isinstance(x["folded"], bool) and len(x["record_sha256"]) == 64 for x in m["members"])
    assert SITES["molecules"]["z3_z13_18_al"]["atlas_id"].startswith("arrowfield.control.")


def test_the_study_says_what_the_cases_say():
    from study import load_results
    r = load_results(RESULTS)
    assert r.site("B_stereo", "C3-H").outcome == "different on average" and r.site("B_stereo", "C3-H").tier == "primary"
    assert r.site("C_position", "C4-H").outcome == "resolved"
    assert r.site("B_stereo", "C5-H").tier == "exploratory"
    h = r.site("A_headgroup (control: Al vs OH)", "head O")
    assert (h.outcome, h.tier, h.x, h.y) == ("unresolved", "secondary", "z3_z13_18_oh", "z3_z13_18_al")
    for k in ("A_headgroup", "A_headgroup (control: Al vs OAc)", "A_headgroup (control: Al vs OH)"):
        assert not [s for s in r.sites if s.contrast == k and re.fullmatch(r"C\d+-H", s.site)], k
    b = r.site("B_stereo", "C3-H")
    assert (b.x, b.y) == ("z3_z13_18_oac", "e3_z13_18_oac")


# --- the card ----------------------------------------------------------------------------
def _run(c):
    from views import card_html
    takes_context = "context" in inspect.signature(card_html).parameters
    try:
        cmp = _cmp(c)
        page = card_html(cmp, context=_context(c)) if takes_context else card_html(cmp)
    except Refused:
        return "REFUSED", None
    except Exception:  # noqa: BLE001 - any other error is a crash, graded as such
        return "ERROR", None
    return "RENDERED", page


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_card_output_is_safe(c, request):
    _xfail(request, f"test_card_output_is_safe[{c['id']}]")
    out = _run(c)[0]
    assert out in c["allowed"], f"{out}: {c['outcomes'][out]}"


@pytest.mark.parametrize("c", [c for c in CASES if c["prefer"]], ids=lambda c: c["id"])
def test_card_output_is_preferred(c, request):
    _xfail(request, f"test_card_output_is_preferred[{c['id']}]")
    assert _run(c)[0] == c["prefer"]


REQ = [(c, i) for c in CASES for i in range(len(c["requires"]))]


@pytest.mark.parametrize("c,i", REQ, ids=[f"{c['id']}-{i}" for c, i in REQ])
def test_a_rendered_card_meets_its_requirements(c, i, request):
    _xfail(request, f"test_a_rendered_card_meets_its_requirements[{c['id']}-{i}]")
    out, page = _run(c)
    if out != "RENDERED":
        return
    q = c["requires"][i]
    text = re.sub(r"<[^>]+>", " ", page)
    if "must" in q:
        assert re.search(q["must"], text, re.I), f"{q['why']} ({q['grade']} if missing)"
    else:
        assert not re.search(q["must_not"], text, re.I), f"{q['why']} ({q['grade']} if present)"
