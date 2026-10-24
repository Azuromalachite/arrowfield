"""The hostile pedagogy set (tests/pedagogy/hostile_set.json) against the
views as they are.

Per case: safe (the output is in `allowed`), preferred (it equals `prefer`),
and what the card must show (`requires`); per card: the language rules.
Everything the views get wrong today is listed below as a strict xfail, so
CI is green now and a change that fixes a case must remove it from its list:
the lists only shrink. Expectations live in make_set.py, are reviewed before
any view changes to meet them, and never change in the same commit as a view.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
PED = Path(__file__).with_name("pedagogy")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PED))
from model import (Atom, ComparisonStatus, MoleculeView, Provenance, Scale,  # noqa: E402
                   Refused, SiteMeasurement, compare)
from views import card_html  # noqa: E402

SET = json.loads((PED / "hostile_set.json").read_text(encoding="utf-8"))
CASES = SET["cases"]
BUILT = {"compare"}                       # components that exist; the rest skip
NOW = [c for c in CASES if c["component"] in BUILT]

# --- what the views get wrong today (measured, not intended) ---------------
# Empty since the views PR that met the set (#2, with the probe cases of #3);
# a case added later that the views fail goes here, with one line saying why.
UNSAFE_NOW = set()
NOT_PREFERRED_NOW = set()
REQUIRES_MISSING_NOW = set()
SCOPE_CAVEAT_MISSING_NOW = set()


HAS_FOLDED = "folded" in MoleculeView.__dataclass_fields__


def views(members, name):
    out = []
    for k, m in enumerate(members):
        part = m.get("partition", None if m["nn"] is None else abs(m["nn"] - m["rs"]))
        sites = {m["site"]: {rule: SiteMeasurement(m["site"], (0,), rule, mean, 10.0, 5.0, 5.0, part)
                             for rule, mean in (("nearest_nucleus", m["nn"]), ("radius_scaled", m["rs"]))}}
        out.append(MoleculeView(
            name=f"{name}{k}", symbols=("C",), xyz=np.zeros((1, 3)), bonds=None, charge=0,
            atoms=(Atom(0, "C", 4),), sites=sites, sites_tile_surface=False, genus=m["genus"],
            scale=Scale("diverging", -20.0, 20.0, "bwr_r"),
            provenance=Provenance(f"{name}{k}", m["key"], "render0"),
            # A13: the producer's folded flag; a member that does not state it
            # is folded iff genus > 0 (the meaning of every case before A10)
            **({"folded": m.get("folded", m["genus"] > 0)} if HAS_FOLDED else {})))
    return out


def run(c):
    """(output label, Comparison or None, card html or None)."""
    label = c.get("label") or (c["a"] or c["b"])[0]["site"]
    try:
        cmp = compare(views(c["a"], "A"), views(c["b"], "B"), label, "A", "B")
    except (KeyError, Refused):          # a refusal that names its reason
        return "REFUSED", None, None
    except Exception:  # noqa: BLE001 - any other error is a crash, graded as such
        return "ERROR", None, None
    out = "NOT_COMPARED" if cmp.status is ComparisonStatus.NOT_COMPARABLE else cmp.verdict.name
    return out, cmp, card_html(cmp)


def shows(req, cmp, page, c=None):
    if req == "extended_n":
        # the conformers left out are exactly the folded ones: None when none is
        if cmp is None:
            return True
        want = c["extended_n"]
        got = None if cmp.extended is None else list(cmp.extended["n"])
        return got == want
    if page is None:
        return True                       # a refusal shows nothing to misread
    text = re.sub(r"<[^>]+>", " ", page)
    if req == "no_numbers":
        return cmp.status is not ComparisonStatus.NOT_COMPARABLE or not re.search(r"\d\.\d\d", text)
    if req == "no_single_range":
        return not re.search(r"(^|[^\d])0\.00 /|/ 0\.00", text)
    if req == "partition_shown":
        return "partition sensitivity" in text
    if req == "extended_row":
        return "extended conformers only" in text
    if req == "extended_verdict":
        return bool(re.search(r"extended conformers[^.]*:\s*(separated|overlapping|not attributable|cannot judge)",
                              text, re.I))
    if req == "no_identity_claim":
        return not re.search(r"identical|no difference|the same", text, re.I)
    if req == "robustness_wording":
        # A15: names what cannot be judged, and why (one geometry on a side)
        return bool(re.search(r"cannot judge conformational robustness", text, re.I)
                    and re.search(r"one geometry", text, re.I))
    raise ValueError(req)


def mark(ids, cid):
    return [pytest.mark.xfail(strict=True, reason="listed as failing now")] if cid in ids else []


def test_the_committed_set_is_what_the_generator_writes():
    import make_set
    assert make_set.document() == (PED / "hostile_set.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_allowed_is_exactly_the_safe_outcomes(c):
    safe = sorted(v for v, g in c["outcomes"].items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert c["allowed"] == safe and (c["prefer"] is None or c["prefer"] in c["allowed"])
    assert c["why"]


@pytest.mark.parametrize("c", [pytest.param(c, id=c["id"], marks=mark(UNSAFE_NOW, c["id"])) for c in NOW])
def test_output_is_safe(c):
    out, _, _ = run(c)
    assert out in c["allowed"], f"{out}: {c['outcomes'][out]}"


@pytest.mark.parametrize("c", [pytest.param(c, id=c["id"], marks=mark(NOT_PREFERRED_NOW, c["id"]))
                               for c in NOW if c["prefer"]])
def test_output_is_preferred(c):
    assert run(c)[0] == c["prefer"]


REQ = [(c, r) for c in NOW for r in c["requires"]]


@pytest.mark.parametrize("c,req", [pytest.param(c, r, id=f"{c['id']}-{r}",
                                                marks=mark(REQUIRES_MISSING_NOW, (c["id"], r))) for c, r in REQ])
def test_card_shows(c, req):
    _, cmp, page = run(c)
    assert shows(req, cmp, page, c)


@pytest.mark.parametrize("c", NOW, ids=[c["id"] for c in NOW])
def test_card_language(c):
    _, _, page = run(c)
    if page is None:
        return
    text = re.sub(r"<[^>]+>", " ", page).lower()
    bad = [w for w in SET["language"]["forbidden"] if re.search(r"\b" + re.escape(w), text)]
    assert not bad, bad


COMPARABLE = [c for c in NOW if (r := run(c))[1] is not None and r[1].status is ComparisonStatus.COMPARABLE]


@pytest.mark.parametrize("c", [pytest.param(c, id=c["id"], marks=mark(SCOPE_CAVEAT_MISSING_NOW, c["id"]))
                               for c in COMPARABLE])
def test_comparable_card_says_what_the_number_is_of(c):
    assert SET["language"]["scope_caveat"] in run(c)[2].lower()


@pytest.mark.parametrize("c", [c for c in CASES if c["component"] not in BUILT
                               and not c.get("superseded_by")], ids=lambda c: c["id"])
def test_component_to_come(c):
    pytest.skip(f"{c['component']} is not built yet; this case's expectations are fixed in hostile_set.json")


def test_superseded_cases_point_at_their_replacement():
    """The abstract evidence cases (B1-B5) are retired, not edited: each names
    the set that replaced it, and that set exists and covers the evidence card."""
    sup = [c for c in CASES if c.get("superseded_by")]
    assert len(sup) == 9 and all(c["component"] == "evidence_card" for c in sup)
    for c in sup:
        target = HERE / c["superseded_by"]
        assert target.is_file()
        comps = {x["component"] for x in json.loads(target.read_text(encoding="utf-8"))["cases"]}
        assert "evidence_card" in comps
