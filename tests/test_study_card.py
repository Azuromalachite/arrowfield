"""The study-verdict card's hostile set (tests/pedagogy/study_card_set.json),
written before the card exists.

Now: the set is what its generator writes, every case's `allowed` is exactly
its safe outputs, all 12 states are covered with two or more variants, and
the set's own premises hold (each row meant to be valid follows from its own
conditions under its rule; each S10 row does not).
When views.study_card_html exists: every case's output is safe and
preferred, and a rendered card meets every requirement. Until then those
tests skip. The card PR may only add a strict-xfail list here, never edit
an expectation.
"""
import json
import re
import sys
from collections import Counter
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
PED = Path(__file__).with_name("pedagogy")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PED))

SET = json.loads((PED / "study_card_set.json").read_text(encoding="utf-8"))
CASES = SET["cases"]
RULE = "espmap-sesiidae-001/range_partition_sign_v1"

# --- what the card gets wrong today (measured, not intended) -------------------
UNSAFE_NOW = set()      # empty since the card met the probe rounds of #11; a new failing case goes here
NOT_PREFERRED_NOW = set(UNSAFE_NOW)
REQUIRES_MISSING_NOW = set()


def _mark(ids, cid):
    return [pytest.mark.xfail(strict=True, reason="listed as failing now")] if cid in ids else []


def test_the_committed_set_is_what_the_generator_writes():
    import make_study_card_set
    assert make_study_card_set.document() == (PED / "study_card_set.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_allowed_is_exactly_the_safe_outputs(c):
    assert c["allowed"] == sorted(k for k, g in c["outcomes"].items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert c["prefer"] is None or c["prefer"] in c["allowed"]
    assert c["why"] and c["outcomes"]["ERROR"] == "CRASH"


def test_all_twelve_states_have_two_or_more_variants():
    per = Counter(c["cls"].split()[0] for c in CASES)
    assert sorted(per) == sorted(f"S{i}" for i in range(1, 14))
    assert min(per.values()) >= 2, per


def _follows(row):
    c1, c2, c3 = (row["conditions"][k] for k in ("1_ranges_disjoint", "2_beyond_partition", "3_same_sign"))
    lo, hi = row["bootstrap_interval"]
    return ("resolved" if c1 and c2 and c3 else
            "different on average" if c2 and c3 and (lo > 0 or hi < 0) else "unresolved")


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_the_set_premises_hold(c):
    """A card case tests the card, not a corrupt input, unless it says so."""
    r = c["row"]
    if c["cls"].startswith("S13"):
        return                                         # malformed numbers: no premise to check
    if c["cls"].startswith(("S9", "S12")) and r["verdict_rule"] != RULE:
        return                                         # the case is about the rule itself
    if r["outcome"] == "not evaluable":
        return
    if c["cls"].startswith("S10"):
        assert _follows(r) != r["outcome"], "an S10 row must contradict itself"
    else:
        assert _follows(r) == r["outcome"], "a row meant to be valid must be self-consistent"


# ---------------------------------------------------- the card (to come)
def _card():
    import views
    if not hasattr(views, "study_card_html"):
        pytest.skip("views.study_card_html is not built yet; its expectations are fixed in study_card_set.json")
    return views.study_card_html


def _run(c):
    card = _card()                                     # skips until the card exists
    from model import Refused, StudySite               # StudySite arrives with study.py (#8)
    row = {k: (tuple(v) if isinstance(v, list) else v) for k, v in c["row"].items()}
    if c.get("row_n_numpy"):                           # JSON cannot carry numpy integers
        import numpy as np
        row["n"] = (np.int64(row["n"][0]), np.int32(row["n"][1]))
    try:
        page = card(StudySite(**row), c["context"])
    except Refused:
        return "REFUSED", None
    except Exception:  # noqa: BLE001 - any other error is a crash, graded as such
        return "ERROR", None
    return "RENDERED", page


@pytest.mark.parametrize("c", [pytest.param(c, id=c["id"], marks=_mark(UNSAFE_NOW, c["id"])) for c in CASES])
def test_card_output_is_safe(c):
    out, _ = _run(c)
    assert out in c["allowed"], f"{out}: {c['outcomes'][out]}"


@pytest.mark.parametrize("c", [pytest.param(c, id=c["id"], marks=_mark(NOT_PREFERRED_NOW, c["id"]))
                               for c in CASES if c["prefer"]])
def test_card_output_is_preferred(c):
    assert _run(c)[0] == c["prefer"]


REQ = [(c, i) for c in CASES for i in range(len(c["requires"]))]


@pytest.mark.parametrize("c,i", [pytest.param(c, i, id=f"{c['id']}-{i}",
                                              marks=_mark(REQUIRES_MISSING_NOW, (c["id"], i)))
                                 for c, i in REQ])
def test_a_rendered_card_meets_its_requirements(c, i):
    out, page = _run(c)
    if out != "RENDERED":
        return                                         # a refusal shows nothing to misread
    q = c["requires"][i]
    target = page if q["on"] == "html" else re.sub(r"<[^>]+>", " ", page)
    if "must" in q:
        assert re.search(q["must"], target, re.I), f"{q['why']} ({q['grade']} if missing)"
    else:
        assert not re.search(q["must_not"], target, re.I), f"{q['why']} ({q['grade']} if present)"
