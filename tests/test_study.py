"""study.py on the real espmap-sesiidae-001 results (tests/data, copied
unedited from the Atlas at 89f5e22) and the real aldehyde control records
(controls/). The adapter shows the study's own verdicts, never re-decides
them, and passes each member's folded flag through as the index states it.
"""
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from model import compare  # noqa: E402
from study import InvalidStudy, load_ensemble, load_results  # noqa: E402

RESULTS = Path(__file__).with_name("data") / "espmap-sesiidae-001.results.json"
AL = next((HERE / "controls" / "espmap" / "z3_z13_18_al").glob("*/ensemble.json"))
RULE = "espmap-sesiidae-001/range_partition_sign_v1"


@pytest.fixture(scope="module")
def r():
    return load_results(RESULTS)


@pytest.fixture(scope="module")
def al():
    return load_ensemble(AL)


def test_every_verdict_carries_the_rule_that_decided_it(r):
    assert r.verdict_rule == RULE and all(s.verdict_rule == RULE for s in r.sites)
    assert len(r.sha256) == 64


def test_the_study_verdicts_are_shown_as_the_study_decided_them(r):
    c4 = r.site("C_position", "C4-H")
    assert c4.outcome == "resolved" and c4.delta == pytest.approx(5.22, abs=0.005)
    assert c4.range_x[1] < c4.range_y[0]                      # the ranges do not overlap
    assert all(c4.conditions.values())
    b3 = r.site("B_stereo", "C3-H")
    assert b3.outcome == "different on average" and not b3.conditions["1_ranges_disjoint"]
    assert r.site("A_headgroup", "global V_max").outcome == "resolved"
    ext = r.site("C_position", "C4-H", population="extended")
    assert ext.outcome == "resolved" and ext.n != c4.n


def test_the_bootstrap_interval_is_kept_under_its_honest_name(r):
    b3 = r.site("B_stereo", "C3-H")
    assert b3.bootstrap_interval[0] < b3.delta < b3.bootstrap_interval[1]
    assert not hasattr(b3, "ci95")


def test_members_come_from_the_index_with_its_folded_flags(r, al):
    source, views = al
    idx = json.loads(AL.read_text(encoding="utf-8"))
    assert source.atlas_id == "arrowfield.control.z3_z13_18_al" and len(views) == 20
    assert [v.folded for v in views] == [m["folded"] for m in idx["members"]]
    assert sum(v.folded for v in views) == r.molecules["z3_z13_18_al"]["n_folded"]
    assert {v.analysis_key for v in views} == {source.analysis_key}


def test_the_members_reproduce_the_studys_numbers(r, al):
    """The control's 20 head-O means, read from its records, span exactly the
    range the study reports for it: the adapter and the study read one data set."""
    _, views = al
    means = [v.site("head O").mean for v in views]
    s = r.site("A_headgroup (control: Al vs OAc)", "head O")
    assert (min(means), max(means)) == pytest.approx(s.range_y, abs=1e-9)


def test_arrowfields_teaching_verdict_is_not_the_studys(r, al):
    """compare() is Arrowfield's own teaching verdict; the study's verdict
    comes from the study. Comparing a molecule with itself shows the two are
    separate things: compare() overlaps, and no study row says anything."""
    _, views = al
    card = compare(views, views, "head O")
    assert card.verdict.name == "OVERLAPPING"
    with pytest.raises(KeyError):
        r.site("Al vs Al", "head O")


def _write(tmp_path, mutate):
    j = json.loads(RESULTS.read_text(encoding="utf-8"))
    mutate(j)
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(j), encoding="utf-8")
    return p


@pytest.mark.parametrize("mutate, match", [
    (lambda j: j.update(verdict_rule="espmap.compare/1"), "never show"),
    (lambda j: j.update(partial=True), "partial"),
    (lambda j: j["contrasts"][0]["sites"]["primary:C3-H"].update(outcome="significant"), "not one of"),
    (lambda j: j["contrasts"][0]["sites"]["primary:C3-H"].update(p_value=1e-11), "not modelled"),
])
def test_a_study_file_arrowfield_cannot_show_faithfully_is_refused(tmp_path, mutate, match):
    with pytest.raises(InvalidStudy, match=match):
        load_results(_write(tmp_path, mutate))


def test_a_member_changed_after_indexing_is_refused(tmp_path):
    import shutil
    d = tmp_path / "z3_z13_18_al" / AL.parent.name
    shutil.copytree(AL.parent, d)
    rec = d / "conformers" / "conf_00.json"
    rec.write_text(rec.read_text(encoding="utf-8").replace('"genus": 0', '"genus": 1', 1), encoding="utf-8")
    with pytest.raises(InvalidStudy, match="sha256"):
        load_ensemble(d / "ensemble.json")
