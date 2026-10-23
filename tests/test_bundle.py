"""bundle.py against a real Atlas teaching bundle (tests/data, unedited).

The adapter keeps every field the Atlas exports, refuses one it does not
model, and derives a card's broad category from the raw fields without
replacing them: a gland negative, a lure with no effect and an effect whose
significance was never tested stay different statements.
"""
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from bundle import InvalidBundle, load  # noqa: E402
from model import EvidenceState  # noqa: E402

DATA = Path(__file__).with_name("data") / "sesiidae-homework-1.bundle.json"
RAW = json.loads(DATA.read_text(encoding="utf-8"))
OAC, OH, E3 = "spa.cmpd.z3_z13_18_oac", "spa.cmpd.z3_z13_18_oh", "spa.cmpd.e3_z13_18_oac"
NAKA_GLAND = "spa.obs.synanthedon_tenuis.naka_2013.gland"
YANG_GLAND = "spa.obs.synanthedon_tenuis.yang_2012.gland"
NAKA_TRAP = "spa.obs.synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1"
TAS_1987 = "spa.obs.synanthedon_tipuliformis.szocs_1990.trap_e3.tasmania_1987"
NZ = "spa.obs.synanthedon_tipuliformis.suckling_2005.trap_e3.new_zealand"


@pytest.fixture(scope="module")
def b():
    return load(DATA)


def test_the_bundle_says_where_it_came_from(b):
    assert b.case_id == "sesiidae-homework-1" and b.atlas_version == RAW["atlas_version"]
    assert len(b.sha256) == 64
    assert set(b.observations) == {o["observation_id"] for o in RAW["observations"]}
    assert len(b.contrasts) == len(RAW["comparisons"])


def test_every_raw_field_is_kept_as_recorded(b):
    o = b.observations[NAKA_TRAP]
    raw = next(x for x in RAW["observations"] if x["observation_id"] == NAKA_TRAP)
    assert o.statement == raw["statement"]                       # the Atlas's sentence, verbatim
    assert (o.effect, o.significance) == ("decreased", "significant")
    assert (o.response.treatment_value, o.response.comparator_value) == (4.6, 17.0)
    assert o.response.qualifier == "as_reported" and o.dose == raw["dose"]
    assert o.author_statement == raw["evidence"]["author_statement"]
    assert o.context == raw["context"] and o.reference["doi"] == raw["reference"]["doi"]
    g = b.observations[NAKA_GLAND]
    assert g.unresolved == ("detection limit for the negatives", "dispersion of the 2 ng titre")
    assert next(d for d in g.detections if d.compound == OAC).measurements[0]["value"] == 2


def test_a_gland_negative_says_what_was_absent_and_at_what_limit(b):
    s, label = b.observations[YANG_GLAND].display_state(OAC)
    assert s is EvidenceState.NEGATIVE and "100 picogram" in label
    s, label = b.observations[NAKA_GLAND].display_state(OH)
    assert s is EvidenceState.NEGATIVE and "limit not stated" in label
    assert b.observations[NAKA_GLAND].display_state(OAC)[0] is EvidenceState.OBSERVED


def test_no_effect_is_negative_only_when_a_test_says_so(b):
    s, label = b.observations[TAS_1987].display_state(E3)
    assert s is EvidenceState.NEGATIVE and "not significant" in label
    s, label = b.observations[NZ].display_state(E3)           # no_effect, not_reported
    assert s is EvidenceState.OBSERVED and "significance not reported" in label


def test_a_compound_outside_an_observation_is_not_tested_there(b):
    assert b.observations[TAS_1987].display_state(OH)[0] is EvidenceState.NOT_TESTED
    assert b.observations[NAKA_GLAND].display_state("spa.cmpd.e2_z13_18_oac_x")[0] is EvidenceState.NOT_TESTED


def _bad(tmp_path, mutate):
    j = json.loads(DATA.read_text(encoding="utf-8"))
    mutate(j)
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(j), encoding="utf-8")
    return p


@pytest.mark.parametrize("mutate, match", [
    (lambda j: j.update(bundle="spa.teaching/1"), "not a spa.teaching/0"),
    (lambda j: j.update(schema_version="0.2"), "schema_version"),
    (lambda j: j.update(atlas_version=""), "atlas_version"),
    (lambda j: j["observations"][0].update(new_field=1), "not modelled"),
    (lambda j: j["observations"][2]["outcome"].update(confidence="high"), "not modelled"),
    (lambda j: j["comparisons"][0]["observation_ids"].append("spa.obs.nowhere"), "does not carry"),
    (lambda j: j["compounds"][0]["occurrences"].append({"observation_id": "spa.obs.nowhere", "roles": ["x"]}),
     "does not carry"),
    (lambda j: j["observations"].append(dict(j["observations"][0])), "repeats"),
])
def test_a_bundle_arrowfield_cannot_show_faithfully_is_refused(tmp_path, mutate, match):
    with pytest.raises(InvalidBundle, match=match):
        load(_bad(tmp_path, mutate))


def test_nan_is_refused(tmp_path):
    p = tmp_path / "nan.json"
    p.write_text(DATA.read_text(encoding="utf-8").replace('"treatment_value": 4.6', '"treatment_value": NaN'),
                 encoding="utf-8")
    with pytest.raises(InvalidBundle, match="NaN"):
        load(p)
