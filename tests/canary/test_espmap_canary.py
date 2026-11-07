"""The non-blocking canary: does Arrowfield read records from an espmap newer
than the pinned v0.2.1, and do its comparisons keep the two apart?

Runs only when CANARY_RECORDS names a directory written by
tests/canary/compute_canary_records.py (the CI job "canary"); otherwise every
test skips, so the required check never depends on it. A failure here is
drift to look at, not a broken build: v0.2.1 stays the pinned instrument.

What it expects of v0.3.0, measured 2026-10-25 on ethanol and acetaldehyde:
the same espmap.record/1 schema and top-level fields; an analysis key that
differs from the pinned records' only in margin_bohr (review 7 R7-6 put the
margin in the key); site means and scales that reproduce the pinned ones.
v0.3.1 (the canary since 2026-11-07) meets the same expectations unchanged:
it writes the same analysis keys as v0.3.0 on all three records (both are
espmap.numeric/1).

Release numbers do not decide comparability; the protocol does. v0.3.0 and
v0.3.1 are two releases with one protocol, and they write equal keys. So the
canary tells two cases apart by the key itself:
- the one decided transition: a pinned record from before v0.3, whose
  analysis key has no margin_bohr, against a canary record that has it. Only
  margin_bohr may move, the keys differ, and the two never compare;
- everything else: a matched canary and pinned record under one protocol.
  No key field may move, the keys are equal, and the records compare. After
  the pin moves to v0.3.1 (tests/espmap_pin.py) and records/ are recomputed,
  this is the case, and it stays so for any later release that keeps the
  protocol. A release that moves a key field fails here until that change
  is decided and named.
"""
import json
import math
import os
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
import model  # noqa: E402
import record  # noqa: E402

CANARY = os.environ.get("CANARY_RECORDS")
STEMS = ("ethanol__s12648430", "ethanol__s1", "acetaldehyde__s12648430")
# key fields allowed to differ from the pinned records' in the one decided transition, legacy
# (before v0.3) to v0.3: review 7 R7-6 put margin_bohr in the analysis key
EXPECTED_KEY_CHANGES = {"margin_bohr"}
REPRODUCE = 0.01                                     # kcal/mol: site means and scale limits

pytestmark = pytest.mark.skipif(not CANARY, reason="the canary runs only in its CI job (CANARY_RECORDS)")


def _new(stem):
    return record.load(Path(CANARY) / stem)


def _pinned(stem):
    return record.load(HERE / "records" / stem)


def _analysis_fields(directory, stem):
    return json.loads((directory / f"{stem}.json").read_text(encoding="utf-8"))["key_fields"]["analysis"]


def _legacy_transition(stem):
    """True only for the decided transition: the pinned record's analysis key
    predates margin_bohr and the canary record's carries it. Read from the
    keys, never from release numbers."""
    new, old = _analysis_fields(Path(CANARY), stem), _analysis_fields(HERE / "records", stem)
    return "margin_bohr" in new and "margin_bohr" not in old


@pytest.mark.parametrize("stem", STEMS)
def test_arrowfield_reads_the_new_record(stem):
    v = _new(stem)
    assert v.provenance.record_schema == "espmap.record/1"
    assert v.provenance.espmap_version == os.environ.get("CANARY_ESPMAP_VERSION", v.provenance.espmap_version)
    assert v.sites and v.surface_loader is not None
    assert len(v.surface("nearest_nucleus").faces)


@pytest.mark.parametrize("stem", STEMS)
def test_the_record_shape_is_unchanged(stem):
    new = json.loads((Path(CANARY) / f"{stem}.json").read_text(encoding="utf-8"))
    old = json.loads((HERE / "records" / f"{stem}.json").read_text(encoding="utf-8"))
    assert set(new) == set(old), {"new only": sorted(set(new) - set(old)), "pinned only": sorted(set(old) - set(new))}


@pytest.mark.parametrize("stem", STEMS)
def test_the_key_moved_only_where_decided(stem):
    new, old = _analysis_fields(Path(CANARY), stem), _analysis_fields(HERE / "records", stem)
    moved = {k for k in set(new) | set(old) if new.get(k) != old.get(k)}
    if _legacy_transition(stem):
        assert moved <= EXPECTED_KEY_CHANGES, f"key fields moved beyond the decided ones: {sorted(moved)}"
        assert _new(stem).analysis_key != _pinned(stem).analysis_key
    else:
        assert not moved, f"one protocol, yet key fields moved: {sorted(moved)}"
        assert _new(stem).analysis_key == _pinned(stem).analysis_key


def test_the_pinned_and_new_records_are_never_compared():
    """Across the legacy-to-v0.3 transition a difference card says NOT
    COMPARABLE: the keys keep a v0.3 number from being read against a legacy
    one. Under one protocol the same molecule's records compare."""
    new = [_new("ethanol__s12648430"), _new("ethanol__s1")]
    if not _legacy_transition("ethanol__s12648430"):
        assert model.comparable(new[0], _pinned("ethanol__s12648430"))
        return
    old = [_pinned("acetaldehyde__s12648430"), _pinned("acetaldehyde__s1")]
    label = new[0].site_labels()[0]
    c = model.compare(new, old, label)
    assert c.status is model.ComparisonStatus.NOT_COMPARABLE


def test_records_of_one_release_compare():
    new = [_new("ethanol__s12648430"), _new("ethanol__s1")]
    assert model.comparable(new[0], new[1])


@pytest.mark.parametrize("stem", STEMS)
def test_the_numbers_reproduce(stem):
    new, old = _new(stem), _pinned(stem)
    assert math.isclose(new.scale.vmin, old.scale.vmin, abs_tol=REPRODUCE)
    assert math.isclose(new.scale.vmax, old.scale.vmax, abs_tol=REPRODUCE)
    for label in old.site_labels():
        for rule in model.RULES:
            a, b = new.site(label, rule).mean, old.site(label, rule).mean
            assert (a is None) == (b is None) and (a is None or math.isclose(a, b, abs_tol=REPRODUCE)), (label, rule)
