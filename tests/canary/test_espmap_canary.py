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

When the pin moves to v0.3.1 (tests/espmap_pin.py), records/ are recomputed
with it and the canary and the pin are one release until a newer one exists.
Then nothing crosses releases: the keys must be equal, no key field may move,
and the canary's records compare with the pinned ones. Which case holds is
read from the two records' own espmap_version, never assumed.
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
# key fields allowed to differ from the pinned records' (each one a decided change)
EXPECTED_KEY_CHANGES = {"margin_bohr"}
REPRODUCE = 0.01                                     # kcal/mol: site means and scale limits

pytestmark = pytest.mark.skipif(not CANARY, reason="the canary runs only in its CI job (CANARY_RECORDS)")


def _new(stem):
    return record.load(Path(CANARY) / stem)


def _pinned(stem):
    return record.load(HERE / "records" / stem)


def _one_release(stem):
    """True when the canary's record and the pinned one come from the same
    espmap release (after the pin move, until a newer release exists)."""
    a, b = _new(stem).provenance.espmap_version, _pinned(stem).provenance.espmap_version
    assert a and b, (a, b)
    return a == b


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
    new = json.loads((Path(CANARY) / f"{stem}.json").read_text(encoding="utf-8"))["key_fields"]["analysis"]
    old = json.loads((HERE / "records" / f"{stem}.json").read_text(encoding="utf-8"))["key_fields"]["analysis"]
    moved = {k for k in set(new) | set(old) if new.get(k) != old.get(k)}
    if _one_release(stem):
        assert not moved, f"one release, yet key fields moved: {sorted(moved)}"
        assert _new(stem).analysis_key == _pinned(stem).analysis_key
    else:
        assert moved <= EXPECTED_KEY_CHANGES, f"key fields moved beyond the decided ones: {sorted(moved)}"
        assert _new(stem).analysis_key != _pinned(stem).analysis_key


def test_the_pinned_and_new_records_are_never_compared():
    """A difference card across two releases says NOT COMPARABLE: the keys
    keep a newer release's number from being read against the pinned one.
    Within one release (canary = pin) the same molecule's records compare."""
    new = [_new("ethanol__s12648430"), _new("ethanol__s1")]
    if _one_release("ethanol__s12648430"):
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
