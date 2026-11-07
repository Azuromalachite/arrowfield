"""What moving the espmap pin from v0.2.1 to v0.3.1 must deliver (expectations
first: the move itself is a later change).

Strict xfails mark what the repository does not do yet; each comes off in the
change that makes it pass, never by editing the expectation. The guards
without a mark hold today and must keep holding through the move.
"""
import json
import re
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from espmap_pin import PIN, STUDY_PRODUCER  # noqa: E402

PIN_MOVE = "the pin move has not happened yet (espmap v0.2.1 -> v0.3.1)"
RECORDS = sorted(p.name for p in (HERE / "records").glob("*.json") if not p.name.endswith(".sites.json"))


def test_the_demo_records_are_listed():
    """The parametrised check below must have something to check."""
    assert len(RECORDS) >= 10, RECORDS


@pytest.mark.xfail(strict=True, reason="the demo records say espmap 0.1.0, code 'unknown': " + PIN_MOVE)
@pytest.mark.parametrize("name", RECORDS)
def test_a_demo_record_names_the_pinned_build(name):
    """records/ are Arrowfield's own demo records: computed with the pinned
    release, from an identifiable build. A record whose build cannot be named
    is one espmap v0.3.1's own compare refuses (COMPARE_CONTRACT_UNKNOWN)."""
    prov = json.loads((HERE / "records" / name).read_text(encoding="utf-8"))["provenance"]
    assert prov.get("espmap_version") == PIN["version"], prov.get("espmap_version")
    assert prov.get("code") in (PIN["tag"], PIN["commit"]), prov.get("code")


@pytest.mark.xfail(strict=True, reason="README still names v0.2.1 as the pin: " + PIN_MOVE)
def test_the_readme_names_the_pin():
    text = (HERE / "README.md").read_text(encoding="utf-8")
    assert f"pinned to a release (`{PIN['tag']}`" in text


@pytest.mark.xfail(strict=True, reason="CONTRIBUTING installs v0.2.1: " + PIN_MOVE)
def test_contributing_installs_the_pin():
    text = (HERE / "CONTRIBUTING.md").read_text(encoding="utf-8")
    assert f"espmap@{PIN['tag']}\"" in text


@pytest.mark.xfail(strict=True, reason="notebooks.md pins v0.2.1: " + PIN_MOVE)
def test_the_notebooks_pin_the_release():
    text = (HERE / "notebooks.md").read_text(encoding="utf-8")
    assert f'ESPMAP_REF = "{PIN["tag"]}"' in text


def _tag_tuple(tag):
    m = re.fullmatch(r"v(\d+)\.(\d+)\.(\d+)", tag)
    assert m, f"not a release tag: {tag!r}"
    return tuple(int(x) for x in m.groups())


def test_the_canary_is_not_behind_the_pin():
    """The canary tracks a release at or after the pin, never before it and
    never a branch."""
    ci = (HERE / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    tags = re.findall(r"git clone --depth 1 --branch (\S+) https://github.com/Azuromalachite/espmap", ci)
    assert len(tags) == 1, tags
    assert _tag_tuple(tags[0]) >= _tag_tuple(PIN["tag"]), tags[0]


def test_the_study_control_keeps_its_producer():
    """The aldehyde control is the study's data (key d0f5ed13, computed with
    v0.2.1): the pin move never recomputes it. Recomputed under v0.3.1 its key
    would change and the homework page would refuse it."""
    ensembles = sorted((HERE / "controls" / "espmap").glob("*/*/ensemble.json"))
    assert ensembles, "no control ensemble found"
    for path in ensembles:
        ens = json.loads(path.read_text(encoding="utf-8"))
        assert ens["espmap"]["version"] == STUDY_PRODUCER["version"], path
        assert ens["espmap"]["commit"] == STUDY_PRODUCER["commit"], path
        assert path.parent.name.endswith(f"espmap-v{STUDY_PRODUCER['version']}"), path


# espmap-sesiidae-001's frozen protocol (B3LYP/def2-SVP + RI, gas phase), as its key: written out
# in full here, not read from model.STUDY_KEYS, so a change to either is caught by the other
STUDY_KEY = "d0f5ed136cda12ada8c6324a79ec9f7e79134d2f"
ALDEHYDE_CONTROL = HERE / "controls" / "espmap" / "z3_z13_18_al" / "b3lyp_def2-svp_gas_espmap-v0.2.1"


def test_the_aldehyde_control_keeps_the_study_protocol():
    """The named control stays under the study's frozen key, in its index and
    in every member: the pin protects the study's protocol, not only its
    producer. A v0.3.1 rerun of the study is a decision of its own."""
    import model
    assert model.STUDY_KEYS["espmap-sesiidae-001/range_partition_sign_v1"] == STUDY_KEY
    ens = json.loads((ALDEHYDE_CONTROL / "ensemble.json").read_text(encoding="utf-8"))
    assert ens["analysis_key"] == STUDY_KEY
    members = sorted((ALDEHYDE_CONTROL / "conformers").glob("conf_*.json"))
    assert len(members) == 20, len(members)
    for m in members:
        rec = json.loads(m.read_text(encoding="utf-8"))
        assert rec["analysis_key"] == STUDY_KEY, m.name
