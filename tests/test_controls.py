"""The teaching controls' espmap records (controls/), checked as the Atlas
gate checks its own computed/: the index's checksums, one molecule, one
protocol, one espmap version, every record valid and verified by espmap,
and loadable by record.py. A control is never an Atlas compound."""
import hashlib
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from record import Record  # noqa: E402

ROOT = HERE / "controls" / "espmap"
ENSEMBLES = sorted(ROOT.glob("*/*/ensemble.json"))


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def canonical(smiles):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_there_are_controls():
    assert ENSEMBLES, "no controls/espmap/<slug>/<protocol>/ensemble.json"


@pytest.mark.parametrize("path", ENSEMBLES, ids=lambda p: p.parent.parent.name)
def test_a_control_index_holds(path):
    d = path.parent
    ens = json.loads(path.read_text(encoding="utf-8"))
    assert ens["schema"] == "spa.espmap_ensemble/0"
    assert ens["atlas_id"].startswith("arrowfield.control.")       # never an Atlas id
    assert ens["atlas_id"].split(".")[-1] == d.parent.name and ens["protocol"] == d.name
    sums = dict(reversed(line.split("  ", 1)) for line in (d / "SHA256SUMS").read_text().splitlines())
    on_disk = {"ensemble.json"} | {f"conformers/{p.name}" for p in (d / "conformers").glob("*.json")}
    assert set(sums) == on_disk
    assert all(sha256(d / n) == h for n, h in sums.items())
    assert len({m["geometry_sha1"] for m in ens["members"]}) == len(ens["members"])
    for m in ens["members"]:
        rec_path = d / m["record"]
        assert sha256(rec_path) == m["record_sha256"] and m["npz"] == f"{m['conformer']}.npz"
        rec = json.loads(rec_path.read_text(encoding="utf-8"))
        assert rec["analysis_key"] == ens["analysis_key"]
        assert rec["key_fields"]["analysis"] == ens["analysis"]
        assert rec["provenance"]["espmap_version"] == ens["espmap"]["version"]
        assert rec["input"]["structure_verified"] is True and rec["validity"]["input"] is True
        assert rec["geometry"]["sha1"] == m["geometry_sha1"]
        assert canonical(rec["input"]["computed_smiles"]) == canonical(ens["smiles"])
        Record(rec_path.with_suffix(""))                           # Arrowfield's own checks
