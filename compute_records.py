"""Compute espmap records (JSON + NPZ) for manifest molecules, for the view
prototypes. Small molecules only: seconds each at the preview basis.

Runs espmap's own CLI, so the records are exactly what a notebook would get:
mapped SMILES (or the manifest's DFT geometry for anions), the manifest's
sites as --sites groups, --json and --npz. Several seeds per molecule give a
small conformer set, which is what the difference card compares.

    python compute_records.py /path/to/espmap [id ...]     # writes records/<id>__s<seed>.*
"""
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ESPMAP = Path(sys.argv[1]).resolve()
OUT = HERE / "records"
SEEDS = (0xC0FFEE, 1, 2)
DEFAULT = ["ethanol", "acetaldehyde", "e2_butene", "z2_butene", "acetate"]


def manifest():
    import yaml
    return {m["id"]: m for m in yaml.safe_load((HERE / "molecules.yaml").read_text(encoding="utf-8"))["molecules"]}


def sites_file(groups):
    """espmap's --sites file (#36): it must say how atoms are named. The
    manifest names every site atom by its atom map, never by index, so map 3
    is C3 in every molecule."""
    return {"by": "map", "groups": groups}


def run(m, seed):
    stem = OUT / f"{m['id']}__s{seed}"
    if stem.with_suffix(".json").exists():
        return
    cmd = [sys.executable, str(ESPMAP / "esp_map.py"), "--basis", "def2-svp", "--yes",
           "--seed", str(seed), "--out", str(stem.with_suffix(".png")),
           "--json", str(stem.with_suffix(".json")), "--npz", str(stem.with_suffix(".npz"))]
    geom = m.get("geometry", "")
    if geom.startswith("geometries/"):          # DFT geometry for an anion
        cmd += ["--file", str(HERE / geom.split()[0])]
    else:
        cmd += ["--smiles", m["mapped_smiles"]]
    if m.get("solvent", "vacuum") != "vacuum":
        cmd += ["--solvent", m["solvent"]]
    if m.get("sites"):
        g = stem.with_suffix(".sites.json")
        g.write_text(json.dumps(sites_file(m["sites"])))
        cmd += ["--sites", str(g)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f"{m['id']} seed {seed} failed:\n{p.stdout[-1500:]}\n{p.stderr[-1500:]}")
    print(f"{m['id']:14s} seed {seed}: ok", flush=True)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    mans = manifest()
    for mid in (sys.argv[2:] or DEFAULT):
        m = mans[mid]
        # a file geometry is one geometry: no seeds to vary
        for seed in (SEEDS[:1] if m.get("geometry", "").startswith("geometries/") else SEEDS):
            run(m, seed)
