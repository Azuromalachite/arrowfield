"""Compute a few small records with an espmap other than the pinned one, for
the non-blocking canary (tests/canary/test_espmap_canary.py).

    python tests/canary/compute_canary_records.py /path/to/espmap /out/dir

Runs that espmap's own CLI, as compute_records.py does for the pinned one
(the manifest's mapped SMILES and sites, def2-SVP, gas): ethanol at two
seeds, and acetaldehyde. Never writes to records/: the committed records stay
the pinned release's.
"""
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE))
import compute_records as cr  # noqa: E402

RUNS = (("ethanol", 0xC0FFEE), ("ethanol", 1), ("acetaldehyde", 0xC0FFEE))


def main(espmap, out):
    out.mkdir(parents=True, exist_ok=True)
    mans = cr.manifest()
    for mid, seed in RUNS:
        m = mans[mid]
        stem = out / f"{mid}__s{seed}"
        g = stem.with_suffix(".sites.json")
        g.write_text(json.dumps(cr.sites_file(m["sites"])))
        cmd = [sys.executable, str(espmap / "esp_map.py"), "--basis", "def2-svp", "--yes", "--seed", str(seed),
               "--smiles", m["mapped_smiles"], "--sites", str(g), "--out", str(stem.with_suffix(".png")),
               "--json", str(stem.with_suffix(".json")), "--npz", str(stem.with_suffix(".npz"))]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode:
            raise SystemExit(f"{mid} seed {seed} failed:\n{p.stdout[-1500:]}\n{p.stderr[-1500:]}")
        print(f"{mid} seed {seed}: ok", flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
