"""Rebuild tests/data/card_label_sites.json from the records it is copied from.

    python tests/pedagogy/extract_card_label_sites.py <sesiidae-pheromone-atlas checkout>

For each conformer of the four Atlas ensembles (computed/espmap/<slug>/
b3lyp_def2-svp_gas_espmap-v0.2.1) and of the aldehyde control (Arrowfield's
controls/espmap/z3_z13_18_al), it copies the named sites' means under both
partition rules from the record's reaction_sites block (nn, rs), the
partition sensitivity |rs - nn| (partition), the index's folded flag, the
record's analysis_key and its sha256. Nothing is recomputed. SITES names
which sites each molecule carries: the ones some hostile set or the homework
plan shows.
"""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
OUT = HERE / "tests" / "data" / "card_label_sites.json"
RUN = "b3lyp_def2-svp_gas_espmap-v0.2.1"
SOURCE = ("Atlas 35112e1 computed/espmap/<slug>/b3lyp_def2-svp_gas_espmap-v0.2.1 (Atlas compounds); "
          "Arrowfield controls/espmap/z3_z13_18_al (the teaching control). Values copied from each record's "
          "sites block: nn and rs are the site means under the nearest-nucleus and radius-scaled rules, "
          "partition the record's partition sensitivity, folded the index's flag.")
# the sites each molecule carries; the homework plan's near-to-far series (act A) and every site
# acts B and C ask about need the ones after the first group
SITES = {
    "z3_z13_18_oac": ["C1-H", "C2-H", "C3-H", "C4-H", "C5-H", "head O", "C6-H", "C10-H", "C16-H"],
    "e3_z13_18_oac": ["C3-H", "C4-H", "C5-H", "C2-H", "head O"],
    "e2_z13_18_oac": ["C4-H", "C2-H", "C3-H"],
    "z3_z13_18_oh": ["C2-H", "C3-H", "head O", "C1-H", "C4-H", "C6-H", "C10-H", "C16-H"],
    "z3_z13_18_al": ["C1-H", "C2-H", "head O"],
}


def _sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def molecule(run_dir, slug, sites):
    index = json.loads((run_dir / "ensemble.json").read_text(encoding="utf-8"))
    members = []
    for m in index["members"]:
        conf = Path(m["record"]).stem
        path = run_dir / "conformers" / f"{conf}.json"
        rec = json.loads(path.read_text(encoding="utf-8"))
        by_label = {s["label"]: s for s in rec["sites"]["reaction_sites"]}
        out = {}
        for site in sites:
            s = by_label[site]
            nn = s["nearest_nucleus"]["V_area_mean_kcal_mol"]
            rs = s["radius_scaled"]["V_area_mean_kcal_mol"]
            out[site] = {"nn": nn, "rs": rs, "partition": abs(rs - nn)}
        members.append({"conformer": conf, "record_sha256": _sha256(path), "folded": m["folded"],
                        "analysis_key": rec["analysis_key"], "sites": out})
    return {"atlas_id": index["atlas_id"], "index_sha256": _sha256(run_dir / "ensemble.json"), "members": members}


def build(atlas):
    mols = {}
    for slug, sites in SITES.items():
        root = (HERE / "controls" / "espmap" if slug == "z3_z13_18_al" else Path(atlas) / "computed" / "espmap")
        mols[slug] = molecule(root / slug / RUN, slug, sites)
    return {"source": SOURCE, "molecules": mols}


if __name__ == "__main__":
    OUT.write_text(json.dumps(build(sys.argv[1]), indent=1) + "\n", encoding="utf-8")
    print(f"-> {OUT}")
