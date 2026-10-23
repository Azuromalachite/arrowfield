"""Sesiidae pilot: do stereo, positional and head-group changes move per-site
ESP by more than the conformers of one molecule do?

Molecules (C1 = the carbon bearing O; maps 1-18 = chain carbons, 31 = the O on
C1, 32 = carbonyl O, 21/22 = acetate methyl/carbonyl C):
  Z3Z13-OAc, E3Z13-OAc   stereo at C3=C4
  E2Z13-OAc, E3Z13-OAc   double bond moved C3=C4 -> C2=C3
  Z3Z13-OH,  Z3Z13-OAc   head group
Per molecule: 100 ETKDG conformers (esp_map's embedder, seed 0xC0FFEE), then the
lowest-MMFF one plus 4 spread over the radius of gyration (bench/conformers.py
select), so folded shapes are included. B3LYP/def2-SVP + density fitting, gas
phase (the molecules are neutral): a SPREAD measurement, like bench/conformers.
Sites: C<n>-H (carbon n with its hydrogens) and the O sites, via caller groups,
so the same label means the same chain position in every molecule.

Scratch only; uses esp_map from the caller-sites worktree (P2 + P1 + --sites).
Appends to sesiidae_pilot.jsonl; resumable.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

SP = Path("/mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/"
          "669a3716-c2cd-43c5-9905-c4ad082b6f87/scratchpad")
sys.path.insert(0, str(SP / "espmap-p1"))
sys.path.insert(0, str(SP / "espmap-p1" / "bench"))
import esp_map as em  # noqa: E402
from conformers import select  # noqa: E402

OUT = SP / "sesiidae_pilot.jsonl"
CH2 = lambda a, b: "".join(f"[CH2:{n}]" for n in range(a, b + 1))  # noqa: E731
TAIL = f"/[CH:13]=[CH:14]\\{CH2(15, 17)}[CH3:18]"                 # Z13 everywhere
ACE = "[CH3:21][C:22](=[O:32])[O:31]"
MOLS = {
    "Z3Z13-OAc": f"{ACE}[CH2:1][CH2:2]/[CH:3]=[CH:4]\\{CH2(5, 12)}{TAIL}",
    "E3Z13-OAc": f"{ACE}[CH2:1][CH2:2]/[CH:3]=[CH:4]/{CH2(5, 12)}{TAIL}",
    "E2Z13-OAc": f"{ACE}[CH2:1]/[CH:2]=[CH:3]/{CH2(4, 12)}{TAIL}",
    "Z3Z13-OH": f"[OH:31][CH2:1][CH2:2]/[CH:3]=[CH:4]\\{CH2(5, 12)}{TAIL}",
}
EXPECT = {  # unmapped canonical SMILES the conformers must reproduce
    "Z3Z13-OAc": "CCCC/C=C\\CCCCCCCC/C=C\\CCOC(C)=O",
    "E3Z13-OAc": "CCCC/C=C\\CCCCCCCC/C=C/CCOC(C)=O",
    "E2Z13-OAc": "CCCC/C=C\\CCCCCCCCC/C=C/COC(C)=O",
    "Z3Z13-OH": "CCCC/C=C\\CCCCCCCC/C=C\\CCO",
}
K = 5


def groups(smiles):
    maps = [int(t.split(":")[1].rstrip("]")) for t in smiles.replace("]", "]|").split("|") if ":" in t]
    g = [{"label": f"C{n}-H", "atoms": [f"{n}+H"]} for n in range(1, 19)]
    g.append({"label": "O31", "atoms": ["31+H"]})                  # ester O / hydroxyl
    if 32 in maps:
        g += [{"label": "O32", "atoms": [32]}, {"label": "acetyl", "atoms": [21, "21+H", 22, 32]}]
    return g


def done_keys():
    return {json.loads(l)["key"] for l in OUT.read_text().splitlines()} if OUT.exists() else set()


def run(name, smiles, check_only=False):
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors3D
    _, info = em.from_smiles(smiles, n_confs=100)
    mol = info["mol"]
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    cids = [c.GetId() for c in mol.GetConformers()]
    energy = np.array([AllChem.MMFFGetMoleculeForceField(mol, props, confId=c).CalcEnergy()
                       for c in cids])
    rg = np.array([Descriptors3D.RadiusOfGyration(mol, confId=c) for c in cids])
    picks = select(rg, energy, K)
    approved = Chem.MolFromSmiles(smiles)
    maps = {a.GetIdx(): a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    done = done_keys()
    for i in picks:
        cid = cids[i]
        key = f"{name}|{cid}"
        spec = em._mol_to_spec(mol, cid)
        ok, computed, why = em.verify_computed(spec, approved)
        assert ok and computed == EXPECT[name], (name, cid, ok, computed, why)
        assert em.verify_mapping(spec, approved)[0], (name, cid)
        rows = [t.split() for t in spec.split(";")]
        sym = [r[0] for r in rows]
        X = np.array([[float(v) for v in r[1:4]] for r in rows])
        grp = em.resolve_site_groups(groups(smiles), sym, X, maps)
        if check_only:
            print(name, "conf", cid, "stereo ok:", computed, "| groups", len(grp))
            return
        if key in done:
            continue
        t = time.perf_counter()
        m, mf = em.run_scf(spec, basis="def2-svp", density_fit=True)
        dm = em.total_dm(mf)
        verts, faces, diag = em.extract_surface(m, dm, 0.26)
        V = em.esp_at_points(m, dm, verts) * em.HA2KCAL
        rep = em.site_report(verts * em.BOHR, faces, V, m.atom_coords() * em.BOHR, sym,
                             atom_maps=maps, groups=grp)
        sites = {s["label"]: {r: {k: s[r][k] for k in ("area_A2", "V_area_mean_kcal_mol",
                                                        "V_p05_kcal_mol", "V_min_kcal_mol")}
                              for r in em.PARTITION_RULES}
                 for s in rep["reaction_sites"]}
        row = dict(key=key, molecule=name, smiles=smiles, conf_id=cid,
                   is_lowest_E=bool(i == int(np.argmin(energy))),
                   mmff_dE=float(energy[i] - energy.min()), rg=float(rg[i]),
                   genus=diag["genus"], V_min=float(V.min()), V_max=float(V.max()),
                   area_A2=em.area_stats(verts * em.BOHR, faces, V)["area_A2"],
                   sites=sites, t_wall=time.perf_counter() - t,
                   method="B3LYP/def2-SVP RI, gas, 0.26 bohr")
        with OUT.open("a") as f:
            f.write(json.dumps(row) + "\n")
        print(f"{name:10s} conf {cid:3d} rg {rg[i]:.2f} genus {diag['genus']} "
              f"V {V.min():+.1f}..{V.max():+.1f}  {row['t_wall']:.0f}s", flush=True)


if __name__ == "__main__":
    check = "--check" in sys.argv
    for name, smi in MOLS.items():
        run(name, smi, check_only=check)
