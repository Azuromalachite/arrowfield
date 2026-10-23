"""DFT geometries for manifest anions that espmap has no fixture for.

Same method as espmap's bench/fixture_geometries.py (whose functions this
reuses): B3LYP/def2-TZVP optimised in PCM water, Cartesian BFGS on analytic
gradients, gtol 3e-4 Ha/A, then checked with espmap's verify_computed so the
optimised geometry is still the molecule asked for. Small anions must not use
SMILES + MMFF: H8 measured up to 6.5 kcal/mol from the geometry alone.

    python make_geometries.py /path/to/espmap        # writes geometries/<id>.sdf
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ESPMAP = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ESPMAP))
sys.path.insert(0, str(ESPMAP / "bench"))
import esp_map as em  # noqa: E402
from fixture_geometries import bonds, build, mapped_mol, optimise, write_sdf  # noqa: E402

WANTED = {"acetate": "[CH3:1][C:2](=[O:3])[O-:4]"}   # mapped SMILES, as in molecules.yaml


def main():
    import numpy as np
    from rdkit import Chem
    out = HERE / "geometries"
    out.mkdir(exist_ok=True)
    record = {}
    for name, smi in WANTED.items():
        if (out / f"{name}.sdf").exists():
            continue
        spec, info = em.from_smiles(smi)
        sym = [t.split()[0] for t in spec.split(";")]
        X0 = np.array([[float(v) for v in t.split()[1:4]] for t in spec.split(";")])
        X, r = optimise(sym, X0, info["charge"])
        ok, computed, why = em.verify_computed(build(sym, X), Chem.MolFromSmiles(smi))
        assert ok, (name, why)
        props = dict(arrowfield_geometry=name, smiles=smi,
                     geometry="B3LYP/def2-TZVP, PCM water, Cartesian BFGS gtol 3e-4 Ha/A",
                     converged=r.success, iterations=r.nit,
                     espmap=em.code_identity()["code"])
        write_sdf(out / f"{name}.sdf", mapped_mol(info, X), props)
        record[name] = dict(props, bonds_mmff=bonds(sym, X0), bonds_dft=bonds(sym, X),
                            computed_smiles=computed)
        print(name, "converged" if r.success else "NOT CONVERGED", r.nit, bonds(sym, X), flush=True)
    if record:
        js = out / "geometries.json"
        old = json.loads(js.read_text()) if js.exists() else {}
        js.write_text(json.dumps({**old, **record}, indent=1) + "\n")


if __name__ == "__main__":
    main()
