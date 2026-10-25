"""The ONLY module that knows espmap's record layout (JSON + NPZ).

espmap's record is a versioned contract, "schema": "espmap.record/1"
(docs/record-schema.json, merged in #37). Every key path lives here, and it
is translated into Arrowfield's own model (model.py) at the door: views and
lessons receive a MoleculeView and never see the record. A change in
espmap's layout is a one-file edit. A record of any other version, or none,
is refused: the key paths below are only known to hold for this one. The
records in records/ were computed on espmap master 1ef8e28 (2026-10-23);
recompute them on the v0.2.1 tag before class.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from model import (RULES, Atom, MoleculeView, Provenance, Scale, SiteMeasurement, Surface,
                   analysis_digest, comparable, same_scale)  # noqa: F401  (re-exported for callers)

BOHR = 0.529177210903          # NPZ verts/coords are bohr; geometry.xyz_angstrom is A
RECORD_SCHEMA = "espmap.record/1"


class InvalidRecord(ValueError):
    pass


class Record:
    """One espmap run: one molecule, one geometry. Validates on load; .view()
    gives the model object everything else uses."""

    def __init__(self, stem):
        stem = Path(stem)
        self.stem = stem
        raw = stem.with_suffix(".json").read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        self.j = json.loads(raw)
        npz = stem.with_suffix(".npz")
        self.z = np.load(npz) if npz.exists() else None
        if self.j.get("schema") != RECORD_SCHEMA:
            raise InvalidRecord(f"{stem.name}: schema {self.j.get('schema')!r}, this reader "
                                f"knows {RECORD_SCHEMA!r}; recompute the record")
        v = self.j.get("validity", {})
        if not (v.get("input") and v.get("surface") and v.get("sites")):
            # the record says whether it can be used; never infer it
            raise InvalidRecord(f"{stem.name}: validity {v}, sites_reason "
                                f"{self.j.get('sites_reason')!r}")
        # the canonical protocol, verified: analysis_key is espmap's digest of key_fields.analysis
        fields = (self.j.get("key_fields") or {}).get("analysis")
        if not isinstance(fields, dict):
            raise InvalidRecord(f"{stem.name}: no canonical analysis fields (key_fields.analysis)")
        if analysis_digest(fields) != self.j.get("analysis_key"):
            raise InvalidRecord(f"{stem.name}: analysis_key {self.j.get('analysis_key')} is not the digest of "
                                "its key_fields.analysis: the key does not describe this record")

    def view(self, name=None):
        """The record as Arrowfield's MoleculeView."""
        j = self.j
        g, s, res = j["geometry"], j["sites"], j["results"]
        atoms = tuple(Atom(index=a["index"], element=a["element"], atom_map=a["atom_map"],
                           parent=a.get("parent"),
                           by_rule={r: (None if a.get(r) is None else
                                        (a[r]["V_area_mean_kcal_mol"], a[r]["area_A2"]))
                                    for r in RULES})
                      for a in s["atoms"])
        sites = {}
        for x in s["reaction_sites"]:
            d = x.get("delta") or {}
            part = (None if d.get("V_area_mean_kcal_mol") is None
                    else abs(d["V_area_mean_kcal_mol"]))
            sites[x["label"]] = {
                r: SiteMeasurement(
                    label=x["label"], atoms=tuple(x["atoms"]), rule=r,
                    mean=x[r]["V_area_mean_kcal_mol"], area=x[r]["area_A2"],
                    area_pos=x[r]["area_pos_A2"], area_neg=x[r]["area_neg_A2"],
                    partition=part, p05=x[r]["V_p05_kcal_mol"], p95=x[r]["V_p95_kcal_mol"],
                    vmin=x[r]["V_min_kcal_mol"], vmax=x[r]["V_max_kcal_mol"],
                    h_generated_polar=x.get("h_generated_polar", False),
                    overlaps=tuple(x["overlaps"]))
                for r in RULES if x.get(r) is not None}
        sc, st, pv = j["scale"], j.get("settings", {}), j.get("provenance", {})
        return MoleculeView(
            name=name or self.stem.name,
            symbols=tuple(g["symbols"]),
            xyz=np.array(g["xyz_angstrom"], float),
            bonds=None if g["bonds"] is None else tuple(tuple(b) for b in g["bonds"]),
            charge=g["charge"],
            atoms=atoms,
            sites=sites,
            sites_tile_surface=s["sites_tile_surface"],
            genus=res["surface"]["genus"],
            scale=Scale(kind=sc["kind"], vmin=sc["vmin"], vmax=sc["vmax"], cmap=sc["cmap"]),
            provenance=Provenance(
                source=self.stem.name, analysis_key=j["analysis_key"], render_key=j["render_key"],
                record_schema=j.get("schema"), espmap_version=pv.get("espmap_version"),
                code=pv.get("code"), sha256=self.sha256, xc=st.get("xc"),
                basis=st.get("basis"), solvent=st.get("solvent"),
                analysis_fields=json.dumps(j["key_fields"]["analysis"], sort_keys=True)),
            mapped_smiles=s["identity"]["mapped_smiles"],
            surface_loader=None if self.z is None else self._surface)

    def _surface(self, rule):
        return Surface(verts=self.z["verts"] * BOHR, faces=self.z["faces"], V=self.z["V"],
                       owner=self.z[f"tri_atom_{rule}"])


def load(stem, name=None):
    """One record, as a MoleculeView."""
    return Record(stem).view(name)


def load_all(pattern_dir, mol_id):
    """Every record of one molecule (its conformers / seeds), as MoleculeViews."""
    return [Record(p.with_suffix("")).view()
            for p in sorted(Path(pattern_dir).glob(f"{mol_id}__*.json"))
            if not p.name.endswith(".sites.json")]          # caller-group files sit beside them
