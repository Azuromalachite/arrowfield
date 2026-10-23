"""The ONLY module that knows espmap's record layout (JSON + NPZ).

espmap's record is a versioned contract, "schema": "espmap.record/1"
(docs/record-schema.json, merged in #37). Every key path lives here, so a
change is a one-file edit. A record of any other version, or none, is
refused: the key paths below are only known to hold for this one. The
records in records/ were computed on espmap master 1ef8e28 (2026-10-23);
recompute them on the v0.2.0 tag before class.
"""
import json
from pathlib import Path

import numpy as np

BOHR = 0.529177210903          # NPZ verts/coords are bohr; geometry.xyz_angstrom is A
RECORD_SCHEMA = "espmap.record/1"
RULES = ("nearest_nucleus", "radius_scaled")


class InvalidRecord(ValueError):
    pass


class Record:
    """One espmap run: one molecule, one geometry."""

    def __init__(self, stem):
        stem = Path(stem)
        self.stem = stem
        self.j = json.loads(stem.with_suffix(".json").read_text())
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

    # --- identity and comparability -------------------------------------
    @property
    def analysis_key(self):
        return self.j["analysis_key"]

    @property
    def render_key(self):
        return self.j["render_key"]

    @property
    def symbols(self):
        return self.j["geometry"]["symbols"]

    @property
    def xyz(self):
        return np.array(self.j["geometry"]["xyz_angstrom"], float)

    @property
    def schema(self):
        """None for a record written before espmap.record/1."""
        return self.j.get("schema")

    @property
    def mapped_smiles(self):
        """The approved input with its atom maps: enough to draw the 2D view and
        bond the 3D one from the record alone. None for an unmapped input."""
        return self.j["sites"]["identity"]["mapped_smiles"]

    @property
    def bonds(self):
        """[[i, j, order], ...] in record atom order (espmap #47, A1), or None
        when espmap could not confirm the atom order."""
        return self.j["geometry"]["bonds"]

    @property
    def sites_tile_surface(self):
        """False when caller sites overlap or leave atoms out, so their areas
        do not sum to the surface (espmap #47)."""
        return self.j["sites"]["sites_tile_surface"]

    @property
    def genus(self):
        """Surface genus: > 0 means the molecule folds onto itself."""
        return self.j["results"]["surface"]["genus"]

    @property
    def atom_maps(self):
        return {a["index"]: a["atom_map"] for a in self.j["sites"]["atoms"] if a["atom_map"]}

    # --- scale -----------------------------------------------------------
    @property
    def scale(self):
        s = self.j["scale"]
        return dict(kind=s["kind"], vmin=s["vmin"], vmax=s["vmax"], cmap=s["cmap"])

    # --- sites -----------------------------------------------------------
    def site_labels(self):
        return [s["label"] for s in self.j["sites"]["reaction_sites"]]

    def site(self, label, rule="nearest_nucleus"):
        """A site's numbers in reading order, kcal/mol and A^2."""
        s = next((x for x in self.j["sites"]["reaction_sites"] if x["label"] == label), None)
        if s is None:
            raise KeyError(f"{self.stem.name}: no site {label!r}; sites are {self.site_labels()}")
        r, d = s[rule], s.get("delta") or {}
        return dict(label=label, atoms=s["atoms"],
                    mean=r["V_area_mean_kcal_mol"], area=r["area_A2"],
                    area_pos=r["area_pos_A2"], area_neg=r["area_neg_A2"],
                    partition=(None if d.get("V_area_mean_kcal_mol") is None
                               else abs(d["V_area_mean_kcal_mol"])),
                    p05=r["V_p05_kcal_mol"], p95=r["V_p95_kcal_mol"],
                    vmin=r["V_min_kcal_mol"], vmax=r["V_max_kcal_mol"],
                    h_generated_polar=s.get("h_generated_polar", False),
                    overlaps=s["overlaps"])

    def site_by_map(self, atom_map, rule="nearest_nucleus"):
        """The site whose heavy atom carries this atom map. Labels follow the
        maps on a mapped input (espmap #47, R5), but a page should not have
        to know that: it asks for 'the site of C3'."""
        idx = next((i for i, m in self.atom_maps.items() if m == atom_map), None)
        if idx is None:
            raise KeyError(f"{self.stem.name}: no atom with map {atom_map}")
        hits = [x["label"] for x in self.j["sites"]["reaction_sites"] if idx in x["atoms"]]
        if not hits:
            raise KeyError(f"{self.stem.name}: map {atom_map} is in no site")
        return self.site(hits[0], rule)

    def site_atoms(self, label):
        """The single atoms of a site under both rules. DIAGNOSTIC ONLY: the
        rule moves one atom by up to 11.9 kcal/mol and its area 5.2x (espmap
        #42), so these are shown folded away, never as the site's value."""
        members = set(self.site(label)["atoms"])
        out = []
        for a in self.j["sites"]["atoms"]:
            if a["index"] not in members:
                continue
            row = dict(index=a["index"], element=a["element"], atom_map=a["atom_map"])
            for r in RULES:
                st = a.get(r)
                row[r] = None if st is None else dict(mean=st["V_area_mean_kcal_mol"],
                                                      area=st["area_A2"])
            out.append(row)
        return out

    # --- surface (NPZ) ---------------------------------------------------
    def surface(self, rule="nearest_nucleus"):
        """verts (A), faces, V (kcal/mol) and the atom owning each triangle."""
        if self.z is None:
            raise InvalidRecord(f"{self.stem.name}: no NPZ beside the JSON")
        return dict(verts=self.z["verts"] * BOHR, faces=self.z["faces"], V=self.z["V"],
                    owner=self.z[f"tri_atom_{rule}"])


def comparable(a, b):
    """Numbers may be compared only between records with one analysis key."""
    return a.analysis_key == b.analysis_key


def same_scale(a, b):
    """Pictures may be compared only between records with one render key."""
    return a.render_key == b.render_key


def load_all(pattern_dir, mol_id):
    """Every record of one molecule (its conformers / seeds)."""
    return [Record(p.with_suffix("")) for p in sorted(Path(pattern_dir).glob(f"{mol_id}__*.json"))
            if not p.name.endswith(".sites.json")]          # caller-group files sit beside them
