"""Checks for the view prototype, against the records in records/.

    ESPMAP_PATH=/path/to/espmap pytest tests/
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from record import InvalidRecord, Record, comparable, load_all  # noqa: E402
from views import difference_card, site_faces, svg_2d, viewer_js  # noqa: E402

REC = HERE / "records"
ETHANOL_SMI = "[CH3:1][CH2:2][OH:3]"


def test_a_record_that_says_it_is_invalid_is_refused(tmp_path):
    j = json.loads((REC / "ethanol__s1.json").read_text())
    j["validity"]["surface"] = False
    (tmp_path / "bad.json").write_text(json.dumps(j))
    with pytest.raises(InvalidRecord):
        Record(tmp_path / "bad")


def test_a_record_of_another_schema_version_is_refused(tmp_path):
    j = json.loads((REC / "ethanol__s1.json").read_text())
    assert j["schema"] == "espmap.record/1"
    for bad in ("espmap.record/2", None):
        j["schema"] = bad
        (tmp_path / "old.json").write_text(json.dumps(j))
        with pytest.raises(InvalidRecord, match="schema"):
            Record(tmp_path / "old")


def test_site_numbers_come_straight_from_the_record():
    r = Record(REC / "ethanol__s1")
    j = r.j["sites"]["reaction_sites"]
    s = r.site("hydroxyl")
    raw = next(x for x in j if x["label"] == "hydroxyl")["nearest_nucleus"]
    assert s["mean"] == raw["V_area_mean_kcal_mol"] and s["area"] == raw["area_A2"]
    with pytest.raises(KeyError, match="sites are"):
        r.site("no such site")


def test_the_2d_view_marks_the_site_and_only_when_asked():
    """RDKit draws text as paths, so the map numbers themselves are checked by
    eye on the demo pages; what a test can check is the site highlight."""
    lilac = "#%02X%02X%02X" % tuple(int(c * 255) for c in (0.85, 0.65, 0.81))
    assert lilac in svg_2d(ETHANOL_SMI, site_maps=[3]).upper()
    assert lilac not in svg_2d(ETHANOL_SMI).upper()


def test_the_highlighted_patch_is_exactly_the_site_triangles():
    r = Record(REC / "ethanol__s1")
    s, site = r.surface(), r.site("hydroxyl")
    patch = site_faces(s, site["atoms"])
    assert patch.sum() == np.isin(s["owner"], site["atoms"]).sum() > 0
    # and the patch's area is the site's area in the record
    P, f = s["verts"], s["faces"][patch]
    area = 0.5 * np.linalg.norm(np.cross(P[f[:, 1]] - P[f[:, 0]], P[f[:, 2]] - P[f[:, 0]]), axis=1).sum()
    assert area == pytest.approx(site["area"], rel=1e-6)


def test_the_3d_view_uses_the_house_style():
    js = viewer_js("v", Record(REC / "ethanol__s1"), ETHANOL_SMI, [2], show_surface=False)
    assert "style:'outline', color:'black', width:0.03" in js
    assert "'O': '#d62728'" in js.replace('"', "'")
    assert "addSphere" in js                                  # the O site is visibly marked


def test_difference_card_verdicts():
    e, z = load_all(REC, "e2_butene"), load_all(REC, "z2_butene")
    alkene = e[0].site_by_map(2)["label"]            # the alkene C, by map not label
    assert alkene == "C2"                            # map-based labels (espmap #47)
    card = difference_card(e, z, alkene)
    assert card["n"] == (3, 3) and card["verdict"].startswith("separated")
    assert card["diff"] == pytest.approx(-1.29, abs=0.01)   # the value before #47
    one = difference_card(e[:1], z[:1], alkene)
    assert one["verdict"] == "single conformer: cannot judge"


def test_site_table_leads_with_the_site_and_folds_atoms_away():
    """Agreed policy (espmap #42): the site mean first, then partition and
    conformer sensitivity; single atoms, tails and extrema only inside
    <details>, never in the first table a student sees."""
    from views import site_table
    recs = load_all(REC, "ethanol")
    t = site_table(recs, ["hydroxyl"])
    primary, advanced = t.split("<details>")
    assert "conformer range" in primary and "partition" in primary
    assert "p05" not in primary and "V<sub>min</sub>" not in primary and "O2" not in primary
    assert "diagnostic" in advanced and "O2" in advanced     # the hydroxyl O, by index
    assert "1 conformer" in site_table(recs[0], ["hydroxyl"])


def test_card_shows_the_extended_only_range_when_a_conformer_folds():
    e, z = load_all(REC, "e2_butene"), load_all(REC, "z2_butene")
    assert "extended" not in difference_card(e, z, "C1")      # nothing folds here
    e[0].j["results"]["surface"]["genus"] = 1                  # pretend one folded
    card = difference_card(e, z, "C1")
    assert card["extended"]["n"] == (2, 3)
    from views import card_html
    assert "extended conformers only" in card_html(card)


def test_sites_file_says_it_names_atoms_by_map():
    sys.argv = sys.argv[:1] + ["unused"]                       # compute_records reads argv[1]
    import compute_records
    assert compute_records.sites_file([{"label": "x", "atoms": [3]}]) == {
        "by": "map", "groups": [{"label": "x", "atoms": [3]}]}


def test_head_o_is_one_definition_in_every_c18_molecule():
    """The pilot compared O31 across OH (O+H) and OAc (bare O). 'head O' is now
    the O bonded to C1, atom alone, in every head group."""
    import yaml
    mols = yaml.safe_load((HERE / "molecules.yaml").read_text(encoding="utf-8"))["molecules"]
    for m in mols:
        if m.get("headgroup"):
            head = next(s for s in m["sites"] if s["label"] == "head O")
            assert len(head["atoms"]) == 1 and isinstance(head["atoms"][0], int), m["id"]


def test_the_record_draws_itself_from_its_own_bonds():
    """espmap #47 (A1): geometry.bonds needs no atom maps, so an unmapped
    record draws too. The bonds route and the mapped-SMILES route agree."""
    from rdkit import Chem
    from views import mol_from_record
    r = Record(REC / "ethanol__s1")
    assert r.bonds is not None
    a = mol_from_record(r)
    b = Chem.Mol(mol_from_record(r, ETHANOL_SMI))
    r.j["geometry"]["bonds"] = None                       # force the SMILES route
    b = mol_from_record(r, ETHANOL_SMI)
    key = lambda m: sorted((x.GetBeginAtomIdx(), x.GetEndAtomIdx(), x.GetBondTypeAsDouble())
                           for x in m.GetBonds())          # noqa: E731
    assert key(a) == key(b)


def test_acetate_charge_comes_back_from_its_bonds():
    from views import mol_from_record
    m = mol_from_record(Record(REC / f"acetate__s{12648430}"))
    assert sum(a.GetFormalCharge() for a in m.GetAtoms()) == -1


def test_sites_are_found_by_atom_map_not_label():
    r = Record(REC / "ethanol__s1")
    assert r.site_by_map(3)["label"] == "hydroxyl"        # the caller group holding O3


def test_overlapping_sites_are_marked_in_the_table():
    from views import site_table
    r = Record(REC / "ethanol__s1")
    r.j["sites"]["reaction_sites"][0]["overlaps"] = ["C-O carbon"]
    r.j["sites"]["sites_tile_surface"] = False
    t = site_table(r, r.site_labels()).split("<details>")[0]
    assert "†" in t and "do not add up to the surface" in t


def test_the_2d_view_draws_from_the_record_alone():
    """No SMILES and no atom maps: an unmapped record still draws its
    formula from geometry.bonds and highlights its site by record index."""
    from views import svg_2d_record
    lilac = "#%02X%02X%02X" % tuple(int(c * 255) for c in (0.85, 0.65, 0.81))
    r = Record(REC / "ethanol__s1")
    site = r.site("hydroxyl")["atoms"]
    assert lilac in svg_2d_record(r, site).upper()
    r.j["sites"]["identity"]["mapped_smiles"] = None
    for a in r.j["sites"]["atoms"]:
        a["atom_map"] = None
    svg = svg_2d_record(r, site)
    assert lilac in svg.upper() and "<svg" in svg


def test_records_that_are_not_comparable_are_refused():
    e = load_all(REC, "e2_butene")
    acetate = Record(REC / f"acetate__s{12648430}")             # PCM water: another analysis key
    assert not comparable(e[0], acetate)
    with pytest.raises(ValueError, match="analysis keys"):
        difference_card(e, [acetate], "C1")
