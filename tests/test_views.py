"""Checks for the view prototype, against the records in records/.

    ESPMAP_PATH=/path/to/espmap pytest tests/

Views take model objects (model.py). A state a test needs (a folded
conformer, overlapping sites, an unmapped record) is made by replacing
fields of a MoleculeView, never by editing an espmap record's JSON: the
tests check Arrowfield's semantics, not espmap's layout.
"""
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from model import ComparisonStatus, Verdict, comparable  # noqa: E402
from record import InvalidRecord, Record, load, load_all  # noqa: E402
from views import (card_html, difference_card, mol_from_view, site_faces, site_table,  # noqa: E402
                   svg_2d, svg_2d_view, viewer_js)

REC = HERE / "records"
ETHANOL_SMI = "[CH3:1][CH2:2][OH:3]"
LILAC = "#%02X%02X%02X" % tuple(int(c * 255) for c in (0.85, 0.65, 0.81))


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
    raw = next(x for x in r.j["sites"]["reaction_sites"] if x["label"] == "hydroxyl")["nearest_nucleus"]
    s = r.view().site("hydroxyl")
    assert s.mean == raw["V_area_mean_kcal_mol"] and s.area == raw["area_A2"]
    with pytest.raises(KeyError, match="sites are"):
        r.view().site("no such site")


def test_the_view_carries_its_provenance():
    v = load(REC / "ethanol__s1")
    p = v.provenance
    assert p.analysis_key and p.render_key and p.record_schema == "espmap.record/1"
    assert len(p.sha256) == 64 and p.source == "ethanol__s1"


def test_the_2d_view_marks_the_site_and_only_when_asked():
    """RDKit draws text as paths, so the map numbers themselves are checked by
    eye on the demo pages; what a test can check is the site highlight."""
    assert LILAC in svg_2d(ETHANOL_SMI, site_maps=[3]).upper()
    assert LILAC not in svg_2d(ETHANOL_SMI).upper()


def test_the_highlighted_patch_is_exactly_the_site_triangles():
    v = load(REC / "ethanol__s1")
    s, site = v.surface(), v.site("hydroxyl")
    patch = site_faces(s, site.atoms)
    assert patch.sum() == np.isin(s.owner, site.atoms).sum() > 0
    # and the patch's area is the site's area in the record
    P, f = s.verts, s.faces[patch]
    area = 0.5 * np.linalg.norm(np.cross(P[f[:, 1]] - P[f[:, 0]], P[f[:, 2]] - P[f[:, 0]]), axis=1).sum()
    assert area == pytest.approx(site.area, rel=1e-6)


def test_the_3d_view_uses_the_house_style():
    js = viewer_js("v", load(REC / "ethanol__s1"), ETHANOL_SMI, [2], show_surface=False)
    assert "style:'outline', color:'black', width:0.03" in js
    assert "'O': '#d62728'" in js.replace('"', "'")
    assert "addSphere" in js                                  # the O site is visibly marked


def test_difference_card_verdicts():
    e, z = load_all(REC, "e2_butene"), load_all(REC, "z2_butene")
    alkene = e[0].site_by_map(2).label              # the alkene C, by map not label
    assert alkene == "C2"                            # map-based labels (espmap #47)
    card = difference_card(e, z, alkene)
    assert card.n == (3, 3) and card.status is ComparisonStatus.COMPARABLE
    assert card.verdict in (Verdict.SEPARATED, Verdict.SEPARATED_WITHIN_PARTITION)
    assert card.diff == pytest.approx(-1.29, abs=0.01)   # the value before #47
    one = difference_card(e[:1], z[:1], alkene)
    assert one.status is ComparisonStatus.SINGLE_CONFORMER
    assert one.verdict is Verdict.CANNOT_JUDGE
    assert one.verdict.value == "single conformer: cannot judge"


def test_site_table_leads_with_the_site_and_folds_atoms_away():
    """Agreed policy (espmap #42): the site mean first, then partition and
    conformer sensitivity; single atoms, tails and extrema only inside
    <details>, never in the first table a student sees."""
    recs = load_all(REC, "ethanol")
    t = site_table(recs, ["hydroxyl"])
    primary, advanced = t.split("<details>")
    assert "conformer range" in primary and "partition" in primary
    assert "p05" not in primary and "V<sub>min</sub>" not in primary and "O2" not in primary
    assert "diagnostic" in advanced and "O2" in advanced     # the hydroxyl O, by index
    assert "1 conformer" in site_table(recs[0], ["hydroxyl"])


def test_card_shows_the_extended_only_range_when_a_conformer_folds():
    e, z = load_all(REC, "e2_butene"), load_all(REC, "z2_butene")
    assert difference_card(e, z, "C1").extended is None        # nothing folds here
    e[0] = replace(e[0], genus=1)                               # one conformer folds
    card = difference_card(e, z, "C1")
    assert card.extended["n"] == (2, 3)
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
    v = load(REC / "ethanol__s1")
    assert v.bonds is not None
    a = mol_from_view(v)
    b = mol_from_view(replace(v, bonds=None), ETHANOL_SMI)    # force the SMILES route
    key = lambda m: sorted((x.GetBeginAtomIdx(), x.GetEndAtomIdx(), x.GetBondTypeAsDouble())
                           for x in m.GetBonds())          # noqa: E731
    assert key(a) == key(b)


def test_acetate_charge_comes_back_from_its_bonds():
    m = mol_from_view(load(REC / f"acetate__s{12648430}"))
    assert sum(a.GetFormalCharge() for a in m.GetAtoms()) == -1


def test_sites_are_found_by_atom_map_not_label():
    assert load(REC / "ethanol__s1").site_by_map(3).label == "hydroxyl"   # the group holding O3


def test_overlapping_sites_are_marked_in_the_table():
    v = load(REC / "ethanol__s1")
    first = v.site_labels()[0]
    sites = dict(v.sites)
    sites[first] = {r: replace(m, overlaps=("C-O carbon",)) for r, m in sites[first].items()}
    v = replace(v, sites=sites, sites_tile_surface=False)
    t = site_table(v, v.site_labels()).split("<details>")[0]
    assert "†" in t and "do not add up to the surface" in t


def test_the_2d_view_draws_from_the_record_alone():
    """No SMILES and no atom maps: an unmapped record still draws its
    formula from its bonds and highlights its site by record index."""
    v = load(REC / "ethanol__s1")
    site = v.site("hydroxyl").atoms
    assert LILAC in svg_2d_view(v, site).upper()
    unmapped = replace(v, mapped_smiles=None,
                       atoms=tuple(replace(a, atom_map=None) for a in v.atoms))
    svg = svg_2d_view(unmapped, site)
    assert LILAC in svg.upper() and "<svg" in svg


def test_records_that_are_not_comparable_are_shown_as_not_comparable():
    """A comparison across analysis keys is never computed: the card says
    NOT COMPARABLE and shows no numbers, so a page cannot quietly drop it."""
    e = load_all(REC, "e2_butene")
    acetate = load(REC / f"acetate__s{12648430}")              # PCM water: another analysis key
    assert not comparable(e[0], acetate)
    card = difference_card(e, [acetate], "C1")
    assert card.status is ComparisonStatus.NOT_COMPARABLE and card.diff is None
    h = card_html(card)
    assert "NOT COMPARABLE" in h and "difference of conformer means" not in h
