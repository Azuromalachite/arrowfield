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


def test_a_record_beside_another_records_surface_is_refused(tmp_path):
    """The NPZ states the geometry it was computed on (geometry_sha1), the record its own: a record
    beside another's surface would draw that surface under its name (#21's probe on dfeffcc)."""
    import shutil
    pair = HERE / "tests" / "data" / "pair" / "e3_z13_18_oac"
    shutil.copy(pair / "conf_00.json", tmp_path / "own.json")
    shutil.copy(pair / "conf_00.npz", tmp_path / "own.npz")
    assert Record(tmp_path / "own").view().provenance.source == "own"
    shutil.copy(pair / "conf_00.json", tmp_path / "swapped.json")
    shutil.copy(pair / "conf_02.npz", tmp_path / "swapped.npz")
    with pytest.raises(InvalidRecord):
        Record(tmp_path / "swapped")
    shutil.copy(pair / "conf_00.json", tmp_path / "rekeyed.json")   # own geometry, another protocol's key
    z = dict(np.load(pair / "conf_00.npz"))
    np.savez(tmp_path / "rekeyed.npz", **{**z, "analysis_key": np.array("2b7bee37" + "0" * 32)})
    with pytest.raises(InvalidRecord):
        Record(tmp_path / "rekeyed")
    for rec in sorted(REC.glob("*.npz")):                      # every record here carries its own
        Record(rec.with_suffix(""))


def test_a_record_whose_key_does_not_digest_its_protocol_is_refused(tmp_path):
    """analysis_key is espmap's sha1 of key_fields.analysis (record_key): a
    record whose fields changed under its key, or that has none, does not say
    what protocol its numbers are under (the user's decision on #21)."""
    j = json.loads((REC / "ethanol__s1.json").read_text())
    assert Record(REC / "ethanol__s1").view().provenance.analysis_fields is not None
    for case in ("changed", "missing", "not_a_mapping"):
        k = json.loads(json.dumps(j))
        if case == "changed":
            k["key_fields"]["analysis"]["solvent"] = "water"
        elif case == "missing":
            del k["key_fields"]
        else:
            k["key_fields"]["analysis"] = "b3lyp/def2svp"
        (tmp_path / f"{case}.json").write_text(json.dumps(k))
        with pytest.raises(InvalidRecord):
            Record(tmp_path / case)
    k = json.loads(json.dumps(j))
    k["settings"]["xc"] = "B3LYP"                     # raw spelling is provenance text, not the protocol
    (tmp_path / "spelled.json").write_text(json.dumps(k))
    assert Record(tmp_path / "spelled").view().provenance.xc == "B3LYP"


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
    assert one.verdict.value == "cannot judge conformational robustness: one geometry on a side"


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
    """folded is the producer's call, used as given: genus is never read as
    folding. Records state no folded flag, so the card says so instead of
    guessing; a caller (a study adapter) states it."""
    e, z = load_all(REC, "e2_butene"), load_all(REC, "z2_butene")
    card = difference_card(e, z, "C1")
    assert card.extended["unstated"] == 6                       # records do not say
    assert "is not stated" in card_html(card)
    e = [replace(v, folded=False) for v in e]
    z = [replace(v, folded=False) for v in z]
    assert difference_card(e, z, "C1").extended is None        # stated: nothing folds
    e[0] = replace(e[0], folded=True, genus=0)                  # folded by contact, genus 0
    card = difference_card(e, z, "C1")
    assert card.extended["n"] == (2, 3)
    assert "extended conformers only" in card_html(card)
    e = [replace(v, folded=False, genus=1) for v in e]          # a hole, nothing folded
    assert difference_card(e, z, "C1").extended is None


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


def test_the_card_names_the_method_of_its_own_measurements():
    """The scope line takes the method from the measurements' canonical
    protocol (key_fields.analysis, which the key digests), never from the raw
    settings' spelling: a hard-coded 'gas phase' would be false on a card
    comparing PCM records."""
    from model import analysis_digest
    e, z = load_all(REC, "e2_butene"), load_all(REC, "z2_butene")
    h = card_html(difference_card(e, z, "C2"))
    assert e[0].provenance.xc == "b3lyp" and "free molecule (b3lyp/def2svp, gas)" in h

    def pcm(v):
        f = json.loads(v.provenance.analysis_fields) | {"solvent": "water"}
        return replace(v, provenance=replace(v.provenance, solvent="water", analysis_key=analysis_digest(f),
                                             analysis_fields=json.dumps(f, sort_keys=True)))
    w = [pcm(v) for v in e + z]
    h = card_html(difference_card(w[:3], w[3:], "C2"))
    assert "PCM water" in h and "gas" not in h.split('class="scope"')[1]


def test_every_real_study_row_renders_as_a_study_card():
    """All rows of the real espmap-sesiidae-001 results (tests/data) render;
    an extended-only row renders beside its own all-member row, and is
    refused beside another site's."""
    import json
    from model import Refused
    from study import load_results
    from views import study_card_html
    r = load_results(HERE / "tests" / "data" / "espmap-sesiidae-001.results.json")
    roles = {s: "atlas_observed" for s in r.molecules} | {"z3_z13_18_al": "teaching_control"}
    ctx = {"provenance": {"results_sha256": r.sha256, "plan_commit": "7e3465859bc2b8950078238f0bffb0ea6b4dd512"},
           "molecules": roles, "limitations": {"basis_sensitivity": "not run", "solvent": "not studied",
                                               "phase": "gas", "method": "B3LYP/def2-SVP"}}
    n = 0
    for s in r.sites:
        c = dict(ctx)
        if s.population != "all":
            c["all_member_row"] = r.site(s.contrast, s.site, "all")
        page = study_card_html(s, c)
        assert "range_partition_sign_v1" in page and "biological discrimination" in page
        n += 1
    assert n == len(r.sites) > 100
    ext = r.site("C_position", "C4-H", "extended")
    with pytest.raises(Refused, match="not this site"):
        study_card_html(ext, dict(ctx, all_member_row=r.site("C_position", "C2-H", "all")))


def test_every_real_record_renders_in_the_evidence_views():
    """The whole homework bundle (Atlas 35112e1): every observation for every
    bundle compound, every contrast, and every compound's classes in both
    taxa render, and none says what the records do not."""
    import re
    from bundle import load
    from views import contrast_card_html, evidence_card_html, evidence_class_panel_html
    b = load(HERE / "tests" / "data" / "sesiidae-homework-1.bundle.json")
    ctx = {"bundle": {"atlas_version": b.atlas_version, "sha256": b.sha256}, "bundle_object": b}
    bad = re.compile(r"inactive|not (made|produced)|repel|inhibitor|antagonist|\boverall\b|on balance|"
                     r"\bESP\b|electrostatic|\bresolved\b|because", re.I)
    pages = [evidence_card_html(st, cid, ctx) for st in b.observations.values() for cid in b.compounds]
    pages += [contrast_card_html(k, [b.observations[o] for o in k.observation_ids], ctx) for k in b.contrasts]
    pages += [evidence_class_panel_html(cid, tax, ctx) for cid in b.compounds
              for tax in ("spa.tax.synanthedon_tenuis", "spa.tax.synanthedon_tipuliformis")]
    assert len(pages) == len(b.observations) * len(b.compounds) + len(b.contrasts) + 2 * len(b.compounds)
    for p in pages:
        assert not bad.search(re.sub(r"<[^>]+>", " ", p)), re.sub(r"<[^>]+>", " ", p)[:300]
