"""Generate the hostile set for the pair view: tests/pedagogy/pair_view_set.json.

    python tests/pedagogy/make_pair_view_set.py     # byte-identical on re-run

Written before the view exists (expectations first). The pair view puts two
ESP pictures side by side: two conformers of one compound (the homework's
Act B, item 2: "the extended and folded conformers of one compound side by
side, on the same scale"), or two compounds (lessons.md: "both butenes on
one shared scale"; notebooks.md: "the control molecule on the same scale").
It re-expresses the two abstract C1 cases of hostile_set.json
(pair_different_render_keys, pair_one_folded) against real records; the
view's PR retires them, each naming its replacements here (RETIRES below).

The view, in views.py, raising model.Refused when the pair cannot be shown
faithfully:

    pair_view_html(left, right, context) -> str

left, right: MoleculeView, one geometry each, with its surface. context:
optional {"labels": {"left": str, "right": str}} (display text only, never
evidence) and "study_rows" (verdicts the page also shows; never linked).

The rendered page carries, for the harness and for any reader's tools:
- a root element with data-relation: "conformers_of_one_compound",
  "same_geometry" or "two_compounds", decided from the structures
  (constitution, charge, and stereo assigned from 3D; the same geometry up
  to rotation and symmetry-equivalent atoms), never from names;
- exactly one element with class "pair-scale" and data-kind, data-vmin,
  data-vmax, data-cmap: the one scale both pictures use;
- exactly two elements with class "pair-panel", each with data-vmin and
  data-vmax (the scale its colours are computed on) and data-folded: "true",
  "false" or "unstated", the producer's flag as given (MoleculeView.folded,
  never inferred from genus).

The shared scale contains both records' own scales (vmin no higher than
either's, vmax no lower), is finite with vmin < vmax, and is diverging, and
then symmetric, when either record's scale is diverging.

Inputs are real records: the committed records/ (espmap v0.1.0, gas phase
and one PCM anion), and four Atlas conformers in tests/data/pair/ (espmap
v0.2.1; folded as the Atlas index states, members.json). A side may carry a
`patch` ((field, value), applied with dataclasses.replace; "scale.vmin" and
"provenance.x" reach one level down) and a `surface` op (none, nan, inf,
empty) for synthetic variants.

Grades as in the other sets (SAFE, SAFE_REFUSAL, CRASH, MISLEADING, FALSE);
`allowed` is exactly the safe outputs. A rendered page must meet every
`requires` entry: must / must_not (a case-insensitive pattern on the text
with tags removed, or on the HTML), or check (a machine check on the data
attributes: shared_scale, relation=<value>, folded=<left>,<right>).
"""
import json
import re
from pathlib import Path

OUT = Path(__file__).with_name("pair_view_set.json")
HERE = Path(__file__).resolve().parents[1]
GRADES = ("SAFE", "SAFE_REFUSAL", "CRASH", "MISLEADING", "FALSE")
MEMBERS = json.loads((HERE / "data" / "pair" / "members.json").read_text(encoding="utf-8"))["members"]
FOLDED = {m["record"][:-5]: m["folded"] for m in MEMBERS}
CASES = []


def key(stem, which):
    return json.loads((HERE.parent / f"{stem}.json").read_text(encoding="utf-8"))[which]

# records
ETOH, ETOH2, ACAL, ACAL2 = "records/ethanol__s1", "records/ethanol__s2", "records/acetaldehyde__s1", "records/acetaldehyde__s2"
E2B, E2B2, Z2B = "records/e2_butene__s1", "records/e2_butene__s2", "records/z2_butene__s1"
ACO = "records/acetate__s12648430"
E3F, E3X = "tests/data/pair/e3_z13_18_oac/conf_00", "tests/data/pair/e3_z13_18_oac/conf_02"   # folded, extended
Z3X, Z3F = "tests/data/pair/z3_z13_18_oac/conf_00", "tests/data/pair/z3_z13_18_oac/conf_01"   # extended, folded


def side(stem, patch=(), surface=None):
    folded = FOLDED.get(stem.removeprefix("tests/data/pair/"))
    p = [list(x) for x in patch]
    if folded is not None and not any(x[0] == "folded" for x in p):
        p = [["folded", folded]] + p          # the Atlas index's flag, as study.load_ensemble sets it
    return {"record": stem, "patch": p, "surface": surface}


def case(id, cls, why, left, right, outcomes, prefer, requires=(), context=None, same_object=False):
    outcomes = dict(outcomes)
    outcomes.setdefault("ERROR", "CRASH")
    allowed = sorted(k for k, g in outcomes.items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert prefer is None or prefer in allowed, id
    CASES.append(dict(id=id, cls=cls, why=why, left=left, right=right, same_object=same_object,
                      context=context or {}, outcomes=outcomes, allowed=allowed, prefer=prefer,
                      requires=list(requires)))


def must(pat, grade, why, on="text"):
    return {"must": pat, "grade": grade, "why": why, "on": on}


def must_not(pat, grade, why, on="text"):
    return {"must_not": pat, "grade": grade, "why": why, "on": on}


def check(what, grade, why):
    return {"check": what, "grade": grade, "why": why}


SHOW = {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}
REFUSE_M = {"RENDERED": "MISLEADING", "REFUSED": "SAFE_REFUSAL"}
REFUSE_F = {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"}
SHARED = check("shared_scale", "MISLEADING", "one scale, containing both records' own, for both pictures")
UNITS = must(r"kcal/mol", "MISLEADING", "the scale's unit")
RANGE = must(r"[−-]\s?\d+(\.\d+)?\s*(to|…|–|\.\.)\s*\+?\d+(\.\d+)?", "MISLEADING", "the shared range, in numbers")
CONVENTION = must(r"(red|warm)[^.]{0,60}negative|negative[^.]{0,60}(red|warm)", "MISLEADING",
                  "which colour is negative: the chemistry convention, said (CLAUDE.md conventions)")
METHOD = must(r"b3lyp", "MISLEADING", "the method both pictures were computed with")
# a picture pair is not evidence about biology, and says nothing about why
NO_CAUSE = must_not(r"because|explain|due to|caus|more reactive|reactivity|attract|binds?\b|binding|"
                    r"pheromone|receptor|\bresolved\b|specific", "FALSE",
                    "two pictures show where the potential differs; they test no mechanism")
CAPTION = [SHARED, UNITS, RANGE, CONVENTION, METHOD, NO_CAUSE]

# --- P1 one scale for both pictures ------------------------------------------------------
case("butene_e_vs_z", "P1 one scale",
     "E- and Z-2-butene, same method, auto-scaled to ±17.3 and ±18.1 kcal/mol (different render keys). "
     "Each on its own scale, the same colour means different kcal/mol and the eye compares the wrong "
     "things: both go on one scale.",
     side(E2B), side(Z2B), SHOW, "RENDERED", CAPTION + [check("relation=two_compounds", "FALSE", "two compounds")])
case("ethanol_vs_acetaldehyde", "P1 one scale",
     "Ethanol (±50.0) and acetaldehyde (±39.3), same method: very different own scales; one shared scale.",
     side(ETOH), side(ACAL), SHOW, "RENDERED", CAPTION)
case("near_equal_scales", "P1 one scale",
     "Two conformers of the E3 acetate at ±40.33 and ±40.77 kcal/mol: close enough to tempt 'the same "
     "scale', different enough that one colour is two values. Still one shared scale.",
     side(E3F), side(E3X), SHOW, "RENDERED", CAPTION)
case("identical_render_keys", "P1 one scale",
     "Two records with one render key (the same scale): the shared scale is that scale.",
     side(E2B), side(E2B2), SHOW, "RENDERED", [SHARED, UNITS, CONVENTION])
case("forged_render_key", "P1 one scale",
     "Ethanol relabelled with acetaldehyde's render key: the key claims one scale, the scales are ±50.0 and "
     "±39.3. A view that trusts the key draws each on its own scale. The scales decide, not the label.",
     side(ETOH, patch=[("provenance.render_key", key(ACAL, "render_key"))]), side(ACAL), SHOW, "RENDERED",
     [SHARED])
case("sequential_scales", "P1 one scale",
     "Two anion pictures never change sign, so their scales are sequential (acetate, −183.4 to −70.4; the "
     "other side synthetic, −150 to −60): one sequential scale over both, never centred on zero, and the "
     "caption says darker is more negative.",
     side(ACO), side(ACO, patch=[("scale.vmin", -150.0), ("scale.vmax", -60.0)]), SHOW, "RENDERED",
     [SHARED, UNITS, must(r"darker[^.]{0,60}negative|negative[^.]{0,60}darker", "MISLEADING",
                          "which end is negative, on a sequential scale")])

# --- P2 one method: analysis keys and the provenance behind them -----------------------------
case("pcm_anion_vs_gas_neutral", "P2 one method",
     "Acetate in PCM water beside ethanol in the gas phase: different analysis keys. One colour would mean "
     "two methods' kcal/mol; refuse, and say to compute both under one protocol.",
     side(ACO), side(ETOH), REFUSE_M, "REFUSED")
case("basis_differs", "P2 one method",
     "Ethanol beside acetaldehyde stated as recomputed at def2-TZVP (its key and its basis both say so): "
     "not one method.",
     side(ETOH), side(ACAL, patch=[("provenance.analysis_key", "f" * 64), ("provenance.basis", "def2-tzvp")]),
     REFUSE_M, "REFUSED")
case("stale_key_solvent", "P2 one method",
     "The keys match, but one side's provenance says PCM water: the key does not describe that record. "
     "Captioned as one gas-phase method, the page would be false.",
     side(ETOH), side(ACAL, patch=[("provenance.solvent", "water")]), REFUSE_F, "REFUSED")
case("stale_key_functional", "P2 one method",
     "The keys match, but one side's provenance names PBE0: a key that disagrees with its own record.",
     side(E2B), side(Z2B, patch=[("provenance.xc", "pbe0")]), REFUSE_F, "REFUSED")
case("ensemble_vs_pilot_versions", "P2 one method",
     "A v0.2.1 Atlas conformer beside a v0.1.0 pilot record, both B3LYP/def2-SVP in the gas phase: the "
     "captions would read alike, but the analysis keys differ. The key decides; refuse.",
     side(E3X), side(ETOH), REFUSE_M, "REFUSED")

# --- P3 what the pair is, from the structures ---------------------------------------------
case("e3_vs_z3_named_alike", "P3 relation from structure",
     "E3 and Z3 acetate, both extended, both named 'conf_0x' by the loader. They differ at one double "
     "bond: two compounds, whatever the names say.",
     side(E3X), side(Z3X), SHOW, "RENDERED",
     CAPTION + [check("relation=two_compounds", "FALSE", "two compounds"),
                must_not(r"conformers of (one|the same) compound|same compound", "FALSE", "not one compound")])
case("e3_vs_z3_labelled_one_compound", "P3 relation from structure",
     "The same pair with display labels 'E3 acetate, conformer 1' and 'E3 acetate, conformer 2': labels are "
     "display text; the structures still say two compounds.",
     side(E3X), side(Z3X), SHOW, "RENDERED",
     [check("relation=two_compounds", "FALSE", "two compounds, from the structures")],
     context={"labels": {"left": "E3 acetate, conformer 1", "right": "E3 acetate, conformer 2"}})
case("true_conformers", "P3 relation from structure",
     "Two genuinely different geometries of the E3 acetate: two conformers of one compound.",
     side(E3F), side(E3X), SHOW, "RENDERED",
     [check("relation=conformers_of_one_compound", "MISLEADING", "one compound, two geometries"),
      must(r"conformers? of (one|the same) compound|same compound", "MISLEADING", "said in words")])
case("identical_coordinates", "P3 relation from structure",
     "Two seed records of E-2-butene with identical coordinates: one geometry twice. Shown as 'two "
     "conformers' it would look like evidence that the picture is robust to conformation.",
     side(E2B), side(E2B2), SHOW, "RENDERED",
     [check("relation=same_geometry", "MISLEADING", "the same geometry"),
      must(r"same geometry|one geometry|identical geometr", "MISLEADING", "said in words")])
case("same_geometry_rotated", "P3 relation from structure",
     "Two acetaldehyde seed records: the same geometry (distances agree to 1.5e-6 Å), rotated by up to 2.2 Å "
     "in place. An exact-coordinate test calls them different.",
     side(ACAL), side(ACAL2), SHOW, "RENDERED", [check("relation=same_geometry", "MISLEADING", "the same geometry")])
case("same_geometry_relabelled_h", "P3 relation from structure",
     "Two ethanol seed records: the same geometry with the methyl hydrogens' labels permuted (interatomic "
     "distances differ by 0.70 Å in atom order, 0 after symmetry). An atom-order test calls them conformers.",
     side(ETOH), side(ETOH2), SHOW, "RENDERED", [check("relation=same_geometry", "MISLEADING", "the same geometry")])
case("same_object_twice", "P3 relation from structure",
     "The same MoleculeView passed as both sides: nothing to compare. Refuse, or say it is one geometry.",
     side(ETOH), side(ETOH), SHOW, "REFUSED", [check("relation=same_geometry", "MISLEADING", "the same geometry")],
     same_object=True)

# --- P4 folded and extended -------------------------------------------------------------------
case("one_compound_folded_vs_extended", "P4 folded",
     "Act B, item 2: the E3 acetate folded (conf_00) beside extended (conf_02), on one scale, each labelled "
     "as the producer states it. conf_00 has genus 0: a view inferring 'folded' from genus calls it extended.",
     side(E3F), side(E3X), SHOW, "RENDERED",
     CAPTION + [check("folded=true,false", "FALSE", "the producer's flags, as given"),
                must(r"\bfolded\b", "MISLEADING", "the folded side named"),
                must(r"\bextended\b", "MISLEADING", "the extended side named")])
case("one_compound_folded_genus_3", "P4 folded",
     "The Z3 acetate extended (conf_00) beside folded (conf_01, genus 3): the same view, the other compound.",
     side(Z3X), side(Z3F), SHOW, "RENDERED",
     [SHARED, check("folded=false,true", "FALSE", "the producer's flags, as given")])
case("two_compounds_folding_differs", "P4 folded",
     "E3 extended beside Z3 folded: the shape difference from folding would read as the E/Z difference. "
     "Refuse, and say to pick conformers with the same folded state (or read the difference card's ranges).",
     side(E3X), side(Z3F), REFUSE_M, "REFUSED")
case("two_compounds_folding_differs_reversed", "P4 folded",
     "The same with the sides swapped and the other conformers: E3 folded beside Z3 extended.",
     side(E3F), side(Z3X), REFUSE_M, "REFUSED")
case("folded_not_stated", "P4 folded",
     "Pilot records state no folded flag (None): each panel says the state is not stated, never 'extended'.",
     side(E2B), side(Z2B), SHOW, "RENDERED",
     [check("folded=unstated,unstated", "FALSE", "unstated stays unstated"),
      must(r"not stated", "MISLEADING", "said in words"), must_not(r"\bextended\b", "FALSE", "not invented")])
case("folded_flag_not_boolean", "P4 folded",
     "A folded flag of 'yes' (a string): outside the vocabulary. Refuse rather than guess.",
     side(E3X, patch=[("folded", "yes")]), side(E3F), REFUSE_M, "REFUSED")

# --- P5 values a picture cannot be drawn from ------------------------------------------------
case("no_surface", "P5 values",
     "A record without its NPZ: no surface to draw. Refuse, naming the missing file.",
     side(ETOH, surface="none"), side(ACAL), REFUSE_M, "REFUSED")
case("nan_in_potential", "P5 values",
     "One vertex's potential is NaN: its colour is undefined, and min/max would skip it. Refuse.",
     side(ETOH), side(ACAL, surface="nan"), REFUSE_M, "REFUSED")
case("inf_in_potential", "P5 values",
     "One vertex's potential is +inf: a shared scale computed from it is infinite. Refuse.",
     side(ETOH, surface="inf"), side(ACAL), REFUSE_M, "REFUSED")
case("empty_surface", "P5 values",
     "A surface with no triangles: nothing to draw. Refuse.",
     side(ETOH), side(ACAL, surface="empty"), REFUSE_M, "REFUSED")
case("scale_not_finite", "P5 values",
     "A record whose scale vmin is NaN: the shared scale cannot contain it. Refuse.",
     side(ETOH, patch=[("scale.vmin", float("nan"))]), side(ACAL), REFUSE_M, "REFUSED")
case("degenerate_scale", "P5 values",
     "Both records' scales are 0 to 0: no range to colour on. Refuse; never divide by zero.",
     side(E2B, patch=[("scale.vmin", 0.0), ("scale.vmax", 0.0)]),
     side(Z2B, patch=[("scale.vmin", 0.0), ("scale.vmax", 0.0)]), REFUSE_M, "REFUSED")

# --- P6 the caption's boundary ----------------------------------------------------------------
case("study_verdict_on_the_page", "P6 caption",
     "The page also shows the study's verdict for E3 vs Z3. The picture pair must not borrow it: no "
     "'resolved', no 'because', nothing about attraction.",
     side(E3X), side(Z3X), SHOW, "RENDERED", [NO_CAUSE, SHARED],
     context={"study_rows": [{"verdict_rule": "espmap-sesiidae-001/range_partition_sign_v1",
                              "contrast": "B_ez", "site": "C3-H", "outcome": "different_on_average"}]})
case("label_is_markup", "P6 caption",
     "A display label containing markup is text, escaped, never executed.",
     side(E2B), side(Z2B), SHOW, "RENDERED",
     [must_not(r"<script>alert", "FALSE", "escaped", on="html"), must(r"&lt;script&gt;", "MISLEADING",
                                                                       "shown as text", on="html")],
     context={"labels": {"left": "<script>alert(1)</script>", "right": "Z-2-butene"}})

# Retired by the view's PR (not edited now): each abstract C1 case names its replacements.
RETIRES = {"pair_different_render_keys": ["butene_e_vs_z", "ethanol_vs_acetaldehyde", "near_equal_scales",
                                          "forged_render_key"],
           "pair_one_folded": ["two_compounds_folding_differs", "two_compounds_folding_differs_reversed",
                               "one_compound_folded_vs_extended"]}


def document():
    ids = [c["id"] for c in CASES]
    assert len(ids) == len(set(ids))
    for c in CASES:
        assert set(c["outcomes"].values()) <= set(GRADES), c["id"]
        for q in c["requires"]:
            assert q["grade"] in GRADES and sum(k in q for k in ("must", "must_not", "check")) == 1, c["id"]
            if "check" in q:
                assert re.fullmatch(r"shared_scale|relation=(conformers_of_one_compound|same_geometry|two_compounds)"
                                    r"|folded=(true|false|unstated),(true|false|unstated)", q["check"]), q
    assert all(r in ids for rs in RETIRES.values() for r in rs)
    return json.dumps({"generator": "tests/pedagogy/make_pair_view_set.py", "grades": list(GRADES),
                       "retires": RETIRES, "cases": CASES}, indent=1, ensure_ascii=False, allow_nan=True) + "\n"


if __name__ == "__main__":
    OUT.write_text(document(), encoding="utf-8")
    print(f"{len(CASES)} cases -> {OUT}")
