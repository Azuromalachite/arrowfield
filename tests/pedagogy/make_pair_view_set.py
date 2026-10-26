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

Two components, both raising model.Refused when they cannot judge or show
the pair faithfully (component "relation" and "pair_view" in each case):

    model.relation(a, b) -> model.Relation
    views.pair_view_html(left, right, context) -> str

model.relation owns the judgement of what the pair is, as model.py owns
`compare`; the view only displays it. Relation values: "same_geometry",
"conformers_of_one_compound", "two_compounds", decided from the structures,
never from names:
- one compound: the same constitution (the record's bonds, or its mapped
  SMILES), the same formal charge, and the same stereo assigned from 3D. If
  the graph or the 3D stereo cannot be established, refuse: never default
  to "two_compounds";
- the same geometry: one compound, and all-atom RMSD after optimal
  superposition, minimised over symmetry-equivalent atoms, at most
  SAME_GEOMETRY_RMSD. Measured before freezing (scratch probe, 2026-10-24):
  the pilot's 12 equivalent seed pairs (identical, rotated, methyl H
  relabelled) are at most 3.4e-6 A apart; the 760 pairs of distinct
  conformers in the Atlas's four ensembles are at least 0.98 A on heavy
  atoms, and the closest five at least 1.40 A on all atoms. 0.01 A is
  ~2900x the first and ~140x below the second; R1 pins it from both sides.

pair_view_html: left, right are MoleculeViews, one geometry each, with
surfaces. context: optional {"labels": {"left", "right"}} (display text
only, never evidence) and "study_rows" (verdicts the page also shows; never
linked). Pair view v1 refuses two different formal charges (P7): mixed-
charge comparison is a later, deliberate decision.

Amended before the view was built (the user, 2026-10-24, carrying #13's
lessons into the frozen contract):
- verbatim: text the view did not write is shown verbatim inside
  class="verbatim", one string per span, and every span equals one of the
  two molecules' names or record sources (MoleculeView.name,
  provenance.source) or one of the caller's display labels, exactly. Word
  bans read only the text outside those spans (P8);
- word families: the bans cover families (explain*, caus*, attract*,
  bind*, pheromon*, receptor*, reactiv*, resol(v|ut)*, specific*), with
  near-misses on the other side tested;
- no roles in v1: a context carrying "roles" is refused; whether a pair
  holds a teaching control is the page builder's to reconcile (the
  user's decision on cross-view roles). A context entry that is present
  but malformed (labels not a mapping of text, study_rows not a list) is
  refused, never read as absent (P9).

The rendered page carries, for the harness and for any reader's tools:
- a root element with data-relation, model.relation's value;
- exactly one element with class "pair-scale": data-kind, data-vmin,
  data-vmax, data-cmap, data-unit ("kcal/mol");
- exactly two elements with class "pair-panel", each with the effective
  scale its colours are computed on (data-vmin, data-vmax, data-cmap),
  data-scale="shared", and data-folded: "true", "false" or "unstated", the
  producer's flag as given (MoleculeView.folded, never inferred from genus).

The shared scale, exactly (lo, hi from the two records' own scales):

    lo = min(left.vmin, right.vmin); hi = max(left.vmax, right.vmax)
    if either record's scale is diverging, or lo < 0 < hi:
        diverging, [-M, +M] with M = max(|lo|, |hi|), cmap CMAP["diverging"]
    else:
        sequential, [lo, hi], cmap CMAP["sequential"]

CMAP is espmap's DEFAULT_CMAP (vik: warm = negative; lajolla: V_min
darkest). The caption states the colour-to-value direction of the scale
actually drawn: on a diverging one, which colour is negative; on a
sequential one, that darker is lower (more negative, or less positive).

Inputs are real records: the committed records/ (espmap v0.1.0, gas phase
and one PCM anion), and four Atlas conformers in tests/data/pair/ (espmap
v0.2.1; folded as the Atlas index states, members.json). A side may carry a
`patch` ((field, value), applied with dataclasses.replace; "scale.vmin" and
"provenance.x" reach one level down), a `surface` op (none, nan, inf,
empty) and a `geometry` op ("zeros", or ["rmsd", target A, seed]: seeded
displacements of exactly that RMSD before superposition) for synthetic
variants.

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
SAME_GEOMETRY_RMSD = 0.01                      # A; measured, see the docstring
CMAP = {"diverging": "vik", "sequential": "lajolla"}
CASES = []


def key(stem, which):
    return json.loads((HERE.parent / f"{stem}.json").read_text(encoding="utf-8"))[which]

# records
ETOH, ETOH2, ACAL, ACAL2 = "records/ethanol__s1", "records/ethanol__s2", "records/acetaldehyde__s1", "records/acetaldehyde__s2"
E2B, E2B2, Z2B = "records/e2_butene__s1", "records/e2_butene__s2", "records/z2_butene__s1"
ACO = "records/acetate__s12648430"
E3F, E3X = "tests/data/pair/e3_z13_18_oac/conf_00", "tests/data/pair/e3_z13_18_oac/conf_02"   # folded, extended
Z3X, Z3F = "tests/data/pair/z3_z13_18_oac/conf_00", "tests/data/pair/z3_z13_18_oac/conf_01"   # extended, folded


def side(stem, patch=(), surface=None, geometry=None):
    folded = FOLDED.get(stem.removeprefix("tests/data/pair/"))
    p = [list(x) for x in patch]
    if folded is not None and not any(x[0] == "folded" for x in p):
        p = [["folded", folded]] + p          # the Atlas index's flag, as study.load_ensemble sets it
    return {"record": stem, "patch": p, "surface": surface, "geometry": geometry}


def case(id, cls, why, left, right, outcomes, prefer, requires=(), context=None, same_object=False,
         component="pair_view"):
    outcomes = dict(outcomes)
    outcomes.setdefault("ERROR", "CRASH")
    allowed = sorted(k for k, g in outcomes.items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert prefer is None or prefer in allowed, id
    CASES.append(dict(id=id, cls=cls, component=component, why=why, left=left, right=right,
                      same_object=same_object, context=context or {}, outcomes=outcomes, allowed=allowed,
                      prefer=prefer, requires=list(requires)))


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
RANGE = must(r"[−+-]?\s?\d+(\.\d+)?\s*(to|…|–|\.\.)\s*[−+-]?\s?\d+(\.\d+)?", "MISLEADING",
             "the shared range, in numbers")
CONVENTION = must(r"(red|warm)[^.]{0,60}negative|negative[^.]{0,60}(red|warm)", "MISLEADING",
                  "on a diverging scale, which colour is negative: the chemistry convention, said "
                  "(CLAUDE.md conventions)")
SEQ_DIRECTION = must(r"darker[^.]{0,60}(more negative|lower|less positive)", "MISLEADING",
                     "on a sequential scale (lajolla, V_min darkest): darker is lower")
# equal folded flags remove one confound; they do not make two geometries equivalent
NO_SHAPE_MATCH = must_not(r"conformation(ally)?[- ](controlled|matched)|shape[- ](matched|controlled)|"
                          r"same shape|matched (shape|conformation)", "MISLEADING",
                          "equal folded flags are not matched shapes")
NOT_ESTABLISHED = must(r"not established", "MISLEADING",
                       "with a folded state unstated, shape comparability is not established, and the page says so")
METHOD = must(r"b3lyp", "MISLEADING", "the method both pictures were computed with")
# a picture pair is not evidence about biology, and says nothing about why
# word families, on the view's own words only (the amendment)
NO_CAUSE = must_not(r"because|explain\w*|due to|caus\w*|reactiv\w*|attract\w*|\bbind\w*|pheromon\w*|"
                    r"receptor\w*|\bresol(v|ut)\w*|specific\w*", "FALSE",
                    "two pictures show where the potential differs; they test no mechanism")
CAPTION = [SHARED, UNITS, RANGE, CONVENTION, METHOD, NO_CAUSE]

# --- P1 one scale for both pictures ------------------------------------------------------
case("butene_e_vs_z", "P1 one scale",
     "E- and Z-2-butene, same method, auto-scaled to ±17.3 and ±18.1 kcal/mol (different render keys). "
     "Each on its own scale, the same colour means different kcal/mol and the eye compares the wrong "
     "things: both go on one scale.",
     side(E2B), side(Z2B), SHOW, "RENDERED",
     CAPTION + [check("relation=two_compounds", "FALSE", "two compounds"), NOT_ESTABLISHED])
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
     "other side synthetic, −150 to −60): one sequential scale, −183.4 to −60, never centred on zero, and "
     "the caption says darker is more negative.",
     side(ACO), side(ACO, patch=[("scale.vmin", -150.0), ("scale.vmax", -60.0)]), SHOW, "RENDERED",
     [SHARED, UNITS, RANGE, SEQ_DIRECTION])
case("positive_only_sequential", "P1 one scale",
     "Synthetic (UI contract, not an observation): two all-positive sequential scales, +10 to +40 and +20 to "
     "+60. One sequential scale, +10 to +60; darker is lower (less positive), never 'more negative'.",
     side(ETOH, patch=[("scale.kind", "sequential"), ("scale.cmap", "lajolla"), ("scale.vmin", 10.0),
                       ("scale.vmax", 40.0)]),
     side(ACAL, patch=[("scale.kind", "sequential"), ("scale.cmap", "lajolla"), ("scale.vmin", 20.0),
                       ("scale.vmax", 60.0)]), SHOW, "RENDERED",
     [SHARED, UNITS, RANGE, SEQ_DIRECTION,
      must_not(r"darker[^.]{0,60}more negative", "FALSE", "nothing on this scale is negative")])
case("opposite_sign_sequential", "P1 one scale",
     "Synthetic (UI contract): one all-negative sequential scale (acetate, −183.4 to −70.4) and one "
     "all-positive (+20 to +80). Neither side changes sign, but the pair does (lo < 0 < hi): diverging, "
     "±183.4, never a sequential scale that runs through zero.",
     side(ACO), side(ACO, patch=[("scale.vmin", 20.0), ("scale.vmax", 80.0)]), SHOW, "RENDERED",
     [SHARED, UNITS, CONVENTION])

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
# The user's decision after #21's round 11: a record's protocol is its canonical key_fields.analysis,
# digested by its key; raw settings.xc is provenance text, which Arrowfield cannot check against
# the canonical functional without PySCF. One verified protocol with a raw label changed renders,
# captioned from the protocol, never from the label. (A raw solvent contradicting the protocol is
# checkable, and stays refused: stale_key_solvent.)
case("stale_key_functional", "P2 one method",
     "The keys match and digest one canonical protocol (B3LYP), but one side's raw settings name PBE0: a "
     "corrupted provenance label. The protocol decides; the caption names it, never the label.",
     side(E2B), side(Z2B, patch=[("provenance.xc", "pbe0")]), SHOW, "RENDERED",
     [METHOD, must_not(r"pbe0", "MISLEADING", "a raw label is not the protocol")])
# #21's round 12 (b6f7f66): one side carries its verified canonical protocol, the other none (a
# MoleculeView built in memory; record.load refuses a file without it). Compared on raw strings,
# both stale "pbe0" labels agreed and the caption named pbe0 for a record whose protocol is B3LYP
case("one_side_without_protocol", "P2 one method",
     "One side's protocol is verified B3LYP with a stale raw 'pbe0'; the other carries no canonical fields and "
     "says 'pbe0': the caption cannot name one method for both. Refuse.",
     side(E2B, patch=[("provenance.xc", "pbe0")]),
     side(Z2B, patch=[("provenance.analysis_fields", None), ("provenance.xc", "pbe0")]), REFUSE_F, "REFUSED")
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
                must_not(r"conformers of (one|the same) compound|same compound", "FALSE", "not one compound"),
                NO_SHAPE_MATCH])
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
case("identity_not_established", "P3 relation from structure",
     "A record with no bonds and no mapped SMILES: its constitution cannot be established, so neither can "
     "the relation. Refuse; never default to 'two compounds'.",
     side(ETOH, patch=[("bonds", None), ("mapped_smiles", None)]), side(ACAL), REFUSE_M, "REFUSED")
case("stereo_not_established", "P3 relation from structure",
     "The E3 acetate with every atom at the origin beside the Z3: stereo cannot be assigned from 3D, so "
     "'two compounds' would be a guess that happens to be right. Refuse.",
     side(E3X, geometry="zeros"), side(Z3X), REFUSE_M, "REFUSED")

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
      must(r"not stated", "MISLEADING", "said in words"), must_not(r"\bextended\b", "FALSE", "not invented"),
      NOT_ESTABLISHED])
case("both_extended_not_shape_matched", "P4 folded",
     "E3 and Z3 both extended: equal flags remove one confound, they do not match the shapes. Render, "
     "never 'conformation-controlled' or 'shape-matched'.",
     side(E3X), side(Z3X), SHOW, "RENDERED", [check("folded=false,false", "FALSE", "as given"), NO_SHAPE_MATCH])
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

# --- P8 verbatim text and the view's own words (the amendment) ------------------------------------
case("label_with_banned_words", "P8 verbatim and own words",
     "Display labels (synthetic) reading 'attractive E isomer' and 'explains the Z response': the caller's "
     "text, shown verbatim in spans; the view's own words stay clean.",
     side(E2B), side(Z2B), SHOW, "RENDERED",
     [must(r"attractive E isomer", "MISLEADING", "the label, unaltered"),
      must(r"explains the Z response", "MISLEADING", "the label, unaltered"), NO_CAUSE, SHARED],
     context={"labels": {"left": "attractive E isomer", "right": "explains the Z response"}})
case("name_with_banned_word", "P8 verbatim and own words",
     "A molecule name (synthetic) 'inhibitor_receptor_conf': the record's name, verbatim in a span, never "
     "the view's words.",
     side(E3F, patch=[("name", "inhibitor_receptor_conf")]), side(E3X), SHOW, "RENDERED",
     [must(r"inhibitor_receptor_conf", "MISLEADING", "the name, unaltered"), NO_CAUSE])

# --- P9 the context: no roles in v1; a present but malformed entry is refused ----------------
for vid, ctx_, why in (
        ("roles_given", {"roles": {"z3_z13_18_al": "teaching_control"}},
         "roles on a pair: v1 takes none (reconciling a control across views is the page builder's)"),
        ("roles_empty", {"roles": {}}, "an empty roles mapping: still roles, which v1 takes none of"),
        ("labels_not_mapping", {"labels": "E vs Z"}, "labels given as a string: not read as absent"),
        ("label_not_text", {"labels": {"left": 3, "right": "Z"}}, "a label that is not text"),
        ("study_rows_not_list", {"study_rows": "resolved"}, "study_rows given as a string")):
    case(f"context_{vid}", "P9 context", f"{why[0].upper()}{why[1:]}: refuse.",
         side(E2B), side(Z2B), REFUSE_M, "REFUSED", context=ctx_)

# --- P7 formal charge: pair view v1 compares one charge only ----------------------------------
case("charge_mismatch", "P7 formal charge",
     "Synthetic (UI contract): ethanol beside acetaldehyde stated as charge −1, one method. Mixed-charge "
     "pictures are out of scope for pair view v1: refuse, and say so.",
     side(ETOH), side(ACAL, patch=[("charge", -1)]), REFUSE_M, "REFUSED")
case("charge_mismatch_same_geometry", "P7 formal charge",
     "Synthetic: one E-2-butene geometry stated once as neutral and once as +1. Not one species: refuse.",
     side(E2B), side(E2B2, patch=[("charge", 1)]), REFUSE_M, "REFUSED")

# --- R1, R2 model.relation: what the pair is, judged in the model -------------------------------
REL = {"SAME_GEOMETRY": "MISLEADING", "CONFORMERS_OF_ONE_COMPOUND": "MISLEADING", "TWO_COMPOUNDS": "FALSE",
       "REFUSED": "SAFE_REFUSAL"}


def rel(id_, cls, why, left, right, answer, wrong=None):
    outcomes = dict(REL, **(wrong or {}), **{answer: "SAFE"})
    case(id_, cls, why, left, right, outcomes, answer, component="relation")


ONE = {"TWO_COMPOUNDS": "FALSE"}               # one compound called two
TWO = {"SAME_GEOMETRY": "FALSE", "CONFORMERS_OF_ONE_COMPOUND": "FALSE"}
rel("rel_identical", "R1 relation", "E-2-butene seeds, identical coordinates.", side(E2B), side(E2B2),
    "SAME_GEOMETRY", ONE)
rel("rel_rotated", "R1 relation", "Acetaldehyde seeds: one geometry, rotated (RMSD 3.4e-6 A at most).",
    side(ACAL), side(ACAL2), "SAME_GEOMETRY", ONE)
rel("rel_h_relabelled", "R1 relation", "Ethanol seeds: one geometry, methyl H labels permuted.",
    side(ETOH), side(ETOH2), "SAME_GEOMETRY", ONE)
rel("rel_below_threshold", "R1 relation",
    f"The E3 acetate against itself displaced by RMSD 0.005 A (half of SAME_GEOMETRY_RMSD, "
    f"{SAME_GEOMETRY_RMSD} A): the same geometry.",
    side(E3X), side(E3X, geometry=["rmsd", 0.005, 1]), "SAME_GEOMETRY", ONE)
rel("rel_above_threshold", "R1 relation",
    "The E3 acetate displaced by RMSD 0.04 A (4x the threshold): no longer the same geometry, one compound.",
    side(E3X), side(E3X, geometry=["rmsd", 0.04, 2]), "CONFORMERS_OF_ONE_COMPOUND", ONE)
rel("rel_true_conformers", "R1 relation", "E3 conf_00 and conf_02: symmetric RMSD 2.37 A.",
    side(E3F), side(E3X), "CONFORMERS_OF_ONE_COMPOUND", ONE)
rel("rel_stereoisomers", "R1 relation", "E3 and Z3 acetate: one double bond's stereo differs.",
    side(E3X), side(Z3X), "TWO_COMPOUNDS", TWO)
rel("rel_constitution", "R1 relation", "Ethanol and acetaldehyde.", side(ETOH), side(ACAL), "TWO_COMPOUNDS", TWO)
rel("rel_no_graph", "R2 identity not established",
    "No bonds and no mapped SMILES on one side: constitution unknown. Refuse.",
    side(ETOH, patch=[("bonds", None), ("mapped_smiles", None)]), side(ETOH2), "REFUSED",
    {"SAME_GEOMETRY": "FALSE"})
rel("rel_flat_one_side", "R2 identity not established",
    "E3 with every atom at the origin, against Z3: no stereo from 3D. 'Two compounds' would be a lucky guess.",
    side(E3X, geometry="zeros"), side(Z3X), "REFUSED", {"SAME_GEOMETRY": "FALSE"})
rel("rel_flat_both_sides", "R2 identity not established",
    "E3 and E3 conf_00, both collapsed to the origin: RMSD 0, but no geometry and no stereo. 'Same "
    "geometry' would be false.", side(E3X, geometry="zeros"), side(E3F, geometry="zeros"), "REFUSED",
    {"SAME_GEOMETRY": "FALSE"})

# probe on #18 (fa0e676): the relation depended on the bond list's order and missed a methyl-H
# relabelling (atom 0 is a CH3 with H22, H23, H24 in E3 conf_02). The minimisation over
# symmetry-equivalent atoms must not depend on how the input lists its atoms or bonds.
rel("rel_methyl_h_swapped", "R1 relation",
    "E3 conf_02 against itself with methyl H22 and H23 exchanging coordinates: one geometry, relabelled "
    "(naive RMSD 0.33 A, 0 over symmetry).", side(E3X), side(E3X, geometry=["swap", 22, 23]), "SAME_GEOMETRY", ONE)
rel("rel_bonds_shuffled_same", "R1 relation",
    "The same geometry with its bond list reversed and each bond's direction flipped: the answer cannot "
    "depend on the list's order.", side(E3X), side(E3X, geometry="bonds_shuffled"), "SAME_GEOMETRY", ONE)
rel("rel_swap_and_shuffle", "R1 relation",
    "Both at once: methyl H swapped on one side, bonds shuffled on the other.",
    side(E3X, geometry=["swap", 22, 23]), side(E3X, geometry="bonds_shuffled"), "SAME_GEOMETRY", ONE)
rel("rel_bonds_shuffled_conformers", "R1 relation",
    "A true conformer pair (conf_00, conf_02) with one side's bonds shuffled: still two conformers.",
    side(E3F), side(E3X, geometry="bonds_shuffled"), "CONFORMERS_OF_ONE_COMPOUND", ONE)
rel("rel_bond_out_of_range", "R2 identity not established",
    "A bond naming an atom the record does not have: the graph is not established. Refuse, with a plain "
    "reason (not a C++ invariant).", side(E3X, geometry="bond_out_of_range"), side(E3F), "REFUSED",
    {"SAME_GEOMETRY": "FALSE"})

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
        assert c["component"] in ("pair_view", "relation")
        if c["component"] == "relation":
            assert not c["requires"] and set(c["outcomes"]) == set(REL) | {"ERROR"}, c["id"]
        for q in c["requires"]:
            assert q["grade"] in GRADES and sum(k in q for k in ("must", "must_not", "check")) == 1, c["id"]
            if "check" in q:
                assert re.fullmatch(r"shared_scale|relation=(conformers_of_one_compound|same_geometry|two_compounds)"
                                    r"|folded=(true|false|unstated),(true|false|unstated)", q["check"]), q
    assert all(r in ids for rs in RETIRES.values() for r in rs)
    return json.dumps({"generator": "tests/pedagogy/make_pair_view_set.py", "grades": list(GRADES),
                       "same_geometry_rmsd": SAME_GEOMETRY_RMSD, "cmap": CMAP, "retires": RETIRES,
                       "cases": CASES}, indent=1, ensure_ascii=False, allow_nan=True) + "\n"


if __name__ == "__main__":
    OUT.write_text(document(), encoding="utf-8")
    print(f"{len(CASES)} cases -> {OUT}")
