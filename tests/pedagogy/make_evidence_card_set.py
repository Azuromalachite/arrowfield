"""Generate the hostile set for the evidence views:
tests/pedagogy/evidence_card_set.json.

    python tests/pedagogy/make_evidence_card_set.py     # byte-identical on re-run

Written before the views exist (expectations first). It re-expresses the
nine abstract evidence cases of hostile_set.json (B1-B5, which predate
bundle.py) against the real Atlas teaching bundle, as bundle.py reads it, and
extends them to the Atlas's evidence classes and population contexts. Inputs
are real records: tests/data/sesiidae-homework-1.bundle.json (Atlas 35112e1),
loaded with bundle.load(), with an optional `bundle_patch` (path, value) for
synthetic variants. The three views, all in views.py, all raising
model.Refused when they cannot show the evidence faithfully:

- evidence_card_html(statement, compound, context): one Atlas observation,
  as it bears on one compound;
- contrast_card_html(contrast, statements, context): one ObservedContrast,
  its two observations side by side;
- evidence_class_panel_html(compound, taxon_id, context): the Atlas's own
  evidence classes for one compound in one taxon (docs/evidence-classes.md
  in the Atlas).

context: {"bundle": {"atlas_version", "sha256"}, "roles": {compound id:
atlas_observed | teaching_control | computational_counterfactual},
optional "study_rows": study verdicts the page also shows}. "sha256":
"from bundle.load" stands for the loaded bundle's own hash, which the
harness fills in (the file's bytes differ with line endings).

The bundle object is the authority (E8, E9): a statement, a contrast or a
provenance claim that is not the bundle's own is refused, never shown as
Atlas evidence. `statement_patch` (path, value) builds the inputs the view
is handed from a patched copy of the bundle, while the view is given the
unpatched bundle: the same id with altered content. A path ending in "+"
appends.

Grades as in the other sets (SAFE, SAFE_REFUSAL, CRASH, MISLEADING, FALSE);
`allowed` is exactly the safe outputs; a rendered view must meet every
`requires` entry (must / must_not, a case-insensitive pattern on the text
with tags removed, or on the HTML). The Atlas's rules this set holds the
views to: no class is netted out; "not detected" is a negative for that
method, sample and population only; absence of a class row is not evidence
of absence; female_pheromone_component is the Atlas's join of two lines of
evidence, not an identification, and its population_match must be read.
"""
import copy
import json
import re
from pathlib import Path

OUT = Path(__file__).with_name("evidence_card_set.json")
BUNDLE = Path(__file__).resolve().parents[1] / "data" / "sesiidae-homework-1.bundle.json"
GRADES = ("SAFE", "SAFE_REFUSAL", "CRASH", "MISLEADING", "FALSE")
RAW = json.loads(BUNDLE.read_text(encoding="utf-8"))
P = "spa.obs."
OAC, OH, E3, E2 = ("spa.cmpd.z3_z13_18_oac", "spa.cmpd.z3_z13_18_oh", "spa.cmpd.e3_z13_18_oac",
                   "spa.cmpd.e2_z13_18_oac")
AL = "arrowfield.control.z3_z13_18_al"
TENUIS, TIPU = "spa.tax.synanthedon_tenuis", "spa.tax.synanthedon_tipuliformis"
ROLES = {OAC: "atlas_observed", OH: "atlas_observed", E3: "atlas_observed", E2: "atlas_observed",
         AL: "teaching_control"}
CASES = []


def ctx(**kw):
    c = {"bundle": {"atlas_version": RAW["atlas_version"], "sha256": "from bundle.load"}, "roles": dict(ROLES)}
    c.update(kw)
    return c


def case(id, cls, component, why, inputs, outcomes, prefer, requires=(), context=None, patch=(), spatch=()):
    outcomes = dict(outcomes)
    outcomes.setdefault("ERROR", "CRASH")
    allowed = sorted(k for k, g in outcomes.items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert prefer is None or prefer in allowed, id
    CASES.append(dict(id=id, cls=cls, component=component, why=why, inputs=inputs,
                      context=context if context is not None else ctx(), bundle_patch=[list(p) for p in patch],
                      statement_patch=[list(p) for p in spatch],
                      outcomes=outcomes, allowed=allowed, prefer=prefer, requires=list(requires)))


def must(pat, grade, why, on="text"):
    return {"must": pat, "grade": grade, "why": why, "on": on}


def must_not(pat, grade, why, on="text"):
    return {"must_not": pat, "grade": grade, "why": why, "on": on}


SHOW = {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}
REFUSE = {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"}
# never on an evidence view: interpretation the records do not carry
# Word bans read the view's own words only (the user, 2026-10-24, after probe round 1):
# text copied from a record (the Atlas statement, the authors' quotation, unresolved
# items, citation, labels, places, names) is shown verbatim inside class="verbatim",
# attributed, and is exempt; the harness checks every verbatim span is record text.
# The bans cover word families (inhibit*, repel*, antagonis*, inactiv*, resol(v|ut)*).
NO_OVERREACH = must_not(r"inactiv\w*|does not (make|produce|attract)|not (made|produced)|absent from|"
                        r"repel\w*|inhibit\w*|antagonis\w*", "FALSE",
                        "an evidence view reports what was measured, never a role the record does not state")
NO_NETTING = must_not(r"\boverall\b|in general|on balance|\bnet\b|taken together|consensus", "MISLEADING",
                      "the Atlas never nets opposite results out; neither does the view")
# \bresol(v|ut)\w*: the study's verdict word and its family (resolved, resolves,
# resolvable, resolution; probe round 1), never the Atlas's own 'unresolved'
# items, which E2 requires the card to show (no word boundary inside them)
NO_ESP = must_not(r"\bESP\b|electrostatic|kcal|\bresol(v|ut)\w*|because|explain|due to|caus", "FALSE",
                  "an evidence view never links a computed difference to a biological response")


def ev(id_, cls, oid, comp, why, requires, outcomes=SHOW, prefer="RENDERED", **kw):
    case(id_, cls, "evidence_card", why, {"observation_id": P + oid, "compound": comp}, outcomes, prefer,
         requires, **kw)


# --- E1 absence is not a negative, and a negative is only as wide as its method ---------
ev("not_detected_with_limit", "E1 absence and negatives", "synanthedon_tenuis.yang_2012.gland", OAC,
   "Z3,Z13-18:OAc not detected in Korean glands, limit 100 pg: a negative for that method, sample "
   "and population only.",
   [must(r"not detected", "FALSE", "the record's word"), must(r"100 ?(pg|picogram)", "MISLEADING", "the limit"),
    must(r"korean population", "MISLEADING", "whose glands"), NO_OVERREACH, NO_ESP])
ev("not_detected_limit_not_stated", "E1 absence and negatives", "synanthedon_tenuis.naka_2013.gland", OH,
   "Z3,Z13-18:OH not detected in Japanese glands, no limit given: the card says the limit is not stated.",
   [must(r"not detected", "FALSE", "the record's word"),
    must(r"limit (is )?not (stated|given|reported)", "MISLEADING", "an unknown limit stays visible"),
    must(r"japanese population", "MISLEADING", "whose glands"), NO_OVERREACH])
ev("compound_not_in_observation", "E1 absence and negatives", "synanthedon_tipuliformis.szocs_1990.trap_e3.tasmania_1987",
   OH, "Z3,Z13-18:OH is not part of this trapping test: not tested here, never 'no effect' or 'negative'.",
   [must(r"not (part of|tested in|included in) this", "FALSE", "absence of a test, said as such"),
    must_not(r"no effect|not detected|negative", "FALSE", "no result is invented")])
ev("no_effect_not_significant", "E1 absence and negatives",
   "synanthedon_tipuliformis.szocs_1990.trap_e3.tasmania_1987", E3,
   "Adding E3,Z13-18:OAc had no effect, not significant, at Slateford 1987: a tested negative for that "
   "lure and place, not 'inactive'.",
   [must(r"no effect", "FALSE", "the record's effect"), must(r"not significant", "MISLEADING", "its test"),
    must(r"slateford", "MISLEADING", "where"), NO_OVERREACH])
ev("no_effect_significance_not_reported", "E1 absence and negatives",
   "synanthedon_tipuliformis.suckling_2005.trap_e3.new_zealand", E3,
   "No effect with significance not reported: the card must not upgrade it to 'not significant'.",
   [must(r"significance not reported|not reported", "MISLEADING", "the record's significance"),
    must_not(r"(?<!significance )\bnot significant\b", "MISLEADING", "no test is invented")])
ev("effect_with_values", "E1 absence and negatives",
   "synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1", OH,
   "A significant decrease when Z3,Z13-18:OH was added (Japan, test IV): the card shows the effect, its "
   "significance and the raw catches as reported.",
   [must(r"decreased", "FALSE", "the effect"), must(r"\bsignificant\b", "MISLEADING", "its test"),
    must(r"4\.6", "MISLEADING", "treatment catch"), must(r"17\.0", "MISLEADING", "comparator catch"),
    must(r"japanese population", "MISLEADING", "which population"), NO_OVERREACH, NO_ESP])

# --- E2 the Atlas's words and the source's words, with provenance ---------------------------
gland = next(o for o in RAW["observations"] if o["observation_id"] == P + "synanthedon_tenuis.naka_2013.gland")
ev("statement_verbatim_and_citation", "E2 provenance", "synanthedon_tenuis.naka_2013.gland", OAC,
   "The Atlas's own statement is shown verbatim, the author's statement as a quotation attributed to the "
   "paper, and the citation with its DOI.",
   [must(re.escape(gland["statement"][:80]), "MISLEADING", "the Atlas's sentence, unaltered"),
    must(re.escape(gland["reference"]["doi"]), "MISLEADING", "the DOI"),
    must(r"[\"“]" + re.escape(gland["evidence"]["author_statement"][:40]), "MISLEADING",
         "the author's words, quoted")])
ev("unresolved_items_shown", "E2 provenance", "synanthedon_tenuis.naka_2013.gland", OH,
   "The record lists what it could not resolve (the detection limit for the negatives): the card shows it.",
   [must(r"unresolved", "MISLEADING", "labelled"), must(r"detection limit for the negatives", "MISLEADING",
                                                        "the item itself")])
# Re-expressed in the trust-boundary commit: the bundle, not the caller's
# context, names the commit (E9), so "missing" is now a bundle whose
# atlas_version is not a commit. The context-absent case is E9's.
ev("atlas_version_missing", "E2 provenance", "synanthedon_tenuis.naka_2013.gland", OAC,
   "Without the Atlas commit the card cannot be traced: a visible warning, or refusal.",
   [must(r"provenance (incomplete|missing)|atlas (version|commit) (not recorded|missing|unknown)",
         "MISLEADING", "a visible warning")],
   context=ctx(bundle={"atlas_version": None, "sha256": None}), patch=[(["atlas_version"], "unknown")])

# --- E3 roles: constructed molecules never become observed --------------------------------------
case("teaching_control_classes", "E3 roles", "class_panel",
     "The aldehyde is a teaching control, not an Atlas compound: its panel says so and claims no evidence.",
     {"compound": AL, "taxon_id": TENUIS}, SHOW, "RENDERED",
     [must(r"teaching control", "FALSE", "labelled"), must(r"no atlas evidence|not an atlas compound", "FALSE",
                                                         "and why it has none"),
      must_not(r"observed|detected|pheromone|attract", "FALSE", "no evidence invented")])
case("counterfactual_classes", "E3 roles", "class_panel",
     "A computational counterfactual is a model, never observed.",
     {"compound": "arrowfield.counterfactual.example", "taxon_id": TENUIS}, SHOW, "RENDERED",
     [must(r"computational (model|counterfactual)", "FALSE", "labelled"),
      must_not(r"observed|detected|pheromone|attract", "FALSE", "no evidence invented")],
     context=ctx(roles=dict(ROLES, **{"arrowfield.counterfactual.example": "computational_counterfactual"})))
ev("teaching_control_in_observation", "E3 roles", "synanthedon_tenuis.naka_2013.gland", AL,
   "Asking an observation about a compound that is not in the bundle (the teaching control) cannot be "
   "answered from the record: refuse.", [], outcomes=REFUSE, prefer="REFUSED")


# --- E4 contexts kept apart: populations, places, years, assays ------------------------------
def ct(id_, label, why, requires):
    case(id_, "E4 contexts kept apart", "contrast_card", why, {"contrast_label": label}, SHOW, "RENDERED",
         list(requires) + [NO_NETTING, NO_ESP])


ct("glands_two_populations", "S. tenuis gland, Japan and Korea",
   "Japanese glands: Z3,Z13-18:OAc, no alcohol; Korean glands: the alcohol, no acetate. Two populations, "
   "shown as two statements, never one 'S. tenuis pheromone'.",
   [must(r"japanese population", "FALSE", "first population"), must(r"korean population", "FALSE",
                                                                    "second population")])
ct("opposite_effects_two_places", "S. tipuliformis, E3,Z13-18:OAc added, Hungary and Tasmania 1986",
   "E3,Z13-18:OAc raised catches in Hungary and lowered them in Tasmania: both shown with their places, "
   "never averaged or called a contradiction to resolve.",
   [must(r"hungary", "FALSE", "first place"), must(r"tasmania", "FALSE", "second place"),
    must(r"increased", "FALSE", "first effect"), must(r"decreased", "FALSE", "second effect"),
    must_not(r"contradict|inconsisten|one of (these|them) is wrong", "MISLEADING",
             "different places are not an error to resolve")])
ct("same_region_two_studies", "S. tipuliformis, E3,Z13-18:OAc added in Tasmania, 1986 and 2005",
   "Two studies in Tasmania, 1986 and 2005: both citations visible, so a reader sees two experiments.",
   [must(r"sz[öo]cs", "MISLEADING", "first study"), must(r"suckling", "MISLEADING", "second study"),
    must(r"1990|1986", "MISLEADING", "first study's date"), must(r"2005", "MISLEADING", "second study's date")])
ct("gland_vs_lure", "S. tipuliformis gland content and the six-compound lure",
   "A gland analysis and a field lure test of the same compounds are different evidence: both assays named.",
   [must(r"gland", "FALSE", "the chemical measurement"), must(r"trap|lure|field", "FALSE",
                                                              "the behavioural test")])

# --- E5 the explanation boundary ----------------------------------------------------------------
ev("study_verdict_on_the_same_page", "E5 explanation boundary",
   "synanthedon_tipuliformis.szocs_1990.trap_e3.hungary", E3,
   "The page also shows a resolved ESP verdict for E3 vs E2 (C4-H). The evidence card must not join them: "
   "no 'because', no ESP, no mechanism.",
   [NO_ESP, NO_OVERREACH],
   context=ctx(study_rows=[{"verdict_rule": "espmap-sesiidae-001/range_partition_sign_v1",
                            "contrast": "C_position", "site": "C4-H", "outcome": "resolved"}]))
case("study_verdict_beside_classes", "E5 explanation boundary", "class_panel",
     "The same on the class panel: E3,Z13-18:OAc's classes in S. tipuliformis, beside a resolved ESP verdict, "
     "with no 'because' and no mechanism.",
     {"compound": E3, "taxon_id": TIPU}, SHOW, "RENDERED", [NO_ESP, NO_OVERREACH],
     context=ctx(study_rows=[{"verdict_rule": "espmap-sesiidae-001/range_partition_sign_v1",
                              "contrast": "C_position", "site": "C4-H", "outcome": "resolved"}]))


# --- E6 the Atlas's evidence classes ---------------------------------------------------------
def cp(id_, comp, tax, why, requires, patch=()):
    case(id_, "E6 evidence classes", "class_panel", why, {"compound": comp, "taxon_id": tax}, SHOW, "RENDERED",
         list(requires) + [NO_NETTING], patch=patch)


CLS = r"female[_ ]pheromone[_ ]component"
cp("population_specific_roles", OAC, TENUIS,
   "Z3,Z13-18:OAc in S. tenuis is a female pheromone component in the Japanese population and not "
   "detected (and inhibitory when added) in the Korean one. Every class row stays, with its population.",
   [must(CLS, "FALSE", "the derived class"), must(r"female[_ ]not[_ ]detected", "FALSE", "the negative row"),
    must(r"addition[_ ]decreased", "FALSE", "the inhibition row"),
    must(r"japanese population", "FALSE", "first population"), must(r"korean population", "FALSE",
                                                                    "second population")])
cp("increase_and_decrease_both_kept", E3, TIPU,
   "E3,Z13-18:OAc in S. tipuliformis: addition_increased (Hungary, Lithuania, NZ, Germany) and "
   "addition_decreased (Tasmania, Bavaria). Both rows, never netted.",
   [must(r"addition[_ ]increased", "FALSE", "one row"), must(r"addition[_ ]decreased", "FALSE", "the other"),
    must(r"tasmania", "MISLEADING", "where it decreased"), must(r"hungary", "MISLEADING",
                                                               "where it increased")])
cp("derived_class_is_a_join", OAC, TENUIS,
   "female_pheromone_component is the Atlas's join of detection and attraction, not an identification; "
   "the panel says so and shows population_match.",
   [must(r"join|combines|derived", "MISLEADING", "what the class is"),
    must(r"not an identification|not (a|an) (identification|claim)", "MISLEADING", "what it is not"),
    must(r"population match:? ?yes", "MISLEADING", "the match")])
cp("outside_this_case", OAC, TENUIS,
   "Most rows rest on observations outside this homework's bundle (7 of 8 for the derived class): the "
   "panel says how many, so a reader does not take the bundle for the whole evidence.",
   [must(r"outside (this|the) (case|homework|bundle)", "MISLEADING", "said"),
    must(r"\b7\b", "MISLEADING", "how many, for the derived class")])
cp("no_class_is_not_absence", OH, TIPU,
   "Z3,Z13-18:OH has no class row for S. tipuliformis. Absence of a row is not evidence of absence.",
   [must(r"no (atlas )?evidence class", "FALSE", "said plainly"),
    must(r"not evidence of absence|does not mean", "MISLEADING", "and what that means"), NO_OVERREACH])
idx = next(i for i, c in enumerate(RAW["compounds"]) if c["atlas_id"] == OAC)
jdx = next(j for j, e in enumerate(RAW["compounds"][idx]["evidence_classes"])
           if e["evidence_class"] == "female_pheromone_component" and e["taxon_id"] == TENUIS)
for vid, val, req in (("no", "no", r"population match:? ?no"), ("unknown", "unknown", r"population match:? ?unknown")):
    cp(f"population_match_{vid}", OAC, TENUIS,
       f"population_match '{val}' (synthetic): the join may cross populations; the panel shows the value and "
       "sends the reader to the records.",
       [must(req, "FALSE", "the value"), must(r"read the records|check the (records|populations)", "MISLEADING",
                                              "what to do about it")],
       patch=[(["compounds", idx, "evidence_classes", jdx, "population_match"], val)])

# --- E7 what cannot be shown --------------------------------------------------------------
def oi(oid):
    return next(i for i, o in enumerate(RAW["observations"]) if o["observation_id"] == P + oid)


def ci(oid, comp):
    return next(j for j, x in enumerate(RAW["observations"][oi(oid)]["components"]) if x["atlas_id"] == comp)


# Since the trust-boundary commit the view itself is handed these (a
# fabricated statement or contrast, built by statement_patch), so the
# refusal is the view's, not the harness's.
fake_obs = dict(copy.deepcopy(RAW["observations"][oi("synanthedon_tenuis.naka_2013.gland")]),
                observation_id=P + "synanthedon_tenuis.not_in_this_bundle")
ev("observation_not_in_bundle", "E7 refusals", "synanthedon_tenuis.not_in_this_bundle", OAC,
   "An observation id the bundle does not carry: refuse, never guess.", [], outcomes=REFUSE, prefer="REFUSED",
   spatch=[(["observations", "+"], fake_obs)])
case("contrast_not_in_bundle", "E7 refusals", "contrast_card",
     "A contrast label the bundle does not carry: refuse.", {"contrast_label": "no such contrast"},
     REFUSE, "REFUSED",
     spatch=[(["comparisons", "+"], dict(copy.deepcopy(RAW["comparisons"][0]), case_label="no such contrast"))])
case("class_panel_unknown_compound", "E7 refusals", "class_panel",
     "A compound neither in the bundle nor given a role: refuse.",
     {"compound": "spa.cmpd.not_here", "taxon_id": TENUIS}, REFUSE, "REFUSED")

# --- trust boundary (review of #13, at 64082a4): cases added before the views change ---------
# E1: the wording of a missing effect, and the method a negative is bounded by
def method_pat(oid, comp):
    m = RAW["observations"][oi(oid)]["components"][ci(oid, comp)]["method"]
    return re.escape(m).replace("_", "[_ ]")


ev("negative_names_its_method", "E1 absence and negatives", "synanthedon_tenuis.yang_2012.gland", OAC,
   "A negative is bounded by its method: the card names the record's method, not just 'this method'.",
   [must(method_pat("synanthedon_tenuis.yang_2012.gland", OAC), "MISLEADING", "the method itself")])
ev("limit_not_stated_names_its_method", "E1 absence and negatives", "synanthedon_tenuis.naka_2013.gland", OH,
   "The same when the limit is not stated: the method is still named.",
   [must(method_pat("synanthedon_tenuis.naka_2013.gland", OH), "MISLEADING", "the method itself")])
ev("effect_not_recorded", "E1 absence and negatives", "synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1", OH,
   "A record with no effect field (synthetic) does not state an effect: 'no effect recorded' reads as a "
   "tested zero effect.",
   [must(r"effect not recorded|does not state an effect", "MISLEADING", "said as missing"),
    must_not(r"no effect", "MISLEADING", "never a tested zero"), NO_OVERREACH],
   patch=[(["observations", oi("synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1"), "outcome", "effect"],
           None)])

# E3: bundle membership wins over a caller's role
for vid, comp, tax, role in (("control", OAC, TENUIS, "teaching_control"),
                             ("counterfactual", E3, TIPU, "computational_counterfactual")):
    case(f"atlas_compound_labelled_{vid}", "E3 roles", "class_panel",
         f"The caller's roles call an Atlas compound a {role}: the bundle carries it, so the roles are wrong. "
         "Saying 'not an Atlas compound' would erase its evidence: refuse the inconsistent input.",
         {"compound": comp, "taxon_id": tax}, REFUSE, "REFUSED", context=ctx(roles=dict(ROLES, **{comp: role})))

# E8: the statement and contrast shown are the bundle's own
GLAND_N, GLAND_Y = "synanthedon_tenuis.naka_2013.gland", "synanthedon_tenuis.yang_2012.gland"
ev("altered_statement_same_id", "E8 bound to the bundle", GLAND_N, OAC,
   "A statement carrying a real observation id with altered text is not the Atlas's: shown, it puts words "
   "in the Atlas's mouth.", [], outcomes=REFUSE, prefer="REFUSED",
   spatch=[(["observations", oi(GLAND_N), "statement"], "Z3,Z13-18:OAc is the sex pheromone of Synanthedon tenuis.")])
ev("altered_negative_same_id", "E8 bound to the bundle", GLAND_Y, OAC,
   "The same id with the detection limit changed from 100 pg to 1 pg: a stronger negative than the record's.",
   [], outcomes=REFUSE, prefer="REFUSED",
   spatch=[(["observations", oi(GLAND_Y), "components", ci(GLAND_Y, OAC), "detection_limit"],
            {"unit": "picogram", "value": 1})])
case("altered_contrast_statement", "E8 bound to the bundle", "contrast_card",
     "The right contrast and ids, but the Korean gland statement altered: refuse.",
     {"contrast_label": "S. tenuis gland, Japan and Korea"}, REFUSE, "REFUSED",
     spatch=[(["observations", oi(GLAND_Y), "statement"], "Korean females make the same blend as Japanese ones.")])
case("altered_contrast_text", "E8 bound to the bundle", "contrast_card",
     "The contrast's own observed text altered under its real label and ids: not the bundle's contrast.",
     {"contrast_label": "S. tenuis gland, Japan and Korea"}, REFUSE, "REFUSED",
     spatch=[(["comparisons", 0, "observed_contrast", "a"], "Both populations use Z3,Z13-18:OAc.")])

# E9: provenance is the bundle's; a caller's claim that disagrees is refused
FULL = "35112e1bc6ecd89525b74cc5c39a77cdbf87018f"
assert FULL.startswith(RAW["atlas_version"]), "update FULL with the fixture's commit"
ev("context_version_disagrees", "E9 provenance from the bundle", GLAND_N, OAC,
   "The page's context says Atlas 30c4a11 (the fixture's old, false label) for a bundle built at another "
   "commit: one of them is wrong, so show neither.",
   [], outcomes=REFUSE, prefer="REFUSED", context=ctx(bundle={"atlas_version": "30c4a11", "sha256": "from bundle.load"}))
case("context_version_disagrees_contrast", "E9 provenance from the bundle", "contrast_card",
     "The same disagreement on the contrast card.", {"contrast_label": "S. tenuis gland, Japan and Korea"},
     REFUSE, "REFUSED", context=ctx(bundle={"atlas_version": "30c4a11", "sha256": "from bundle.load"}))
case("context_sha_disagrees", "E9 provenance from the bundle", "class_panel",
     "The context names a bundle file hash that is not this bundle's: a stale page.",
     {"compound": OAC, "taxon_id": TENUIS}, REFUSE, "REFUSED",
     context=ctx(bundle={"atlas_version": RAW["atlas_version"], "sha256": "0" * 64}))
ev("context_provenance_absent", "E9 provenance from the bundle", GLAND_N, OAC,
   "No provenance in the context: the bundle names its commit, and the card shows that one.",
   [must(r"atlas commit " + re.escape(RAW["atlas_version"]), "MISLEADING", "the bundle's commit")],
   context={"roles": dict(ROLES)})
case("full_commit_in_context_agrees", "E9 provenance from the bundle", "class_panel",
     "The context gives the full hash of the bundle's commit: that agrees, so the panel renders.",
     {"compound": OAC, "taxon_id": TENUIS}, SHOW, "RENDERED",
     [must(r"atlas commit " + re.escape(RAW["atlas_version"]), "MISLEADING", "the bundle's commit")],
     context=ctx(bundle={"atlas_version": FULL, "sha256": "from bundle.load"}))

# --- held-out probe round 1 (Publishing, on 64082a4): its findings, as cases ----------------
# (the author-quote findings wait on the user's decision between verbatim quotes and the word ban)
T86 = "synanthedon_tipuliformis.szocs_1990.trap_e3.tasmania_1986"
ev("altered_effect_same_id", "E8 bound to the bundle", T86, E3,
   "Tasmania 1986 with its effect and statement turned from decreased to increased, same id: the card "
   "would show the opposite of the record.", [], outcomes=REFUSE, prefer="REFUSED",
   spatch=[(["observations", oi(T86), "outcome", "effect"], "increased"),
           (["observations", oi(T86), "statement"], "Adding E3,Z13-18:OAc increased catches at Huonville.")])
case("contrast_invented_label", "E8 bound to the bundle", "contrast_card",
     "An invented contrast carrying the real Hungary and Tasmania 1986 ids and texts: not one of the bundle's.",
     {"contrast_label": "S. tipuliformis, probe-invented pair"}, REFUSE, "REFUSED",
     spatch=[(["comparisons", "+"], dict(copy.deepcopy(RAW["comparisons"][2]),
                                         case_label="S. tipuliformis, probe-invented pair"))])
ev("card_role_outside_vocabulary", "E3 roles", "synanthedon_tipuliformis.szocs_1990.trap_e3.hungary", E3,
   "The page's roles call E3 'confirmed_pheromone': not one of the three roles. Refuse rather than ignore it.",
   [], outcomes=REFUSE, prefer="REFUSED", context=ctx(roles=dict(ROLES, **{E3: "confirmed_pheromone"})))
ev("card_atlas_compound_as_counterfactual", "E3 roles", "synanthedon_tipuliformis.szocs_1990.trap_e3.hungary",
   E3, "The page calls an Atlas compound a computational counterfactual: the evidence card refuses the "
   "inconsistent roles, as the class panel does.", [], outcomes=REFUSE, prefer="REFUSED",
   context=ctx(roles=dict(ROLES, **{E3: "computational_counterfactual"})))
ev("doi_not_recorded", "E2 provenance", "synanthedon_tenuis.yang_2012.gland", OAC,
   "The reference without its DOI (synthetic): the authors' words still shown, and the missing DOI said, "
   "so a reader does not think it was left off by the card.",
   [must(r"doi not (recorded|given|stated)", "MISLEADING", "the gap, said")],
   patch=[(["observations", oi("synanthedon_tenuis.yang_2012.gland"), "reference", "doi"], None)])
ev("comparator_value_not_recorded", "E1 absence and negatives",
   "synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1", OH,
   "A response with no comparator value (synthetic): 'not recorded', never 'None'.",
   [must(r"not recorded", "MISLEADING", "the gap, said"), must_not(r"\bNone\b", "MISLEADING", "no raw null")],
   patch=[(["observations", oi("synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1"), "outcome", "response",
            "comparator_value"], None)])
cp("tentative_counts_shown", OAC, TENUIS,
   "A class row with tentative observations (synthetic: n_tentative 2): the panel shows them, as it shows "
   "the significant ones.", [must(r"\b2 tentative|tentative:? 2", "MISLEADING", "the count")],
   patch=[(["compounds", idx, "evidence_classes", jdx, "n_tentative"], 2)])

# --- E10 verbatim record text and the view's own words (probe round 1, finding C) ------------
YT, HU = "synanthedon_tenuis.yang_2012.trap_acetate_added", "synanthedon_tipuliformis.szocs_1990.trap_e3.hungary"
ev("real_quote_inhibited", "E10 verbatim and own words", YT, OAC,
   "The real Yang 2012 quotation says the acetate 'strongly inhibited attraction of males'. It is the "
   "measured result in the authors' words: shown verbatim and attributed, while the view's own words "
   "still never add a role.",
   [must(r"strongly inhibited attraction", "MISLEADING", "the authors' words, unaltered"),
    must(r"the authors'? words", "MISLEADING", "attributed"), NO_OVERREACH, NO_ESP])
for vid, oid, comp, quote in (
        ("inactive", HU, E3, "E3,Z13-18:OAc was inactive in this test."),
        ("repel", YT, OAC, "the blend appeared to repel males at high dose."),
        ("antagonist", YT, OAC, "Z3,Z13-18:OAc acted as a behavioural antagonist."),
        ("because_esp", HU, E3, "catches rose because the electrostatic potential of the E3 isomer explains its fit.")):
    ev(f"quote_says_{vid}", "E10 verbatim and own words", oid, comp,
       f"A quotation (synthetic) containing a banned word ('{quote}'): the authors' words are shown as they "
       "are, attributed; the bans apply to the view's own words.",
       [must(re.escape(quote), "MISLEADING", "the quotation, unaltered"),
        must(r"the authors'? words", "MISLEADING", "attributed"), NO_OVERREACH, NO_ESP],
       patch=[(["observations", oi(oid), "evidence", "author_statement"], quote)])
# probe round 2 (on 48b27f5): record fields interpolated into the view's own sentence must be
# spanned too, one field each, or a banned word in the record lands in the view's words
GY, GN, T41 = ("synanthedon_tenuis.yang_2012.gland", "synanthedon_tenuis.naka_2013.gland",
               "synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1")
ev("population_in_negative_sentence", "E10 verbatim and own words", GY, OAC,
   "A population (synthetic) named 'inactive-strain population', repeated in the not-detected sentence: it "
   "stays the record's words, in a span.",
   [must(r"inactive-strain population", "MISLEADING", "the record's population"), NO_OVERREACH],
   patch=[(["observations", oi(GY), "context", "population"], "inactive-strain population")])
ev("population_resolved_in_negative_sentence", "E10 verbatim and own words", GN, E3,
   "A population (synthetic) named 'resolved-lineage population', in the same sentence.",
   [must(r"resolved-lineage population", "MISLEADING", "the record's population"), NO_ESP],
   patch=[(["observations", oi(GN), "context", "population"], "resolved-lineage population")])
ev("method_in_negative_sentence", "E10 verbatim and own words", GY, OAC,
   "A detection method (synthetic) reading 'gland_extract_inhibition_assay': the record's method string, "
   "verbatim in a span, never the view's words.",
   [must(r"gland_extract_inhibition_assay", "MISLEADING", "the record's method, verbatim"), NO_OVERREACH],
   patch=[(["observations", oi(GY), "components", ci(GY, OAC), "method"], "gland_extract_inhibition_assay")])
ev("response_unit_in_effect_sentence", "E10 verbatim and own words", T41, OH,
   "A response unit (synthetic) reading 'catches_explained_per_trap', interpolated into the effect sentence.",
   [must(r"catches_explained_per_trap", "MISLEADING", "the record's unit"), NO_ESP],
   patch=[(["observations", oi(T41), "outcome", "response", "unit"], "catches_explained_per_trap")])
ev("atlas_statement_says_repel", "E10 verbatim and own words", HU, E3,
   "The Atlas statement itself (synthetic) containing 'repellent': shown verbatim as the Atlas's sentence.",
   [must(r"no repellent effect at Budapest", "MISLEADING", "the Atlas's sentence"), NO_OVERREACH],
   patch=[(["observations", oi(HU), "statement"], "Adding E3,Z13-18:OAc showed no repellent effect at Budapest.")])

# --- E11 controlled words (probe round 3, on 3795473) ---------------------------------------------
# effect, significance and occurrence roles are the Atlas's enums (observation.schema.json,
# teaching-bundle.schema.json). A value outside them is refused (bundle.load or the view; the
# harness reads an adapter refusal as REFUSED); every value inside them is the view's own,
# explicit phrase, never the raw string passed through.
T86 = "synanthedon_tipuliformis.szocs_1990.trap_e3.tasmania_1986"


def occ(comp, oid):
    c = next(i for i, x in enumerate(RAW["compounds"]) if x["atlas_id"] == comp)
    k = next(j for j, o in enumerate(RAW["compounds"][c]["occurrences"]) if o["observation_id"] == P + oid)
    return ["compounds", c, "occurrences", k, "roles"]


for vid, oid, path, value in (
        ("effect_inhibited", HU, ["outcome", "effect"], "inhibited"),
        ("significance_resolved", T86, ["outcome", "significance"], "resolved"),
        ("significance_explains", T86, ["outcome", "significance"], "significant because of ESP")):
    ev(f"vocabulary_{vid}", "E11 controlled words", oid, E3,
       f"An {path[-1]} of {value!r} (synthetic): not one of the Atlas's words. Printed as the view's own "
       "phrase it would be the view's claim; refuse.", [], outcomes=REFUSE, prefer="REFUSED",
       patch=[(["observations", oi(oid)] + path, value)])
for vid, value in (("role_antagonist", ["antagonist"]), ("role_explains", ["inhibitor because ESP"])):
    ev(f"vocabulary_{vid}", "E11 controlled words", HU, E3,
       f"Occurrence roles {value} (synthetic): not the Atlas's role words. Refuse.", [], outcomes=REFUSE,
       prefer="REFUSED", patch=[(occ(E3, HU), value)])
ev("vocabulary_antennal_inactive_phrased", "E11 controlled words", HU, E3,
   "'antennal_inactive' is a real Atlas role word (no antennal response recorded). The view phrases it in "
   "its own words, which never say 'inactive'.",
   [must(r"no antennal response", "MISLEADING", "the role, phrased"), NO_OVERREACH],
   patch=[(occ(E3, HU), ["antennal_inactive"])])
ev("vocabulary_caught_without_control_phrased", "E11 controlled words", HU, E3,
   "'caught_without_control' is a real Atlas effect word: phrased as males caught with no comparator run.",
   [must(r"no comparator", "MISLEADING", "the effect, phrased"), must_not(r"caught_without_control", "MISLEADING",
                                                                       "not the raw enum")],
   patch=[(["observations", oi(HU), "outcome", "effect"], "caught_without_control")])

# --- E12 every printed record value is a controlled word, a number, or verbatim (probe round 4) --
# observation_type and evidence_class are the Atlas's enums; detection-limit values, response
# values and class counts are numbers (bool excluded); detected is a bool; roles are a list.
cls_idx = next(i for i, c in enumerate(RAW["compounds"]) if c["atlas_id"] == E3)
cls_row = next(j for j, e in enumerate(RAW["compounds"][cls_idx]["evidence_classes"])
               if e["evidence_class"] == "addition_decreased" and e["taxon_id"] == TIPU)
CLS = ["compounds", cls_idx, "evidence_classes", cls_row]
for vid, comp, oid, path, value, why in (
        ("observation_type", E3, HU, ["observations", oi(HU), "observation_type"], "repellency_assay",
         "an observation type outside the Atlas's three"),
        ("limit_value_text", OAC, GY, ["observations", oi(GY), "components", ci(GY, OAC), "detection_limit"],
         {"value": "inhibited", "unit": "picogram"}, "a detection limit whose value is text"),
        ("detected_as_text", OAC, GY, ["observations", oi(GY), "components", ci(GY, OAC), "detected"], "false",
         "detected given as the string 'false' (read as true by bool())"),
        ("response_value_bool", OH, T41, ["observations", oi(T41), "outcome", "response", "treatment_value"],
         True, "a response value of true (an int to Python)"),
        ("roles_as_string", E3, HU, occ(E3, HU), "treatment", "roles given as a string, not a list")):
    ev(f"value_{vid}", "E12 printed values", oid, comp, f"{why[0].upper()}{why[1:]} (synthetic): refuse.",
       [], outcomes=REFUSE, prefer="REFUSED", patch=[(path, value)])
for vid, path, value, why in (
        ("class_name_unknown", CLS + ["evidence_class"], "antagonist_effect", "a class name outside the Atlas's"),
        ("class_name_explains", CLS + ["evidence_class"], "resolved_because_esp", "a class name with explanation words"),
        ("class_count_text", CLS + ["n_observations"], "many because ESP", "an observation count given as text"),
        ("class_count_bool", CLS + ["n_significant"], True, "a significant count of true")):
    case(f"panel_{vid}", "E12 printed values", "class_panel", f"{why[0].upper()}{why[1:]} (synthetic): refuse.",
         {"compound": E3, "taxon_id": TIPU}, REFUSE, "REFUSED", patch=[(path, value)])
case("panel_antennal_inactive_class", "E12 printed values", "class_panel",
     "'antennal_inactive' is a real Atlas class: shown as the Atlas's own term, in a span; the panel's own "
     "words never say 'inactive'.",
     {"compound": E3, "taxon_id": TIPU}, SHOW, "RENDERED",
     [must(r"antennal_inactive", "MISLEADING", "the Atlas's class name"), NO_OVERREACH],
     patch=[(CLS + ["evidence_class"], "antennal_inactive")])
case("panel_taxon_not_in_bundle", "E12 printed values", "class_panel",
     "A taxon id the bundle carries in no class row (and named with explanation words): refuse, rather than "
     "print the caller's id as the panel's title.",
     {"compound": OAC, "taxon_id": "spa.tax.because_esp"}, REFUSE, "REFUSED")


def document():
    obs = {o["observation_id"] for o in RAW["observations"]}
    labels = {k["case_label"] for k in RAW["comparisons"]}
    for c in CASES:
        assert set(c["outcomes"].values()) <= set(GRADES), c["id"]
        for q in c["requires"]:
            assert q["grade"] in GRADES and ("must" in q) != ("must_not" in q) and q["on"] in ("text", "html")
        i = c["inputs"]
        if c["cls"].startswith("E7"):
            continue
        if "observation_id" in i and c["prefer"] == "RENDERED":
            assert i["observation_id"] in obs, c["id"]
        if "contrast_label" in i and not c["statement_patch"]:      # an appended fabrication is not
            assert i["contrast_label"] in labels, c["id"]
    return json.dumps({"generator": "tests/pedagogy/make_evidence_card_set.py",
                       "bundle": "tests/data/sesiidae-homework-1.bundle.json", "atlas_version": RAW["atlas_version"],
                       "grades": list(GRADES), "cases": CASES}, indent=1, ensure_ascii=False) + "\n"


if __name__ == "__main__":
    OUT.write_text(document(), encoding="utf-8")
    print(f"{len(CASES)} cases -> {OUT}")
