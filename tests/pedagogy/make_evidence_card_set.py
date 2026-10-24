"""Generate the hostile set for the evidence views:
tests/pedagogy/evidence_card_set.json.

    python tests/pedagogy/make_evidence_card_set.py     # byte-identical on re-run

Written before the views exist (expectations first). It re-expresses the
nine abstract evidence cases of hostile_set.json (B1-B5, which predate
bundle.py) against the real Atlas teaching bundle, as bundle.py reads it, and
extends them to the Atlas's evidence classes and population contexts. Inputs
are real records: tests/data/sesiidae-homework-1.bundle.json (Atlas 30c4a11),
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
optional "study_rows": study verdicts the page also shows}.

Grades as in the other sets (SAFE, SAFE_REFUSAL, CRASH, MISLEADING, FALSE);
`allowed` is exactly the safe outputs; a rendered view must meet every
`requires` entry (must / must_not, a case-insensitive pattern on the text
with tags removed, or on the HTML). The Atlas's rules this set holds the
views to: no class is netted out; "not detected" is a negative for that
method, sample and population only; absence of a class row is not evidence
of absence; female_pheromone_component is the Atlas's join of two lines of
evidence, not an identification, and its population_match must be read.
"""
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


def case(id, cls, component, why, inputs, outcomes, prefer, requires=(), context=None, patch=()):
    outcomes = dict(outcomes)
    outcomes.setdefault("ERROR", "CRASH")
    allowed = sorted(k for k, g in outcomes.items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert prefer is None or prefer in allowed, id
    CASES.append(dict(id=id, cls=cls, component=component, why=why, inputs=inputs,
                      context=context if context is not None else ctx(), bundle_patch=[list(p) for p in patch],
                      outcomes=outcomes, allowed=allowed, prefer=prefer, requires=list(requires)))


def must(pat, grade, why, on="text"):
    return {"must": pat, "grade": grade, "why": why, "on": on}


def must_not(pat, grade, why, on="text"):
    return {"must_not": pat, "grade": grade, "why": why, "on": on}


SHOW = {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}
REFUSE = {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"}
# never on an evidence view: interpretation the records do not carry
NO_OVERREACH = must_not(r"inactive|does not (make|produce|attract)|not (made|produced)|absent from|"
                        r"repel|inhibitor|antagonist", "FALSE",
                        "an evidence view reports what was measured, never a role the record does not state")
NO_NETTING = must_not(r"\boverall\b|in general|on balance|\bnet\b|taken together|consensus", "MISLEADING",
                      "the Atlas never nets opposite results out; neither does the view")
NO_ESP = must_not(r"\bESP\b|electrostatic|kcal|resolved|because|explain|due to|caus", "FALSE",
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
ev("atlas_version_missing", "E2 provenance", "synanthedon_tenuis.naka_2013.gland", OAC,
   "Without the Atlas commit the card cannot be traced: a visible warning, or refusal.",
   [must(r"provenance (incomplete|missing)|atlas (version|commit) (not recorded|missing|unknown)",
         "MISLEADING", "a visible warning")],
   context=ctx(bundle={"atlas_version": None, "sha256": None}))

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
ev("observation_not_in_bundle", "E7 refusals", "synanthedon_tenuis.not_in_this_bundle", OAC,
   "An observation id the bundle does not carry: refuse, never guess.", [], outcomes=REFUSE, prefer="REFUSED")
case("contrast_not_in_bundle", "E7 refusals", "contrast_card",
     "A contrast label the bundle does not carry: refuse.", {"contrast_label": "no such contrast"},
     REFUSE, "REFUSED")
case("class_panel_unknown_compound", "E7 refusals", "class_panel",
     "A compound neither in the bundle nor given a role: refuse.",
     {"compound": "spa.cmpd.not_here", "taxon_id": TENUIS}, REFUSE, "REFUSED")


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
        if "contrast_label" in i:
            assert i["contrast_label"] in labels, c["id"]
    return json.dumps({"generator": "tests/pedagogy/make_evidence_card_set.py",
                       "bundle": "tests/data/sesiidae-homework-1.bundle.json", "atlas_version": RAW["atlas_version"],
                       "grades": list(GRADES), "cases": CASES}, indent=1, ensure_ascii=False) + "\n"


if __name__ == "__main__":
    OUT.write_text(document(), encoding="utf-8")
    print(f"{len(CASES)} cases -> {OUT}")
