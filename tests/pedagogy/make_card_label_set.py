"""Generate the hostile set for the difference card's label:
tests/pedagogy/card_label_set.json.

    python tests/pedagogy/make_card_label_set.py     # byte-identical on re-run

Written before the fix (expectations first). The difference card
(views.card_html over model.compare) gives a teaching verdict under its own
rule: "separated" when the conformer ranges do not overlap. The frozen study
(espmap-sesiidae-001) decides under another (range_partition_sign_v1, with
its own partition sensitivity) and words its outcomes differently
(resolved, different on average, unresolved). Measured on the real
ensembles (every site the study reports, plus Act A's chain sites, which it
does not), the two rules never disagree where the study has a row, but:

- at Act A's chain sites the study has no row, and the card says
  "separated: every conformer differs in the same direction" for the
  acetate and the alcohol against the aldehyde (a constructed teaching
  control) at C1-H and C2-H: the sites the homework sends students to;
- the study's "different on average" and "unresolved" both read on the
  card as "not attributable to the change", so B's C3-H looks like a
  contradiction when the two are side by side.

The card's contract after the fix:

    card_html(card, context=None) -> str

context, optional: {"study": StudyResults (study.load_results), "contrast":
the study's contrast id, "x", "y": the slugs of the card's first and second
molecule, "roles": {slug: "atlas_observed" | "teaching_control" |
"computational_counterfactual"}}. The card looks its row up in the study
itself (contrast, its own site, population "all"); it never takes a tier or
an outcome from the caller. A study whose contrast is not (x, y), in that
order, or has no such contrast, is refused (model.Refused): the card's
difference would carry the other sign, or describe another pair. A study
that is not a StudyResults is refused too.

Every rendered card, with or without context:
- calls itself a teaching comparison, a look and not a prespecified test,
  and attributes its verdict to the card's own rule;
- never uses the study's outcome words (resolved, unresolved, different on
  average) for anything;
- with a study row: points to the study card for the study's verdict, and
  says the study decides under a different rule; an exploratory row is
  called exploratory;
- with a study but no row for this site: says the study gives no verdict
  here;
- with a teaching control on either side: says it is constructed, and the
  comparison descriptive only;
- with no study in the context: names no study.

Inputs are real: per-conformer site values extracted from the Atlas
ensembles at 35112e1 and the aldehyde control (tests/data/card_label_sites.json,
each member with its record's sha256), the pilot butenes (records/), and the
study's results file (tests/data/espmap-sesiidae-001.results.json).

Grades as in the other sets; outputs RENDERED, REFUSED, ERROR; `allowed` is
exactly the safe outputs; a rendered card must meet every `requires` entry
(a case-insensitive pattern on the text with tags removed).
"""
import json
from pathlib import Path

OUT = Path(__file__).with_name("card_label_set.json")
GRADES = ("SAFE", "SAFE_REFUSAL", "CRASH", "MISLEADING", "FALSE")
OAC, E3, E2, OH, AL = "z3_z13_18_oac", "e3_z13_18_oac", "e2_z13_18_oac", "z3_z13_18_oh", "z3_z13_18_al"
A, A_AL_OAC, A_AL_OH = "A_headgroup", "A_headgroup (control: Al vs OAc)", "A_headgroup (control: Al vs OH)"
ROLES = {OAC: "atlas_observed", E3: "atlas_observed", E2: "atlas_observed", OH: "atlas_observed",
         AL: "teaching_control"}
CASES = []


def ctx(contrast, x, y, study="results", **kw):
    c = {"study": study, "contrast": contrast, "x": x, "y": y, "roles": dict(ROLES)}
    c.update(kw)
    return c


def case(id, cls, why, x, y, site, context, outcomes, prefer, requires=(), measured=None, card_patch=(),
         study_patch=()):
    outcomes = dict(outcomes)
    outcomes.setdefault("ERROR", "CRASH")
    allowed = sorted(k for k, g in outcomes.items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert prefer is None or prefer in allowed, id
    CASES.append(dict(id=id, cls=cls, why=why, x=x, y=y, site=site, context=context, outcomes=outcomes,
                      allowed=allowed, prefer=prefer, requires=list(requires), measured=measured,
                      card_patch=[list(p) for p in card_patch], study_patch=[list(p) for p in study_patch]))


def must(pat, grade, why):
    return {"must": pat, "grade": grade, "why": why}


def must_not(pat, grade, why):
    return {"must_not": pat, "grade": grade, "why": why}


SHOW = {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}
RENDER = {"RENDERED": "SAFE", "REFUSED": "MISLEADING"}       # refusing hides a legitimate look
REFUSE = {"RENDERED": "MISLEADING", "REFUSED": "SAFE_REFUSAL"}

LABEL = must(r"teaching comparison", "MISLEADING", "the card says what it is")
NOT_A_TEST = must(r"not a (prespecified )?test|a look,? not a test", "MISLEADING", "a look, not a test")
CARD_RULE = must(r"(the |this )?card'?s (own )?rule", "MISLEADING", "the verdict is the card's rule's, not the study's")
# Bans read the card's own words (the user, 2026-10-24): record text (the site label, the
# study's plan reference) is shown verbatim inside class="verbatim" and is exempt; the harness
# checks each verbatim span is record text. Word families: resol(v|ut)*, unresol(v|ut)*.
NO_STUDY_WORDS = must_not(r"\bresol(v|ut)\w*|\bunresol(v|ut)\w*|different on average", "FALSE",
                          "the study's outcome words belong to the study's rule and the study card")
ALWAYS = [LABEL, NOT_A_TEST, CARD_RULE, NO_STUDY_WORDS]
TO_STUDY_CARD = must(r"study card", "MISLEADING", "where the study's verdict is")
OTHER_RULE = must(r"(a |its own |another |a different )rule|rules? differ", "MISLEADING",
                  "the study decides under a different rule")
NO_ROW = must(r"(gives|has|makes) no (verdict|test|comparison)|no (study )?verdict (here|at this site)|"
              r"does not test this site", "MISLEADING", "the study is silent here, and the card says so")
CONTROL = must(r"constructed|teaching control", "MISLEADING", "the aldehyde is not an Atlas compound")
DESCRIPTIVE = must(r"descriptive", "MISLEADING", "a control's comparison is descriptive chemistry only")
NO_CONTROL = must_not(r"teaching control|constructed", "FALSE", "no control on this card")

# --- L1 sites the study does not test (Act A's chain) ------------------------------------
case("control_c1_separated", "L1 no study row",
     "Acetate vs the aldehyde control at C1-H: the card measures 'separated'; the study has no row at any "
     "Act A chain site. Unlabelled, 'every conformer differs in the same direction' reads as a result.",
     OAC, AL, "C1-H", ctx(A_AL_OAC, OAC, AL), RENDER, "RENDERED",
     ALWAYS + [NO_ROW, CONTROL, DESCRIPTIVE], measured="SEPARATED")
case("control_c2_separated", "L1 no study row",
     "The alcohol vs the aldehyde control at C2-H: 'separated', no study row.",
     OH, AL, "C2-H", ctx(A_AL_OH, OH, AL), RENDER, "RENDERED",
     ALWAYS + [NO_ROW, CONTROL, DESCRIPTIVE], measured="SEPARATED")
case("chain_no_row", "L1 no study row",
     "Acetate vs alcohol at C3-H (Act A item 2's chain comparison): the ranges overlap, the study has no row.",
     OAC, OH, "C3-H", ctx(A, OAC, OH), RENDER, "RENDERED", ALWAYS + [NO_ROW, NO_CONTROL], measured="OVERLAPPING")

# --- L2 sites the study decides, under its own rule --------------------------------------
case("b_c3_primary_vocabulary", "L2 study row",
     "B C3-H, primary: the study says 'different on average'; the card's rule says 'not attributable to the "
     "change'. Side by side they read as a contradiction unless the card says the study decides under "
     "another rule, on the study card.",
     OAC, E3, "C3-H", ctx("B_stereo", OAC, E3), RENDER, "RENDERED",
     ALWAYS + [TO_STUDY_CARD, OTHER_RULE, must(r"prespecified", "MISLEADING", "a primary site is a test")],
     measured="OVERLAPPING")
case("c_c4_primary_agrees", "L2 study row",
     "C C4-H, primary: the study resolves it and the card says 'separated'. Agreement does not make the card "
     "the test: it still points to the study card and never says 'resolved'.",
     E3, E2, "C4-H", ctx("C_position", E3, E2), RENDER, "RENDERED", ALWAYS + [TO_STUDY_CARD, OTHER_RULE],
     measured="SEPARATED")
case("b_c5_exploratory", "L2 study row",
     "B C5-H is an exploratory row of the study: the card calls it exploratory.",
     OAC, E3, "C5-H", ctx("B_stereo", OAC, E3), RENDER, "RENDERED",
     ALWAYS + [must(r"exploratory", "MISLEADING", "the row's tier")], measured="OVERLAPPING")
case("head_o_control_secondary", "L2 study row",
     "The alcohol vs the aldehyde at head O, a secondary row: the card's rule says the direction depends on "
     "the partition rule; the study, 'unresolved'. Control and study card both named.",
     OH, AL, "head O", ctx(A_AL_OH, OH, AL), RENDER, "RENDERED",
     ALWAYS + [TO_STUDY_CARD, CONTROL, DESCRIPTIVE], measured="RULE_DEPENDENT")

# --- L3 a study context that is not this card's ------------------------------------------
case("reversed_order", "L3 study mismatch",
     "E3 → Z3 on the card, but B is Z3 → E3: the card's difference has the opposite sign to the study's. "
     "Pointing to the study card as if they matched invites a sign error. Refuse.",
     E3, OAC, "C3-H", ctx("B_stereo", E3, OAC), REFUSE, "REFUSED", measured="OVERLAPPING")
case("wrong_contrast", "L3 study mismatch",
     "Acetate vs aldehyde at C1-H handed B's contrast (Z3 vs E3): another pair. Refuse.",
     OAC, AL, "C1-H", ctx("B_stereo", OAC, AL), REFUSE, "REFUSED", measured="SEPARATED")
case("unknown_contrast", "L3 study mismatch",
     "A contrast id the study does not have: refuse, never guess a tier.",
     OAC, E3, "C3-H", ctx("D_new", OAC, E3), REFUSE, "REFUSED", measured="OVERLAPPING")
case("study_not_results", "L3 study mismatch",
     "The study given as the raw results JSON, not a StudyResults read by study.load_results (which refuses "
     "corrupt and partial files): refuse.",
     OAC, E3, "C3-H", ctx("B_stereo", OAC, E3, study="raw"), REFUSE, "REFUSED", measured="OVERLAPPING")

# --- L4 no study on the page, or a context that cannot be read ---------------------------
case("butene_no_study", "L4 no study",
     "The organic notebook's E/Z-butene card (alkene C2, 'separated'): no study on the page. Labelled a "
     "teaching comparison, naming no study.",
     "pilot:e2_butene", "pilot:z2_butene", "C2", None, RENDER, "RENDERED",
     ALWAYS + [must_not(r"study card|frozen study|espmap-sesiidae", "MISLEADING", "no study here")],
     measured="SEPARATED")
case("context_not_mapping", "L4 no study",
     "A context that is a string: read as no context. The card renders its label and claims nothing about "
     "a study.",
     OAC, AL, "C1-H", "A_headgroup (control: Al vs OAc)", SHOW, "RENDERED",
     ALWAYS + [must_not(r"study card|prespecified", "MISLEADING", "no study read")], measured="SEPARATED")
case("role_outside_vocabulary", "L4 no study",
     "The aldehyde's role given as 'control-ish': outside the vocabulary. Refuse rather than drop the "
     "control label.",
     OAC, AL, "C1-H", ctx(A_AL_OAC, OAC, AL, roles=dict(ROLES, **{AL: "control-ish"})), REFUSE, "REFUSED",
     measured="SEPARATED")

# --- held-out probe round 1 (Publishing, on ab58c1a): its findings, as cases ---------------
# card_patch / study_patch: (field, value) applied with dataclasses.replace to the Comparison or
# the StudyResults, for inputs compare() and load_results() never produce themselves.
case("card_pair_not_context_pair", "L3 study mismatch",
     "The card compares Z3 → E2 at C4-H; the context names B_stereo's own pair (Z3 → E3). The study's "
     "primary-test status would be attached to a pair it never tested. Refuse unless the card's molecules "
     "are the context's.",
     OAC, E2, "C4-H", ctx("B_stereo", OAC, E3), REFUSE, "REFUSED", measured="SEPARATED")
case("card_pair_swapped_against_context", "L3 study mismatch",
     "The card compares E3 → Z3; the context declares Z3 → E3, the study's order. The card's sign is the "
     "reverse of the study's direction. Refuse.",
     E3, OAC, "C3-H", ctx("B_stereo", OAC, E3), REFUSE, "REFUSED", measured="OVERLAPPING")
case("study_unknown_rule", "L3 study mismatch",
     "A StudyResults whose verdict rule is not one Arrowfield knows (model.KNOWN_VERDICT_RULES): never "
     "describe a study's tiers under a rule it did not use. Refuse, as the study card does.",
     OAC, E3, "C3-H", ctx("B_stereo", OAC, E3), REFUSE, "REFUSED", measured="OVERLAPPING",
     study_patch=[("verdict_rule", "someone-else/v9")])

# --- L6 a record label that carries an outcome word (probe round 1, item 4) ---------------------
case("outcome_word_in_site_label", "L6 verbatim labels",
     "A site label (synthetic) reading 'resolved site': the record's label, shown verbatim; the card's "
     "own words still never use an outcome word.",
     OAC, E3, "C3-H", None, RENDER, "RENDERED",
     ALWAYS + [must(r"resolved site", "MISLEADING", "the record's label, unaltered")], measured="OVERLAPPING",
     card_patch=[("site", "resolved site")])
case("outcome_word_in_site_label_with_study", "L6 verbatim labels",
     "The same label with the study on the page: the site is not one of the study's, so the card says the "
     "study gives no verdict here, in its own words.",
     OAC, E3, "C3-H", ctx("B_stereo", OAC, E3), RENDER, "RENDERED",
     ALWAYS + [must(r"resolved site", "MISLEADING", "the record's label"), NO_ROW], measured="OVERLAPPING",
     card_patch=[("site", "resolved site")])

# --- L5 numbers the card cannot print ----------------------------------------------------
case("nan_difference", "L5 values",
     "A Comparison whose difference is NaN (compare() refuses these; a hand-built one need not): never "
     "print 'nan' as a difference.",
     OAC, E3, "C3-H", None, SHOW, "REFUSED", [must_not(r"\bnan\b", "MISLEADING", "no raw NaN")],
     measured="OVERLAPPING", card_patch=[("diff", float("nan"))])
case("infinite_partition", "L5 values",
     "A partition sensitivity of +inf: not a number a reader can weigh a difference against.",
     OAC, E3, "C3-H", None, SHOW, "REFUSED", [must_not(r"\binf\b", "MISLEADING", "no raw inf")],
     measured="OVERLAPPING", card_patch=[("partition", float("inf"))])


def document():
    ids = [c["id"] for c in CASES]
    assert len(ids) == len(set(ids))
    for c in CASES:
        assert set(c["outcomes"].values()) <= set(GRADES), c["id"]
        for q in c["requires"]:
            assert q["grade"] in GRADES and ("must" in q) != ("must_not" in q)
    return json.dumps({"generator": "tests/pedagogy/make_card_label_set.py", "grades": list(GRADES),
                       "cases": CASES}, indent=1, ensure_ascii=False) + "\n"


if __name__ == "__main__":
    OUT.write_text(document(), encoding="utf-8")
    print(f"{len(CASES)} cases -> {OUT}")
