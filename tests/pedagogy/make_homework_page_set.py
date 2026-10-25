"""Generate the hostile set for the homework page builder:
tests/pedagogy/homework_page_set.json.

    python tests/pedagogy/make_homework_page_set.py     # byte-identical on re-run

Written before the builder exists (expectations first). The builder turns
the homework's plan (sesiidae-homework.md, as data) and its sources into
one page that assembles the views already built: study cards, labelled
difference cards, the pair view, evidence cards, contrast cards and class
panels. It is where four things the views cannot settle alone are settled:

- ROLES. A molecule has one role on a page, from the plan alone, checked
  against the data: an Atlas compound must be a compound of the bundle, a
  teaching control must not be, and an ensemble's own source id must agree
  (arrowfield.control.* for a control). Every block shows the role of each
  molecule it holds ("Every page and card must show which of the two kinds a
  compound is", sesiidae-homework.md). A block never carries its own role.
- ACTS. Exactly A, B, C, D, in that order, ids unique. An act's study cards
  and difference cards belong to its own contrast (act A: A_headgroup and its
  control contrasts; B: B_stereo; C: C_position); act D holds no view (it is
  the act about what free-molecule ESP cannot show). In each act the
  prediction comes before any view. The "expected" answers appear only on
  the instructor's page, labelled for the instructor, never on a student's.
- PROVENANCE. One page, one bundle, one results file, one analysis key (the
  homework: "one analysis_key throughout"). A footer states the Atlas commit,
  the results file's sha256, the espmap version(s) and the analysis key(s).
  A results file under a rule Arrowfield does not know refuses the page; a
  bundle that names no commit renders with the provenance marked incomplete.
- STUDY / TEACHING. The study's outcome words (resolved, unresolved,
  different on average, and their families) appear only inside study cards:
  never in the page's own prose, a difference card or a pair view. The
  evidence views keep their own frozen vocabulary (#12-#13), which uses the
  Atlas's term "Unresolved in the record" for open items: the check does not
  read inside them (a premise conflict found before freezing). A difference card at a site the act's contrast tests (a primary or
  secondary row) needs that site's study card in the same act, so its
  pointer "on the study card" resolves on the page; an exploratory site needs
  none. Nothing the page itself writes joins an ESP difference to a
  biological result.

Contract:

    homework.page_html(plan, sources, audience="student") -> str, or Refused

plan: {"case_id", "title", "molecules": {slug: {"name", "role":
atlas_observed | teaching_control, "atlas_id"}}, "acts": [{"id", "contrasts":
[study contrast ids], "title", "predict", "tasks": [text], "expected", "blocks":
[block]}]}, block: {"view": study_card | difference_card | pair_view |
evidence_card | contrast_card | class_panel, ...its arguments}. sources:
{"bundle": TeachingBundle, "study": StudyResults, "ensembles": {slug:
(EnsembleSource, [MoleculeView])}, "conformers": {"slug/conf": MoleculeView}}.

Plan text (titles, prompts, tasks, expected answers) is authored teaching
text, reviewed by a person: shown verbatim in class="verbatim" spans and not
linted. The page's own words are under the bans.

The rendered page carries, for the harness:
- one element per act, class "act", data-act; inside it, one class "predict"
  element before every class "block" element;
- one element per block, class "block", data-view, data-molecules (slugs,
  comma-joined) and data-roles (their roles, same order);
- one element class "page-provenance" with data-atlas, data-results-sha256,
  data-espmap and data-analysis-keys (comma-joined);
- study cards wrapped in class "study-card" (as study_card_html writes them).

Inputs are real: the plan transcribed from sesiidae-homework.md; the Atlas
bundle (35112e1); the study's results file; per-conformer site values for the
five ensembles (tests/data/card_label_sites.json, one analysis key d0f5ed13,
the control's source id arrowfield.control.z3_z13_18_al); the E3 conformers
in tests/data/pair. A case patches the plan ((path, value); a path ending in
"+" appends, DROP removes) or the sources (a named op) or sets the audience.
"""
import copy
import json
import re
from pathlib import Path

OUT = Path(__file__).with_name("homework_page_set.json")
GRADES = ("SAFE", "SAFE_REFUSAL", "CRASH", "MISLEADING", "FALSE")
DROP = "<drop>"
OAC, OH, AL, E3, E2 = "z3_z13_18_oac", "z3_z13_18_oh", "z3_z13_18_al", "e3_z13_18_oac", "e2_z13_18_oac"
P = "spa.obs."
CASES = []

PLAN = {
    "case_id": "sesiidae-homework-1",
    "title": "Sesiidae homework: when does molecular electronics explain biological specificity?",
    "molecules": {
        OH: {"name": "(3Z,13Z)-octadeca-3,13-dien-1-ol", "role": "atlas_observed", "atlas_id": "spa.cmpd.z3_z13_18_oh"},
        AL: {"name": "(3Z,13Z)-octadeca-3,13-dienal", "role": "teaching_control", "atlas_id": None},
        OAC: {"name": "(3Z,13Z)-octadeca-3,13-dien-1-yl acetate", "role": "atlas_observed",
              "atlas_id": "spa.cmpd.z3_z13_18_oac"},
        E3: {"name": "(3E,13Z)-octadeca-3,13-dien-1-yl acetate", "role": "atlas_observed",
             "atlas_id": "spa.cmpd.e3_z13_18_oac"},
        E2: {"name": "(2E,13Z)-octadeca-2,13-dien-1-yl acetate", "role": "atlas_observed",
             "atlas_id": "spa.cmpd.e2_z13_18_oac"},
    },
    "acts": [
        {"id": "A", "contrasts": ["A_headgroup", "A_headgroup (control: Al vs OAc)", "A_headgroup (control: Al vs OH)"],
         "title": "What the head group does",
         "predict": "Where is each molecule most electron-rich? Which has an H-bond donor? Does adding the acetyl "
                    "make \"everything more negative\"?",
         "tasks": ["Identify the largest electrostatic differences between the three head groups and explain them "
                   "with ordinary organic chemistry: electronegativity, resonance in the ester, donor and acceptor.",
                   "Explore how far the head-group perturbation appears to extend, with the difference card. These "
                   "chain-site comparisons are descriptive and exploratory, not a prespecified test."],
         "expected": "This act is the \"ESP works\" act. Global V_max: resolved, the prespecified result.",
         "blocks": [{"view": "study_card", "contrast": "A_headgroup", "site": "global V_max"},
                    {"view": "study_card", "contrast": "A_headgroup", "site": "head O"},
                    {"view": "difference_card", "contrast": "A_headgroup", "x": OAC, "y": OH, "site": "head O"},
                    {"view": "difference_card", "contrast": "A_headgroup (control: Al vs OAc)", "x": OAC, "y": AL,
                     "site": "C1-H"},
                    {"view": "contrast_card", "label": "S. tenuis gland, Japan and Korea"},
                    {"view": "class_panel", "molecule": OAC, "taxon_id": "spa.tax.synanthedon_tenuis"}]},
        {"id": "B", "contrasts": ["B_stereo"], "title": "What E/Z does (E3 vs Z3 acetate)",
         "predict": "Does changing E to Z at C3 mainly change the charge distribution, or the shape?",
         "tasks": ["Read the difference card at C3-H, C4-H and head O. Is each difference larger than that site's "
                   "conformer range? Than its partition sensitivity?",
                   "Look at the extended and folded conformers of one compound side by side, on the same scale."],
         "expected": "E/Z is not resolved at any site. C3-H: different on average.",
         "blocks": [{"view": "study_card", "contrast": "B_stereo", "site": "C3-H"},
                    {"view": "difference_card", "contrast": "B_stereo", "x": OAC, "y": E3, "site": "C3-H"},
                    {"view": "pair_view", "left": {"slug": E3, "conformer": "conf_00"},
                     "right": {"slug": E3, "conformer": "conf_02"}}]},
        {"id": "C", "contrasts": ["C_position"], "title": "Where the double bond is (E2 vs E3 acetate)",
         "predict": "Moving C=C by one carbon. Which sites should change, and by how much compared with act B?",
         "tasks": ["Compare C2-H, C3-H and C4-H between the two compounds, with the difference card. Which "
                   "differences exceed the conformer range?",
                   "Does free-molecule ESP give a convincing reason an organism or enzyme would treat these two "
                   "compounds differently?"],
         "expected": "Resolved at C4 only. Question 2's honest answer stays no.",
         "blocks": [{"view": "study_card", "contrast": "C_position", "site": "C4-H"},
                    {"view": "difference_card", "contrast": "C_position", "x": E3, "y": E2, "site": "C4-H"},
                    {"view": "evidence_card", "observation_id": P + "synanthedon_tipuliformis.szocs_1990.trap_e3.hungary",
                     "molecule": E3}]},
        {"id": "D", "contrasts": [], "title": "What is missing",
         "predict": "The organism makes, or responds to, one particular compound. Suppose free-molecule ESP does "
                    "not explain the distinction. What information is missing, and how would you get it?",
         "tasks": [], "expected": "The receptor's or enzyme's geometry, active-site residues, orientation, "
                                  "dynamics, and experimental evidence.",
         "blocks": []},
    ],
}


def case(id, cls, why, outcomes, prefer, requires=(), plan_patch=(), source_ops=(), audience="student"):
    outcomes = dict(outcomes)
    outcomes.setdefault("ERROR", "CRASH")
    allowed = sorted(k for k, g in outcomes.items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert prefer is None or prefer in allowed, id
    CASES.append(dict(id=id, cls=cls, why=why, plan_patch=[list(p) for p in plan_patch],
                      source_ops=list(source_ops), audience=audience, outcomes=outcomes, allowed=allowed,
                      prefer=prefer, requires=list(requires)))


def must(pat, grade, why):
    return {"must": pat, "grade": grade, "why": why}


def must_not(pat, grade, why):
    return {"must_not": pat, "grade": grade, "why": why}


def check(what, grade, why):
    return {"check": what, "grade": grade, "why": why}


RENDER = {"RENDERED": "SAFE", "REFUSED": "MISLEADING"}       # a sound plan must render
REFUSE = {"RENDERED": "MISLEADING", "REFUSED": "SAFE_REFUSAL"}
REFUSE_F = {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"}
ROLES_SHOWN = check("block_roles", "FALSE", "every block names each molecule's one role, the plan's")
ACT_ORDER = check("act_order", "MISLEADING", "acts A-D in order; each act's prediction before its views")
FOOTER = check("footer", "MISLEADING", "one provenance footer naming the commit, results sha, espmap, keys")
STUDY_SPLIT = check("study_words_in_study_cards", "FALSE", "the study's outcome words only inside study cards")
NO_BRIDGE = check("no_bridge", "FALSE", "nothing the page writes joins an ESP difference to a biological result")
SOUND = [ROLES_SHOWN, ACT_ORDER, FOOTER, STUDY_SPLIT, NO_BRIDGE]


def blk(act, k):
    """Path to block k of act index act."""
    return ["acts", act, "blocks", k]


# --- H1 the real plan --------------------------------------------------------------------------
case("real_plan_student", "H1 the real plan",
     "The homework as written, for a student: renders, with every rule below holding at once, and no "
     "expected answer on the page.",
     RENDER, "RENDERED", SOUND + [
         must(r"teaching control", "FALSE", "the aldehyde named as a teaching control"),
         must_not(r"ESP works|honest answer stays no|for the instructor", "MISLEADING",
                  "no expected answers on a student's page")])
case("real_plan_instructor", "H1 the real plan",
     "The same plan for the instructor: the expected answers appear, labelled for the instructor.",
     RENDER, "RENDERED", SOUND + [must(r"for the instructor", "MISLEADING", "labelled"),
                                  must(r"ESP works", "MISLEADING", "the expected answer, verbatim")],
     audience="instructor")

# --- H2 roles: one per molecule, from the plan, checked against the data ------------------------
case("role_conflict_in_block", "H2 roles",
     "A block carrying its own role for the aldehyde ('atlas_observed'): the page's one role would be "
     "contradicted on one card. A block never carries a role: refuse.",
     REFUSE_F, "REFUSED", plan_patch=[(blk(0, 3) + ["roles"], {AL: "atlas_observed"})])
case("control_given_an_atlas_id", "H2 roles",
     "The aldehyde stays a teaching control but is given the acetate's Atlas id: a control the bundle "
     "carries is not a control.", REFUSE_F, "REFUSED",
     plan_patch=[(["molecules", AL, "atlas_id"], "spa.cmpd.z3_z13_18_oac")])
case("observed_not_in_bundle", "H2 roles",
     "The alcohol called an Atlas compound under an id the bundle does not carry.", REFUSE_F, "REFUSED",
     plan_patch=[(["molecules", OH, "atlas_id"], "spa.cmpd.not_in_this_bundle")])
case("control_called_observed", "H2 roles",
     "The aldehyde called an Atlas compound: its ensemble's own source id is arrowfield.control.*. "
     "The data contradicts the plan; refuse.", REFUSE_F, "REFUSED",
     plan_patch=[(["molecules", AL, "role"], "atlas_observed"),
                 (["molecules", AL, "atlas_id"], "spa.cmpd.z3_z13_18_al")])
case("role_outside_vocabulary", "H2 roles",
     "A role of 'control-ish': outside the vocabulary.", REFUSE, "REFUSED",
     plan_patch=[(["molecules", AL, "role"], "control-ish")])
case("molecule_not_in_plan", "H2 roles",
     "A difference card naming a slug the plan's molecules do not list: no role to show. Refuse.",
     REFUSE, "REFUSED", plan_patch=[(blk(1, 1) + ["y"], "z3_z13_18_unknown")])

# --- H3 acts ------------------------------------------------------------------------------------
case("acts_out_of_order", "H3 acts", "Act B before act A: the homework's order is the argument. Refuse.",
     REFUSE, "REFUSED", plan_patch=[(["acts", 0, "id"], "B"), (["acts", 1, "id"], "A")])
case("act_missing", "H3 acts", "No act C: the plan is not the homework. Refuse.", REFUSE, "REFUSED",
     plan_patch=[(["acts", 2], DROP)])
case("act_duplicated", "H3 acts", "Two acts B.", REFUSE, "REFUSED", plan_patch=[(["acts", 2, "id"], "B")])
case("study_card_from_another_act", "H3 acts",
     "Act B showing C_position's C4-H study card: another act's verdict. Refuse.", REFUSE_F, "REFUSED",
     plan_patch=[(blk(1, 0), {"view": "study_card", "contrast": "C_position", "site": "C4-H"})])
case("difference_card_other_molecules", "H3 acts",
     "Act B's card comparing E3 with E2 under B_stereo: not the act's pair. Refuse.", REFUSE_F, "REFUSED",
     plan_patch=[(blk(1, 1) + ["x"], E2)])
case("view_in_act_d", "H3 acts",
     "A difference card in act D, the act about what free-molecule ESP cannot show. Refuse.", REFUSE, "REFUSED",
     plan_patch=[(["acts", 3, "blocks", "+"], {"view": "difference_card", "contrast": "C_position", "x": E3,
                                               "y": E2, "site": "C4-H"})])
case("prediction_missing", "H3 acts", "Act C without its prediction: the reveal would come first. Refuse.",
     REFUSE, "REFUSED", plan_patch=[(["acts", 2, "predict"], "")])

# --- H4 provenance ------------------------------------------------------------------------------
case("analysis_keys_differ", "H4 provenance",
     "One ensemble under another analysis key: the homework requires one key throughout. Refuse.",
     REFUSE_F, "REFUSED", source_ops=["e2_other_key"])
case("results_unknown_rule", "H4 provenance",
     "A results file under a verdict rule Arrowfield does not know: no study card can be shown. Refuse "
     "the page.", REFUSE_F, "REFUSED", source_ops=["results_unknown_rule"])
case("bundle_without_commit", "H4 provenance",
     "A bundle whose atlas_version is not a commit: the page renders, its footer saying the provenance is "
     "incomplete.", RENDER, "RENDERED",
     [must(r"provenance incomplete", "MISLEADING", "said, in the footer"), FOOTER], source_ops=["bundle_unknown"])
case("conformer_not_in_ensemble", "H4 provenance",
     "The pair view's conformer from outside the molecule's sources: no record to show. Refuse.",
     REFUSE, "REFUSED", plan_patch=[(blk(1, 2) + ["right", "conformer"], "conf_99")])

# --- H5 the study / teaching split --------------------------------------------------------------
case("tested_site_without_study_card", "H5 study and teaching",
     "Act B's difference card at C3-H with no C3-H study card in the act: the card points to a study card "
     "the page does not hold. Refuse.", REFUSE, "REFUSED", plan_patch=[(blk(1, 0), DROP)])
case("exploratory_site_without_study_card", "H5 study and teaching",
     "Act A's card at the aldehyde's C1-H, a site no study row tests: it needs no study card. Renders.",
     RENDER, "RENDERED", SOUND, plan_patch=[(blk(0, 0), DROP), (blk(0, 0), DROP), (blk(0, 0), DROP)])
case("study_card_site_not_in_study", "H5 study and teaching",
     "A study card for a site the study has no row for (B_stereo, C9-H is exploratory; C20-H does not "
     "exist). Refuse.", REFUSE, "REFUSED",
     plan_patch=[(blk(1, 0), {"view": "study_card", "contrast": "B_stereo", "site": "C20-H"})])

# --- H6 nothing bridges ESP and biology ----------------------------------------------------------
case("evidence_beside_resolved_card", "H6 no bridge",
     "Act C as written: a resolved C4-H study card and difference card beside an evidence card recording "
     "a trap response. Both render; the page's own words add no 'because', no 'explains'.",
     RENDER, "RENDERED", [NO_BRIDGE, STUDY_SPLIT])
case("plan_text_is_verbatim", "H6 no bridge",
     "A plan task carrying 'because' (authored text): shown verbatim in a span, not linted, and not the "
     "page's own words.", RENDER, "RENDERED",
     [must(r"because the ester is a weaker donor", "MISLEADING", "the task, verbatim"), NO_BRIDGE],
     plan_patch=[(["acts", 0, "tasks", "+"], "Say whether the acetate looks less negative because the ester "
                                            "is a weaker donor.")])

# --- probe on #21 (fe2fe08): its findings, as cases -------------------------------------------
# H7: plan text is authored text: a field that is not text is refused, never shown as str(value)
# (the user's "a validator, not a sanitizer").
for vid, path, value in (("task_not_text", ["acts", 0, "tasks", "+"], 7),
                         ("predict_not_text", ["acts", 1, "predict"], ["Does changing E to Z", "matter?"]),
                         ("title_not_text", ["title"], 42),
                         ("act_title_not_text", ["acts", 2, "title"], 3.5)):
    case(f"plan_{vid}", "H7 plan text",
         f"The plan's {path[-1] if path[-1] != '+' else 'task'} given as {type(value).__name__} {value!r}: not text. "
         "Refuse, naming where.", REFUSE, "REFUSED", plan_patch=[(path, value)])
case("plan_expected_not_text", "H7 plan text",
     "An expected answer given as a mapping, on the instructor's page: not text. Refuse.", REFUSE, "REFUSED",
     plan_patch=[(["acts", 0, "expected"], {"answer": "resolved"})], audience="instructor")
# the rerun on #21 (09097e4): an optional text field of None reads as absent (as #13's context
# entries do) and prints nothing; it is never the text "None"
for vid, path, aud in ((("act_title_none", ["acts", 0, "title"], "student")),
                       (("page_title_none", ["title"], "student")),
                       (("expected_none", ["acts", 1, "expected"], "instructor"))):
    case(f"plan_{vid}", "H7 plan text",
         f"The plan's {'/'.join(map(str, path))} given as None (synthetic): read as absent. The page renders "
         "without it and never prints 'None'.", RENDER, "RENDERED",
         # case-sensitive: the class panel's own "none is weighed against another" is not this
         [must_not(r"(?-i:\bNone\b)", "MISLEADING", "a missing field is not the text 'None'"), ACT_ORDER],
         plan_patch=[(path, None)], audience=aud)
# the rerun on #21 (318da98): text with no visible character (spaces, tabs, newlines, NBSP,
# zero-width) reads like "": an optional field is absent; a task or a prediction with nothing to
# read is refused. Never an empty heading, an "Act A:" with nothing after it, or an empty bullet.
BLANKS = (("spaces", "   "), ("tab_newline", "\n\t "), ("nbsp", " "), ("zwsp", "​"))
for bid, blank in BLANKS:
    for vid, path, aud in (("page_title", ["title"], "student"),
                           ("act_title", ["acts", 0, "title"], "student"),
                           ("expected", ["acts", 2, "expected"], "instructor")):
        case(f"plan_{vid}_blank_{bid}", "H7 plan text",
             f"The plan's {'/'.join(map(str, path))} is {blank!r}, no visible character: read as absent, like \"\". "
             "No empty heading, no 'Act A:' with nothing after it, no empty instructor box.", RENDER, "RENDERED",
             [check("no_blank_text", "MISLEADING", "every shown piece of plan text has something to read"),
              ACT_ORDER], plan_patch=[(path, blank)], audience=aud)
for bid, blank in (("empty", ""),) + BLANKS:
    case(f"plan_task_blank_{bid}", "H7 plan text",
         f"A task {blank!r}: nothing to do, and an empty bullet reads as a task left out. Refuse, naming where.",
         REFUSE, "REFUSED", plan_patch=[(["acts", 0, "tasks", "+"], blank)])
case("plan_predict_blank_zwsp", "H7 plan text",
     "A prediction that is only a zero-width space: no prediction, so the reveal would come first. Refuse.",
     REFUSE_F, "REFUSED", plan_patch=[(["acts", 1, "predict"], "​")])
# the rerun on #21 (2a67bdd): "something to read" is a letter, number, punctuation mark or symbol
# that is not default-ignorable and not the braille blank. Fillers, variation selectors, a lone
# combining mark and a lone private-use character have nothing to read.
INVISIBLE = (("cgj", "͏"), ("vs16", "️"), ("vs17", "\U000e0100"), ("mongolian_fvs1", "᠋"),
             ("khmer_aq", "឴"), ("hangul_filler", "ㅤ"), ("hangul_choseong_filler", "ᅟ"),
             ("hangul_jungseong_filler", "ᅠ"), ("halfwidth_hangul_filler", "ﾠ"),
             ("braille_blank", "⠀"), ("mix", "ㅤ⠀​  "),
             ("combining_acute", "́"), ("private_use", "\U000f0000"))
for bid, blank in INVISIBLE:
    case(f"plan_predict_invisible_{bid}", "H7 plan text",
         f"A prediction that is only {blank!r}: nothing to read, so the reveal would come first. Refuse.",
         REFUSE_F, "REFUSED", plan_patch=[(["acts", 1, "predict"], blank)])
for bid, blank in INVISIBLE[5:]:  # the fillers, the braille blank, the mix, a lone mark, a lone private-use
    for vid, path, aud in (("page_title", ["title"], "student"),
                           ("act_title", ["acts", 0, "title"], "student"),
                           ("expected", ["acts", 2, "expected"], "instructor")):
        case(f"plan_{vid}_invisible_{bid}", "H7 plan text",
             f"The plan's {'/'.join(map(str, path))} is {blank!r}, nothing to read: absent.", RENDER, "RENDERED",
             [check("no_blank_text", "MISLEADING", "every shown piece of plan text has something to read"),
              ACT_ORDER], plan_patch=[(path, blank)], audience=aud)
    case(f"plan_task_invisible_{bid}", "H7 plan text",
         f"A task {blank!r}: an empty bullet. Refuse, naming where.",
         REFUSE, "REFUSED", plan_patch=[(["acts", 0, "tasks", "+"], blank)])
# text that cannot be written as UTF-8 (a lone surrogate, which json.loads makes from "\ud800"), or
# that carries a control character other than tab, newline and carriage return, is not text
for vid, path, value in (("title_lone_surrogate", ["title"], "\ud800"),
                         ("title_inner_surrogate", ["title"], "A\ud800B"),
                         ("predict_inner_surrogate", ["acts", 1, "predict"], "Which is \udfff more negative?"),
                         ("title_nul", ["title"], "Homework\x00one"),
                         ("title_esc", ["title"], "Homework\x1b[31m one"),
                         ("task_vt", ["acts", 0, "tasks", "+"], "Compare\x0bthe two maps.")):
    case(f"plan_{vid}", "H7 plan text",
         f"The plan's {path[-1] if path[-1] != '+' else 'task'} is {value!r}: not writable text. Refuse, naming where.",
         REFUSE, "REFUSED", plan_patch=[(path, value)])
# the rerun on #21 (9eae791): a noncharacter (U+FDD0-FDEF, U+xFFFE/xFFFF) is an HTML parse error,
# so text holding one is not text a page can carry
for vid, path, value in (("title_nonchar_fffe", ["title"], "Real￾text"),
                         ("title_nonchar_fdd0", ["title"], "Real﷐text"),
                         ("predict_nonchar_10ffff", ["acts", 1, "predict"], "Which is\U0010ffff more negative?")):
    case(f"plan_{vid}", "H7 plan text",
         f"The plan's {path[-1]} is {value!r}: a noncharacter inside text. Refuse, naming where.",
         REFUSE, "REFUSED", plan_patch=[(path, value)])
# U+1D159 MUSICAL SYMBOL NULL NOTEHEAD is a symbol drawn as nothing, by design: nothing to read
for bid, blank in (("null_notehead", "\U0001d159"), ("null_notehead_mix", "\U0001d159ㅤ⠀")):
    case(f"plan_predict_invisible_{bid}", "H7 plan text",
         f"A prediction that is only {blank!r}: nothing to read. Refuse.",
         REFUSE_F, "REFUSED", plan_patch=[(["acts", 1, "predict"], blank)])
    case(f"plan_page_title_invisible_{bid}", "H7 plan text",
         f"The page title is {blank!r}, nothing to read: absent.", RENDER, "RENDERED",
         [check("no_blank_text", "MISLEADING", "every shown piece of plan text has something to read"), ACT_ORDER],
         plan_patch=[(["title"], blank)])
    case(f"plan_task_invisible_{bid}", "H7 plan text", f"A task {blank!r}: an empty bullet. Refuse.",
         REFUSE, "REFUSED", plan_patch=[(["acts", 0, "tasks", "+"], blank)])
# plan.case_id names the page (its <title>): plan text like any other. None, "" or nothing to read
# is absent (the page is named "homework"); not text, or not writable, is refused
for vid, value in (("none", None), ("empty", ""), ("filler", "ㅤ")):
    case(f"plan_case_id_{vid}", "H7 plan text",
         f"The plan's case_id is {value!r}: absent, so the page is named 'homework', never 'None' or blank.",
         RENDER, "RENDERED",
         [check("page_name", "MISLEADING", "the page's <title> has something to read and is never 'None'"), ACT_ORDER],
         plan_patch=[(["case_id"], value)])
for vid, value in (("not_text", 3), ("lone_surrogate", "\ud800"), ("nul", "Case\x00id"), ("nel", "Case\x85id"),
                   ("nonchar", "Case￾id")):
    case(f"plan_case_id_{vid}", "H7 plan text",
         f"The plan's case_id is {value!r}: not text a page can carry. Refuse, naming where.",
         REFUSE, "REFUSED", plan_patch=[(["case_id"], value)])

# --- the round-7 probe on #21 (2b87a0a): every path from plan or source text to the page --------
# H8 source text. A molecule's slug names it in the role legend, in data-molecules (split on ",")
# and in the views: text with something to read, carriable, and no comma. Renamed everywhere it
# is named (plan, ensembles, conformers, the study), so only the slug itself is under test.
case("slug_renamed_plainly", "H8 source text",
     "The acetate's slug renamed 'oac_renamed' everywhere: the harness is consistent, so the page renders.",
     RENDER, "RENDERED", [ROLES_SHOWN, ACT_ORDER], source_ops=[["rename_molecule", OAC, "oac_renamed"]])
for vid, value in (("lone_surrogate", "x\ud800"), ("nul", "x\x00y"), ("nel", "x\x85y"), ("nonchar", "x￾y"),
                   ("filler", "ㅤ"), ("empty", ""), ("comma", "a,b")):
    case(f"slug_{vid}", "H8 source text",
         f"The acetate's slug is {value!r} everywhere: no name a reader can read, not carriable, or two names in "
         "data-molecules for one molecule. Refuse, naming the slug.", REFUSE, "REFUSED",
         source_ops=[["rename_molecule", OAC, value]])
# a difference card at a site the molecules do not have is refused, not a crash
for vid, value in (("unknown", "zz_no_such_site"), ("none", None), ("not_text", 31337)):
    case(f"difference_card_site_{vid}", "H8 source text",
         f"A difference card at site {value!r}, which the acetate and alcohol records do not carry. Refuse.",
         REFUSE, "REFUSED", plan_patch=[(["acts", 0, "blocks", 2, "site"], value)])
# the study's plan reference, a record's espmap version and the analysis key are shown on the page
for vid, value in (("none", None), ("not_text", 31337), ("nul", "analysis-plan.md @ 7e34658\x00"),
                   ("nonchar", "analysis-plan.md￾ @ 7e34658"), ("lone_surrogate", "analysis-plan.md\ud800 @ 7e34658"),
                   ("filler", "ㅤ")):
    case(f"study_plan_{vid}", "H8 source text",
         f"The study's plan reference is {value!r}: no plan a reader can read, or not carriable. Refuse.",
         REFUSE, "REFUSED", source_ops=[["study_plan", value]])
for vid, value in (("not_text", 31337), ("nul", "0.3.0\x00"), ("c1", "0.3.0\x9b"), ("lone_surrogate", "0.3.0\ud800")):
    case(f"espmap_version_{vid}", "H8 source text",
         f"One acetate record's espmap version is {value!r}: not text a page can carry. Refuse.",
         REFUSE, "REFUSED", source_ops=[["espmap_version", OAC, value]])
for vid, value in (("not_text", 31337), ("nul_late", "d0f5ed13" + "\x00" + "6cda12ada8c6324a79ec9f7e79134d2f"),
                   ("nul_early", "d0f5\x00d136cda12ada8c6324a79ec9f7e79134d2f"),
                   ("lone_surrogate", "d0f5\ud800d136cda12ada8c6324a79ec9f7e79134d2f")):
    case(f"analysis_key_{vid}", "H8 source text",
         f"Every record's analysis key is {value!r}: not a key a page can carry. Refuse.",
         REFUSE, "REFUSED", source_ops=[["analysis_key", value]])
for vid, value in (("nul", "B\x00stereo"), ("lone_surrogate", "B\ud800")):
    case(f"contrast_id_{vid}", "H8 source text",
         f"Act B's contrast renamed {value!r} in plan and study: not text a page can carry. Refuse.",
         REFUSE, "REFUSED", source_ops=[["rename_contrast", "B_stereo", value]])

# H4: a plan reference naming two commits is ambiguous: never resolved by position
case("plan_reference_two_commits", "H4 provenance",
     "The study's plan reference names two commits ('@ 7e34658 @ deadbee'): the page cannot tell which froze "
     "the plan. It renders with the plan commit marked not recorded, never 'Plan frozen at' either one.",
     {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}, "RENDERED",
     [must_not(r"Plan frozen at", "MISLEADING", "no commit picked by position"),
      must(r"provenance incomplete", "MISLEADING", "the gap, said")], source_ops=["plan_two_commits"])
case("plan_reference_no_commit", "H4 provenance",
     "A plan reference naming no commit: the plan commit is marked not recorded.",
     {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}, "RENDERED",
     [must_not(r"Plan frozen at", "MISLEADING", "no commit invented"),
      must(r"provenance incomplete", "MISLEADING", "the gap, said")], source_ops=["plan_no_commit"])

CHECKS = ("block_roles", "act_order", "footer", "study_words_in_study_cards", "no_bridge", "no_blank_text", "page_name")


def document():
    ids = [c["id"] for c in CASES]
    assert len(ids) == len(set(ids))
    for c in CASES:
        assert set(c["outcomes"].values()) <= set(GRADES), c["id"]
        for q in c["requires"]:
            assert sum(k in q for k in ("must", "must_not", "check")) == 1 and q["grade"] in GRADES
            assert "check" not in q or q["check"] in CHECKS
        assert c["audience"] in ("student", "instructor")
    s = json.dumps({"generator": "tests/pedagogy/make_homework_page_set.py", "grades": list(GRADES),
                    "drop": DROP, "plan": PLAN, "cases": CASES}, indent=1, ensure_ascii=False) + "\n"
    # a lone surrogate cannot be written as UTF-8: keep it as a JSON escape, which json.load turns
    # back into the lone surrogate, the way a plan read from JSON would carry one
    return re.sub("[\ud800-\udfff]", lambda m: f"\\u{ord(m.group()):04x}", s)


if __name__ == "__main__":
    OUT.write_text(document(), encoding="utf-8")
    print(f"{len(CASES)} cases -> {OUT}")
