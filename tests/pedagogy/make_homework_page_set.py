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
# act A task 2, "how far the head-group perturbation appears to extend": a near-to-far series of
# chain sites, OAc -> OH (both carry these sites; the fixture copies them from the records)
A_SERIES = ("C1-H", "C2-H", "C3-H", "C4-H", "C6-H", "C10-H", "C16-H")
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
                    # task 2's near-to-far series (amendment): chain sites the study does not test, so
                    # each card is exploratory and needs no study card
                    *({"view": "difference_card", "contrast": "A_headgroup", "x": OAC, "y": OH, "site": s}
                      for s in A_SERIES),
                    {"view": "contrast_card", "label": "S. tenuis gland, Japan and Korea"},
                    {"view": "class_panel", "molecule": OAC, "taxon_id": "spa.tax.synanthedon_tenuis"}]},
        {"id": "B", "contrasts": ["B_stereo"], "title": "What E/Z does (E3 vs Z3 acetate)",
         "predict": "Does changing E to Z at C3 mainly change the charge distribution, or the shape?",
         "tasks": ["Read the difference card at C3-H, C4-H and head O. Is each difference larger than that site's "
                   "conformer range? Than its partition sensitivity?",
                   "Look at the extended and folded conformers of one compound side by side, on the same scale."],
         "expected": "E/Z is not resolved at any site. C3-H: different on average.",
         # every site task 1 names (amendment): each is tested, so its study card comes with it
         "blocks": [*(b for s in ("C3-H", "C4-H", "head O")
                      for b in ({"view": "study_card", "contrast": "B_stereo", "site": s},
                                {"view": "difference_card", "contrast": "B_stereo", "x": OAC, "y": E3, "site": s})),
                    {"view": "pair_view", "left": {"slug": E3, "conformer": "conf_00"},
                     "right": {"slug": E3, "conformer": "conf_02"}}]},
        {"id": "C", "contrasts": ["C_position"], "title": "Where the double bond is (E2 vs E3 acetate)",
         "predict": "Moving C=C by one carbon. Which sites should change, and by how much compared with act B?",
         "tasks": ["Compare C2-H, C3-H and C4-H between the two compounds, with the difference card. Which "
                   "differences exceed the conformer range?",
                   "Does free-molecule ESP give a convincing reason an organism or enzyme would treat these two "
                   "compounds differently?"],
         "expected": "Resolved at C4 only. Question 2's honest answer stays no.",
         # every site task 1 names (amendment), each tested, with its study card
         "blocks": [*(b for s in ("C2-H", "C3-H", "C4-H")
                      for b in ({"view": "study_card", "contrast": "C_position", "site": s},
                                {"view": "difference_card", "contrast": "C_position", "x": E3, "y": E2, "site": s})),
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
# the fresh probe on #21 (565e303): every study card printed "solvent: not stated" beside "gas
# phase", though the study's protocol is the gas phase and every record says so
STUDY_SCOPE = [must_not(r"solvent: not stated", "MISLEADING", "the study's protocol states its phase"),
               must(r"solvent: [^.]*gas phase", "MISLEADING", "the study's solvent: none, the gas phase"),
               # the fresh probe on e17af42: "basis sensitivity: not stated", though the study states its
               # def2-TZVP check was not run (plan section 6; results.md)
               must_not(r"basis sensitivity: not stated", "MISLEADING", "the study states its basis check"),
               must(r"basis sensitivity: not run", "MISLEADING", "the study's basis check: not run")]


def blk(act, k):
    """Path to block k of act index act."""
    return ["acts", act, "blocks", k]


def blk_of(act, view):
    """Path to the one block of a kind in act index act, wherever the plan puts it."""
    ks = [k for k, b in enumerate(PLAN["acts"][act]["blocks"]) if b["view"] == view]
    assert len(ks) == 1, (act, view, ks)
    return blk(act, ks[0])


# --- H1 the real plan --------------------------------------------------------------------------
case("real_plan_student", "H1 the real plan",
     "The homework as written, for a student: renders, with every rule below holding at once, and no "
     "expected answer on the page.",
     RENDER, "RENDERED", SOUND + STUDY_SCOPE + [
         must(r"teaching control", "FALSE", "the aldehyde named as a teaching control"),
         must_not(r"ESP works|honest answer stays no|for the instructor", "MISLEADING",
                  "no expected answers on a student's page")])
case("real_plan_instructor", "H1 the real plan",
     "The same plan for the instructor: the expected answers appear, labelled for the instructor.",
     RENDER, "RENDERED", SOUND + STUDY_SCOPE + [must(r"for the instructor", "MISLEADING", "labelled"),
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
     REFUSE, "REFUSED", plan_patch=[(blk_of(1, "pair_view") + ["right", "conformer"], "conf_99")])

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

# --- the round-8 probe on #21 (23d5ba2) ------------------------------------------------------------
# H9 display. A bidi override, embedding or isolate (U+202A-202E, U+2066-2069) in any text reorders
# the page's own words after it (nothing on the page isolates bidi): not text a page carries.
for vid, op in (("study_plan", ["study_plan", "‮AB @ 7e34658"]),
                ("espmap_version", ["espmap_version", OAC, "‮0.3.0"]),
                ("contrast_id", ["rename_contrast", "B_stereo", "B‮stereo"]),
                ("shorthand", ["bundle_shorthand", "spa.cmpd.z3_z13_18_oac", "‮AB"])):
    case(f"bidi_override_{vid}", "H9 display",
         f"An unterminated RLO in the {vid.replace('_', ' ')}: the page's own words after it read reversed. Refuse.",
         REFUSE_F, "REFUSED", source_ops=[op])
case("bidi_override_task", "H9 display", "An unterminated RLO in a task: the page's words after it reverse. Refuse.",
     REFUSE_F, "REFUSED", plan_patch=[(["acts", 0, "tasks", "+"], "‮AB")])
case("bidi_isolate_balanced_task", "H9 display",
     "A balanced RLI...PDI in a task reorders nothing outside it; refusing it is also safe (the rule is one "
     "class of characters, not a balance count).", {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}, "REFUSED",
     plan_patch=[(["acts", 0, "tasks", "+"], "⁧AB⁩")])
# identifiers have a grammar: a slug is lowercase ASCII (so no two molecules look alike), an
# analysis key is hex, an espmap version is a version. None of them holds the page's separators
# or its "not recorded" sentinel.
for vid, value in (("zwsp_lookalike", E3 + "​"), ("trailing_space", E3 + " "), ("vs_lookalike", E3 + "️"),
                   ("word_joiner", "⁠" + E3), ("cyrillic_lookalike", "e3_z13_18_оас"),
                   ("literal_none", "None")):
    case(f"slug_{vid}", "H9 display",
         f"The acetate's slug is {value!r}: it reads as another molecule's name, or as 'None'. Refuse.",
         REFUSE_F, "REFUSED", source_ops=[["rename_molecule", OAC, value]])
for vid, value in (("sentinel", "not recorded"), ("comma_space", "0.3.0, 0.3.1"), ("comma", "0.3.0,0.3.1"),
                   ("literal_none", "None")):
    case(f"espmap_version_{vid}", "H9 display",
         f"One acetate record's espmap version is {value!r}: it reads as another value, or as two. Refuse.",
         REFUSE_F, "REFUSED", source_ops=[["espmap_version", OAC, value]])
for vid, value in (("comma", "aaaaaaaa,bbbbbbbb"), ("literal_none", "None")):
    case(f"analysis_key_{vid}", "H9 display",
         f"Every record's analysis key is {value!r}: two keys on a page that has one, or 'None'. Refuse.",
         REFUSE_F, "REFUSED", source_ops=[["analysis_key", value]])
# shown source text with nothing to read, or the literal "None" (which reads as the None-as-absent
# defect), is refused wherever it reaches the page
for vid, op in (("study_plan_literal_none", ["study_plan", "None"]),
                ("shorthand_literal_none", ["bundle_shorthand", "spa.cmpd.z3_z13_18_oac", "None"]),
                ("shorthand_filler", ["bundle_shorthand", "spa.cmpd.z3_z13_18_oac", "ㅤ"])):
    case(vid, "H9 display", f"{op[0]} {op[-1]!r}: shown, but nothing a reader can tell from a missing value. Refuse.",
         REFUSE, "REFUSED", source_ops=[op])

# --- round 9 on #21 (8fe0749): the record settings, shown outside verbatim spans ------------------
# A record's xc and basis appear in the pair view's caption and the study card's scope line; its
# solvent decides "gas phase" or "PCM <solvent>". espmap always writes xc and basis, and solvent as
# None or a lowercase name; anything else came from another producer or a hand edit. Refuse it
# rather than read it: never "b3lyp/None", never a gas-phase run called condensed.
for vid, ops in (("basis_none", [["settings", "basis", None]]),
                 ("xc_none", [["settings", "xc", None]]),
                 ("both_none", [["settings", "xc", None], ["settings", "basis", None]]),
                 ("basis_literal_none", [["settings", "basis", "None"]]),
                 ("xc_literal_none", [["settings", "xc", "None"]]),
                 ("xc_filler", [["settings", "xc", "ㅤ"]]),
                 ("xc_empty", [["settings", "xc", ""]])):
    case(f"record_settings_{vid}", "H9 display",
         f"Every record's settings: {ops}. The method would read with 'None' or nothing in it. Refuse.",
         REFUSE_F, "REFUSED", source_ops=ops)
for vid, value in (("Gas", "Gas"), ("GAS", "GAS"), ("gas_space", "gas "), ("empty", ""), ("literal_none", "None"),
                   ("Vacuum", "Vacuum")):
    case(f"record_solvent_{vid}", "H9 display",
         f"Every record's solvent is {value!r}: plainly gas phase, but not a token the page knows, so it would "
         "read as 'PCM ...', condensed phase. Refuse, never interpret.", REFUSE_F, "REFUSED",
         source_ops=[["settings_ensembles", "xc", "b3lyp"], ["settings_ensembles", "basis", "def2-svp"]] + [["settings", "solvent", value]])
case("record_solvent_water", "H11 the study's protocol",
     "Every record computed in PCM water (canonical protocol and settings alike, re-keyed): not the study's "
     "protocol, so its study cards would scope the frozen verdict under water (the fresh probe on 565e303). Refuse.",
     REFUSE_F, "REFUSED",
     source_ops=[["settings_ensembles", "xc", "b3lyp"], ["settings_ensembles", "basis", "def2-svp"]] + [["canonical_rekey", "solvent", "water"], ["settings", "solvent", "water"]])

# --- round 10 on #21 (64b5fa2) -------------------------------------------------------------------
# L: a solvent is one espmap can write (esp_map.SOLVENTS, v0.2.1 and v0.3.0 alike, None for gas)
# or "gas"; a gas-phase synonym or a placeholder inside the grammar is not, and would be called
# "PCM <word>", condensed phase
for value in ("vacuo", "in_vacuo", "gas_phase", "gas-phase", "vac", "no_solvent", "isolated", "null", "unknown"):
    case(f"record_solvent_{value.replace('-', '_hyphen_')}", "H9 display",
         f"Every record's solvent is {value!r}: not a solvent espmap writes; shown, it reads as PCM {value}, "
         "condensed phase. Refuse.", REFUSE, "REFUSED", source_ops=[["settings_ensembles", "xc", "b3lyp"], ["settings_ensembles", "basis", "def2-svp"]] + [["settings", "solvent", value]])
case("record_solvent_dmso", "H11 the study's protocol",
     "Every record computed in PCM dmso (canonical protocol and settings alike, re-keyed): not the study's "
     "protocol. Refuse.", REFUSE_F, "REFUSED",
     source_ops=[["settings_ensembles", "xc", "b3lyp"], ["settings_ensembles", "basis", "def2-svp"]] + [["canonical_rekey", "solvent", "dmso"], ["settings", "solvent", "dmso"]])
case("record_solvent_gas_renders", "H9 display",
     "Every record's settings say solvent 'gas', which the views read as the gas phase, as the protocol says.",
     RENDER, "RENDERED", [must(r"gas phase", "MISLEADING", "the phase, as recorded"), ACT_ORDER],
     source_ops=[["settings_ensembles", "xc", "b3lyp"], ["settings_ensembles", "basis", "def2-svp"]] + [["settings", "solvent", "gas"]])

# --- the user's decision after round 11 on #21 (7fb0bed): method identity is espmap's -------------
# A record's protocol is its canonical key_fields.analysis, and its analysis_key the sha1 of those
# fields' sorted JSON (espmap's record_key, v0.2.1 and v0.3.0 alike). Arrowfield verifies the
# digest, compares keys, and names the method from the canonical fields; it never normalises xc or
# basis names itself. Raw settings are provenance text: spelled differently under one canonical
# protocol, they render; a raw solvent that contradicts the canonical one, a solvent with no method,
# a digest that does not match, or no canonical fields at all are refused.
SAME = [["settings_ensembles", "xc", "b3lyp"], ["settings_ensembles", "basis", "def2-svp"],
        ["settings_ensembles", "solvent", "vacuum"]]
CANON = must(r"b3lyp/def2svp", "MISLEADING", "the method, named from the canonical protocol")
case("one_method_baseline", "H10 canonical protocol",
     "Every record states b3lyp/def2-svp in the gas phase: one protocol, named once, canonically.",
     RENDER, "RENDERED", [CANON, ACT_ORDER], source_ops=SAME)
# N (round 11): one canonical protocol, spelled as each run's command line had it
for vid, ops in (("ensembles_upper", [["settings_ensembles", "xc", "B3LYP"], ["settings_ensembles", "basis", "def2-SVP"],
                                      ["settings_ensembles", "solvent", "vacuum"]]),
                 ("conformers_def2svp", SAME + [["settings_conformers", "basis", "def2svp"]]),
                 ("conformers_def2_svp_upper", SAME + [["settings_conformers", "basis", "DEF2_SVP"]]),
                 ("conformers_libxc_name", SAME + [["settings_conformers", "xc", "HYB_GGA_XC_B3LYP"]]),
                 ("one_ensemble_upper", SAME + [["settings_ensemble", OAC, "xc", "B3LYP"],
                                                ["settings_ensemble", OAC, "basis", "def2-SVP"]]),
                 # O (round 11): a PySCF custom functional is the record's own text, spaces and all
                 ("custom_functional", [["settings", "xc", "0.2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*VWN"],
                                        ["settings", "basis", "def2-svp"], ["settings", "solvent", "vacuum"]])):
    case(f"one_protocol_{vid}", "H10 canonical protocol",
         f"Records under one analysis key, its digest verified, spelling the method differently ({ops[-1]}): one "
         "protocol. Render, naming it from the canonical fields.", RENDER, "RENDERED", [CANON, ACT_ORDER],
         source_ops=ops)
# round 10's M cases: a raw xc or basis changed under one verified canonical protocol is a
# corrupted provenance label, which Arrowfield cannot check without PySCF; the protocol decides.
# The page renders and names the canonical method, never the label (the user's decision)
for vid, ops, label in (("one_ensemble_pbe0", [["settings_ensemble", OAC, "xc", "pbe0"]], r"pbe0"),
                        ("one_ensemble_tzvp", [["settings_ensemble", OAC, "basis", "def2-tzvp"]], r"tzvp"),
                        ("conformers_pbe0", [["settings_conformers", "xc", "pbe0"]], r"pbe0")):
    case(f"one_method_{vid}", "H10 canonical protocol",
         f"One verified canonical protocol, one raw label changed ({ops}): render, naming the protocol's method, "
         "never the label.", RENDER, "RENDERED",
         [CANON, must_not(label, "MISLEADING", "a raw label is not the protocol"), ACT_ORDER], source_ops=SAME + ops)
# a raw solvent that contradicts the canonical one, and P (round 11): a solvent with no method
for vid, ops in (("one_ensemble_water", SAME + [["settings_ensemble", OAC, "solvent", "water"]]),
                 ("conformers_water", SAME + [["settings_conformers", "solvent", "water"]]),
                 ("ensembles_water", SAME + [["settings_ensembles", "solvent", "water"]]),
                 ("solvent_without_method", [["settings_ensembles", "solvent", "water"]]),
                 ("canonical_water_raw_vacuum", [["canonical_rekey", "solvent", "water"]])):
    case(f"one_method_{vid}", "H10 canonical protocol",
         f"A record's raw solvent contradicts its canonical protocol, or stands with no method ({ops[-1]}): the "
         "phase the page names would be wrong for some record. Refuse.", REFUSE_F, "REFUSED", source_ops=ops)
# the digest and the fields themselves
for vid, ops in (("digest_mismatch", [["canonical", "solvent", "water"]]),
                 ("digest_mismatch_xc", [["canonical", "xc", "pbe0"]]),
                 ("no_canonical_fields", [["drop_canonical"]])):
    case(f"canonical_{vid}", "H10 canonical protocol",
         f"Every record's canonical fields {ops[0]}: the key no longer digests what the record says its protocol "
         "is, or there is none to name. Refuse.", REFUSE_F, "REFUSED", source_ops=ops)
# round 12 on #21 (b6f7f66): canonical fields that digest correctly but that espmap never writes.
# The page prints the canonical functional and basis into "xc/basis, phase": a canonical name is
# readable text with no whitespace or "/" (espmap's keys have neither, its libxc form included),
# never "None"; and the fields name their solvent (None for gas), never leave it out
for vid, field, value in (("basis_with_phase", "basis", "def2svp, PCM water"), ("xc_with_slash", "xc", "b3lyp/def2tzvp"),
                          ("xc_literal_none", "xc", "None"), ("basis_literal_none", "basis", "None"),
                          ("xc_space", "xc", " "), ("xc_filler", "xc", "ㅤ"), ("no_solvent_key", "solvent", DROP)):
    case(f"canonical_{vid}", "H10 canonical protocol",
         f"Every record's canonical {field} is {value!r} (the key re-digested): the caption would read another "
         "method or phase, or nothing. Refuse.", REFUSE_F, "REFUSED",
         source_ops=SAME + [["canonical_rekey", field, value]])
case("canonical_libxc_form", "H11 the study's protocol",
     "Every record's canonical xc is espmap's key for a functional without a display name "
     "(libxc:hyb(0.2,0.08,0)+..., re-keyed): a protocol espmap writes, but not the study's. Refuse. (The name "
     "rule itself is unit-tested in tests/test_protocol.py.)", REFUSE_F, "REFUSED",
     source_ops=SAME + [["canonical_rekey", "xc", "libxc:hyb(0.2,0.08,0)+1*0.08"]])
# round 13 on #21 (f082152): commas and phase words inside the old name grammar made the caption
# name another method or phase. espmap's keys (xc_key, basis_key; 50 common inputs checked against
# v0.2.1) put a comma only inside parentheses and never a phase or solvent word
for vid, field, value, solvent in (("xc_two_names", "xc", "b3lyp,def2tzvp", None),
                                   ("xc_two_functionals", "xc", "b3lyp,pbe0", None),
                                   ("basis_comma_pcm", "basis", "def2svp,PCM", None),
                                   ("basis_comma_pcm_water", "basis", "def2svp,PCM-water", None),
                                   ("basis_comma_gas_on_water", "basis", "def2svp,gas", "water"),
                                   ("basis_paren_gas_on_water", "basis", "def2svp(gas)", "water")):
    phase = [["canonical_rekey", "solvent", solvent], ["settings", "solvent", solvent]] if solvent else []
    case(f"canonical_{vid}", "H10 canonical protocol",
         f"Every record's canonical {field} is {value!r} (re-digested{', in ' + solvent if solvent else ''}): the "
         "caption would name a second method or another phase. Refuse.", REFUSE_F, "REFUSED",
         source_ops=SAME + phase + [["canonical_rekey", field, value]])
case("canonical_basis_with_parenthesised_comma", "H11 the study's protocol",
     "Every record's canonical basis is espmap's key for 6-311++G(2d,2p) (re-keyed): not the study's protocol. "
     "Refuse. (The name rule is unit-tested in tests/test_protocol.py.)", REFUSE_F, "REFUSED",
     source_ops=SAME + [["canonical_rekey", "basis", "6311++g(2d,2p)"]])
# the fresh probe on #21 (565e303): records rewritten as espmap writes a --solvent water, --xc pbe0
# or --basis def2-tzvp run, all under one (new) key. The page rendered, and its study cards scoped
# the frozen study's verdicts under that protocol; the study was decided under one key, d0f5ed13
# (B3LYP/def2-SVP + RI, gas phase; its plan's section 3 at 7e34658). Records under another key are
# not the study's: refuse
for vid, field, value, raw in (("water", "solvent", "water", ["settings", "solvent", "water"]),
                               ("pbe0", "xc", "pbe0", ["settings", "xc", "pbe0"]),
                               ("def2tzvp", "basis", "def2tzvp", ["settings", "basis", "def2-tzvp"])):
    case(f"study_protocol_{vid}", "H11 the study's protocol",
         f"Every record is a {value} run, written as espmap writes one (canonical and settings, one new key): the "
         "page's study cards would say the frozen study decided under it. Refuse.", REFUSE_F, "REFUSED",
         source_ops=SAME + [["canonical_rekey", field, value], raw])
# the fresh probe on #21 (e17af42): the key names a protocol, not a member set. Each difference card
# beside a study card sends the student to the study's verdict, decided on the study's members; a
# card computed on other members (10 of 20, 19 of 20) is not that evidence. A count that differs
# from the study's n_ok proves the members are not the study's: refuse (the user's decision). An
# equal count does not prove they are; the results carry no member ids. The teaching control is
# not one of the study's Atlas molecules, and its count is not checked
for vid, ops in (("first_ten_each", [["keep_members", slug, 10] for slug in (OAC, OH, E3, E2)]),
                 ("e3_nineteen", [["keep_members", E3, 19]])):
    case(f"study_members_{vid}", "H11 the study's protocol",
         f"The page's ensembles hold other member counts than the study's n_ok ({ops}): its cards are not the "
         "evidence the study cards beside them decided on. Refuse, naming the molecule and both counts.",
         REFUSE, "REFUSED", source_ops=ops)
case("study_members_control_nineteen", "H11 the study's protocol",
     "The aldehyde teaching control with 19 members: not an Atlas molecule of the study, so its count is not the "
     "study's to check. Renders.", RENDER, "RENDERED", [ACT_ORDER], source_ops=[["keep_members", AL, 19]])
# the fresh probe on #21 (9e2a397): a molecule's ensemble holding another molecule's record (one
# swapped in, or all of them) passes the key and the count, and its cards compare the wrong
# molecules beside the study's verdict. Every record carries its structure
# (sites.identity.mapped_smiles); every real ensemble holds exactly one: an ensemble of more than
# one structure, or two molecules of the page with one structure, are not the plan's molecules.
# This is chemical identity, which the records state; member identity stays a follow-up
case("ensemble_with_another_molecules_record", "H11 the study's protocol",
     "E3 acetate's 20 members include one Z3 acetate record (count and key pass; Z3's record carries every "
     "site E3's cards show): its cards mix two molecules. Refuse.", REFUSE, "REFUSED",
     source_ops=[["swap_member", E3, 19, OAC, 19]])
case("ensemble_of_another_molecule", "H11 the study's protocol",
     "E3 acetate's ensemble holds Z3 acetate's 20 records: act B compares Z3 with itself beside the study's "
     "verdict on E3. Refuse.", REFUSE, "REFUSED", source_ops=[["swap_all", E3, OAC]])
# a count is a number: True is not 1 member
case("study_n_ok_not_a_count", "H11 the study's protocol",
     "The results file gives E2's n_ok as true, and E2's ensemble has one member: true is not a count. Refuse.",
     REFUSE, "REFUSED", source_ops=[["study_n_ok", E2, True], ["keep_members", E2, 1]])

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
