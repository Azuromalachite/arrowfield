"""Generate the hostile set for the study-verdict card:
tests/pedagogy/study_card_set.json.

    python tests/pedagogy/make_study_card_set.py     # byte-identical on re-run

Written before the card exists (expectations first). The card is a view
that shows ONE study row, the study's own verdict at one site of one
contrast, with what the study allows a reader to say:

    views.study_card_html(row, context) -> HTML, or raises model.Refused

- row: a frozen study's row, as study.py reads it (StudySite fields, given
  here as a dict). Arrowfield never re-decides it.
- context: what the card needs that the row does not carry.
  - provenance: results file sha256 and the study's frozen-plan commit;
  - molecules: role per slug, atlas_observed or teaching_control;
  - limitations: basis sensitivity run or not, solvent studied or not;
  - the all-member row, when the row shown is the extended-only sensitivity.

Each case grades the card's two possible outputs (RENDERED, REFUSED; a crash
is ERROR, graded CRASH) on the comparison set's scale: SAFE, SAFE_REFUSAL,
CRASH, MISLEADING, FALSE. `allowed` is exactly the SAFE and SAFE_REFUSAL
ones. A RENDERED card must also meet every entry of `requires`:
- must: a pattern the card's text has to contain;
- must_not: a pattern it must not contain.
Each requirement carries the grade a violation earns. Patterns are
case-insensitive regular expressions on the card's text with tags removed
(on: "text"), or on its HTML (on: "html") where styling is the point.
The 12 states come from the cross-repo review of 2026-10-24; each has two or
more variants, so a card cannot pass by special-casing one row.
"""
import json
from pathlib import Path

OUT = Path(__file__).with_name("study_card_set.json")
RULE = "espmap-sesiidae-001/range_partition_sign_v1"
GRADES = ("SAFE", "SAFE_REFUSAL", "CRASH", "MISLEADING", "FALSE")
C_ALL = {"1_ranges_disjoint": True, "2_beyond_partition": True, "3_same_sign": True}
C_AVG = {"1_ranges_disjoint": False, "2_beyond_partition": True, "3_same_sign": True}
C_NONE = {"1_ranges_disjoint": False, "2_beyond_partition": False, "3_same_sign": True}


def row(contrast, x, y, site, outcome, cond, delta, rx, ry, part, interval, tier="primary",
        population="all", rule=RULE, n=(20, 20)):
    return dict(verdict_rule=rule, contrast=contrast, x=x, y=y, population=population, tier=tier,
                site=site, outcome=outcome, n=list(n), delta=delta, delta_other_rule=delta,
                partition=part, partition_max_member=part, range_x=list(rx), range_y=list(ry),
                mean_x=None, mean_y=None, bootstrap_interval=list(interval), conditions=dict(cond))


OBSERVED = {"z3_z13_18_oac": "atlas_observed", "e3_z13_18_oac": "atlas_observed",
            "e2_z13_18_oac": "atlas_observed", "z3_z13_18_oh": "atlas_observed",
            "z3_z13_18_al": "teaching_control"}
PROV = {"results_sha256": "f" * 64, "plan_commit": "7e3465859bc2b8950078238f0bffb0ea6b4dd512"}
LIMITS = {"basis_sensitivity": "not run", "solvent": "not studied", "phase": "gas",
          "method": "B3LYP/def2-SVP"}


def ctx(**kw):
    c = {"provenance": dict(PROV), "molecules": dict(OBSERVED), "limitations": dict(LIMITS)}
    c.update(kw)
    return c


# the same rows the real study reports (Atlas f286bc1), rounded
C4 = row("C_position", "e3_z13_18_oac", "e2_z13_18_oac", "C4-H", "resolved", C_ALL, 5.22,
         (0.22, 5.54), (6.66, 9.53), 0.69, (4.57, 5.88))
VMAX = row("A_headgroup", "z3_z13_18_oac", "z3_z13_18_oh", "global V_max", "resolved", C_ALL, 26.50,
           (19.57, 20.98), (34.00, 51.68), 0.0, (24.02, 28.80))
B3 = row("B_stereo", "z3_z13_18_oac", "e3_z13_18_oac", "C3-H", "different on average", C_AVG, 2.31,
         (-2.78, 5.02), (1.04, 6.60), 0.93, (1.20, 3.42))
C2 = row("C_position", "e3_z13_18_oac", "e2_z13_18_oac", "C2-H", "different on average", C_AVG, -2.50,
         (3.30, 10.06), (-0.67, 6.89), 0.94, (-3.89, -1.17))
B4 = row("B_stereo", "z3_z13_18_oac", "e3_z13_18_oac", "C4-H", "unresolved", C_NONE, -0.04,
         (0.27, 5.08), (0.22, 5.54), 0.69, (-0.84, 0.74))
C3 = row("C_position", "e3_z13_18_oac", "e2_z13_18_oac", "C3-H", "unresolved",
         {"1_ranges_disjoint": False, "2_beyond_partition": False, "3_same_sign": False}, -0.35,
         (1.04, 6.60), (-1.79, 8.80), 0.93, (-1.95, 1.27))

CASES = []


def case(id, cls, why, r, c, outcomes, prefer, requires=(), extra=None):
    outcomes = dict(outcomes)
    outcomes.setdefault("ERROR", "CRASH")
    allowed = sorted(k for k, g in outcomes.items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert prefer is None or prefer in allowed, id
    CASES.append(dict(id=id, cls=cls, component="study_card", why=why, row=r, context=c,
                      outcomes=outcomes, allowed=allowed, prefer=prefer, requires=list(requires),
                      **(extra or {})))


def must(pat, grade, why, on="text"):
    return {"must": pat, "grade": grade, "why": why, "on": on}


def must_not(pat, grade, why, on="text"):
    return {"must_not": pat, "grade": grade, "why": why, "on": on}


SHOW = {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}      # a refusal is safe but hides a valid result
RULE_NAMED = must(r"range_partition_sign_v1", "MISLEADING",
                  "every verdict names the rule that decided it")
NO_BIOLOGY = must_not(r"biologically|significan|explain|caus|attract|receptor|specificity|active\b|prove",
                      "FALSE", "a free-molecule ESP verdict says nothing about biology")
NO_CI = must_not(r"confidence|\bci\b|95 ?%|p[- ]?value|\bp ?[<=]", "FALSE",
                 "the bootstrap interval describes the selected conformers, not a population")
BOOT = must(r"bootstrap.{0,60}selected conformers|selected conformers.{0,60}bootstrap", "MISLEADING",
            "the interval says what it resamples")

# 1. resolved: under this study and protocol, never biology
for vid, r in (("C4H", C4), ("Vmax", VMAX)):
    case(f"resolved_{vid}", "S1 resolved", "A resolved verdict is a statement about this protocol's "
         "free-molecule ESP: the card says so, names the rule, and claims nothing biological.",
         r, ctx(), SHOW, "RENDERED",
         [RULE_NAMED, NO_BIOLOGY, must(r"under (this|the) (study|protocol)", "MISLEADING",
                                      "scope words stay beside the verdict"),
          must(r"resolved", "FALSE", "the study's own word is shown")])

# 2. different on average: never success semantics
for vid, r in (("B_C3H", B3), ("C_C2H", C2)):
    case(f"different_on_average_{vid}", "S2 different on average",
         "Means differ, the observed conformer ranges overlap. The card must not look like a success "
         "(no 'resolved' except 'not resolved', no check-mark verdict, no success styling) and must "
         "say the ranges overlap.",
         r, ctx(), SHOW, "RENDERED",
         [RULE_NAMED, NO_BIOLOGY,
          must(r"different on average", "FALSE", "the study's own words"),
          must(r"ranges overlap", "MISLEADING", "why it is not resolved"),
          must_not(r"(?<!not )\bresolved\b", "MISLEADING", "only 'not resolved' may appear"),
          must_not(r"✔|✅", "MISLEADING", "no check-mark verdict"),
          must_not(r"class=\"[^\"]*(verdict-resolved|verdict-success)", "MISLEADING",
                   "no success styling on the verdict", on="html")])

# 3. unresolved: never "no difference"
for vid, r in (("B_C4H", B4), ("C_C3H", C3)):
    case(f"unresolved_{vid}", "S3 unresolved",
         "Unresolved means this protocol does not resolve the two at this site, never that they are "
         "the same.",
         r, ctx(), SHOW, "RENDERED",
         [RULE_NAMED, must(r"not resolved|unresolved", "FALSE", "the study's word"),
          must_not(r"no difference|identical|the same|equivalent|no effect", "FALSE",
                   "absence of resolution is not evidence of sameness")])

# 4. not evaluable: no verdict, whatever numbers are present
for vid, why in (("too_many_excluded", "3 of 20 members excluded for one molecule"),
                 ("excluded_folded", "an excluded member is folded")):
    r = dict(B3, outcome="not evaluable")
    case(f"not_evaluable_{vid}", "S4 not evaluable",
         f"The plan makes the contrast not evaluable ({why}). The row still carries numbers; the card "
         "must not present them as a finding.",
         r, ctx(), SHOW, "RENDERED",
         [must(r"not evaluable", "FALSE", "the study's word"),
          must_not(r"Δ|difference of (the )?means", "MISLEADING", "no difference value shown as a result"),
          must_not(r"(?<!not )\bresolved\b|different on average", "FALSE", "no verdict word")])

# 5. exploratory: visibly not a prespecified test
for vid, out, cond in (("resolved", "resolved", C_ALL), ("unresolved", "unresolved", C_NONE)):
    r = row("B_stereo", "z3_z13_18_oac", "e3_z13_18_oac", "C11-H", out, cond, 1.4 if out == "resolved" else 0.2,
            (0.1, 1.0), (1.5, 2.6) if out == "resolved" else (0.0, 2.0), 0.3, (1.1, 1.7) if out == "resolved"
            else (-0.4, 0.8), tier="exploratory")
    case(f"exploratory_{vid}", "S5 exploratory",
         "A non-primary chain site meeting (or not) the rule is exploratory: never interpreted, never "
         "presented like a prespecified result (the pilot's C11-H was a chance hit at K = 5).",
         r, ctx(), SHOW, "RENDERED",
         [must(r"exploratory", "MISLEADING", "labelled"),
          must(r"not (a )?prespecified|never interpreted", "MISLEADING", "and what that means")])

# 6. extended-only sensitivity: never replaces the all-member result
for vid, with_all in (("alone", False), ("with_primary", True)):
    r = dict(C4, population="extended", n=[17, 14])
    c = ctx(all_member_row=C4) if with_all else ctx()
    case(f"extended_only_{vid}", "S6 extended-only sensitivity",
         "The extended-only row is a sensitivity analysis, shown beside the all-member result and "
         "never instead of it. Without the all-member row the card cannot show it honestly.",
         r, c, ({"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"} if with_all else
                {"RENDERED": "MISLEADING", "REFUSED": "SAFE_REFUSAL"}),
         "RENDERED" if with_all else "REFUSED",
         [must(r"sensitivity", "MISLEADING", "labelled as a sensitivity"),
          must(r"extended", "MISLEADING", "which population"),
          must(r"all (members|conformers)", "MISLEADING", "the primary result it sits beside")])

# 7. bootstrap interval: never a confidence interval
for vid, r in (("on_resolved", C4), ("on_different_on_average", B3)):
    case(f"bootstrap_named_{vid}", "S7 bootstrap interval",
         "The interval resamples the 20 selected conformers (Atlas 30c4a11): not a confidence "
         "interval, not a p-value.",
         r, ctx(), SHOW, "RENDERED", [BOOT, NO_CI])

# 8. teaching control: constructed, no biological claim
for vid, r in (("Al_vs_OAc_Vmax", row("A_headgroup (control: Al vs OAc)", "z3_z13_18_oac", "z3_z13_18_al",
                                        "global V_max", "resolved", C_ALL, 5.76, (19.57, 20.98),
                                        (24.89, 26.85), 0.0, (5.47, 6.05))),
               ("Al_vs_OH_headO", row("A_headgroup (control: Al vs OH)", "z3_z13_18_oh", "z3_z13_18_al",
                                      "head O", "unresolved",
                                      {"1_ranges_disjoint": True, "2_beyond_partition": False,
                                       "3_same_sign": False}, 4.24, (-31.95, -29.03), (-26.79, -25.81),
                                      5.97, (3.81, 4.65), tier="secondary"))):
    case(f"teaching_control_{vid}", "S8 teaching control",
         "The aldehyde is a constructed teaching control, not an Atlas compound: descriptive chemistry, "
         "never 'observed', never a pheromone.",
         r, ctx(), SHOW, "RENDERED",
         [must(r"teaching control", "FALSE", "labelled as constructed"),
          must_not(r"observed|pheromone|natural|attract", "FALSE", "no biological status")])

# 9. unknown verdict rule: refuse
for vid, rule in (("generic_compare", "espmap.compare/1"), ("missing", None)):
    case(f"unknown_rule_{vid}", "S9 unknown verdict rule",
         "A verdict under a rule Arrowfield does not know (espmap's generic resolved_under), or under "
         "none, cannot be shown as this study's.",
         dict(C4, verdict_rule=rule), ctx(), {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"}, "REFUSED")

# 10. internally inconsistent: refuse
for vid, r in (("resolved_ranges_overlap", dict(C4, conditions=dict(C_ALL, **{"1_ranges_disjoint": False}))),
               ("average_interval_spans_0", dict(B3, bootstrap_interval=[-0.5, 3.4]))):
    case(f"inconsistent_{vid}", "S10 internally inconsistent",
         "The row's outcome does not follow from its own conditions under its rule: corrupt, refused, "
         "never re-decided.",
         r, ctx(), {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"}, "REFUSED")

# 11. limitations stay visible
for vid, r in (("resolved", C4), ("different_on_average", B3)):
    case(f"limitations_{vid}", "S11 limitations",
         "Basis sensitivity was not run and solvent not studied: a verdict on this protocol carries "
         "both, with the scope (free molecule, gas phase, method).",
         r, ctx(), SHOW, "RENDERED",
         [must(r"basis.{0,40}not run", "MISLEADING", "unmeasured basis sensitivity named"),
          must(r"solvent.{0,40}not studied", "MISLEADING", "solvent named"),
          must(r"free molecule", "MISLEADING", "scope"), must(r"gas", "MISLEADING", "phase"),
          must(r"def2-svp", "MISLEADING", "method")])
case("limitations_basis_run", "S11 limitations",
     "When the context says basis sensitivity was run, the card must not claim it was not.",
     C4, ctx(limitations=dict(LIMITS, basis_sensitivity="run")), SHOW, "RENDERED",
     [must_not(r"basis.{0,40}not run", "FALSE", "do not understate what was done")])

# 12. provenance
case("provenance_missing_rule_and_sha", "S12 provenance",
     "No results sha256 and no verdict rule: nothing ties the card to a study. Refuse.",
     dict(C4, verdict_rule=None), ctx(provenance={}), {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"},
     "REFUSED")
case("provenance_missing_sha", "S12 provenance",
     "The rule is known but the results file's sha256 is missing: the card may show the verdict with "
     "a visible provenance warning, or refuse; it may not look fully sourced.",
     C4, ctx(provenance={"plan_commit": PROV["plan_commit"]}), SHOW, "RENDERED",
     [must(r"provenance (incomplete|missing)|unverified source|source not (recorded|verified)",
           "MISLEADING", "a visible warning")])
case("provenance_shown", "S12 provenance",
     "With full provenance, the card names the rule and the frozen plan's commit.",
     C4, ctx(), SHOW, "RENDERED",
     [RULE_NAMED, must(r"7e34658", "MISLEADING", "the frozen plan's commit")])

# 13. malformed numbers in a row that claims a verdict (held-out probe of #11)
# A reversed interval read naively 'excludes 0' and lets 'different on
# average' through; a NaN difference renders as '+nan'. A verdict-bearing row
# must carry finite numbers and ordered intervals, or be refused.
NAN, INF = float("nan"), float("inf")
REFUSE = {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"}
for vid, r, grade_r in (
        ("interval_reversed", dict(B3, bootstrap_interval=[2.2, -0.6]), "FALSE"),
        ("interval_nan", dict(B3, bootstrap_interval=[NAN, 3.4]), "FALSE"),
        ("delta_nan", dict(C4, delta=NAN), "MISLEADING"),
        ("range_inf", dict(C4, range_y=[6.66, INF]), "MISLEADING"),
        ("range_reversed", dict(C4, range_x=[5.54, 0.22]), "MISLEADING"),
        ("partition_nan", dict(C4, partition=NAN), "MISLEADING")):
    case(f"malformed_{vid}", "S13 malformed numbers",
         "A row that claims a verdict must carry finite numbers and ordered (low, high) intervals; "
         "anything else is a corrupt record, refused rather than rendered as if measured.",
         r, ctx(), {"RENDERED": grade_r, "REFUSED": "SAFE_REFUSAL"}, "REFUSED")

# Probe round 2 of #11: Python's True and False are numbers to isinstance(.., int),
# so a boolean slipped through as a delta, a partition or an interval bound,
# and reached the card through the all-member row too.
for vid, r, c in (("bool_delta", dict(C4, delta=True), ctx()),
                  ("bool_partition", dict(C4, partition=False), ctx()),
                  ("bool_interval", dict(B3, bootstrap_interval=[True, 2.0]), ctx()),
                  ("bool_in_all_member_row", dict(C4, population="extended", n=[17, 14]),
                   ctx(all_member_row=dict(C4, partition=True)))):
    case(f"malformed_{vid}", "S13 malformed numbers",
         "A boolean is not a measurement: a row (or the all-member row beside it) carrying True or "
         "False where a number belongs is corrupt and refused.",
         r, c, {"RENDERED": "FALSE", "REFUSED": "SAFE_REFUSAL"}, "REFUSED")
# ...and a not-evaluable row's ranges were formatted without any check: a
# string crashed the card, an infinity was shown. The card may still show a
# not-evaluable row, but only well-formed ranges, and never crash.
for vid, rx in (("range_garbage", ["a", None]), ("range_infinite", [-INF, 1.0])):
    case(f"not_evaluable_{vid}", "S4 not evaluable",
         "A not-evaluable row carries no verdict, so stray numbers do not make it corrupt; but a "
         "range that is not a finite (low, high) pair must not be shown, and must not crash the card.",
         dict(B3, outcome="not evaluable", range_x=rx, delta=NAN, partition=-4.0,
              bootstrap_interval=[5.0, 1.0]), ctx(), SHOW, "RENDERED",
         [must(r"not evaluable", "FALSE", "the study's word"),
          must_not(r"\binf\b|infinity|\bnan\b|∞", "MISLEADING", "no infinite or NaN value shown"),
          must_not(r"Δ|difference of (the )?means", "MISLEADING", "no difference value shown")])

# Probe round 3 of #11: provenance the card cannot verify. A plan commit that
# is not a hex string (an int crashed the card) or is missing, or a results
# sha256 that is not 64 hex characters, must not look fully sourced.
WARN = must(r"provenance (incomplete|missing)|unverified source|source not (recorded|verified)",
            "MISLEADING", "a visible warning")
for vid, prov in (("plan_commit_int", dict(PROV, plan_commit=7346580)),
                  ("plan_commit_missing", dict(PROV, plan_commit=None)),
                  ("sha_not_hex", dict(PROV, results_sha256="not-a-sha"))):
    case(f"provenance_{vid}", "S12 provenance",
         "Provenance the card cannot verify: show the verdict with a visible warning (or refuse); "
         "never crash, never present it as fully sourced.",
         C4, ctx(provenance=prov), SHOW, "RENDERED", [WARN])
# ...and a verdict row whose n is not two non-negative integers.
for vid, n in (("n_garbage", ["x", -3]), ("n_bool", [True, 20])):
    case(f"malformed_{vid}", "S13 malformed numbers",
         "A verdict rests on n members per side; n that is not two non-negative integers is a "
         "corrupt record.",
         dict(C4, n=n), ctx(), {"RENDERED": "MISLEADING", "REFUSED": "SAFE_REFUSAL"}, "REFUSED")


# Probe round 4 of #11: context entries that are not mappings (a provenance
# given as a list crashed the card), and n given as numpy integers, which a
# study built from array counts produces and which are perfectly good counts.
for vid, c in (("provenance_list", ctx(provenance=["f" * 64, PROV["plan_commit"]])),
               ("molecules_list", ctx(molecules=["z3_z13_18_al"])),
               ("limitations_string", ctx(limitations="gas, def2-SVP"))):
    case(f"context_not_a_mapping_{vid}", "S12 provenance",
         "A context entry that is not a mapping is treated as missing: the card never crashes, and "
         "never claims what it could not read (a missing provenance shows the warning).",
         C4, c, SHOW, "RENDERED",
         [WARN] if vid == "provenance_list" else [must_not(r"basis.{0,40}\brun\b(?!.{0,5}not)", "MISLEADING",
                                                          "nothing claimed as run")])
case("n_numpy_integers", "S13 malformed numbers",
     "n given as numpy integers (np.int64, np.int32) is a valid count; the card must not refuse it. "
     "The harness converts n to numpy integers for this case (row_n_numpy).",
     C4, ctx(), {"RENDERED": "SAFE", "REFUSED": "SAFE_REFUSAL"}, "RENDERED", extra={"row_n_numpy": True})


def document():
    for c in CASES:
        assert set(c["outcomes"].values()) <= set(GRADES), c["id"]
        for q in c["requires"]:
            assert q["grade"] in GRADES and ("must" in q) != ("must_not" in q), c["id"]
            assert q["on"] in ("text", "html"), c["id"]
    return json.dumps({"generator": "tests/pedagogy/make_study_card_set.py", "grades": list(GRADES),
                       "component": "study_card", "cases": CASES}, indent=1, ensure_ascii=False) + "\n"


if __name__ == "__main__":
    OUT.write_text(document(), encoding="utf-8")
    print(f"{len(CASES)} cases -> {OUT}")
