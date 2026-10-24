"""Generate the hostile pedagogy set: tests/pedagogy/hostile_set.json.

    python tests/pedagogy/make_set.py        # rewrites the JSON; byte-identical on re-run

Every case is a situation in which a teaching view could leave a student
believing something the evidence does not support. Each case says, before any
view is changed to handle it, what every possible output would do to the
student, graded by how visible the damage is:

  SAFE          correct, and nothing in it invites a wrong reading
  SAFE_REFUSAL  the view declines to judge and says why
  CRASH         an error with no reason a student can act on
  MISLEADING    true numbers arranged so that they invite a false conclusion
  FALSE         a statement that is wrong, shown as a finding

`allowed` is exactly the outputs graded SAFE or SAFE_REFUSAL (a test checks
that). `prefer` is what a good view says; null when two answers are equally
defensible. `requires` lists what the rendered card must show whatever the
verdict. `component` names the code under test. A component that does not
exist yet (evidence_card, pair_view) has its expectations fixed here now, and
its tests skip until it lands.

Synthetic, not real: the numbers are chosen to sit where a rule breaks, not
taken from a record. Real-record variants are the second half of roadmap
step 2.
"""
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).with_name("hostile_set.json")
SEED = 20261023
KEY, KEY_TZ, KEY_PCM = "k_def2svp_gas", "k_def2tzvp_gas", "k_def2svp_pcm_water"
GRADES = ("SAFE", "SAFE_REFUSAL", "CRASH", "MISLEADING", "FALSE")
VERDICTS = ("SEPARATED", "SEPARATED_WITHIN_PARTITION", "OVERLAPPING", "CANNOT_JUDGE",
            "NOT_COMPARED", "RULE_DEPENDENT", "REFUSED", "ERROR")
rng = np.random.default_rng(SEED)


def ens(n, nn_mean, nn_sd, rs_shift=0.0, rs_sd=0.0, genus=0, key=KEY, site="C4-H", fixed=None,
        rs_exact=None, folded=None, partition=None):
    """n conformers of one molecule at one site: nearest-nucleus means, the
    radius-scaled mean (nn + rs_shift + noise), genus and analysis key each.
    rs_exact: radius-scaled means as given, unrounded (float-noise cases).
    folded: each member's folded flag, stated apart from genus (A12); when
    absent a member is folded iff genus > 0, as in the cases before A9."""
    nn = fixed if fixed is not None else nn_mean + nn_sd * rng.standard_normal(n)
    noise = rs_sd * rng.standard_normal(len(nn))
    rs = [None if a is None else a + rs_shift + e for a, e in zip(nn, noise)]   # None: buried
    g = genus if isinstance(genus, list) else [genus] * len(nn)
    k = key if isinstance(key, list) else [key] * len(nn)
    out = [{"site": site, "nn": None if a is None else round(float(a), 3),
            "rs": None if a is None else round(float(b), 3), "genus": gg, "key": kk}
           for a, b, gg, kk in zip(nn, rs, g, k)]
    for i, m in enumerate(out):
        if rs_exact is not None:
            m["rs"] = rs_exact[i]
        if folded is not None:
            m["folded"] = folded if isinstance(folded, bool) else folded[i]
        if partition is not None and partition[i] is not None:
            # as given, not |nn - rs| (A14); "missing" writes an explicit null
            m["partition"] = None if partition[i] == "missing" else partition[i]
    return out


def grades(**kw):
    """Outcome grades; every verdict not named is FALSE unless given."""
    out = {v: kw.pop(v, "FALSE") for v in VERDICTS}
    assert not kw, kw
    return out


CASES = []


def case(id, cls, component, why, a=None, b=None, outcomes=None, prefer=None, requires=(), extra=None):
    allowed = sorted(v for v, g in (outcomes or {}).items() if g in ("SAFE", "SAFE_REFUSAL"))
    CASES.append(dict(id=id, cls=cls, component=component, a=a, b=b, outcomes=outcomes,
                      allowed=allowed, prefer=prefer, requires=list(requires), why=why, **(extra or {})))


# A refusal that says why is always acceptable when the view cannot judge;
# an unexplained error never is.
REFUSE_OK = dict(NOT_COMPARED="SAFE_REFUSAL", REFUSED="SAFE_REFUSAL", ERROR="CRASH")

# ------------------------------------------------ A1: different protocols
for vid, kb, why in [
    ("basis", KEY_TZ, "def2-TZVP moves C18 V_max by +1.5 to +2.8 kcal/mol, the size of the differences taught"),
    ("solvent", KEY_PCM, "PCM water moves chlorobenzene's V_max by -3.60 under one convention and +0.05 under the other"),
]:
    case(f"key_mismatch_{vid}", "A1 protocol mismatch", "compare",
         f"Same molecule pair, one ensemble computed under another protocol: {why}. Any number compared "
         "across analysis keys is a protocol difference shown as a chemical one.",
         a=ens(10, 3.0, 1.0), b=ens(10, 7.5, 0.6, key=kb),
         outcomes=grades(**REFUSE_OK), prefer="NOT_COMPARED", requires=["no_numbers"])
case("key_mismatch_one_stray_member", "A1 protocol mismatch", "compare",
     "19 conformers share a key and one does not (a re-run added to the folder). Averaging it in "
     "silently mixes protocols; the stray member must stop the comparison.",
     a=ens(20, 3.0, 1.0, key=[KEY] * 19 + [KEY_TZ]), b=ens(20, 7.5, 0.6),
     outcomes=grades(**REFUSE_OK), prefer="NOT_COMPARED", requires=["no_numbers"])
case("key_mismatch_stray_in_b", "A1 protocol mismatch", "compare",
     "As above, the stray member in the second molecule, last in the list.",
     a=ens(20, 3.0, 1.0), b=ens(20, 7.5, 0.6, key=[KEY] * 19 + [KEY_PCM]),
     outcomes=grades(**REFUSE_OK), prefer="NOT_COMPARED", requires=["no_numbers"])
case("key_mismatch_single_conformers", "A1 protocol mismatch", "compare",
     "One conformer each, different keys. Both refusals are safe; 'not compared' names the real "
     "problem, 'single conformer' hides it behind a lesser one.",
     a=ens(1, 3.0, 0.0), b=ens(1, 7.5, 0.0, key=KEY_TZ),
     outcomes=grades(CANNOT_JUDGE="SAFE_REFUSAL", **REFUSE_OK), prefer="NOT_COMPARED",
     requires=["no_numbers"])

# ------------------------------------------------ A2: one conformer
for vid, na, nb in [("1v1", 1, 1), ("1v20", 1, 20), ("20v1", 20, 1)]:
    case(f"single_conformer_{vid}", "A2 single conformer", "compare",
         "With one conformer on either side there is no conformer range to beat; a large "
         "difference between two single conformers is routine for a C18 chain (pilot ranges 2.5-6.2). "
         "A range of 0.00 printed for n = 1 reads as 'no variation', which is false.",
         a=ens(na, 2.0, 1.0), b=ens(nb, 8.0, 1.0),
         outcomes=grades(CANNOT_JUDGE="SAFE_REFUSAL", SEPARATED="MISLEADING", OVERLAPPING="MISLEADING",
                         SEPARATED_WITHIN_PARTITION="MISLEADING", **REFUSE_OK),
         prefer="CANNOT_JUDGE", requires=["no_single_range"])

# ------------------------------------------------ A3: separated, but within the partition
for vid, diff, part in [("well_inside", 0.6, 1.5), ("just_inside", 1.15, 1.2), ("negative", -0.7, 1.4)]:
    a_vals = np.linspace(2.0, 2.4, 10)   # range 0.4 < |diff|: disjoint by construction
    b_vals = a_vals + diff
    case(f"within_partition_{vid}", "A3 within partition", "compare",
         "Every conformer differs in the same direction, but by less than the two partition rules "
         "disagree at this site: the 'difference' is the size of an arbitrary convention.",
         a=ens(10, 0, 0, rs_shift=part, fixed=a_vals), b=ens(10, 0, 0, rs_shift=part, fixed=b_vals),
         outcomes=grades(SEPARATED_WITHIN_PARTITION="SAFE", SEPARATED="MISLEADING", OVERLAPPING="FALSE",
                         RULE_DEPENDENT="SAFE", **REFUSE_OK),
         prefer="SEPARATED_WITHIN_PARTITION", requires=["partition_shown"])

# ------------------------------------------------ A4: sign flips under the other rule
for vid, nn_a, nn_b, rs_a, rs_b in [("clear", 0.0, 3.0, 1.75, 1.25), ("narrow", 1.0, 3.2, 2.6, 2.4),
                                    ("negative_side", 3.0, 0.0, 1.25, 1.75)]:
    case(f"rule_sign_flip_{vid}", "A4 rule-dependent sign", "compare",
         "Nearest-nucleus separates the molecules beyond the (per-member) partition sensitivity, but "
         "the radius-scaled rule puts the difference the other way. The direction of the effect is a "
         "property of the convention, not of the molecules (frozen plan section 4, condition 3).",
         a=ens(10, nn_a, 0.15, rs_shift=rs_a - nn_a, rs_sd=0.05),
         b=ens(10, nn_b, 0.15, rs_shift=rs_b - nn_b, rs_sd=0.05),
         outcomes=grades(RULE_DEPENDENT="SAFE", SEPARATED_WITHIN_PARTITION="SAFE", OVERLAPPING="MISLEADING",
                         SEPARATED="MISLEADING", **REFUSE_OK),
         prefer="RULE_DEPENDENT", requires=["partition_shown"])

# ------------------------------------------------ A5: folded conformers
case("folded_drive_separation", "A5 folded conformers", "compare",
     "All 10 v 10 separate, but only because A's two folded conformers sit at its low end; over "
     "extended conformers the ranges overlap. 'Separated' alone teaches that the molecules differ, "
     "when their extended shapes do not (CLAUDE.md: extrema are robust only while nothing folds).",
     a=ens(10, 0, 0, rs_shift=0.2, genus=[2, 3] + [0] * 8,
           fixed=[0.5, 0.8] + list(np.linspace(3.0, 4.5, 8))),
     b=ens(10, 0, 0, rs_shift=0.2, fixed=list(np.linspace(4.6, 6.0, 10))),
     outcomes=grades(SEPARATED="SAFE", SEPARATED_WITHIN_PARTITION="SAFE", OVERLAPPING="FALSE", **REFUSE_OK),
     prefer="SEPARATED", requires=["extended_row", "extended_verdict"])
case("folded_hide_separation", "A5 folded conformers", "compare",
     "The reverse: extended conformers separate, one folded outlier makes the full ranges overlap. "
     "'Overlapping' is the honest all-member answer; the card must still show the extended picture.",
     a=ens(10, 0, 0, rs_shift=0.2, genus=[0] * 9 + [2], fixed=list(np.linspace(1.0, 2.5, 9)) + [5.5]),
     b=ens(10, 0, 0, rs_shift=0.2, fixed=list(np.linspace(4.0, 6.0, 10))),
     outcomes=grades(OVERLAPPING="SAFE", SEPARATED="FALSE", **REFUSE_OK),
     prefer="OVERLAPPING", requires=["extended_row", "extended_verdict"])
case("folded_extended_single", "A5 folded conformers", "compare",
     "Only one of A's conformers is extended: there is no extended range to show, and the card must "
     "say so, not print a range of 0.00.",
     a=ens(5, 0, 0, rs_shift=0.2, genus=[2, 2, 3, 2, 0], fixed=[1.0, 1.4, 1.8, 2.0, 2.2]),
     b=ens(5, 0, 0, rs_shift=0.2, fixed=[5.0, 5.2, 5.4, 5.6, 6.0]),
     outcomes=grades(SEPARATED="SAFE", SEPARATED_WITHIN_PARTITION="SAFE", OVERLAPPING="FALSE", **REFUSE_OK),
     prefer="SEPARATED", requires=["extended_row", "no_single_range"])

# ------------------------------------------------ A6: buried sites
case("buried_in_one_conformer", "A6 buried site", "compare",
     "In one folded conformer the site has no surface (mean None). Dropping it silently changes n; "
     "averaging fails. The view should say the site is buried in that conformer.",
     a=ens(10, 0, 0, rs_shift=0.3, genus=[2] + [0] * 9, fixed=[None] + list(np.linspace(1.0, 2.0, 9))),
     b=ens(10, 0, 0, rs_shift=0.3, fixed=list(np.linspace(4.0, 5.0, 10))),
     outcomes=grades(CANNOT_JUDGE="SAFE_REFUSAL", SEPARATED="MISLEADING",
                     SEPARATED_WITHIN_PARTITION="MISLEADING", OVERLAPPING="MISLEADING", **REFUSE_OK),
     prefer="REFUSED")
case("buried_everywhere", "A6 buried site", "compare",
     "The site is buried in every conformer of both molecules: there is nothing to compare.",
     a=ens(5, 0, 0, fixed=[None] * 5), b=ens(5, 0, 0, fixed=[None] * 5),
     outcomes=grades(CANNOT_JUDGE="SAFE_REFUSAL", **REFUSE_OK), prefer="REFUSED")

# ------------------------------------------------ A7: different on average, ranges overlap
for vid, mb in [("strong", 3.2), ("weak", 1.2)]:
    case(f"overlap_mean_shift_{vid}", "A7 overlapping but shifted", "compare",
         "The means differ consistently but the ranges overlap: some single conformations of the "
         "two molecules are indistinguishable. 'Overlapping' is right; a card that then reads as "
         "'no difference' or 'identical' is false (frozen plan section 5, 'unresolved').",
         a=ens(20, 0, 0, rs_shift=0.3, fixed=list(np.linspace(0.0, 4.0, 20))),
         b=ens(20, 0, 0, rs_shift=0.3, fixed=list(np.linspace(mb, 4.0 + mb, 20))),   # overlap: mb < 4
         outcomes=grades(OVERLAPPING="SAFE", SEPARATED="FALSE", SEPARATED_WITHIN_PARTITION="FALSE",
                         **REFUSE_OK),
         prefer="OVERLAPPING", requires=["no_identity_claim"])

# ------------------------------------------------ A8: site missing
case("site_missing_in_b", "A8 site missing", "compare",
     "The second molecule has no site of that label (C2-H asked of a molecule whose C2 is sp2 and "
     "labelled differently). The refusal must name the sites that exist.",
     a=ens(5, 3.0, 0.5, site="C2-H"), b=ens(5, 3.0, 0.5, site="C3-H"),
     outcomes=grades(**REFUSE_OK), prefer="REFUSED", extra={"label": "C2-H"})

# ------------------------------------------------ A10: values that are not numbers
# Added from the held-out probe of PR #2 (x_nan_mean): min() and max() skip a
# NaN, so a NaN mean quietly drops out of the ranges and the card claims a
# separation it never measured.
nan, inf = float("nan"), float("inf")
case("nan_in_one_conformer", "A10 non-finite value", "compare",
     "One conformer's site mean is NaN (a failed integration, a bad record). It must stop the "
     "comparison, naming the conformer; any verdict computed around it is unmeasured.",
     a=ens(5, 0, 0, rs_shift=0.3, fixed=[1.0, nan, 1.2, 1.3, 1.4]),
     b=ens(5, 0, 0, rs_shift=0.3, fixed=[4.0, 4.1, 4.2, 4.3, 4.4]),
     outcomes=grades(CANNOT_JUDGE="MISLEADING", SEPARATED_WITHIN_PARTITION="FALSE",
                     OVERLAPPING="MISLEADING", RULE_DEPENDENT="FALSE", **REFUSE_OK), prefer="REFUSED")
case("inf_in_other_molecule", "A10 non-finite value", "compare",
     "An infinite mean in the second molecule: the same refusal.",
     a=ens(5, 0, 0, rs_shift=0.3, fixed=[1.0, 1.1, 1.2, 1.3, 1.4]),
     b=ens(5, 0, 0, rs_shift=0.3, fixed=[4.0, 4.1, inf, 4.3, 4.4]),
     outcomes=grades(CANNOT_JUDGE="MISLEADING", SEPARATED_WITHIN_PARTITION="FALSE",
                     OVERLAPPING="MISLEADING", RULE_DEPENDENT="FALSE", **REFUSE_OK), prefer="REFUSED")
case("nan_under_other_rule_only", "A10 non-finite value", "compare",
     "Finite under nearest-nucleus, NaN under radius-scaled: the sign test and the partition "
     "sensitivity both need the second rule, so the comparison is not measured.",
     a=ens(5, 0, 0, fixed=[1.0, 1.1, 1.2, 1.3, 1.4], rs_exact=[1.2, 1.3, nan, 1.5, 1.6]),
     b=ens(5, 0, 0, fixed=[4.0, 4.1, 4.2, 4.3, 4.4], rs_exact=[4.2, 4.3, 4.4, 4.5, 4.6]),
     outcomes=grades(CANNOT_JUDGE="MISLEADING", SEPARATED_WITHIN_PARTITION="FALSE",
                     OVERLAPPING="MISLEADING", RULE_DEPENDENT="FALSE", **REFUSE_OK), prefer="REFUSED")

# ------------------------------------------------ A11: an empty side
# From the probe (x_empty_side): an empty list reached views_a[0].
for vid, a, b in [("a", [], ens(5, 4.0, 0.2)), ("b", ens(5, 1.0, 0.2), []), ("both", [], [])]:
    case(f"empty_side_{vid}", "A11 empty side", "compare",
         "No conformers on one side (a filter removed them all, a load found none). There is "
         "nothing to compare, and the refusal must say which side is empty.",
         a=a, b=b, outcomes=grades(CANNOT_JUDGE="SAFE_REFUSAL", **REFUSE_OK), prefer="REFUSED",
         extra={"label": "C4-H"})

# ------------------------------------------------ A12: a sign at noise level
# From the probe (x_flip_float_noise): Delta under the radius-scaled rule of
# -2.2e-16 is floating-point rounding, not a direction. The card prints two
# decimals, so the threshold must lie between rounding noise and a flip it
# can show; these cases pin both sides and leave the middle undecided.
_nn_a, _nn_b = [1.0, 1.1], [4.0, 4.1]     # Delta_nn = +3.0, per-member partition ~1.55
for vid, rs_b, outc, pref, why in [
    ("float_noise", [2.55, 2.55],
     grades(SEPARATED="SAFE", RULE_DEPENDENT="MISLEADING", SEPARATED_WITHIN_PARTITION="MISLEADING",
            OVERLAPPING="FALSE", **REFUSE_OK), "SEPARATED",
     "Delta_rs is -2.2e-16: rounding noise in the means, not a direction. Calling it rule-dependent "
     "tells a student the convention flips a +3.0 difference."),
    ("flip_0.01", [2.54, 2.54],
     grades(RULE_DEPENDENT="SAFE", SEPARATED="MISLEADING", SEPARATED_WITHIN_PARTITION="SAFE",
            OVERLAPPING="FALSE", **REFUSE_OK), "RULE_DEPENDENT",
     "Delta_rs is -0.01: small, but a real opposite sign the card can print."),
    # probe round 2 (y_tol_exact): exactly the tolerance, which floats put at
    # 0.004999999999999893; the rule as written (>= SIGN_TOL) counts it
    ("flip_0.005", [2.545, 2.545],
     grades(RULE_DEPENDENT="SAFE", SEPARATED="MISLEADING", SEPARATED_WITHIN_PARTITION="SAFE",
            OVERLAPPING="FALSE", **REFUSE_OK), "RULE_DEPENDENT",
     "Delta_rs is -0.005, exactly the stated tolerance: a rule that says '0.005 or more' must "
     "count it whatever the last bit of the float subtraction is."),
    ("flip_0.004", [2.546, 2.546],
     grades(RULE_DEPENDENT="SAFE", SEPARATED="SAFE", SEPARATED_WITHIN_PARTITION="SAFE",
            OVERLAPPING="FALSE", **REFUSE_OK), None,
     "Delta_rs is -0.004, which the card prints as -0.00: either reading is defensible, so "
     "no preference; only OVERLAPPING is wrong."),
]:
    case(f"sign_noise_{vid}", "A12 sign at noise level", "compare", why,
         a=ens(2, 0, 0, fixed=_nn_a, rs_exact=[2.5500000000000003, 2.5500000000000003]),
         b=ens(2, 0, 0, fixed=_nn_b, rs_exact=rs_b), outcomes=outc, prefer=pref)

# ------------------------------------------------ A14: a partition that is not a number
# Probe round 2 (y_nan_partition*, y_inf_partition, y_negative_partition):
# the partition sensitivity is an absolute difference, so it is finite and
# >= 0. Used as given otherwise, NaN is kept or skipped by max() depending on
# order, inf turns every separation "within", a negative value is corrupt.
for vid, side, pos, val in [("nan_first", "a", 0, nan), ("nan_last", "b", 4, nan),
                            ("inf", "a", 2, inf), ("negative", "b", 1, -5.0)]:
    parts = {"a": [None] * 5, "b": [None] * 5}
    parts[side][pos] = val
    case(f"partition_{vid}", "A14 partition not a number", "compare",
         f"Every mean is finite, but one conformer's partition sensitivity is {val}. A partition "
         "is |difference between the two rules' means|: not finite, or negative, means the record "
         "is corrupt, and the comparison must stop, naming the conformer.",
         a=ens(5, 0, 0, rs_shift=0.3, fixed=[1.0, 1.1, 1.2, 1.3, 1.4], partition=parts["a"]),
         b=ens(5, 0, 0, rs_shift=0.3, fixed=[4.0, 4.1, 4.2, 4.3, 4.4], partition=parts["b"]),
         outcomes=grades(SEPARATED="FALSE", SEPARATED_WITHIN_PARTITION="FALSE", OVERLAPPING="FALSE",
                         RULE_DEPENDENT="FALSE", CANNOT_JUDGE="MISLEADING", **REFUSE_OK),
         prefer="REFUSED")

# Probe round 3 (z_part_none_hides): a partition left unstated (None) while
# both rules' means are finite. Read as 0 it hid an 8 kcal/mol disagreement
# between the rules; the two means fix it, so it must count.
case("partition_missing_means_finite", "A14 partition not a number", "compare",
     "One conformer states no partition sensitivity, but both its means are there (1.0 and -7.0): "
     "the rules disagree by 8 kcal/mol at that conformer, more than the +3.0 difference. Treating "
     "the missing value as 0 shows a clean separation that the convention could undo.",
     a=ens(5, 0, 0, fixed=[1.0, 1.1, 1.2, 1.3, 1.4], rs_exact=[-7.0, 1.1, 1.2, 1.3, 1.4],
           partition=["missing", 0.0, 0.0, 0.0, 0.0]),
     b=ens(5, 0, 0, fixed=[4.0, 4.1, 4.2, 4.3, 4.4], rs_exact=[4.0, 4.1, 4.2, 4.3, 4.4],
           partition=[0.0] * 5),
     outcomes=grades(SEPARATED_WITHIN_PARTITION="SAFE", SEPARATED="MISLEADING", OVERLAPPING="FALSE",
                     RULE_DEPENDENT="SAFE", CANNOT_JUDGE="MISLEADING", **REFUSE_OK),
     prefer="SEPARATED_WITHIN_PARTITION")

# Probe round 4 (w_part_contradicts_means): a stated partition of 0.0 while the
# member's own means differ by 8. The means are the measurement; a partition
# smaller than |nn - rs| contradicts them and must not hide the difference.
for vid, given in (("zero", 0.0), ("small", 0.5)):
    case(f"partition_contradicts_means_{vid}", "A14 partition not a number", "compare",
         f"One conformer states a partition of {given}, but its means (1.0 and -7.0) differ by 8: "
         "the stated number contradicts the record's own measurement. Trusted, it shows a clean "
         "separation that the convention could undo.",
         a=ens(5, 0, 0, fixed=[1.0, 1.1, 1.2, 1.3, 1.4], rs_exact=[-7.0, 1.1, 1.2, 1.3, 1.4],
               partition=[given, 0.0, 0.0, 0.0, 0.0]),
         b=ens(5, 0, 0, fixed=[4.0, 4.1, 4.2, 4.3, 4.4], rs_exact=[4.0, 4.1, 4.2, 4.3, 4.4],
               partition=[0.0] * 5),
         outcomes=grades(SEPARATED_WITHIN_PARTITION="SAFE", SEPARATED="MISLEADING", OVERLAPPING="FALSE",
                         RULE_DEPENDENT="SAFE", CANNOT_JUDGE="MISLEADING", **REFUSE_OK),
         prefer="SEPARATED_WITHIN_PARTITION")

# ------------------------------------------------ A13: what "folded" means
# The cross-repo review (2026-10-23) found three definitions in use: the
# frozen study (contact OR genus > 0), espmap v0.3 (contact; genus only a
# diagnostic) and Arrowfield (genus > 0). A member's folded flag is a fact the
# producer states; the view must use it, never re-derive it from genus.
case("folded_by_contact_genus0", "A13 folded flag", "compare",
     "Two of A's conformers are folded by an intramolecular contact but their surfaces have genus "
     "0. They drive the separation; the extended view must leave them out (n = 8 / 10).",
     a=ens(10, 0, 0, rs_shift=0.2, fixed=[0.5, 0.8] + list(np.linspace(3.0, 4.5, 8)),
           folded=[True, True] + [False] * 8),
     b=ens(10, 0, 0, rs_shift=0.2, fixed=list(np.linspace(4.6, 6.0, 10)), folded=False),
     outcomes=grades(SEPARATED="SAFE", SEPARATED_WITHIN_PARTITION="SAFE", OVERLAPPING="FALSE", **REFUSE_OK),
     prefer="SEPARATED", requires=["extended_n", "extended_verdict"], extra={"extended_n": [8, 10]})
case("genus_without_folding", "A13 folded flag", "compare",
     "Every conformer's surface has genus 1 (a macrocycle's hole), and none is folded. Treating "
     "genus as folding would empty the extended set; there is nothing to leave out.",
     a=ens(6, 0, 0, rs_shift=0.2, genus=1, fixed=list(np.linspace(1.0, 2.0, 6)), folded=False),
     b=ens(6, 0, 0, rs_shift=0.2, genus=1, fixed=list(np.linspace(4.0, 5.0, 6)), folded=False),
     outcomes=grades(SEPARATED="SAFE", SEPARATED_WITHIN_PARTITION="SAFE", OVERLAPPING="FALSE", **REFUSE_OK),
     prefer="SEPARATED", requires=["extended_n"], extra={"extended_n": None})

# ------------------------------------------------ A15: what a one-geometry refusal says
# The user's decision (2026-10-24): espmap's compare allows one geometry on a
# side (flagged, never resolved); refusing it belongs here, and the refusal
# must name what cannot be judged: conformational robustness. "Single
# conformer: cannot judge" leaves a student free to read it as "too little
# data" rather than "shape alone could produce this difference".
for vid, na, nb in [("1v1", 1, 1), ("1v8", 1, 8), ("8v1", 8, 1)]:
    case(f"one_geometry_says_why_{vid}", "A15 one-geometry wording", "compare",
         "With one geometry on a side, the conformational spread is unmeasured. The card must say "
         "it cannot judge conformational robustness, and that it is because one side has one "
         "geometry: the reason a single structure cannot show an effect on a flexible molecule.",
         a=ens(na, 2.0, 1.0), b=ens(nb, 6.0, 1.0),
         outcomes=grades(CANNOT_JUDGE="SAFE_REFUSAL", SEPARATED="MISLEADING", OVERLAPPING="MISLEADING",
                         SEPARATED_WITHIN_PARTITION="MISLEADING", **REFUSE_OK),
         prefer="CANNOT_JUDGE", requires=["robustness_wording", "no_single_range"])

# ------------------------------------------------ A9: language on every card
# Checked on every A case's card, not as a case of its own: see LANGUAGE below.
LANGUAGE = {
    "forbidden": ["explain", "receptor", "specificity", "attract", "binds", "binding", "active",
                  "responsible", "causes", "identical", "no difference"],
    "scope_caveat": "free molecule",
    "why": "A difference card shows free-molecule, gas-phase ESP. Words that suggest a biological "
           "mechanism, or that overlapping ranges mean 'the same', turn a measurement into a claim "
           "nobody tested; and every comparable card should say what the number is of.",
}

# ------------------------------------------------ B: evidence (component to come)
EV = [
    ("not_tested_is_not_negative", "B1 absence of evidence",
     dict(effect=None, significance="not_tested", role="treatment"),
     {"NOT_TESTED": "SAFE", "NEGATIVE": "FALSE", "OBSERVED": "FALSE"}, "NOT_TESTED",
     "A compound nobody tested shows as 'not tested', never as a negative result."),
    ("no_effect_significant_test", "B1 absence of evidence",
     dict(effect="no_effect", significance="not_significant", role="treatment"),
     {"NEGATIVE": "SAFE", "NOT_TESTED": "MISLEADING", "OBSERVED": "FALSE"}, "NEGATIVE",
     "Tested with no detectable effect is a negative result, and must not be shown as untested."),
    ("effect_not_reported", "B1 absence of evidence",
     dict(effect="increased", significance="not_reported", role="treatment"),
     {"OBSERVED": "SAFE", "CONTESTED": "MISLEADING", "NEGATIVE": "FALSE"}, "OBSERVED",
     "An effect reported without a test is still an observation; the card shows 'significance not "
     "reported' beside it, not a stronger or weaker state."),
    ("two_studies_disagree", "B2 contested",
     dict(effects=["increased", "decreased"], same_taxon=True), {"CONTESTED": "SAFE", "OBSERVED": "MISLEADING",
                                                                   "NEGATIVE": "FALSE"}, "CONTESTED",
     "Two studies on one taxon point opposite ways: show both, labelled contested, never the newer one alone."),
    ("teaching_control_not_observed", "B3 teaching control",
     dict(compound="arrowfield.control.z3_z13_18_al", evidence_status="teaching_control"),
     {"TEACHING_CONTROL": "SAFE", "OBSERVED": "FALSE", "COMPUTATIONAL_MODEL": "MISLEADING"}, "TEACHING_CONTROL",
     "The aldehyde control is constructed, not found in any Sesiidae: labelling it observed invents evidence."),
    ("counterfactual_labelled", "B3 teaching control",
     dict(compound="(hypothetical)", evidence_status="computational_counterfactual"),
     {"COMPUTATIONAL_MODEL": "SAFE", "OBSERVED": "FALSE", "TEACHING_CONTROL": "MISLEADING"}, "COMPUTATIONAL_MODEL",
     "A structure computed to ask 'what if' is a model, and is never shown as observed."),
    ("assay_contexts_kept_apart", "B4 incompatible contexts",
     dict(assays=["field_trapping", "gc_ead"], same_compounds=True),
     {"TWO_STATEMENTS": "SAFE", "MERGED": "MISLEADING"}, "TWO_STATEMENTS",
     "A trap catch and an antennal response to the same pair are different evidence; one merged "
     "'active' badge erases which was measured."),
    ("populations_kept_apart", "B4 incompatible contexts",
     dict(assays=["field_trapping", "field_trapping"], localities=["gongju", "jinju"], effects=["increased", "no_effect"]),
     {"TWO_STATEMENTS": "SAFE", "MERGED": "MISLEADING", "CONTESTED": "SAFE"}, "TWO_STATEMENTS",
     "Two populations of one species responding differently is a population difference, not noise to average."),
    ("difference_not_specificity", "B5 explanation",
     dict(comparison="resolved", observation="species_specific_attraction"),
     {"SIDE_BY_SIDE": "SAFE", "LINKED_AS_CAUSE": "FALSE"}, "SIDE_BY_SIDE",
     "A resolved ESP difference shown next to a species-specific response must not be joined to it "
     "by 'because', 'explains' or an arrow: nothing here tests the link."),
]
# Retired (not edited) by the evidence views' PR: these abstract cases were
# re-expressed against the real bundle in tests/pedagogy/evidence_card_set.json
# (#12), which the built views are tested against. Each names the cases there
# that replace it.
REPLACED_BY = {
    "not_tested_is_not_negative": ["compound_not_in_observation"],
    "no_effect_significant_test": ["no_effect_not_significant"],
    "effect_not_reported": ["no_effect_significance_not_reported"],
    "two_studies_disagree": ["opposite_effects_two_places", "increase_and_decrease_both_kept"],
    "teaching_control_not_observed": ["teaching_control_classes", "teaching_control_in_observation"],
    "counterfactual_labelled": ["counterfactual_classes"],
    "assay_contexts_kept_apart": ["gland_vs_lure"],
    "populations_kept_apart": ["glands_two_populations", "population_specific_roles"],
    "difference_not_specificity": ["study_verdict_on_the_same_page", "study_verdict_beside_classes"],
}
for eid, cls, inp, outc, pref, why in EV:
    case(eid, cls, "evidence_card", why, outcomes=outc, prefer=pref,
         extra={"input": inp, "superseded_by": [f"tests/pedagogy/evidence_card_set.json#{r}"
                                                for r in REPLACED_BY[eid]]})

# ------------------------------------------------ C: pictures (component to come)
PAIR_RETIRES = json.loads((Path(__file__).with_name("pair_view_set.json")).read_text(encoding="utf-8"))["retires"]
for eid, why in [
    ("pair_different_render_keys", "Two ESP pictures with different colour scales side by side: the same "
                                   "colour means different kcal/mol, so the eye compares the wrong things."),
    ("pair_one_folded", "A folded and an extended conformer side by side read as two molecules' shapes, "
                        "not as one molecule's conformational range."),
]:
    # Retired (not edited) by the pair view's PR: re-expressed against real
    # records in tests/pedagogy/pair_view_set.json (#14, #17), whose "retires"
    # map names each replacement; pair_one_folded's same-compound half was
    # reversed there by the user's decision (Act B renders folded/extended).
    case(eid, "C1 pictures", "pair_view", why,
         outcomes={"REFUSED": "SAFE_REFUSAL", "SHARED_SCALE": "SAFE", "SIDE_BY_SIDE_AS_IS": "MISLEADING"},
         prefer="SHARED_SCALE" if "render" in eid else "REFUSED",
         extra={"superseded_by": [f"tests/pedagogy/pair_view_set.json#{r}" for r in PAIR_RETIRES[eid]]})


def premises():
    """Each class's construction must hold, or the case tests something else."""
    def nn(ms):
        return [m["nn"] for m in ms if m["nn"] is not None]

    def disjoint(a, b):
        return max(a) < min(b) or max(b) < min(a)
    for c in CASES:
        if c["a"] is None:
            continue
        a, b = nn(c["a"]), nn(c["b"])
        cls = c["cls"][:2]
        if cls == "A3":
            part = max(abs(m["nn"] - m["rs"]) for m in c["a"] + c["b"])
            assert disjoint(a, b) and abs(np.mean(b) - np.mean(a)) < part, c["id"]
        if cls == "A4":
            part = max(abs(m["nn"] - m["rs"]) for m in c["a"] + c["b"])
            rs = [np.mean([m["rs"] for m in x]) for x in (c["a"], c["b"])]
            d = np.mean(b) - np.mean(a)
            assert disjoint(a, b) and abs(d) > part and np.sign(rs[1] - rs[0]) != np.sign(d), c["id"]
        if cls == "A7":
            assert not disjoint(a, b) and np.mean(b) > np.mean(a), c["id"]


def document():
    premises()
    for c in CASES:
        if c["outcomes"]:
            assert set(c["outcomes"].values()) <= set(GRADES), c["id"]
            assert c["prefer"] is None or c["prefer"] in c["allowed"], c["id"]
    doc = {"generator": "tests/pedagogy/make_set.py", "seed": SEED, "grades": list(GRADES),
           "language": LANGUAGE, "cases": CASES}
    return json.dumps(doc, indent=1, ensure_ascii=False) + "\n"


if __name__ == "__main__":
    OUT.write_text(document(), encoding="utf-8")
    print(f"{len(CASES)} cases -> {OUT}")
