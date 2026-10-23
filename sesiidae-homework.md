# Sesiidae homework: when does molecular electronics explain biological specificity?

Status: **draft.** The questions are fixed. The expected findings for acts
A–C come from the frozen study espmap-sesiidae-001 (Atlas results at
f286bc1: 5 C18 molecules × 20 shape-diverse conformers, def2-SVP, gas),
which replicated the pilot (`sesiidae-pilot.md`) where it mattered. Before
class, a person must still check every observation the homework uses
against its paper. The aim is that a student finishes able to make a
chemical argument from an ESP surface, and to say where that argument stops.

## What students already need

Lessons 1, 2, 4 and 5 of `lessons.md`: polar bonds, the carbonyl, esters, and
E/Z. No quantum chemistry vocabulary is needed.

## Molecules

All five from `molecules.yaml`, chosen from a dropdown, never typed. The
manifest guarantees every double bond's E/Z is set, and every compound is the
C18 skeleton with C1 bearing the oxygen:

| id | compound | in the homework | Atlas compound |
| --- | --- | --- | --- |
| `z3z13_18_oh` | (3Z,13Z)-octadeca-3,13-dien-1-ol | act A: alcohol | `spa.cmpd.z3_z13_18_oh` |
| `z3z13_18_al` | (3Z,13Z)-octadeca-3,13-dienal | act A: aldehyde | **teaching control**, none |
| `z3z13_18_ac` | (3Z,13Z)-octadeca-3,13-dien-1-yl acetate | acts A, B: acetate, Z3 | `spa.cmpd.z3_z13_18_oac` |
| `e3z13_18_ac` | (3E,13Z)-octadeca-3,13-dien-1-yl acetate | acts B, C: E3 | `spa.cmpd.e3_z13_18_oac` |
| `e2z13_18_ac` | (2E,13Z)-octadeca-2,13-dien-1-yl acetate | act C: double bond moved to C2 | `spa.cmpd.e2_z13_18_oac` |

Four of the five are compounds the Sesiidae Pheromone Atlas records;
`molecules.yaml` names each by its Atlas ID, and the manifest builder checks
that the structures agree. The aldehyde is different: it is a **teaching
control**, built here so act A changes only the head group on the Z3,Z13
backbone. It is not an observed pheromone. The aldehyde the Atlas records is
(2E,13Z)-octadeca-2,13-dienal (`spa.cmpd.e2_z13_18_al`). Every page and
card must show which of the two kinds a compound is.

The biological context (which species produce or respond to which compound,
in what blend) comes from the Atlas, through a teaching bundle exported from
its records. It is never written here by hand, and the assignment must not
state a biological fact the bundle does not carry.

## Data: precomputed, with the uncertainty attached

**Status:** computed and analysed. The five ensembles are the study
espmap-sesiidae-001: 20 conformers per molecule, espmap v0.2.1, seed
20261024, 100 of 100 members, no exclusions.
- **Frozen before the computation:** which molecules it compares and the
  rule that decides each comparison. They are in the Atlas at
  [7e34658](https://github.com/Azuromalachite/sesiidae-pheromone-atlas/tree/7e3465859bc2b8950078238f0bffb0ea6b4dd512/analyses/espmap-sesiidae-001):
  the case bundle (`case-bundle.yaml`) and the prospective analysis plan
  (`analysis-plan.md`). A dated clarification on reading its statistics is
  at 30c4a11.
- **Results:** at
  [f286bc1](https://github.com/Azuromalachite/sesiidae-pheromone-atlas/blob/f286bc1/analyses/espmap-sesiidae-001/results.md)
  (`results.md`). The "expected" answers below come from it, not from the
  pilot.
- **Where the records live:** the Atlas's `computed/espmap/<Atlas ID>/` for
  the four Atlas compounds; Arrowfield's `controls/` for the aldehyde.

Every verdict below is the study's, under its frozen rule
(`espmap-sesiidae-001/range_partition_sign_v1`). A site is **resolved** only
if all three hold:
1. the two molecules' conformer ranges do not overlap;
2. the difference exceeds the partition sensitivity;
3. the sign is the same under both partition rules.

**Different on average** means the means differ, but the observed conformer
ranges overlap; it is descriptive, not a weaker "resolved".

The 20 conformers are a designed spread of shapes, not random draws. So
non-overlap is a robustness criterion, not a probability.

- Each compound is to be precomputed as espmap records on a **conformer
  ensemble**: the lowest-energy
  conformer plus shapes spread by radius of gyration, so folded ones are
  included. A single conformer is not enough here. For C18 chains, site
  means move a median 5.1 kcal/mol across conformers, and V_min 5.2 once
  folded shapes are included. That is larger than most differences students
  would compare. **The pilot decided it: the ensemble is mandatory.** Act B's
  differences sit inside the conformer range, so a single conformer would let
  a student read a conformational accident as a stereo effect. Use 10 or
  more conformers, folded ones included. In the pilot, 5 showed the spread
  but were too few to measure it.
- Settings: B3LYP, gas phase (neutral molecules), one basis for the whole
  class, one `analysis_key` throughout. Every image in one comparison shares
  one `render_key`.
- Sites, the same labels on every compound: `C1-H` … `C5-H`, `C12-H` …
  `C15-H` (each carbon with its hydrogens) and `head O`, defined once in the
  manifest and passed to espmap as caller groups (`--sites`).
- "Recompute this yourself" stays available as an option. It is never the
  centre of the assignment, because a C18 calculation takes minutes even at
  the preview basis.

## Act A: what the head group does

**Compare** the alcohol, aldehyde and acetate (Z3,Z13 throughout).

**Predict first:** where is each molecule most electron-rich? Which has an
H-bond donor? Does adding the acetyl make "everything more negative"?

**Then:**
1. Identify the largest electrostatic differences between the three head
   groups and explain them with ordinary organic chemistry: electronegativity,
   resonance in the ester, donor and acceptor.
2. Explore how far the head-group perturbation appears to extend: compare
   the `C1-H` … `C5-H` sites with `C12-H` … `C15-H`, using the difference
   card. These chain-site comparisons are **descriptive and exploratory**,
   not a prespecified test: the frozen study (espmap-sesiidae-001) tests
   act A only on the global V_max and V_min, with `head O` as a secondary
   site.

**Expected (study results):** this act is the "ESP works" act.
- **Global V_max: resolved, the prespecified result.** The alcohol's O–H
  hydrogen gives V_max +34.0 to +51.7 kcal/mol across its 20 conformers;
  the acetate's peaks at +19.6 to +21.0. The difference of means is +26.5,
  and the ranges are far apart.
- **Global V_min: different on average** (−1.2). The ranges overlap.
- **`head O` (the O atom alone, secondary):** the alcohol's O is far more
  negative than the acetate's ester O (−29.0 to −32.0 against −8.0 to
  −13.8). The difference is −20.7, resolved.
- **The aldehyde control sits between them:** V_max +24.9 to +26.9, and head
  O −25.8 to −26.8. It is a constructed molecule, so treat these as
  descriptive chemistry only.

The chain sites are for looking, not for verdicts; no chain site met the
rule in any contrast. A student who calls a chain site "resolved" from an
exploratory comparison has missed the difference between a test and a look.
A student who describes how the difference fades with distance, with its
conformer range beside it, has used the card well.

`head O` is the O atom alone in every molecule. In the pilot, the obvious
comparison (`O31`) was not like-for-like: in the alcohol it held O plus its
H, in the acetate an ester O. That mismatch is worth teaching as a point of
its own.

## Act B: what E/Z does (E3 vs Z3 acetate)

**Predict first:** does changing E to Z at C3 mainly change the charge
distribution, or the shape?

**Then:**
1. Read the difference card at `C3-H`, `C4-H` and `head O`. Is each
   difference larger than that site's conformer range? Than its partition
   sensitivity?
2. Look at the extended and folded conformers of one compound side by side,
   on the same scale. What changed because the chain folded back on itself?
3. Are the differences mainly electronic or mainly geometric? How confident
   can you be, given the numbers in (1)?

**Expected (study results):** E/Z is **not resolved** at any site, as the
pilot predicted.
- **`C3-H`: different on average.** E3 is +2.31 kcal/mol above Z3 on
  average, the same direction under both partition rules and beyond the
  partition sensitivity (0.93). But the conformer ranges overlap widely:
  Z3 −2.8 to +5.0, E3 +1.0 to +6.6. The observed conformer ranges overlap
  at C3-H.
- **`C4-H`: unresolved** (−0.04, ranges 0.3–5.1 and 0.2–5.5).
- **`head O` and global V_min: different on average,** both under 1 kcal/mol.
  Global V_max is unresolved.
- No chain site met the rule. The pilot's chance hit at 5 conformers (C11-H)
  did not come back.

This is the act that teaches the difference between **a mean shift and a
resolved difference**. The correct answer is that the average moves, but by
less than the shapes alone move it, so a single structure cannot show an
E/Z effect at C3-H. It sits in the same order as espmap's stereo row
(1.6–1.9 kcal/mol, threonine).

## Act C: where the double bond is (E2 vs E3 acetate)

**Predict first:** moving C=C by one carbon. Which sites should change, and
by how much compared with act B?

**Then:**
1. Compare `C2-H`, `C3-H` and `C4-H` between the two compounds, with the
   difference card. Which differences exceed the conformer range?
2. Does free-molecule ESP give a convincing reason an organism or enzyme
   would treat these two compounds differently?

**Expected (study results):** resolved **at C4 only**, as the pilot found.
The pilot's call held up under the stricter test with twice the conformers.
- **`C4-H`: resolved.** It moves +5.22 kcal/mol. The conformer ranges do not
  overlap (E3 +0.2 to +5.5, E2 +6.7 to +9.5), the difference is well beyond
  the partition sensitivity (0.69), and both rules agree on the sign. C4 is
  the alkene carbon in E3 and allylic in E2. The gap between the
  selected-conformer ranges is about 1.1 kcal/mol: robust under this
  protocol (gas phase, def2-SVP, basis sensitivity not run), but not wide.
- **`C2-H`: different on average** (−2.50). The ranges overlap: E3 +3.3 to
  +10.1, E2 −0.7 to +6.9.
- **`C3-H`: unresolved** (−0.35). The two partition rules do not even agree
  on its sign.
- So moving the bond is more visible than flipping it, but only at one site.

Non-overlap of 20 shape-spread conformers is a robustness criterion, not a
probability. Don't translate it into a p-value: the conformers were chosen
to span shapes, not drawn at random.

Question 2's honest answer stays "no". A 4 kcal/mol change at one CH₂ on a
free molecule does not by itself explain discrimination by an enzyme or
receptor. That is Act D.

## Act D: what is missing

The organism makes, or responds to, one particular compound. Suppose
free-molecule ESP does not explain the distinction. What information is
missing, and how would you get it?

Answers worth credit: the receptor's or enzyme's geometry, active-site
residues and any metal or cofactor, how the substrate is held (orientation,
conformation), dynamics, transition-state energetics, and experimental
evidence (biosynthetic labelling, mutants, blends). This leads to
Arrowfield's later layers, and needs no AF3 or docking in the homework
itself.

## What the homework does not ask

- **"Which carbon will the desaturase act on?"** ESP cannot answer that, and
  asking invites over-reading. The question is always "does the electrostatic
  picture distinguish these sites strongly enough to plausibly matter?".
- A ranking from a single conformer, a single partition rule, an extremum, or
  a 5th percentile.

## Grading notes

- Credit the reasoning and the limits stated, not agreement with a number.
- A student who writes "the difference is smaller than the conformer range,
  so I cannot attribute it to E/Z" has the central point.
- A student who quotes a raw V_min difference with no uncertainty has missed
  it, even if the sign is right.
- A student who calls act B's C3-H "resolved" because the means differ has
  confused a mean shift with a resolved difference. A student who says "the
  E isomer is more positive at C3-H on average, but their conformer ranges
  overlap" has it exactly.

## Before class

- [x] Sesiidae pilot finished; acts B/C "expected" filled in from it
  (`sesiidae-pilot.md`).
- [x] Ensemble size decided: mandatory, 10+ conformers with folded ones.
- [x] Re-check the pilot's calls at 10 conformers (2026-10-23): C11 was
  chance, C4 holds (`sesiidae-pilot.md`).
- [ ] Re-check at the class basis if it is not def2-SVP. Basis moves C18
  V_max by +1.5 to +2.8, the same order as these differences.
- [x] `head O` defined identically for OH, CHO and OAc (act A).
- [x] Four compounds named by their Atlas IDs, structures checked against the
  Atlas; the aldehyde labelled a teaching control (2026-10-23).
- [x] Case bundle and prospective analysis plan frozen in the Atlas
  (`analyses/espmap-sesiidae-001/`, 7e34658, 2026-10-23).
- [x] Ensembles (20 conformers) for the five compounds computed with espmap
  v0.2.1 (2026-10-23/24, 100/100, no exclusions). The four Atlas compounds
  are in the Atlas's `computed/espmap/` (results at f286bc1); the aldehyde is
  in Arrowfield's `controls/`.
- [x] Acts A–C's "expected" rewritten from the study's results (2026-10-24).
- [ ] Notebook pinned to the Atlas teaching bundle and the records; the
  study adapter reads them.
- [ ] Every observation the homework uses checked by a person against its
  paper (a hard gate before class).
- [ ] `--require-stereo` in a tagged espmap. It is merged on master (3746cad)
  but not yet released, so the manifest-only guarantee still holds.
- [ ] Biological context from the Atlas teaching bundle, not written here.
- [ ] McMurry section numbers for the class edition.
