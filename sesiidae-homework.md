# Sesiidae homework: when does molecular electronics explain biological specificity?

Status: **draft.** The questions are fixed. The expected findings for acts B
and C now come from the Sesiidae pilot (`sesiidae-pilot.md`, 2026-10-22:
4 C18 molecules × 5 shape-diverse conformers, def2-SVP, gas). They must be
re-checked on the class ensemble before use. The aim is that a student finishes able to make a
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

**Status:** the five ensembles are being computed (from 2026-10-23) as
the study espmap-sesiidae-001. Which molecules it compares, and the rule that
decides each comparison, were fixed before the computation started. They are
in the Atlas, frozen at
[7e34658](https://github.com/Azuromalachite/sesiidae-pheromone-atlas/tree/7e3465859bc2b8950078238f0bffb0ea6b4dd512/analyses/espmap-sesiidae-001):
the case bundle (`case-bundle.yaml`) and the prospective analysis plan
(`analysis-plan.md`). Settings: K = 20, seed 20261024, espmap v0.2.1.
The plan treats the study as a confirmatory replication of the pilot below:
the homework's "expected" answers came from the pilot. Where the study's
results differ, the homework changes to match the study, never the reverse.
The ensembles' canonical home is the Atlas's `computed/espmap/<Atlas ID>/`;
Arrowfield reads them from the teaching bundle.

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
2. How far along the chain does the head-group change reach? Compare the
   `C1-H` … `C5-H` sites with `C12-H` … `C15-H`, using the difference card.

**Expected:** this act is the "ESP works" act. The head group should differ
locally and clearly, and the distant chain should barely move. Confirm with
the records before stating it.

From the pilot (OH vs OAc only), the clear signal is in the global extremum.
The alcohol's O–H hydrogen gives V_max +38 to +52 kcal/mol; the acetate
peaks at +19 to +21. No site mean separated, but the obvious site
comparison, `O31`, is not like-for-like: in the alcohol it holds O plus its
H, in the acetate an ester O. Define `head O` as the O atom alone, or as one
caller group with the same atoms in all three molecules, and teach the
mismatch as a point of its own.

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

**Expected (pilot, re-checked at 10 conformers):** E/Z is **not resolved**
on site means.
- At 10 conformers, C3-H and C4-H move by +1.3 and +1.0 kcal/mol, inside
  conformer ranges of 2.5–6.2.
- Global V_min and V_max move by under 0.5.
- No site separates. The pilot's one hit at 5 conformers (C11-H) disappeared,
  as a chance hit should.

The correct student answer is the one the grading notes reward: the
difference is smaller than the conformer range, so it cannot be attributed
to E/Z. It sits in the same order as espmap's stereo row (1.6–1.9 kcal/mol,
threonine).

## Act C: where the double bond is (E2 vs E3 acetate)

**Predict first:** moving C=C by one carbon. Which sites should change, and
by how much compared with act B?

**Then:**
1. Compare `C2-H`, `C3-H` and `C4-H` between the two compounds, with the
   difference card. Which differences exceed the conformer range?
2. Does free-molecule ESP give a convincing reason an organism or enzyme
   would treat these two compounds differently?

**Expected (pilot, re-checked at 10 conformers):** resolved **at C4 only**.
- C4-H moves +4.2 kcal/mol, just above both conformer ranges (4.1 and 2.0)
  and above the partition sensitivity (1.4). C4 is the alkene carbon in E3
  and allylic in E2.
- At 10 conformers a chance separation has probability about 1×10⁻⁵ per
  site, so this one is real. The margin over the wider range is small.
- C2-H (−3.1) and C3-H (−0.7) stay inside their ranges (5.9–8.1).
- So moving the bond is more visible than flipping it, but only at one site.

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
- [ ] Ensembles (20 conformers) for the five compounds computed with espmap
  v0.2.1 (started 2026-10-23), stored in the Atlas's `computed/espmap/`, and
  exported as a teaching bundle; notebook pinned to both.
- [ ] Acts B and C's "expected" re-checked against the study's results.
- [ ] `--require-stereo` in a tagged espmap. It is merged on master (3746cad)
  but not yet released, so the manifest-only guarantee still holds.
- [ ] Biological context from the Atlas teaching bundle, not written here.
- [ ] McMurry section numbers for the class edition.
