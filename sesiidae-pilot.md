# Sesiidae pilot: site differences against conformer spread

Run 2026-10-22, finished the same day. Data and scripts in `pilot/`:
`sesiidae_pilot.jsonl` (20 rows), `sesiidae_pilot.py` (the run) and
`sesiidae_analyse.py` (this analysis; its data path still points at the
scratchpad where it ran). **Pilot only:** move them into Arrowfield's repo,
or into espmap `bench/`, before anything here is quoted in class.

## Protocol

- **Molecules:** Z3Z13-18:OAc, E3Z13-18:OAc, E2Z13-18:OAc, Z3Z13-18:OH.
  All are mapped C1–C18, so a site label means the same carbon in every
  molecule.
- **Conformers:** 5 per molecule, chosen from the MMFF94s ensemble: the
  lowest-energy one plus shapes spread by radius of gyration. Folded
  (genus > 0) surfaces are kept: 2, 1, 1 and 1 per molecule.
- **Method:** B3LYP/def2-SVP with RI, **gas phase**, 0.26 bohr.
- **Sites:** each carbon with its hydrogens (`Cn-H`), the two acetate O, and
  the acetyl group. Site means are taken under the nearest-nucleus rule; the
  radius-scaled rule gives the partition column.
- **Test:** a site difference counts as **separated** only if the two
  molecules' 5-conformer ranges do not overlap at all.

**Chance level.** For two sets of 5 values from one distribution, the
probability that they do not overlap is 2/252 ≈ 0.8 % per site. Over about 21
sites that is ≈ 0.17 false "separated" per comparison, so a single hit at an
unexpected site is not evidence on its own.

## Results (site area-weighted means, kcal/mol)

| comparison | sites tested | separated | the sites the change is at |
| --- | --- | --- | --- |
| **E/Z at C3** (E3 − Z3 acetate) | 21 | C11-H (−1.29) | C3-H +1.81, C4-H +2.24: inside conformer ranges of 2.4–6.4 |
| **C=C position** (E2 − E3 acetate) | 21 | **C4-H (+3.57)**, C11-H (+1.72) | C2-H −2.69 and C3-H −1.44 inside ranges of 3.3–8.6; C4-H's ranges are 2.37 and 1.28 |
| **Head group** (OAc − OH, Z3Z13) | 19 | none | O31 −5.06, inside ranges of 9.6 and 4.4 |

Global extrema:
- E/Z changes neither V_min nor V_max (+0.02 / −0.05).
- Position changes them by +0.13 / +0.34.
- The head group changes V_max by **−28.3**: the alcohol's O–H hydrogen
  (V_max +38 to +52) has no counterpart in the acetate (+19 to +21).

The partition sensitivity at these sites is ≤ 1.8, below every separated
difference.

## Reading

1. **E/Z at C3 is not resolved.** The changes at C3 and C4 (+1.8, +2.2) lie
   inside the conformer spread. The one separated site, C11, is eight
   carbons away. A single hit is at chance level, so it gets no
   interpretation.
2. **Moving the double bond is resolved where it moved.** C4-H changes by
   +3.6, above both ranges and the partition sensitivity. In E2, C4 is an
   allylic CH₂; in E3 it is the alkene carbon. C11-H separates here too, and
   both C11 hits involve E3, whose C11 range is narrow (1.43). That looks
   like a feature of E3's 5-conformer sample, not chemistry. Check it with a
   larger ensemble before saying anything about it.
3. **The head group shows in the global V_max, not in a site mean.** O31 is
   not the same site in the two molecules. In the alcohol it is O plus its
   hydroxyl H, which pulls the mean positive. In the acetate it is an ester O
   with no H. Comparing site means across a change in what the site contains
   is not like-for-like; the record's `h_generated_polar` flag marks the
   alcohol's site.
4. **Five conformers are enough to show the spread, not to measure it.** A
   site's range across 5 conformers runs from 0.2 to 10.6 kcal/mol, and the
   sites with folded members are the widest. Class data should use a larger
   ensemble, and every comparison should be shown against it.

## What this decides for the homework

- **The ensemble is mandatory.** Act B's differences sit inside the
  conformer range; a single conformer would let a student read a
  conformational accident as a stereo effect.
- **Act B's expected finding:** E/Z is not resolved on site means at this
  protocol. That is the lesson, and the grading notes already reward saying
  so.
- **Act C's expected finding:** resolved at C4 only, not at C2 or C3.
- **Act A:** compare the head groups by their global extrema and by a site
  defined the same way in all three molecules (the O atom alone, or a caller
  group), not by `O31` across OH and OAc.

## Re-check at 10 conformers (2026-10-23)

The same protocol with K = 10 per molecule, folded conformers still
included. The pilot's 5 were reused and 5 new ones computed on espmap master
1ef8e28. The scripts and data are `pilot/sesiidae_k10.py`,
`pilot/sesiidae_analyse_k10.py` and `pilot/sesiidae_k10.jsonl`.

With 10 vs 10 conformers, two ranges fail to overlap by chance with
probability 2/184,756 ≈ 1×10⁻⁵ per site. That is far stricter than the
pilot's 0.8 %.

| comparison | separated at K = 5 | separated at K = 10 |
| --- | --- | --- |
| E/Z at C3 | C11-H | **none** |
| C=C position | C4-H, C11-H | **C4-H only: +4.20** (ranges 4.07 / 2.03, partition 1.43) |
| head group | none | none (global V_max −27.5) |

- **The C11 calls were chance**, as suspected. Both disappear with more
  conformers.
- **C4 holds.** It is the one site where moving the double bond moves the
  electrostatics by more than the shape does. The margin is small, though:
  +4.20 against a 4.07 range.
- **The ranges grew.** Single sites now span up to 10.2 kcal/mol across 10
  conformers, against at most 10.6 across 5. More sampling did not narrow
  the spread, which is the case for ensembles in class.

## Not settled by the pilot

- Gas phase and def2-SVP only. Basis moves C18 V_max by +1.5 to +2.8
  (espmap README), which is the same order as the act B differences.
- 5 conformers per molecule. The separated calls should be repeated with a
  larger ensemble (10+, folded included) before class.
- The aldehyde (act A) was not in the pilot.
