# Arrowfield teaching layer: draft designs

Status: **draft, 2026-10-23.** Nothing here changes espmap. The repository
opened once espmap v0.2.0 shipped its P3 fixture records and tag; these are
the draft designs and working prototype it was seeded with.

| file | what it is |
| --- | --- |
| `README.md` | this overview: boundary, principles, how to read a result |
| `molecules.yaml` | the curated manifest. Never edit it by hand: `build_manifest.py` generates it and checks every entry with RDKit (full maps, charge, E/Z, no open stereo, site atoms exist) |
| `build_manifest.py` | `python build_manifest.py molecules.yaml` |
| `lessons.md` | the organic-chemistry companion: six lessons following the order of a first course (McMurry) |
| `sesiidae-homework.md` | the class assignment: four acts, data plan, questions, what is still pending |
| `notebooks.md` | cell-by-cell outlines of both notebooks and the three views, plus the prototype's status |
| `record.py`, `views.py`, `build_demo.py`, `compute_records.py`, `make_geometries.py`, `tests/` | the working prototype of the three views; see `notebooks.md` |
| `records/`, `geometries/` | computed records (small molecules), acetate's DFT geometry |
| `demo/` | rendered pages, not committed: `python build_demo.py` writes them |
| `LICENSE` | MIT, as espmap |

## Boundary

- **espmap is the instrument**, pinned to a release (`v0.2.1`, never
  `master`). It computes and records; it knows nothing about lessons.
- **Arrowfield is the teaching layer.** It chooses molecules (from the
  manifest, never from typed SMILES), names reaction sites, asks for
  predictions, and decides how to present uncertainty.
- **The [Sesiidae Pheromone Atlas](https://github.com/Azuromalachite/sesiidae-pheromone-atlas)
  owns the Sesiidae science:** compound identity (`spa.cmpd.*`), the
  biological observations, and the canonical espmap records for those
  compounds (`computed/espmap/`). Arrowfield names an Atlas compound by its ID
  and reads observations and records from a teaching bundle the Atlas
  exports. It never restates them by hand. A molecule Arrowfield constructs
  for a lesson is labelled a teaching control.
- Everything the notebooks need from espmap is in
  [v0.2.1](https://github.com/Azuromalachite/espmap/releases/tag/v0.2.1)
  (41a989d), with one exception: `--require-stereo` (an unspecified E/Z or
  stereocentre becomes a BLOCK). Until it lands, the manifest itself
  guarantees full stereo: `build_manifest.py` refuses an entry with any
  unspecified stereo element.

## Principles

1. **Chemistry before computation.** Every exercise runs
   *predict → reveal → explain → compare a control → state what ESP does not tell you*.
   The prediction is not graded; its job is to make you build an organic
   model before seeing the quantum picture.
2. **Structure beside surface.** A 2D skeletal formula with the same atom
   numbers as the site table sits next to every 3D surface.
3. **One scale, or no comparison.** Two pictures are compared only when they
   share `render_key`; two sets of numbers only when they share
   `analysis_key`. Textbook ESP maps (red = electron-rich, as in espmap) are
   each scaled to their own molecule; learning why that matters is an early
   lesson, not an advanced one.
4. **Significance is shown, not inferred from decimals.** Any difference is
   printed next to the measured partition sensitivity and the conformer range
   at that site (the difference card, `notebooks.md`).
5. **ESP is not reactivity.** A permanent box after every result:

   | you may conclude | you may not conclude from ESP alone |
   | --- | --- |
   | local polarity; the electrostatic environment around atoms | which atom is most nucleophilic or acidic |
   | how related molecules differ, on one scale | which bond breaks, how fast, by which mechanism |
   | a hint for interaction reasoning (H-bonding, ion pairing) | where an enzyme reacts |

## How to read a site result (three tiers, as `views.site_table` shows them)

**Primary:** the **site area-weighted mean**, the number to lead with. A site
is a heavy atom with its hydrogens, or a caller group, never a single atom.

**Secondary:** the two numbers any difference must beat.
- **Partition sensitivity** (`delta`): ≤ 3.0 kcal/mol on small anions,
  ≤ 3.95 on C18 site means.
- **Conformer range:** the median C18 site moves 5.1 kcal/mol across
  conformers with folded ones included, 2.69 over extended ones only. The
  difference card shows both.

**Advanced** (folded away under "Advanced"): positive/negative area, the
5th/95th percentiles, raw extrema, and **single atoms**.
- These are the fragile statistics. The 5th percentile moves up to 8.9
  kcal/mol with the partition rule alone.
- A single atom moves up to 11.9 kcal/mol, and its area up to 5.2×,
  because the two rules trade surface between a heavy atom and its own
  hydrogens (espmap #42).
- Atom values are diagnostic. Assigning surface to an atom is a convention,
  and that is itself a lesson: ESP is a field in space, not a property each
  atom owns.

Source for all figures: espmap `bench/site_partition.jsonl`.

## Sources

- **Chemistry claims** cite keys from espmap's `docs/references.bib`, so both
  projects share one bibliography, e.g. `politzer2002` for the ESP as a guide
  to interactions. A claim with no key there is marked `TODO: cite`.
- **Reaction outcomes** ("alkyl halides alkylate the enolate at C") are cited
  with their conditions (solvent, counter-ion, electrophile) before a lesson
  states them. Ratios depend on all three.
- **McMurry** is cited by chapter and section, `McMurry §x.y` (edition to be
  fixed for the class). The lessons never reproduce its text, figures or
  problems.

## Settings

- Neutral molecules: gas phase is the intrinsic baseline.
- **Anions** (acetate, enolate): PCM water, always. In the gas phase these
  small anions come out unbound (positive HOMO), so their ESP is an artefact
  of the basis set (calibration pilot, 2026-10-22). Compare an acid with its
  conjugate base only in the same solvent setting.
- **Anion geometries are DFT, never SMILES + MMFF.** H8 measured MMFF
  geometries moving small anions by up to 6.5 kcal/mol (espmap README,
  `bench/anion_geometry.jsonl`). The enolate uses espmap's
  `bench/fixtures/enolate.sdf` (B3LYP/def2-TZVP, PCM water). Acetate needs
  one made the same way before lesson 3 runs (marked `TODO` in the manifest).
- Small molecules (≤ 6 heavy atoms) run live in seconds. The C18 Sesiidae
  compounds are to be precomputed as conformer ensembles, kept in the
  Sesiidae Pheromone Atlas (`computed/espmap/`) and loaded from its teaching
  bundle; none is computed yet. A "recompute yourself" option stays
  available.
