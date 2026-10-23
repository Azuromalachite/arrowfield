# Notebooks and views: outline

Two notebooks on top of espmap pinned to a release (`ESPMAP_REF = "v0.2.0"`,
never `master`), and three views they share. Neither notebook reimplements
anything espmap computes.

**Record fields are named loosely here, on purpose.** Adversarial review 4 is
auditing espmap's JSON record against the spec now, and may rename or add
fields before v0.2.0. Below, "the site block", "a site's mean", "the
partition sensitivity" and so on refer to concepts, not key paths. The key
paths are fixed in one small adapter (`record.py`, below) once v0.2.0 is
tagged, and nowhere else.

## Prototype (2026-10-22): built and tested

| file | what it does |
| --- | --- |
| `record.py` | the only module that knows espmap's JSON/NPZ layout; refuses a record whose `validity` is false |
| `views.py` | `svg_2d`, `viewer_js` (3Dmol.js, house style), `difference_card`, `page` |
| `compute_records.py` | runs espmap's CLI on manifest molecules (mapped SMILES or the DFT geometry, sites as `--sites`, 3 seeds) into `records/` |
| `make_geometries.py` | DFT (PCM water) geometries for anions without an espmap fixture: `geometries/acetate.sdf` |
| `build_demo.py` | writes `demo/acetaldehyde.html`, `demo/acetate.html`, `demo/difference_card.html` |
| `tests/test_views.py` | 7 checks: invalid records refused, site numbers from the record, the patch area is the site's area, house style present, card verdicts, incomparable records refused |

**House style** (from the reference figure): white background; thin black
outline (`setViewStyle outline, width 0.03`); O red, N blue, S yellow, P
orange; carbons muted blue, site carbons lilac; hydrogens hidden except the
site's own; a translucent lilac halo on site atoms, so an O or N site, which
keeps its element colour, is still visibly selected. Molecules load as mol
blocks so double bonds show. The camera is fixed to the molecule's principal
axes, so every student starts from the same view.

**Static figures** (handouts, slides): ChimeraX on the desktop, using the
other session's tested recipe (simple lighting, silhouettes width ~3, ortho
camera, role colour then `color byhetero`, fixed `turn`). Labels and layout
come in a 2D step. ChimeraX's stick style draws no bond orders, so for
organic teaching the interactive 3Dmol.js view is the one that shows C=O.

**Seen in the prototype:**
- E- vs Z-but-2-ene differ by 1.29 kcal/mol at the alkene carbon, while
  three seeds of each vary by only 0.01–0.03 (a rigid molecule has one
  conformer). The card says "separated". The C18 homework is the opposite
  case, and the same card makes that visible.
- Surfaces from different auto-scaled runs carry different `render_key`s, and
  the page says so in words above the card: compare numbers, not colours.

## Shared pieces

### `record.py`: the only code that knows espmap's JSON layout

- Load a record; refuse it unless the record says it is **valid** for input,
  surface and sites (one check, read from the record, never inferred).
- Expose: molecule identity, the two comparability keys, the global surface
  statistics, and per site (atoms or caller groups): area, ± areas, mean,
  percentiles, extrema, and the difference between the two partition rules.
- `comparable(a, b)` is true only when the analysis keys match. `same_scale(a, b)`
  is true only when the render keys match.

### View 1: 2D structure beside the 3D surface

- RDKit draws the skeletal formula from the manifest's mapped SMILES with
  atom numbers shown. For the C18s, C1 is the carbon bearing O and numbering
  follows the chain.
- The 3D surface, from espmap's NPZ, sits next to it at the same size.
- Clicking a row of the site table highlights that site in both panels.

### View 2: site highlighting

- espmap's NPZ records which atom owns each triangle, under each partition
  rule. A site is a union of atoms, so its triangles are known with no new
  computation.
- The selected site stays opaque, everything else fades to about 30 %, and
  the site boundary is outlined.
- A toggle switches the partition rule, so the student sees the boundary move.
  This shows the partition sensitivity directly.

### View 3: the difference card

For two records with the **same analysis key** (otherwise it refuses), one
card per site:

```
E3 → Z3 acetate · site C3-H                 site mean, kcal/mol
  difference of conformer means        +0.8
  conformer range (E3 / Z3)             4.6 / 4.1
  partition sensitivity (max)           0.7
  verdict   smaller than the conformer range: not attributable to E/Z
```

(The numbers above are placeholders for layout, not results.)

Verdict rules, deliberately simple:

- **separated** when the two molecules' conformer ranges do not overlap, so
  every conformer of one differs from every conformer of the other in the
  same direction;
- otherwise **"not attributable to the change"**;
- with a single conformer per molecule: **"single conformer: cannot judge"**,
  never a verdict.

Underneath it always shows the reading order: mean, ± area, partition, then
conformers, with the 5th percentile and extrema behind a "more" fold.

## Notebook 1: organic-chemistry companion

One notebook, one lesson per section (`lessons.md`). Each section:

1. **Pick** the lesson's molecules from the manifest. A dropdown shows the
   name, the 2D structure and the SMILES underneath, for learning.
2. **Predict**: multiple choice from the lesson, recorded and not graded. The
   reveal stays hidden until something is picked.
3. **Run** live: small molecules, a few seconds. espmap at the class basis.
   Anions in PCM water on their DFT fixture geometry, which the manifest says.
4. **Reveal**: View 1 + View 2, the site table in the reading order, and the
   scale caption ("darker = more negative" on ions).
5. **Explain**: the lesson's reasoning, with the McMurry section reference.
6. **Control**: the control molecule on the same scale (same render key).
7. **Limits**: the lesson's "what this does not tell you" box, plus the
   permanent table from `README.md`.

## Notebook 2: Sesiidae homework

1. **Setup:** pin espmap, load the precomputed records and NPZs for the five
   compounds and their conformers. No SCF unless the student opts in.
2. **Act A:** dropdown alcohol / aldehyde / acetate; predict; View 1 + 2;
   difference cards for the head-group sites and the chain sites; the
   student writes the explanation.
3. **Act B:** E3 vs Z3; predict; difference cards at C3-H, C4-H and head O;
   the extended and folded conformers of one compound side by side.
4. **Act C:** E2 vs E3; predict; difference cards at C2-H to C4-H.
5. **Act D:** free text: what is missing.
6. **Submit:** export the student's predictions, chosen cards and answers,
   each card carrying its records' analysis keys and espmap version, so a
   grader can see every comparison was legitimate.

## What this needs from espmap

| need | status on master 6c66bbf |
| --- | --- |
| per-site statistics under two partition rules | present (P2) |
| caller-defined sites by atom map | present (`--sites`) |
| triangle → atom ownership in the NPZ | present |
| comparability keys and validity in the record | present (names may change in review 4) |
| **unspecified E/Z refused** | **missing**: `--require-stereo`; until then the manifest guarantees it |
| precomputed records for the five C18 compounds | P3, with the fixture set extended |
| a tagged release to pin | v0.2.0, after P3 |
