# Test data from the Sesiidae Pheromone Atlas

`sesiidae-homework-1.bundle.json` is the Atlas teaching bundle for the case
`sesiidae-homework-1`, exported unchanged by the Atlas's
`scripts/export_teaching_bundle.py` at the commit named in its
`atlas_version`. It is here so `tests/test_bundle.py` reads a real bundle.

The Atlas's records are © 2026 Thanakrit (Austin) Wongsatit, licensed
CC BY 4.0 (the Atlas's `LICENSE-DATA.md`): Sesiidae Pheromone Atlas,
https://github.com/Azuromalachite/sesiidae-pheromone-atlas. Each
observation's `author_statement` is a short quotation from the cited paper.
The Atlas's licence does not cover it, and it stays under its source's terms.

To refresh it, re-export from a newer Atlas commit and replace the file
whole. Never edit it by hand: the tests compare it with what the Atlas
states.

`espmap-sesiidae-001.results.json` is the frozen study's results file,
copied unedited from the Atlas at dabaf15 (main, after PR #1)
(`analyses/espmap-sesiidae-001/results.json`), CC BY 4.0 like the bundle.
`tests/test_study.py` reads it, together with the aldehyde control's records
in `controls/`.

`pair/` holds four conformer records of the Atlas's espmap ensembles, with
their NPZ surfaces, for the pair view's hostile set
(`tests/pedagogy/pair_view_set.json`): the E3 and Z3 acetates, one folded
and one extended conformer each. The records are copied unedited from the
Atlas at 35112e1 (`computed/espmap/<compound>/b3lyp_def2-svp_gas_espmap-v0.2.1/conformers/`);
the NPZs are the same run's surfaces, which the Atlas will publish as
release assets. `members.json` gives each file's sha256 and the folded flag
from the Atlas's `ensemble.json`; `tests/test_pair_view.py` checks both.
CC BY 4.0, as above.
