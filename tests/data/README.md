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
