# Changelog

What each Arrowfield release contains, and what it does not. Arrowfield
teaches from espmap records and Atlas evidence; it computes nothing itself
(`AGENTS.md`).

## v0.1.0: the teaching platform (2026-11-08)

The first tagged release. espmap is pinned to **v0.3.1** (5f6eaad).

**Runs from a checkout:**
- **The data model and views** (`model.py`, `record.py`, `views.py`): 2D
  structure beside the 3D surface, site highlighting under both partition
  rules, the difference card, the study card, the evidence views and the
  pair view, all in the house style.
- **Hostile pedagogy sets** under `tests/pedagogy/`, written before the code
  they test, with held-out probe rounds behind each view. The full suite
  runs with `python -m pytest -q tests`.
- **The Sesiidae homework page builder** (`homework.page_html`), with
  student and instructor pages, tested against real Atlas data (bundle
  35112e1) and the frozen study espmap-sesiidae-001. It ships with the
  study's records and the aldehyde control (`controls/`).
- **Small-molecule demos:** `python build_demo.py` with espmap v0.3.1 on
  `ESPMAP_PATH`. The 13 demo records in `records/` were computed with
  v0.3.1.
- **The curated manifest** (`molecules.yaml`, regenerated and checked by
  `build_manifest.py`).

**Provenance kept:** the Sesiidae study's records and the aldehyde control
were computed with espmap v0.2.1 under the study's key d0f5ed13. They keep
that producer: the pin move does not recompute study data (guarded by
`tests/test_espmap_pin.py`).

**Not in v0.1.0:**
- **No one-command render of the Sesiidae homework for a class.** A fresh
  checkout does not yet fetch the exact Atlas inputs and NPZ surfaces it
  needs. That is v0.2.0.
- **The homework's observations are not all human-checked** against their
  papers, which is a hard gate before class (`sesiidae-homework.md`).
- **The two notebooks are outlines** (`notebooks.md`), not runnable
  notebooks.
- **The CYP450 case is planning only.** No CYP450 code is in this release.

## Next

- **v0.1.x, CYP450 slice 1 preview (planned):** Compound I bookkeeping, with
  the student page showing no reveal content and no espmap. Its gates are
  the slice's own: the decisions, the expectations-first test set, held-out
  probes and a reviewed render.
- **v0.2.0, the Sesiidae homework, reproducible:** a fresh checkout fetches
  the frozen study's inputs and renders both pages, with the observations
  checked by a person.
