# Contributing to Arrowfield

Arrowfield teaches organic chemistry and Sesiidae pheromone reasoning from
computed ESP surfaces. Its rules for contributors (and for AI agents) are in
`AGENTS.md`; this file is the workflow.

## Set up

Python 3.12 or newer (the pinned numpy needs it), with the versions CI uses:

```bash
pip install -r requirements-ci.txt
```

That is enough for the manifest and the tests, which CI runs on every push
and pull request (`.github/workflows/ci.yml`). Computing or rendering
records (`compute_records.py`, `build_p3.py`, `build_demo.py`) also needs
espmap at the pinned tag:

```bash
pip install "espmap @ git+https://github.com/Azuromalachite/espmap@v0.2.1"
```

## Where things go

| path | what it is | edit by hand? |
| --- | --- | --- |
| `build_manifest.py` | the curated molecules, their sites, roles and Atlas IDs | yes |
| `molecules.yaml` | generated from `build_manifest.py` | **no** |
| `lessons.md`, `sesiidae-homework.md` | the lessons and the homework | yes |
| `notebooks.md` | cell-by-cell outlines of the notebooks and views | yes |
| `model.py` | Arrowfield's data model; views consume only these objects | yes |
| `record.py` | the only reader of espmap records; maps them into the model | yes, carefully |
| `views.py`, `build_demo.py` | the three views and their rendered demo | yes |
| `records/`, `geometries/` | computed small-molecule records, a DFT geometry | regenerate, don't edit |
| `demo/` | rendered pages | generated, not committed |
| `pilot/` | the Sesiidae pilot's scripts and raw results | add, don't rewrite |
| `requirements-ci.txt` | the pinned versions CI installs | yes, with a rebuilt manifest |
| `.github/` | CI workflow and the pull request template | yes |

## Definition of done

A change is done when all of these hold:

1. `python build_manifest.py molecules.yaml` prints "all checked", and
   regenerating after your commit changes nothing.
2. `python -m pytest -q tests` passes.
3. If `record.py` or `views.py` changed, `python build_demo.py` ran and you
   looked at the pages.
4. Every new chemistry or biology claim has its source: a key from espmap's
   `docs/references.bib` plus a DOI for chemistry, the Atlas teaching bundle
   for Sesiidae biology. Say whether you read the full text.
5. No expectation was changed to make code pass. If one was wrong, it changes
   in its own commit, with the reason.
6. The pull request template is filled in, AI use included.

## AI-assisted work

AI tools are welcome; much of this repository was written with them. The
person opening the pull request is responsible for every line: reviewed,
understood and checked. The template asks how AI was used and how the
result was checked independently (a second reviewer, a test written first, a
source reread). "The model said so" is not a check.

## Changes that reach another repository

espmap computes, the Sesiidae Pheromone Atlas owns the Sesiidae science, and
Arrowfield teaches. If a lesson needs a new espmap feature or a new Atlas
field, open an issue in that repository; don't work around it here.
