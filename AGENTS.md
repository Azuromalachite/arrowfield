# AGENTS.md — rules for anyone changing Arrowfield, human or AI

Arrowfield is the **teaching layer**. It shows espmap surfaces, asks for
predictions and explains what can and cannot be concluded. It owns pedagogy
only. Read `README.md` for the boundary; this file is the short, binding
version of it.

## Who owns which fact

| fact | owner | how Arrowfield gets it |
| --- | --- | --- |
| an ESP number, surface, site statistic | espmap (pinned tag, never `master`) | a record, read through `record.py` only |
| a Sesiidae compound's identity | Sesiidae Pheromone Atlas (`spa.cmpd.*`) | `atlas_compound` in the manifest, InChIKey checked by `build_manifest.py` |
| who produces or responds to what, under which conditions | the Atlas | its exported teaching bundle (`spa.teaching/0`) |
| a reaction outcome in a lesson | the primary literature | a citation key from espmap's `docs/references.bib`, with its DOI and conditions |
| which molecule, which site, which question, in what order | **Arrowfield** | `build_manifest.py`, `lessons.md`, `sesiidae-homework.md` |

## Never

- Compute or "correct" an espmap quantity here. Read it from a record, or
  don't show it.
- Read espmap's JSON/NPZ layout anywhere except `record.py`.
- Write a biological fact by hand. If the teaching bundle doesn't carry it,
  the lesson doesn't say it.
- State a Sesiidae compound's identity apart from its Atlas ID, or show a
  teaching control (`role: teaching_control`) without saying so.
- Edit `molecules.yaml` by hand. Change `build_manifest.py` and regenerate.
- Compare two records whose `analysis_key` differs, or two images whose
  `render_key` differs.
- Quote a difference without the uncertainty it must beat: the partition
  delta and the conformer range.
- State a reaction outcome without its conditions and a source, or
  reproduce a textbook's text, figures or problems.
- Pin espmap to `master`, or the Atlas to anything but a recorded commit or
  tag.
- Change an expectation (a lesson's "answer", a test's expected value) to
  make new code pass. Argue for the change in its own commit.

## Every lesson has the same shape

**Predict → Reveal → Explain → Control → Limits.** The prediction is
recorded, never graded. The Limits section says what this picture cannot
tell you.

## Before you call a change done

```bash
python build_manifest.py molecules.yaml   # must print "all checked"; commit the result
git diff --exit-code molecules.yaml       # after committing: regenerating changes nothing
python -m pytest -q tests
python build_demo.py                      # when record.py or views.py changed; open the pages
```

The definition of done, in full, is in `CONTRIBUTING.md`.

## Working style

- Respect human edits: never rewrite text you weren't asked to change.
- One owner per change across the three repositories. A change that needs
  espmap or the Atlas is proposed to whoever maintains that repository, not
  made here.
- Commit subjects are imperative and say what changed; the body says why.
  AI-assisted commits carry a `Co-Authored-By` line.
