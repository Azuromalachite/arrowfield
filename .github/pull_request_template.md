## What and why

<!-- What changed, and the problem it solves. Link the issue if there is one. -->

## Teaching or scientific contract

- [ ] No change to what a lesson claims, which molecules or sites it uses, or how records are read
- [ ] Changed — described here:

## Sources

<!-- Every new chemistry claim: espmap references.bib key + DOI, and whether you read the full text.
     Every Sesiidae biology claim: the Atlas teaching bundle (version) it comes from. "None" is a fine answer. -->

## Validation

- [ ] `python build_manifest.py molecules.yaml` prints "all checked"; regenerating changes nothing
- [ ] `python -m pytest -q tests` passes
- [ ] `python build_demo.py` ran and I looked at the pages (only if `record.py` or `views.py` changed)
- [ ] No expectation was changed to make code pass

## Known limits

<!-- What this does not do yet, or cannot conclude. -->

## AI use

<!-- Was AI used? For what? How was the result checked independently
     (a reviewer, a test written first, a source reread)? -->
