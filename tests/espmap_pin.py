"""The espmap release Arrowfield is moving its pin to, for the tests only
(the user, 2026-11-07: the whole repo moves from v0.2.1 to v0.3.1).

One statement, imported by the tests that check the pin: the pin move itself
(docs, vendored tables, regenerated records) is a later change, and these
tests say what it must deliver. The Sesiidae study's records and the
aldehyde control are study data under their own producer (v0.2.1, key
d0f5ed13) and never move with the pin.
"""
PIN = {"tag": "v0.3.1", "version": "0.3.1", "commit": "5f6eaade3c0d60af921808cc08e7225288c96e36"}
# what the study's data were computed with, and stay computed with
STUDY_PRODUCER = {"version": "0.2.1", "commit": "41a989d6240b56e3c17938dac2084dfc5e1423e4"}
