# Teaching controls

Computed molecules a lesson needs that are **not** in the Sesiidae
Pheromone Atlas: constructed so that a comparison changes one thing, and
never presented as observed. The Atlas holds only its own compounds, so
their records live here, in the same layout and index format the Atlas uses
for its `computed/espmap/` (`spa.espmap_ensemble/0`), and
`tests/test_controls.py` checks them as the Atlas gate would.

| control | why | study |
| --- | --- | --- |
| `arrowfield.control.z3_z13_18_al`, (3Z,13Z)-octadeca-3,13-dienal | holds the Z3,Z13 backbone fixed while only the head group changes (act A); the Atlas's aldehyde is E2,Z13-18:Al, which is not in the case | espmap-sesiidae-001 (Atlas `analyses/espmap-sesiidae-001/`, frozen at 7e34658) |

Layout: `espmap/<slug>/<protocol>/{ensemble.json, SHA256SUMS,
conformers/conf_NN.json}`. Twenty members, computed with espmap v0.2.1
(41a989d) under the study's protocol: B3LYP/def2-SVP + RI, gas phase,
0.26 bohr, seed 20261024. NPZ surfaces are not in git; each has its sha256
in `ensemble.json` and ships as a release asset. How the conformers were
chosen and run (selection, run log, folded calls) is recorded with the
study in the Atlas, at `analyses/espmap-sesiidae-001/selection/z3_z13_18_al/`.
