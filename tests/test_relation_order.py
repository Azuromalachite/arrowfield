"""model.relation does not depend on how a record orders its atoms or bonds
(the probe on #18, fa0e676: a methyl-H relabelling and a reordered bond list
changed the answer). A record's atoms are permuted at random, its bonds
remapped and listed in a random order, and the relation must not move."""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import model  # noqa: E402
import record  # noqa: E402

PAIR = HERE / "tests" / "data" / "pair"


def _permuted(v, seed):
    """The same molecule and geometry with its atoms in another order: every
    per-atom field reindexed, bonds remapped, shuffled and flipped."""
    rng = np.random.default_rng(seed)
    n = len(v.symbols)
    perm = rng.permutation(n)                     # new index k holds old atom perm[k]
    new_of = np.empty(n, int)
    new_of[perm] = np.arange(n)
    bonds = [(int(new_of[j]), int(new_of[i]), o) for i, j, o in v.bonds]
    bonds = [bonds[k] for k in rng.permutation(len(bonds))]
    atoms = tuple(replace(v.atoms[perm[k]], index=k,
                          parent=None if v.atoms[perm[k]].parent is None else int(new_of[v.atoms[perm[k]].parent]))
                  for k in range(n))
    return replace(v, symbols=tuple(v.symbols[p] for p in perm), xyz=v.xyz[perm], bonds=tuple(bonds), atoms=atoms,
                   mapped_smiles=None)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_one_geometry_in_any_atom_order_is_the_same_geometry(seed):
    e3 = record.load(PAIR / "e3_z13_18_oac" / "conf_02")
    assert model.relation(e3, _permuted(e3, seed)) is model.Relation.SAME_GEOMETRY


@pytest.mark.parametrize("seed", [1, 2])
def test_two_conformers_in_any_atom_order_stay_two_conformers(seed):
    a = record.load(PAIR / "e3_z13_18_oac" / "conf_00")
    b = record.load(PAIR / "e3_z13_18_oac" / "conf_02")
    assert model.relation(_permuted(a, seed), _permuted(b, seed + 10)) is model.Relation.CONFORMERS_OF_ONE_COMPOUND


@pytest.mark.parametrize("seed", [1, 2])
def test_two_compounds_in_any_atom_order_stay_two_compounds(seed):
    e3 = record.load(PAIR / "e3_z13_18_oac" / "conf_02")
    z3 = record.load(PAIR / "z3_z13_18_oac" / "conf_00")
    assert model.relation(_permuted(e3, seed), _permuted(z3, seed)) is model.Relation.TWO_COMPOUNDS
