"""A record's canonical protocol (model.protocol_problems), checked on its own.

The homework page refuses records under any key but its study's, so its
hostile set no longer reaches a re-keyed protocol's names; the pair view and
the difference card still print them. The names are espmap's xc_key and
basis_key outputs (checked against v0.2.1 on 50 common inputs); the refused
ones are #21's rounds 12 and 13.
"""
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from model import Provenance, analysis_digest, protocol_problems, protocol_words  # noqa: E402

BASE = {"xc": "b3lyp", "basis": "def2svp", "solvent": None, "solvent_esp": None, "spacing_bohr": 0.26,
        "density_fit": True, "geometry_method": "etkdg_mmff94s",
        "surface_definition": {"property": "electron_density", "isovalue_e_bohr3": 0.002}}


def _prov(raw_solvent=None, **fields):
    f = {**BASE, **fields}
    return Provenance("r", analysis_digest(f), "k", xc="b3lyp", basis="def2-svp", solvent=raw_solvent,
                      analysis_fields=json.dumps(f, sort_keys=True))


ESPMAP_NAMES = [("xc", "b3lyp"), ("xc", "m062x"), ("xc", "camb3lyp"), ("xc", "libxc:hyb(0,0,0)+531*1"),
                ("xc", "libxc:hyb(0,0,0)+263*1+267*1"), ("xc", "libxc:hyb(0.2,0.08,0)+1*0.08"),
                ("basis", "def2tzvpp"), ("basis", "631g**"), ("basis", "6311++g(2d,2p)"), ("basis", "631+g*"),
                ("basis", "augccpvtz"), ("basis", "pcseg1")]
REFUSED_NAMES = [("basis", "def2svp, PCM water"), ("xc", "b3lyp/def2tzvp"), ("xc", "None"), ("basis", "None"),
                 ("xc", " "), ("xc", "ㅤ"), ("xc", "b3lyp,def2tzvp"), ("xc", "b3lyp,pbe0"),
                 ("basis", "def2svp,PCM"), ("basis", "def2svp,PCM-water"), ("basis", "def2svp,gas"),
                 ("basis", "def2svp(gas)"), ("xc", ""), ("xc", 3)]


@pytest.mark.parametrize("field,name", ESPMAP_NAMES)
def test_a_name_espmap_writes_is_named(field, name):
    p = _prov(**{field: name})
    assert protocol_problems(p) == []
    assert name in protocol_words(p)[0]


@pytest.mark.parametrize("field,name", REFUSED_NAMES)
def test_a_name_that_would_misread_the_caption_is_refused(field, name):
    assert protocol_problems(_prov(**{field: name}))


def test_the_phase_comes_from_the_protocol_and_must_agree_with_the_settings():
    assert protocol_words(_prov(solvent="water", raw_solvent="water"))[1] == "water"
    assert protocol_problems(_prov(solvent="water", raw_solvent="water")) == []
    assert protocol_problems(_prov(solvent="water", raw_solvent="vacuum"))
    assert protocol_problems(_prov(raw_solvent="water"))
    f = {k: v for k, v in BASE.items() if k != "solvent"}
    assert protocol_problems(Provenance("r", analysis_digest(f), "k", analysis_fields=json.dumps(f, sort_keys=True)))


def test_a_key_that_does_not_digest_its_fields_is_refused():
    p = _prov()
    assert protocol_problems(Provenance("r", "0" * 40, "k", analysis_fields=p.analysis_fields))
    assert protocol_problems(Provenance("r", p.analysis_key, "k", analysis_fields=None))
