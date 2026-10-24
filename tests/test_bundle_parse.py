"""bundle.load's parse boundary: a malformed JSON shape is refused with
InvalidBundle, naming where, never a crash on the way to the validator.

Probe round 5 on #13 (1099b7d), part 2: 48 malformed shapes, 20 of which
crashed in bundle.load's own parsing (a mapping indexed where something
else stood, a required key missing, or the error path itself failing),
and one loaded with a coercion (espmap as a string became its letters).
Semantic checks stay in model.bundle_problems; this file covers only the
shapes the parser must refuse before a model object exists.
"""
import copy
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from bundle import InvalidBundle, load  # noqa: E402

RAW = json.loads((HERE / "tests" / "data" / "sesiidae-homework-1.bundle.json").read_text(encoding="utf-8"))
DROP = object()


def _first(kind, key):
    """Index of the first compound / observation carrying `key`."""
    return next(i for i, x in enumerate(RAW[kind]) if x.get(key))


OBS_R = next(i for i, o in enumerate(RAW["observations"]) if (o.get("outcome") or {}).get("response"))
OBS_C = _first("observations", "components")
CMP_O = _first("compounds", "occurrences")
CMP_E = _first("compounds", "evidence_classes")

# (id, path, value): value DROP removes the key
SHAPES = [
    # (a) a mapping expected, something else given
    ("observations_as_mapping", ["observations"], {"x": {}}),
    ("observation_as_string", ["observations", 0], "spa.obs.x"),
    ("compounds_as_mapping", ["compounds"], {"x": {}}),
    ("compound_as_list", ["compounds", 0], ["spa.cmpd.x"]),
    ("comparison_as_string", ["comparisons", 0], "a contrast"),
    ("comparisons_as_mapping", ["comparisons"], {"a": {"case_label": "x"}}),
    ("comparisons_as_empty_mapping", ["comparisons"], {}),
    ("taxon_as_string", ["observations", 0, "taxon"], "Synanthedon tenuis"),
    ("taxon_as_list", ["observations", 0, "taxon"], ["Synanthedon tenuis"]),
    ("observed_contrast_as_list", ["comparisons", 0, "observed_contrast"], ["a", "b"]),
    # (b) a required key missing
    ("taxon_without_accepted_name", ["observations", 0, "taxon", "accepted_name"], DROP),
    ("response_without_measure", ["observations", OBS_R, "outcome", "response", "measure"], DROP),
    ("component_without_atlas_id", ["observations", OBS_C, "components", 0, "atlas_id"], DROP),
    ("component_without_detected", ["observations", OBS_C, "components", 0, "detected"], DROP),
    ("observation_without_statement", ["observations", 0, "statement"], DROP),
    ("observation_without_type", ["observations", 0, "observation_type"], DROP),
    ("occurrence_without_roles", ["compounds", CMP_O, "occurrences", 0, "roles"], DROP),
    ("class_without_observation_ids", ["compounds", CMP_E, "evidence_classes", 0, "observation_ids"], DROP),
    ("compound_without_atlas_id", ["compounds", 0, "atlas_id"], DROP),
    ("comparison_without_case_label", ["comparisons", 0, "case_label"], DROP),
    # (c) the error path itself must not crash
    ("evidence_as_list", ["observations", 0, "evidence"], [{"author_statement": "x"}]),
    ("component_as_list_with_bool", ["observations", OBS_C, "components", 0], [True, "unknown_key"]),
    # a coercion: a string where a list belongs
    ("espmap_as_string", ["espmap"], "none"),
]


def _variant(path, value, tmp_path):
    raw = copy.deepcopy(RAW)
    node = raw
    for key in path[:-1]:
        node = node[key]
    if value is DROP:
        del node[path[-1]]
    else:
        node[path[-1]] = value
    p = tmp_path / "bundle.json"
    p.write_text(json.dumps(raw), encoding="utf-8")
    return p


# measured on 1099b7d: bundle.load indexes before checking shapes and keys
XFAIL = {
    "observations_as_mapping",
    "observation_as_string",
    "compounds_as_mapping",
    "compound_as_list",
    "comparison_as_string",
    "comparisons_as_mapping",
    "comparisons_as_empty_mapping",
    "taxon_as_string",
    "taxon_as_list",
    "observed_contrast_as_list",
    "taxon_without_accepted_name",
    "response_without_measure",
    "component_without_atlas_id",
    "component_without_detected",
    "observation_without_statement",
    "observation_without_type",
    "occurrence_without_roles",
    "class_without_observation_ids",
    "compound_without_atlas_id",
    "comparison_without_case_label",
    "evidence_as_list",
    "component_as_list_with_bool",
    "espmap_as_string",
}


@pytest.mark.parametrize("case_id,path,value", SHAPES, ids=[s[0] for s in SHAPES])
def test_a_malformed_shape_is_refused_not_crashed(case_id, path, value, tmp_path, request):
    if case_id in XFAIL:
        request.applymarker(pytest.mark.xfail(strict=True, reason="parse boundary on 1099b7d"))
    with pytest.raises(InvalidBundle):
        load(_variant(path, value, tmp_path))


def test_the_real_bundle_still_loads():
    load(HERE / "tests" / "data" / "sesiidae-homework-1.bundle.json")
