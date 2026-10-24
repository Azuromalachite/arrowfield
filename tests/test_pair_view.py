"""The pair view's hostile set (tests/pedagogy/pair_view_set.json), written
before the view exists.

Now: the set is what its generator writes; every case's `allowed` is exactly
its safe outputs; all six classes have two or more variants; the fixtures in
tests/data/pair are the Atlas's records, byte for byte; and the facts each
case's reason cites are true of the records (a case that misdescribes its
input would test the wrong thing).
When views.pair_view_html exists: every output is safe and preferred, and
every requirement holds. Until then those tests skip. The view's PR may add
strict-xfail lists here, never edit an expectation.
"""
import hashlib
import json
import math
import re
import sys
from dataclasses import replace
from html.parser import HTMLParser
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
PED = Path(__file__).with_name("pedagogy")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PED))
import record  # noqa: E402

SET = json.loads((PED / "pair_view_set.json").read_text(encoding="utf-8"))
CASES = SET["cases"]
PAIR = HERE / "tests" / "data" / "pair"


def test_the_committed_set_is_what_the_generator_writes():
    import make_pair_view_set
    assert make_pair_view_set.document() == (PED / "pair_view_set.json").read_text(encoding="utf-8")


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_allowed_is_exactly_the_safe_outputs(c):
    assert c["allowed"] == sorted(k for k, g in c["outcomes"].items() if g in ("SAFE", "SAFE_REFUSAL"))
    assert c["prefer"] is None or c["prefer"] in c["allowed"]
    assert c["why"] and c["outcomes"]["ERROR"] == "CRASH"


def test_every_class_has_two_or_more_variants():
    from collections import Counter
    per = Counter(c["cls"].split()[0] for c in CASES)
    assert sorted(per) == [f"P{i}" for i in range(1, 7)] and min(per.values()) >= 2, per


def test_retired_cases_name_existing_replacements():
    ids = {c["id"] for c in CASES}
    old = {c["id"] for c in json.loads((PED / "hostile_set.json").read_text(encoding="utf-8"))["cases"]
           if c["component"] == "pair_view"}
    assert set(SET["retires"]) == old
    assert all(r in ids for rs in SET["retires"].values() for r in rs)


# --- the fixtures, and the facts the cases rest on -------------------------------------------
def _sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_pair_fixtures_are_the_atlas_records():
    members = json.loads((PAIR / "members.json").read_text(encoding="utf-8"))["members"]
    assert len(members) == 4
    for m in members:
        j = PAIR / m["record"]
        assert _sha(j) == m["record_sha256"] and _sha(j.with_suffix(".npz")) == m["npz_sha256"], m["record"]
        z = np.load(j.with_suffix(".npz"))
        r = json.loads(j.read_text(encoding="utf-8"))
        assert str(z["analysis_key"]) == r["analysis_key"] and str(z["render_key"]) == r["render_key"]


def _load(stem):
    return record.load(HERE / stem)


def _dist(v):
    return np.linalg.norm(v.xyz[:, None] - v.xyz[None], axis=-1)


def test_the_records_say_what_the_cases_say():
    import make_pair_view_set as g
    v = {s: _load(s) for s in (g.ETOH, g.ETOH2, g.ACAL, g.ACAL2, g.E2B, g.E2B2, g.Z2B, g.ACO,
                               g.E3F, g.E3X, g.Z3X, g.Z3F)}
    # P1: same method, different own scales
    assert v[g.E2B].analysis_key == v[g.Z2B].analysis_key and v[g.E2B].render_key != v[g.Z2B].render_key
    assert round(v[g.E2B].scale.vmax, 1) == 17.3 and round(v[g.Z2B].scale.vmax, 1) == 18.1
    assert round(v[g.ETOH].scale.vmax, 1) == 50.0 and round(v[g.ACAL].scale.vmax, 1) == 39.3
    assert v[g.ETOH].analysis_key == v[g.ACAL].analysis_key
    assert (round(v[g.E3F].scale.vmax, 2), round(v[g.E3X].scale.vmax, 2)) == (40.33, 40.77)
    assert v[g.E2B].render_key == v[g.E2B2].render_key
    # P2: the keys that differ, and the ones that do not
    assert v[g.ACO].analysis_key != v[g.ETOH].analysis_key and v[g.ACO].charge == -1
    assert v[g.E3X].analysis_key != v[g.ETOH].analysis_key
    assert v[g.E3X].provenance.method == v[g.ETOH].provenance.method == "b3lyp/def2-svp, gas"
    # P3: same geometry three ways, and true conformers
    assert np.array_equal(v[g.E2B].xyz, v[g.E2B2].xyz)
    assert np.abs(_dist(v[g.ACAL]) - _dist(v[g.ACAL2])).max() < 1e-5
    assert np.abs(v[g.ACAL].xyz - v[g.ACAL2].xyz).max() > 2.0
    assert round(float(np.abs(_dist(v[g.ETOH]) - _dist(v[g.ETOH2])).max()), 2) == 0.70
    assert np.abs(_dist(v[g.E3F]) - _dist(v[g.E3X])).max() > 1.0
    # P4: the producer's flags, and genus that does not follow them
    folded = {m["record"][:-5]: m["folded"] for m in
              json.loads((PAIR / "members.json").read_text(encoding="utf-8"))["members"]}
    assert folded == {"e3_z13_18_oac/conf_00": True, "e3_z13_18_oac/conf_02": False,
                      "z3_z13_18_oac/conf_00": False, "z3_z13_18_oac/conf_01": True}
    assert v[g.E3F].genus == 0 and v[g.Z3F].genus == 3
    assert all(v[s].folded is None for s in (g.E2B, g.Z2B))
    # the loader names both compounds' conformers alike
    assert v[g.E3X].name.startswith("conf_0") and v[g.Z3X].name.startswith("conf_0")


# --- the view (to come) --------------------------------------------------------------------
def _view():
    import views
    if not hasattr(views, "pair_view_html"):
        pytest.skip("views.pair_view_html is not built yet; its expectations are fixed in pair_view_set.json")
    return views.pair_view_html


def _side(s):
    from model import Surface
    v = _load(s["record"])
    for field, value in s["patch"]:
        if "." in field:
            outer, inner = field.split(".")
            v = replace(v, **{outer: replace(getattr(v, outer), **{inner: value})})
        else:
            v = replace(v, **{field: value})
    op = s["surface"]
    if op == "none":
        v = replace(v, surface_loader=None)
    elif op is not None:
        base = v.surface_loader

        def loader(rule, base=base, op=op):
            x = base(rule)
            V, faces = x.V.copy(), x.faces
            if op == "nan":
                V[0] = math.nan
            elif op == "inf":
                V[0] = math.inf
            elif op == "empty":
                faces = np.zeros((0, 3), dtype=x.faces.dtype)
            return Surface(x.verts, faces, V, x.owner[:len(faces)])
        v = replace(v, surface_loader=loader)
    return v


def _run(c):
    view = _view()
    from model import Refused
    left = _side(c["left"])
    right = left if c["same_object"] else _side(c["right"])
    try:
        page = view(left, right, dict(c["context"]))
    except Refused:
        return "REFUSED", None, (left, right)
    except Exception:  # noqa: BLE001 - any other error is a crash, graded as such
        return "ERROR", None, (left, right)
    return "RENDERED", page, (left, right)


@pytest.mark.parametrize("c", CASES, ids=[c["id"] for c in CASES])
def test_view_output_is_safe(c):
    out = _run(c)[0]
    assert out in c["allowed"], f"{out}: {c['outcomes'][out]}"


@pytest.mark.parametrize("c", [c for c in CASES if c["prefer"]], ids=lambda c: c["id"])
def test_view_output_is_preferred(c):
    assert _run(c)[0] == c["prefer"]


class _Attrs(HTMLParser):
    def __init__(self):
        super().__init__()
        self.found = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if any(k.startswith("data-") for k in a):
            self.found.append(a)


def _check(what, page, sides):
    p = _Attrs()
    p.feed(page)
    classes = [(a.get("class") or "").split() for a in p.found]
    scale = [a for a, k in zip(p.found, classes) if "pair-scale" in k]
    panels = [a for a, k in zip(p.found, classes) if "pair-panel" in k]
    if what == "shared_scale":
        assert len(scale) == 1 and len(panels) == 2, "one pair-scale and two pair-panels"
        s = scale[0]
        lo, hi = float(s["data-vmin"]), float(s["data-vmax"])
        assert math.isfinite(lo) and math.isfinite(hi) and lo < hi
        for q in panels:
            assert math.isclose(float(q["data-vmin"]), lo, abs_tol=1e-9)
            assert math.isclose(float(q["data-vmax"]), hi, abs_tol=1e-9)
        assert lo <= min(v.scale.vmin for v in sides) + 1e-6 and hi >= max(v.scale.vmax for v in sides) - 1e-6
        if any(v.scale.kind == "diverging" for v in sides):
            assert s["data-kind"] == "diverging" and math.isclose(lo, -hi, rel_tol=1e-12)
        else:
            assert s["data-kind"] == "sequential", "neither side changes sign: no zero to centre on"
        return
    if what.startswith("relation="):
        rel = [a["data-relation"] for a in p.found if "data-relation" in a]
        assert rel == [what.split("=")[1]], rel
        return
    want = what.split("=")[1].split(",")
    assert [q.get("data-folded") for q in panels] == want


REQ = [(c, i) for c in CASES for i in range(len(c["requires"]))]


@pytest.mark.parametrize("c,i", REQ, ids=[f"{c['id']}-{i}" for c, i in REQ])
def test_a_rendered_view_meets_its_requirements(c, i):
    out, page, sides = _run(c)
    if out != "RENDERED":
        return
    q = c["requires"][i]
    if "check" in q:
        _check(q["check"], page, sides)
        return
    target = page if q["on"] == "html" else re.sub(r"<[^>]+>", " ", page)
    if "must" in q:
        assert re.search(q["must"], target, re.I), f"{q['why']} ({q['grade']} if missing)"
    else:
        assert not re.search(q["must_not"], target, re.I), f"{q['why']} ({q['grade']} if present)"
