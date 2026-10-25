"""The vendored espmap v0.2.1 colour tables and palette's lookup, tested
without espmap (CI), plus an exact parity test where espmap at v0.2.1 is
available (ESPMAP_PATH at commit 41a989d)."""
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import palette  # noqa: E402

T = palette.TABLES
COMMIT = "41a989d6240b56e3c17938dac2084dfc5e1423e4"


def _sha(rgba):
    return hashlib.sha256(json.dumps(rgba, separators=(",", ":")).encode()).hexdigest()


def _lum(c):
    return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]


def test_the_tables_are_espmaps_and_unedited():
    assert T["espmap"] == {"version": "0.2.1", "commit": COMMIT, "function": "esp_map.build_maps()[name]"}
    assert sorted(T["maps"]) == ["lajolla", "vik"]
    for name, m in T["maps"].items():
        assert m["N"] == 256 and len(m["rgba"]) == 256, name
        assert m["sha256"] == _sha(m["rgba"]), f"{name}: the table was edited"
        assert all(len(c) == 4 and all(0.0 <= v <= 1.0 for v in c) for c in m["rgba"]), name
    assert "Fabio Crameri" in T["attribution"] and (HERE / "vendor" / "LICENSE-cmcrameri").is_file()


def test_vik_is_warm_at_the_negative_end():
    vik = T["maps"]["vik"]["rgba"]
    assert vik[0][0] > vik[0][2]                 # t = 0, V at vmin (negative): red over blue
    assert vik[-1][2] > vik[-1][0]               # t = 1, positive: blue over red
    assert all(v > 0.8 for v in vik[128][:3])    # the middle is light: zero is not a colour


def test_lajolla_is_darkest_at_v_min():
    laj = [_lum(c) for c in T["maps"]["lajolla"]["rgba"]]
    assert laj[0] == min(laj) and laj[-1] == max(laj)


def test_the_lookup_is_espmaps_rule():
    for name in ("vik", "lajolla"):
        assert palette.index(name, 0.0) == 0 and palette.index(name, 1.0) == 255
        assert palette.index(name, 0.5) == 128
        assert palette.index(name, -3.0) == 0 and palette.index(name, 7.0) == 255     # clipped
        for k in (0, 1, 127, 254):
            assert palette.index(name, (k + 0.999) / 256) == k                       # truncation, not rounding
        for bad in (math.nan, math.inf, True, "0.5"):
            with pytest.raises(ValueError):
                palette.index(name, bad)
    with pytest.raises(ValueError):
        palette.index("bwr_r", 0.5)                                                 # no frozen table, no guess


def test_equal_values_on_one_scale_get_one_colour():
    a = palette.hex_colours([-12.5, 0.0, 30.0], -40.0, 40.0, "vik")
    b = palette.hex_colours([-12.5, 0.0, 30.0], -40.0, 40.0, "vik")
    assert a == b and len(set(a)) == 3
    assert palette.hex_colours([-50.0], -40.0, 40.0, "vik") == palette.hex_colours([-40.0], -40.0, 40.0, "vik")
    with pytest.raises(ValueError):
        palette.hex_colours([1.0], 0.0, 0.0, "vik")                                # a degenerate scale


def test_palette_never_reorients_a_table():
    """No code in palette.py reverses, flips or reverse-slices: checked on
    its syntax tree, not its prose (the docstring says it never does)."""
    import ast
    tree = ast.parse((HERE / "palette.py").read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        name = (node.id if isinstance(node, ast.Name) else node.attr if isinstance(node, ast.Attribute) else "")
        assert not name.startswith(("reversed", "flip")), name
        if isinstance(node, ast.Slice) and node.step is not None:
            step = node.step
            assert not (isinstance(step, ast.UnaryOp) and isinstance(step.op, ast.USub)), "a negative-step slice"


def _espmap_at_v021():
    root = os.environ.get("ESPMAP_PATH")
    if not root or not (Path(root) / "esp_map.py").is_file():
        return None
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True).stdout.strip()
    return root if head == COMMIT else None


def test_parity_with_espmap_v021_where_available():
    """Exact equality with the pinned source only: a newer espmap is never
    authoritative; a difference there is drift to report, not to follow."""
    root = _espmap_at_v021()
    if root is None:
        pytest.skip("espmap at v0.2.1 (41a989d) is not available here")
    sys.path.insert(0, root)
    import esp_map as em
    built = em.build_maps()
    for name, m in T["maps"].items():
        assert [[float(v) for v in built[name](i / 256)] for i in range(256)] == m["rgba"], name
