"""The pair view's pictures: both panels colour through one function, on the
shared scale, from the vendored espmap v0.2.1 tables (the user, 2026-10-24:
"identical scalar values on the two panels receive identical colours under
the shared scale"). pair_view_set.json checks the page's attributes; this
checks the colours the pictures are drawn with."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import model  # noqa: E402
import palette  # noqa: E402
import record  # noqa: E402
import views  # noqa: E402

PAIR = HERE / "tests" / "data" / "pair"


def _captured(left, right, monkeypatch):
    seen = []
    real = views.viewer_js

    def spy(div_id, rec, mapped_smiles, *a, colour=None, **kw):
        seen.append((rec, colour))
        return real(div_id, rec, mapped_smiles, *a, colour=colour, **kw)
    monkeypatch.setattr(views, "viewer_js", spy)
    views.pair_view_html(left, right, {})
    return seen


def test_both_panels_colour_one_value_alike(monkeypatch):
    e3f = record.load(PAIR / "e3_z13_18_oac" / "conf_00")
    e3x = record.load(PAIR / "e3_z13_18_oac" / "conf_02")
    seen = _captured(e3f, e3x, monkeypatch)
    assert len(seen) == 2 and all(c is not None for _, c in seen)
    probe = np.array([-40.77, -12.0, 0.0, 7.5, 19.8, 40.77])
    a, b = (c(probe) for _, c in seen)
    assert a == b                                      # one value, one colour, on either panel
    s = model.shared_scale(e3f, e3x)
    assert a == palette.hex_colours(probe, s.vmin, s.vmax, "vik")
    assert a[0] != a[-1] and a[2] == palette.hex_colours([0.0], s.vmin, s.vmax, "vik")[0]


def test_a_panel_is_not_coloured_on_its_own_record_scale(monkeypatch):
    """The two records' own scales differ (+-40.33 and +-40.77): a value near
    the end takes the shared scale's colour, not its record's."""
    e3f = record.load(PAIR / "e3_z13_18_oac" / "conf_00")
    e3x = record.load(PAIR / "e3_z13_18_oac" / "conf_02")
    (r0, c0), _ = _captured(e3f, e3x, monkeypatch)
    s = model.shared_scale(e3f, e3x)
    grid = np.linspace(-40.0, 40.0, 801)
    own = palette.hex_colours(grid, r0.scale.vmin, r0.scale.vmax, "vik")
    shared = palette.hex_colours(grid, s.vmin, s.vmax, "vik")
    differ = [i for i in range(len(grid)) if own[i] != shared[i]]
    assert differ, "premise: the two scales must colour some value differently"
    v = grid[differ]
    assert c0(v) == [shared[i] for i in differ]       # the shared scale's colours, never the record's
