"""#13's evidence views and #16's labelled difference card on one page.

Both PRs changed views.py; a clean textual merge proves no line conflict,
not that one view's helpers, wording rules, escaping or context handling
left the other's alone (the user's integration gate). Each view is
rendered into one page, and each keeps its own contract there.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "tests"))
from bundle import load  # noqa: E402

BUNDLE = HERE / "tests" / "data" / "sesiidae-homework-1.bundle.json"
RAW = json.loads(BUNDLE.read_text(encoding="utf-8"))
P = "spa.obs."
OAC, E3 = "spa.cmpd.z3_z13_18_oac", "spa.cmpd.e3_z13_18_oac"
TIPU = "spa.tax.synanthedon_tipuliformis"


def _fields(node, out):
    if isinstance(node, str):
        out.add(node)
    elif isinstance(node, dict):
        for v in node.values():
            _fields(v, out)
    elif isinstance(node, list):
        for v in node:
            _fields(v, out)
    return out


def _render_all():
    import views
    import test_card_label as tl
    from study import load_results
    b = load(BUNDLE)
    ctx = {"bundle_object": b, "roles": {OAC: "atlas_observed", E3: "atlas_observed"}}
    evidence = views.evidence_card_html(b.observations[P + "synanthedon_tenuis.yang_2012.trap_acetate_added"],
                                        OAC, ctx)
    k = b.contrasts[0]
    contrast = views.contrast_card_html(k, [b.observations[o] for o in k.observation_ids], ctx)
    panel = views.evidence_class_panel_html(E3, TIPU, ctx)
    x, y = "z3_z13_18_oac", "e3_z13_18_oac"
    from model import compare
    card = compare(tl._views(x), tl._views(y), "C3-H", x, y)
    diff = views.card_html(card, {"study": load_results(tl.RESULTS), "contrast": "B_stereo", "x": x, "y": y,
                                  "roles": {x: "atlas_observed", y: "atlas_observed"}})
    page = views.page("Together", "evidence views and the difference card", [evidence, contrast, panel, diff])
    return page, evidence, contrast, panel, diff


def test_each_view_keeps_its_contract_on_one_page():
    import test_card_label as tl
    import test_evidence_cards as te
    from study import load_results
    page, evidence, contrast, panel, diff = _render_all()
    for fragment in (evidence, contrast, panel, diff):
        assert fragment in page                       # the page adds nothing inside a view
    # the evidence views: every verbatim span is a record field, exactly; own words clean
    record = _fields(RAW, set())
    for fragment in (evidence, contrast, panel):
        own, spans = te._split(fragment)
        assert spans and all(s in record for s in spans), [s for s in spans if s not in record]
        assert not re.search(r"inhibit\w*|\bresol(v|ut)\w*|because|electrostatic", own, re.I), own
    assert "strongly inhibited attraction" in evidence  # the authors' words, in their span
    # the difference card: its label, the study pointer, no study outcome word in its own words
    own, spans = tl._split(diff)
    assert re.search(r"teaching comparison", own, re.I) and re.search(r"study card", own, re.I)
    assert not re.search(r"\bresol(v|ut)\w*|\bunresol(v|ut)\w*|different on average", own, re.I), own
    assert set(spans) <= {"C3-H", load_results(tl.RESULTS).plan}
    # one escaping rule: no raw apostrophe entity anywhere on the page
    assert "&#x27;" not in page
