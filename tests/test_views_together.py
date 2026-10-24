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


# The integration probe (on 1cbab22): context entries that are present but
# malformed are refused, never read as absent (the user: "a validator, not a
# sanitizer"); an absent entry stays absent, and a whole context that is not
# a mapping still reads as none (#15's choice 6). And card_html refuses an
# object that is not a Comparison, as the evidence views refuse one.
import pytest  # noqa: E402

MALFORMED = ["card_given_a_statement", "card_roles_as_list", "evidence_provenance_as_list", "evidence_roles_as_list"]
XFAIL = set()                # the 4 measured on 1cbab22 are met; only shrinks


def _malformed(case_id):
    import views
    import test_card_label as tl
    from model import compare
    from study import load_results
    b = load(BUNDLE)
    x, y = "z3_z13_18_oac", "z3_z13_18_al"
    if case_id == "card_given_a_statement":
        return lambda: views.card_html(b.observations[P + "synanthedon_tenuis.yang_2012.gland"], None)
    if case_id == "card_roles_as_list":
        card = compare(tl._views(x), tl._views(y), "C1-H", x, y)
        return lambda: views.card_html(card, {"study": load_results(tl.RESULTS), "contrast": "A_headgroup (control: Al vs OAc)",
                                              "x": x, "y": y, "roles": ["teaching_control"]})
    st = b.observations[P + "synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1"]
    bad = {"evidence_roles_as_list": {"roles": ["atlas_observed"]},
           "evidence_provenance_as_list": {"bundle": ["35112e1"]}}[case_id]
    return lambda: views.evidence_card_html(st, "spa.cmpd.z3_z13_18_oh", dict({"bundle_object": b}, **bad))


@pytest.mark.parametrize("case_id", MALFORMED)
def test_a_malformed_context_entry_or_object_is_refused(case_id, request):
    from model import Refused
    if case_id in XFAIL:
        request.applymarker(pytest.mark.xfail(strict=True, reason="integration probe on 1cbab22"))
    with pytest.raises(Refused):
        _malformed(case_id)()


def test_absent_entries_and_a_non_mapping_context_still_render():
    """The other side of the rule: absent roles and provenance render, and a
    whole non-mapping context reads as none for the card (#15)."""
    import views
    import test_card_label as tl
    from model import compare
    b = load(BUNDLE)
    st = b.observations[P + "synanthedon_tenuis.naka_2013.trap_z3z13oh_added.test4_1"]
    assert "evidence-card" in views.evidence_card_html(st, "spa.cmpd.z3_z13_18_oh", {"bundle_object": b})
    x, y = "z3_z13_18_oac", "z3_z13_18_al"
    card = compare(tl._views(x), tl._views(y), "C1-H", x, y)
    assert "Teaching comparison" in views.card_html(card, "not a mapping")
