"""The three teaching views, in the Arrowfield house style.

House style (agreed 2026-10-22): white background, thin black outline,
heteroatoms in their familiar colours (O red, N blue, S yellow, P orange),
carbons by role (muted blue; the selected site lilac), hydrogens hidden
except a selected site's own, sparse black labels. Rendered with 3Dmol.js,
which is what py3Dmol emits in a notebook, so the same calls work there.

Views:
  svg_2d        2D skeletal formula with the atom-map numbers the site table uses
  viewer_js     3D sticks, optionally with the ESP surface, the site highlighted
  difference_card   one site, two molecules, with the uncertainty beside it
  page          the three together, as one HTML page
Views take Arrowfield's model objects (model.py: MoleculeView, Comparison),
never an espmap record or an Atlas record. record.py and bundle.py translate.
"""
import html
import json
import os
import sys

import numpy as np

from model import SIGN_TOL, ComparisonStatus, compare

THREEDMOL = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.0/3Dmol-min.js"
STYLE = dict(background="white", outline="black", outline_width=0.03,
             C="#7ea6d8", C_site="#c98bbf", O="#d62728", N="#3366cc", S="#d9b300", P="#e67e22",
             H="#e8e8e8", stick=0.17)
ELEM = {"O": STYLE["O"], "N": STYLE["N"], "S": STYLE["S"], "P": STYLE["P"], "H": STYLE["H"]}


# ------------------------------------------------------------ geometry utils
def principal_frame(xyz):
    """Centroid and axes so every viewer starts from the same, reproducible view:
    longest axis across the screen, flattest axis toward the viewer."""
    c = xyz.mean(0)
    _, _, vt = np.linalg.svd(xyz - c)
    if np.linalg.det(vt) < 0:
        vt[2] *= -1
    return c, vt


_USUAL_VALENCE = {"H": 1, "C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "F": 1, "Cl": 1,
                  "Br": 1, "I": 1, "B": 3, "Si": 4, "Se": 2}


def mol_from_view(rec, mapped_smiles=None):
    """RDKit mol in the record's atom order, with its coordinates.

    The record's order is the input's (espmap spec section 6): SMILES order
    with hydrogens appended for --smiles, the file's own order for --file. So
    nothing here assumes a canonical order. Bonds come from the record's own
    geometry.bonds (espmap #47, A1), which needs no atom maps and so also
    draws an unmapped input. Formal charges follow from each atom's bond
    valence and must add up to the record's charge. A record without bonds
    (espmap could not confirm its atom order, or a pre-#47 record) falls back
    to rebuilding them from the mapped SMILES."""
    from rdkit import Chem
    if rec.bonds is not None:
        m = Chem.RWMol()
        for s in rec.symbols:
            x = Chem.Atom(s)
            x.SetNoImplicit(True)
            m.AddAtom(x)
        order = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
        for i, j, k in rec.bonds:
            m.AddBond(i, j, order[k])
        for a in m.GetAtoms():
            if a.GetSymbol() in _USUAL_VALENCE:
                a.SetFormalCharge(sum(int(b.GetBondTypeAsDouble()) for b in a.GetBonds())
                                  - _USUAL_VALENCE[a.GetSymbol()])
        for a in m.GetAtoms():
            if a.GetAtomMapNum() == 0 and a.GetIdx() in rec.atom_maps:
                a.SetAtomMapNum(rec.atom_maps[a.GetIdx()])
        m = m.GetMol()
        if sum(a.GetFormalCharge() for a in m.GetAtoms()) != rec.charge:
            raise ValueError("charges read from the record's bonds do not add up to its charge")
        Chem.SanitizeMol(m)
        return _with_coords(m, rec)
    if mapped_smiles is None:
        raise ValueError("record has no bonds and no mapped SMILES was given")
    tpl = Chem.MolFromSmiles(mapped_smiles)
    Chem.Kekulize(tpl, clearAromaticFlags=True)        # explicit single/double bonds
    by_map = {a.GetAtomMapNum(): a for a in tpl.GetAtoms() if a.GetAtomMapNum()}
    atoms = rec.atoms
    if any(a.element != "H" and a.atom_map not in by_map for a in atoms):
        raise ValueError("record has heavy atoms the mapped SMILES does not name")
    m = Chem.RWMol()
    rec_of_tpl = {}
    for a in atoms:
        x = Chem.Atom(a.element)
        if a.element != "H":
            t = by_map[a.atom_map]
            x.SetFormalCharge(t.GetFormalCharge())
            x.SetAtomMapNum(a.atom_map)
            rec_of_tpl[t.GetIdx()] = a.index
        x.SetNoImplicit(True)
        m.AddAtom(x)
    for b in tpl.GetBonds():
        m.AddBond(rec_of_tpl[b.GetBeginAtomIdx()], rec_of_tpl[b.GetEndAtomIdx()], b.GetBondType())
    for a in atoms:
        if a.element == "H":
            m.AddBond(a.parent, a.index, Chem.BondType.SINGLE)
    m = m.GetMol()
    Chem.SanitizeMol(m)
    if tuple(x.GetSymbol() for x in m.GetAtoms()) != tuple(rec.symbols):
        raise ValueError("record atoms and their element list disagree")
    return _with_coords(m, rec)


def _with_coords(m, rec):
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    conf = Chem.Conformer(m.GetNumAtoms())
    for i, p in enumerate(rec.xyz):
        conf.SetAtomPosition(i, Point3D(*map(float, p)))
    m.AddConformer(conf, assignId=True)
    return m


def colours(V, scale):
    """Vertex colours from espmap's own colormap and the record's scale."""
    sys.path.insert(0, os.environ.get("ESPMAP_PATH", ""))
    import esp_map as em
    cmap = em.build_maps()[scale.cmap]
    t = np.clip((V - scale.vmin) / (scale.vmax - scale.vmin), 0, 1)
    return ["#%02x%02x%02x" % tuple(int(255 * c) for c in cmap(float(x))[:3]) for x in t]


def vertex_normals(verts, faces):
    n = np.zeros_like(verts)
    fn = np.cross(verts[faces[:, 1]] - verts[faces[:, 0]], verts[faces[:, 2]] - verts[faces[:, 0]])
    for k in range(3):
        np.add.at(n, faces[:, k], fn)
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)


def site_faces(surface, atoms):
    """Triangles owned by the site's atoms (espmap's NPZ ownership): the patch."""
    return np.isin(surface.owner, atoms)


# ------------------------------------------------------------------- 2D view
def svg_2d(mapped_smiles, site_maps=(), size=(340, 260)):
    """Skeletal formula; every heavy atom shows its map number, the number
    the site table uses. Site atoms are highlighted lilac."""
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    m = Chem.MolFromSmiles(mapped_smiles)
    hl = []
    for a in m.GetAtoms():
        n = a.GetAtomMapNum()
        a.SetProp("atomNote", str(n))
        a.SetAtomMapNum(0)
        if n in site_maps:
            hl.append(a.GetIdx())
    return _draw_2d(m, hl, size)


def svg_2d_view(rec, site_atoms=(), size=(340, 260)):
    """The same skeletal formula drawn from the record alone: its bonds
    (espmap #47), E/Z assigned from its 3D coordinates, no SMILES needed. An
    unmapped record draws too, numbered by atom index instead of map.
    site_atoms are record indices, as in reaction_sites[*].atoms."""
    from rdkit import Chem
    m = Chem.Mol(mol_from_view(rec, rec.mapped_smiles))
    Chem.AssignStereochemistryFrom3D(m)
    for a in m.GetAtoms():
        a.SetIntProp("rec_idx", a.GetIdx())
        a.SetProp("atomNote", str(rec.atom_maps.get(a.GetIdx(), a.GetIdx())))
        a.SetAtomMapNum(0)
    m = Chem.RemoveHs(m)
    m.RemoveAllConformers()
    hl = [a.GetIdx() for a in m.GetAtoms() if a.GetIntProp("rec_idx") in set(site_atoms)]
    return _draw_2d(m, hl, size)


def _draw_2d(m, hl, size):
    from rdkit.Chem.Draw import rdMolDraw2D
    colour = {i: (0.85, 0.65, 0.81) for i in hl}
    d = rdMolDraw2D.MolDraw2DSVG(*size)
    o = d.drawOptions()
    o.annotationFontScale = 0.75
    o.addStereoAnnotation = True
    o.bondLineWidth = 2
    d.DrawMolecule(m, highlightAtoms=hl, highlightAtomColors=colour, highlightBonds=[])
    d.FinishDrawing()
    return d.GetDrawingText()


# ------------------------------------------------------------------- 3D view
def viewer_js(div_id, rec, mapped_smiles, site_atoms=(), show_surface=True, surface_opacity=0.85,
              rule="nearest_nucleus", labels=()):
    """3Dmol.js calls for one viewer. Site atoms (atom indices) and their own
    hydrogens are drawn lilac; with the surface on, the site's patch stays
    opaque and the rest fades, so the table row and the picture are one thing."""
    from rdkit import Chem
    m = mol_from_view(rec, mapped_smiles)
    c, R = principal_frame(rec.xyz)
    conf = m.GetConformer()
    for i, p in enumerate((rec.xyz - c) @ R.T):
        conf.SetAtomPosition(i, p.tolist())
    Chem.Kekulize(m, clearAromaticFlags=True)
    molblock = Chem.MolToMolBlock(m)
    site = sorted(set(site_atoms))
    site_h = [n.GetIdx() for i in site for n in m.GetAtomWithIdx(i).GetNeighbors()
              if n.GetAtomicNum() == 1]
    js = [f"var v = $3Dmol.createViewer(document.getElementById('{div_id}'), "
          f"{{backgroundColor:'{STYLE['background']}'}});",
          f"v.addModel({json.dumps(molblock)}, 'mol');",
          f"v.setViewStyle({{style:'outline', color:'{STYLE['outline']}', width:{STYLE['outline_width']}}});",
          "v.setStyle({}, {});",
          f"v.setStyle({{not:{{elem:'H'}}}}, {{stick:{{radius:{STYLE['stick']}, "
          f"colorscheme:{{prop:'elem', map:{json.dumps(dict(ELEM, C=STYLE['C']))}}}}}}});"]
    if site:
        js.append(f"v.setStyle({{index:{json.dumps(site + site_h)}}}, {{stick:{{radius:{STYLE['stick'] * 1.15}, "
                  f"colorscheme:{{prop:'elem', map:{json.dumps(dict(ELEM, C=STYLE['C_site']))}}}}}}});")
        # a translucent halo, so an O or N site (which keeps its element
        # colour) is still visibly the selected one
        for i in site:
            p = conf.GetAtomPosition(i)
            js.append(f"v.addSphere({{center:{{x:{p.x:.3f},y:{p.y:.3f},z:{p.z:.3f}}}, radius:0.75, "
                      f"color:'{STYLE['C_site']}', opacity:0.5}});")
    if show_surface:
        s = rec.surface(rule)
        verts = (s.verts - c) @ R.T
        cols = colours(s.V, rec.scale)
        norms = vertex_normals(verts, s.faces)
        inside = site_faces(s, site) if site else np.ones(len(s.faces), bool)
        for mask, op in ((inside, surface_opacity), (~inside, 0.22 if site else surface_opacity)):
            f = s.faces[mask]
            if not len(f):
                continue
            used = np.unique(f)
            remap = -np.ones(len(verts), int)
            remap[used] = np.arange(len(used))
            js.append("v.addCustom({vertexArr:%s, normalArr:%s, faceArr:%s, color:%s, opacity:%.2f});" % (
                json.dumps([dict(x=round(float(a), 3), y=round(float(b), 3), z=round(float(cc), 3))
                            for a, b, cc in verts[used]]),
                json.dumps([dict(x=round(float(a), 3), y=round(float(b), 3), z=round(float(cc), 3))
                            for a, b, cc in norms[used]]),
                json.dumps(remap[f].ravel().tolist()),
                json.dumps([cols[i] for i in used]), op))
    for idx, text in labels:
        js.append(f"v.addLabel({json.dumps(text)}, {{position:{{x:{conf.GetAtomPosition(idx).x:.3f}, "
                  f"y:{conf.GetAtomPosition(idx).y + 0.6:.3f}, z:{conf.GetAtomPosition(idx).z:.3f}}}, "
                  "fontSize:16, fontColor:'black', backgroundOpacity:0, inFront:true});")
    # zoomTo leaves a wide margin around a small molecule; tighten it
    js += ["v.zoomTo();", f"v.zoom({1.15 if show_surface else 1.9});", "v.render();"]
    return "\n".join(js)


# ----------------------------------------------------------- difference card
def difference_card(views_a, views_b, label, name_a="A", name_b="B"):
    """One site, two molecules, each a list of conformers: a Comparison
    (model.compare). Measurements that are not comparable come back as a
    NOT_COMPARABLE comparison, which the card shows, not as an error."""
    return compare(views_a, views_b, label, name_a, name_b)


def _scope(card):
    """What the numbers are of, with the method from the measurements' own
    provenance: a hard-coded 'gas phase' would be false on a PCM card."""
    how = f" ({html.escape(card.method)})" if card.method else ""
    return (f"Site means of the free molecule{how}, all under one calculation protocol. "
            "On its own, a difference here says nothing about how an organism responds.")


def _range(x):
    """A conformer range, or a dash when there is one conformer: 0.00 would
    read as 'no variation'."""
    return "–" if x is None else f"{x:.2f}"


def _extended_row(card):
    e = card.extended
    if e is None:
        return ""
    if e.get("unstated"):
        return (f"\n<tr><td>… extended conformers only: not shown, because whether a conformer "
                f"is folded is not stated for {e['unstated']} of them</td><td>–</td></tr>")
    return (f"\n<tr><td>… extended conformers only (n = {e['n'][0]} / {e['n'][1]}; "
            f"folded ones removed)</td><td>{_range(e['range_a'])} / {_range(e['range_b'])}</td></tr>")


def _extended_verdict(card):
    e = card.extended
    if e is None or e.get("unstated") or card.status is not ComparisonStatus.COMPARABLE:
        return ""
    v = ("cannot judge (fewer than two extended conformers on one side)" if e["verdict"] is None
         else e["verdict"].value)
    return (f'<div class="verdict-ext">extended conformers only (n = {e["n"][0]} / {e["n"][1]}): '
            f"{html.escape(v)}</div>")


# The card's verdict is a teaching look under its own rule, never the frozen
# study's test (tests/pedagogy/card_label_set.json, #15): it labels itself,
# never borrows the study's outcome words, and sends the reader to the study
# card for the study's verdict.
_CARD_LABEL = ("Teaching comparison: a look, not a test. The verdict below is the card's own rule "
               "(every conformer of one molecule beyond every conformer of the other, beyond the partition "
               "sensitivity, with one sign under both partition rules).")
_CARD_ROLES = ("atlas_observed", "teaching_control", "computational_counterfactual")
_TIER_WORDS = {"primary": "a prespecified primary test", "secondary": "a prespecified secondary test"}


def _esc(x):
    # text content: quotes need no escaping, and an escaped apostrophe reads as &#x27; in plain text
    return html.escape(str(x), quote=False)


def _verbatim(x):
    """Record text (a site label, a study's plan reference), verbatim and
    marked: the card's word bans read only the text outside these spans."""
    return f'<span class="verbatim">{_esc(x)}</span>'


def _card_notes(card, context):
    """Notes from an optional context: the study's stance on this site, read
    from the study itself, and any teaching control. Refused when the context
    cannot be this card's."""
    from model import Refused, StudyResults
    if not isinstance(context, dict):
        return []                          # not a mapping: read as no context
    roles = context.get("roles") if isinstance(context.get("roles"), dict) else {}
    bad = {k: v for k, v in roles.items() if v not in _CARD_ROLES}
    if bad:
        raise Refused(f"role(s) {bad} are not one of {_CARD_ROLES}: the card cannot say what those molecules are")
    x, y = context.get("x", card.name_a), context.get("y", card.name_b)
    notes = []
    study = context.get("study")
    if study is not None:
        from model import KNOWN_VERDICT_RULES
        if not isinstance(study, StudyResults):
            raise Refused("the study must be a StudyResults read by study.load_results, which checks the file")
        if study.verdict_rule not in KNOWN_VERDICT_RULES:
            raise Refused(f"the study's verdict rule {study.verdict_rule!r} is not one Arrowfield knows "
                          f"({KNOWN_VERDICT_RULES}): its tiers cannot be described")
        if (card.name_a, card.name_b) != (x, y):
            raise Refused(f"this card compares {card.name_a} → {card.name_b}, but the study context names "
                          f"{x} → {y}: the study's statement would be attached to another pair or direction")
        contrast = context.get("contrast")
        rows = [s for s in study.sites if s.contrast == contrast and s.population == "all"]
        if not rows:
            raise Refused(f"the study has no contrast {contrast!r}")
        if (rows[0].x, rows[0].y) != (x, y):
            raise Refused(f"the study's {contrast} compares {rows[0].x} → {rows[0].y}, this card {x} → {y}: "
                          "its difference would carry another sign or describe another pair")
        row = next((s for s in rows if s.site == card.site), None)
        plan = _verbatim(study.plan)
        if row is None:
            notes.append(f"The study ({plan}) gives no verdict at this site: it does not test it.")
        elif row.tier in _TIER_WORDS:
            notes.append(f"The study ({plan}) decides this site as {_TIER_WORDS[row.tier]}, under a "
                         "different rule: its verdict is on the study card, not here.")
        else:
            notes.append(f"The study ({plan}) reports this site as exploratory, under a different rule: "
                         "its row is on the study card.")
    for slug, name in ((x, card.name_a), (y, card.name_b)):
        if roles.get(slug) == "teaching_control":
            notes.append(f"{_esc(name)} is a teaching control: a constructed molecule, not an Atlas compound. "
                         "This comparison is descriptive only.")
        elif roles.get(slug) == "computational_counterfactual":
            notes.append(f"{_esc(name)} is a computational counterfactual, not an observed molecule. "
                         "This comparison is descriptive only.")
    return notes


def card_html(card, context=None):
    """The difference card. context (optional): {"study": StudyResults,
    "contrast", "x", "y": the card's molecules' slugs, "roles": {slug: role}};
    model.Refused when it cannot be this card's study or a role is unknown,
    and when a number it would print is not finite."""
    import math
    from model import Refused
    shown = (card.diff, card.partition, card.diff_other_rule, card.range_a, card.range_b)
    if card.status is not ComparisonStatus.NOT_COMPARABLE and any(
            v is not None and not (isinstance(v, (int, float)) and math.isfinite(v)) for v in shown):
        raise Refused("a difference, range or partition sensitivity on this card is not a finite number: "
                      "recompute the comparison")
    notes = "".join(f'<div class="limits">{n}</div>' for n in _card_notes(card, context))     # HTML, escaped
    label = f'<div class="label">{_esc(_CARD_LABEL)}</div>'
    title = (f'<div class="card-title">{html.escape(card.name_a)} → {html.escape(card.name_b)} '
             f'· site {_verbatim(card.site)}<span class="unit">site mean, kcal/mol</span></div>')
    banner = (f'<div class="status status-{card.status.name.lower()}">'
              f'{card.status.value.upper()}</div>')
    if card.status is ComparisonStatus.NOT_COMPARABLE:
        return (f'<div class="card">{banner}{title}{label}'
                f'<div class="verdict">{html.escape(card.verdict.value)}</div>{notes}</div>')
    return f"""
<div class="card">{banner}{title}{label}
<table>
<tr><td>difference of conformer means</td><td>{card.diff:+.2f}</td></tr>
<tr><td>conformer range ({html.escape(card.name_a)} / {html.escape(card.name_b)}, n = {card.n[0]} / {card.n[1]})</td>
<td>{_range(card.range_a)} / {_range(card.range_b)}</td></tr>{_extended_row(card)}
<tr><td>partition sensitivity (max)</td><td>{card.partition:.2f}</td></tr>
<tr><td>difference of means under the other partition rule (a sign within ±{SIGN_TOL} is not counted)</td>
<td>{card.diff_other_rule:+.2f}</td></tr>
</table>
<div class="verdict">{html.escape(card.verdict.value)}</div>{_extended_verdict(card)}{notes}
<div class="scope">{_scope(card)}</div></div>"""


# --------------------------------------------------------- study-verdict card
_HEADLINE = {
    "resolved": "Resolved under this study's protocol",
    "different on average": "Different on average: not resolved",
    "unresolved": "Not resolved under this study's protocol",
    "not evaluable": "Not evaluable: the study gives no verdict for this contrast",
    "no data": "No data: the study has no value for this site",
}
_SAYS = {
    "resolved": "The free molecule's ESP at this site separates the two by more than their conformational "
                "variation, under this protocol.",
    "different on average": "The means differ, but the observed conformer ranges overlap.",
    "unresolved": "Under this protocol, free-molecule ESP does not separate the two at this site.",
    "not evaluable": "Too many members were excluded, or an excluded member was folded (plan section 3).",
    "no data": "",
}
_CONDITION_WORDS = (("1_ranges_disjoint", "conformer ranges separated"),
                    ("2_beyond_partition", "beyond the partition sensitivity"),
                    ("3_same_sign", "same sign under the radius-scaled rule"))


def _sfmt(x, fmt="+.2f"):
    return "–" if x is None else format(x, fmt)


def _span(r):
    """A (low, high) range, or why it is not shown: a not-evaluable row may
    carry anything, and only a finite ordered pair is displayed (probe round
    2 of #11: a string range crashed the card, an infinite one was shown)."""
    from model import _bad_pair
    if r is None:
        return "–"
    if _bad_pair(r):
        return "not shown (not a finite range)"
    return f"{_sfmt(r[0])} … {_sfmt(r[1])}"


def study_card_html(row, context):
    """One frozen-study row: the study's own verdict at one site of one
    contrast, never re-decided here. row is a model.StudySite; context
    carries provenance ({results_sha256, plan_commit}), each molecule's role
    ({slug: atlas_observed | teaching_control}), the limitations ({basis_
    sensitivity, solvent, phase, method}) and, for the extended-only
    sensitivity, the all-member row it sits beside.

    Refused (model.Refused), never rendered: a row under a rule Arrowfield
    does not know; a row whose outcome contradicts its own conditions; an
    extended-only row without its all-member row. Specified by the hostile
    set tests/pedagogy/study_card_set.json (#10)."""
    from model import KNOWN_VERDICT_RULES, Refused, study_row_inconsistency, study_row_malformed
    if row.verdict_rule not in KNOWN_VERDICT_RULES:
        raise Refused(f"verdict rule {row.verdict_rule!r} is not one Arrowfield can show: a study's "
                      "verdict is shown only under the rule that decided it")
    why = study_row_inconsistency(row.verdict_rule, row.outcome, row.conditions or {},
                                  row.bootstrap_interval) or study_row_malformed(row)
    if why:
        raise Refused(f"this study row is corrupt and is not shown: {why}")
    # A context, or an entry of it, that is not a mapping is read as missing:
    # never a crash, never a claim the card could not read (probe round 4).
    context = context if isinstance(context, dict) else {}
    primary = context.get("all_member_row")
    if isinstance(primary, dict):                      # given as data: read it as a StudySite
        from model import StudySite
        primary = StudySite(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in primary.items()})
    if row.population != "all":
        if primary is None:
            raise Refused("an extended-only sensitivity row is shown only beside the all-member result; "
                          "pass the all-member row as context['all_member_row']")
        if (primary.population, primary.contrast, primary.site, primary.verdict_rule) != \
                ("all", row.contrast, row.site, row.verdict_rule):
            raise Refused("the all-member row given is not this site's all-member result "
                          f"({primary.population} {primary.contrast} {primary.site}); a sensitivity is "
                          "shown only beside its own primary result")
        why = study_row_inconsistency(primary.verdict_rule, primary.outcome, primary.conditions or {},
                                      primary.bootstrap_interval) or study_row_malformed(primary)
        if why:
            raise Refused(f"the all-member row is corrupt and is not shown: {why}")
    prov, roles, lim = (x if isinstance(x, dict) else {}
                        for x in (context.get("provenance"), context.get("molecules"),
                                  context.get("limitations")))
    def esc(t):
        return html.escape(str(t), quote=False)      # text, not attribute values
    title = (f'<div class="card-title">{esc(row.site)} · {esc(row.x)} → {esc(row.y)}'
             f'<span class="unit">{esc(row.contrast)}, kcal/mol</span></div>')
    parts = [f'<div class="card study-card">{title}',
             '<div class="status">STUDY VERDICT</div>',
             f'<div class="verdict verdict-{row.outcome.replace(" ", "-")}">{esc(_HEADLINE[row.outcome])}</div>',
             f'<div class="rule">Decided by {esc(row.verdict_rule)}.</div>']
    if row.tier == "exploratory":
        parts.append('<div class="flag">Exploratory site: not a prespecified test, and never interpreted.</div>')
    if row.population != "all":
        parts.append(f'<div class="flag">Sensitivity analysis: extended (unfolded) conformers only '
                     f'(n = {row.n[0]} / {row.n[1]}), shown beside the all-members result, never instead '
                     f'of it. All members: {esc(_HEADLINE[primary.outcome])}.</div>')
    controls = [s for s in (row.x, row.y) if roles.get(s) == "teaching_control"]
    if controls:
        parts.append(f'<div class="flag">Teaching control: {esc(", ".join(controls))} is a constructed molecule, '
                     'not an Atlas compound; descriptive chemistry only.</div>')
    if _SAYS[row.outcome]:
        parts.append(f"<p>{esc(_SAYS[row.outcome])}</p>")
    rows = [f"<tr><td>{esc(row.x)} conformer range (n = {row.n[0]})</td><td>{_span(row.range_x)}</td></tr>",
            f"<tr><td>{esc(row.y)} conformer range (n = {row.n[1]})</td><td>{_span(row.range_y)}</td></tr>"]
    if row.outcome not in ("not evaluable", "no data"):
        rows = ([f"<tr><td>Δ mean ({esc(row.y)} − {esc(row.x)})</td><td>{_sfmt(row.delta)}</td></tr>"] + rows +
                [f"<tr><td>partition sensitivity</td><td>{_sfmt(row.partition, '.2f')}</td></tr>",
                 f"<tr><td>bootstrap interval over the selected conformers</td>"
                 f"<td>{_span(row.bootstrap_interval)}</td></tr>"])
    parts.append("<table>" + "".join(rows) + "</table>")
    if row.outcome not in ("not evaluable", "no data") and row.conditions:
        parts.append("<ul class='conditions'>" + "".join(
            f"<li>{'✓' if row.conditions.get(k) else '✗'} {w}</li>" for k, w in _CONDITION_WORDS) + "</ul>")
        parts.append('<p class="limits">The bootstrap interval resamples the selected conformers: it describes '
                     "this set of shapes, not a population.</p>")
    method = lim.get("method", "method not stated")
    parts.append(f'<div class="scope">Scope: free molecule · {esc(method)} · {esc(lim.get("phase", "phase not stated"))} '
                 "phase.</div>")
    parts.append("<ul class='limits'>"
                 f"<li>basis sensitivity: {esc(lim.get('basis_sensitivity', 'not stated'))}</li>"
                 f"<li>solvent: {esc(lim.get('solvent', 'not stated'))}</li>"
                 "<li>protein or membrane context: outside this study</li></ul>")
    parts.append('<p class="limits">This does not establish biological discrimination.</p>')
    # Provenance is shown only in a form that can be checked: a 64-hex results
    # sha256 and a 7-40-hex plan commit. Anything else (missing, an int, a
    # malformed string; probe round 3 of #11) is named, with a visible warning.
    import re
    sha, plan = prov.get("results_sha256"), prov.get("plan_commit")
    sha_ok = isinstance(sha, str) and re.fullmatch(r"[0-9a-f]{64}", sha) is not None
    plan_ok = isinstance(plan, str) and re.fullmatch(r"[0-9a-f]{7,40}", plan) is not None
    known = []
    if plan_ok:
        known.append(f"Plan frozen at {plan[:7]}")
    if sha_ok:
        known.append(f"results sha256 {sha[:12]}…")
    if known:
        parts.append(f'<div class="provenance">{esc(" · ".join(known))}</div>')
    missing = [w for ok, w in ((plan_ok, "the frozen plan's commit"), (sha_ok, "the results file's sha256"))
               if not ok]
    if missing:
        parts.append(f'<div class="provenance warn">Provenance incomplete: {esc(" and ".join(missing))} '
                     f'{"is" if len(missing) == 1 else "are"} not recorded in a verifiable form, so this '
                     "card cannot be fully traced to the frozen study.</div>")
    parts.append("</div>")
    return "\n".join(parts)


# ----------------------------------------------------------------------- page
CSS = """
:root { --ink:#111; --muted:#666; --line:#ddd; --card:#fafafa; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --ink:#111; } }
body { background:#fff; color:var(--ink); font:15px/1.45 system-ui, sans-serif; margin:0; padding:16px; }
h1 { font-size:20px; margin:0 0 4px; } .sub { color:var(--muted); margin:0 0 14px; }
.row { display:flex; flex-wrap:wrap; gap:14px; align-items:flex-start; }
.panel { border:1px solid var(--line); border-radius:8px; padding:8px; background:#fff; }
.viewer { width:min(460px, 92vw); height:340px; position:relative; }
table { border-collapse:collapse; font-variant-numeric:tabular-nums; }
td, th { padding:3px 10px; border-bottom:1px solid var(--line); text-align:right; }
td:first-child, th:first-child { text-align:left; }
.card { border:1px solid var(--line); border-radius:8px; padding:10px 12px; background:var(--card); max-width:520px; }
.card-title { font-weight:600; margin-bottom:6px; } .unit { font-weight:400; color:var(--muted); margin-left:8px; font-size:13px; }
.verdict { margin-top:8px; font-weight:600; } .verdict-ext { margin-top:4px; }
.scope { margin-top:8px; color:var(--muted); font-size:13px; }
.study-card .rule, .study-card .provenance { color:var(--muted); font-size:13px; margin-top:4px; }
.study-card .flag { margin-top:6px; padding:4px 8px; border-left:3px solid var(--ink); }
.study-card .warn { color:#a40000; } .study-card ul { margin:6px 0; padding-left:18px; }
.status { display:inline-block; font-size:12px; font-weight:700; letter-spacing:.04em; padding:1px 6px;
          border:1px solid var(--ink); border-radius:4px; margin-bottom:6px; }
.limits { color:var(--muted); font-size:13px; max-width:900px; margin-top:12px; }
"""


def _num(x, fmt):
    return "" if x is None else format(x, fmt)


def site_table(recs, labels):
    """Sites in three tiers, in reading order.

    primary    the site's area-weighted mean (over conformers, when several);
    secondary  partition sensitivity and conformer range, the two numbers a
               difference must beat;
    advanced   folded away: areas, the 5th/95th percentiles, raw extrema, and
               the site's single atoms. Atom values are diagnostic only (the
               rule moves an atom by up to 11.9 kcal/mol, espmap #42); a
               beginner should never meet "the ESP of C3" as a number.
    recs: one Record or a list (one molecule's conformers)."""
    recs = recs if isinstance(recs, list) else [recs]
    head, adv, atoms = [], [], []
    for lab in labels:
        ss = [r.site(lab) for r in recs]
        means = [s.mean for s in ss if s.mean is not None]
        part = max((s.partition for s in ss if s.partition is not None), default=None)
        rng = f"{max(means) - min(means):.2f}" if len(means) > 1 else "1 conformer"
        mean = f"{np.mean(means):+.1f}" if means else "buried"
        mark = " †" if ss[0].overlaps else ""
        head.append(f"<tr><td>{html.escape(lab)}{mark}</td><td><b>{mean}</b></td>"
                    f"<td>{_num(part, '.2f')}</td><td>{rng}</td></tr>")
        s = ss[0]
        adv.append(f"<tr><td>{html.escape(lab)}</td><td>{s.area:.1f}</td>"
                   f"<td>{s.area_neg:.1f} / {s.area_pos:.1f}</td>"
                   f"<td>{_num(s.p05, '+.1f')} / {_num(s.p95, '+.1f')}</td>"
                   f"<td>{_num(s.vmin, '+.1f')} / {_num(s.vmax, '+.1f')}</td></tr>")
        for a in recs[0].site_atoms(lab):
            nn, rs = a.by_rule.get("nearest_nucleus"), a.by_rule.get("radius_scaled")
            mean_nn = _num(nn[0], "+.1f") or "buried"
            mean_rs = "" if rs is None else (_num(rs[0], "+.1f") or "buried")
            area_rs = "" if rs is None else format(rs[1], ".1f")
            atoms.append(f"<tr><td>{html.escape(lab)}</td><td>{a.element}{a.index}</td>"
                         f"<td>{mean_nn} / {mean_rs}</td>"
                         f"<td>{nn[1]:.1f} / {area_rs}</td></tr>")
    first = "" if len(recs) == 1 else " (first conformer)"
    shown = {lab: recs[0].site(lab).overlaps for lab in labels}
    shared = {lab: ov for lab, ov in shown.items() if ov}
    tiling = ""
    if shared or not recs[0].sites_tile_surface:
        pairs = "; ".join(f"{html.escape(a)} with {', '.join(map(html.escape, b))}"
                          for a, b in shared.items())
        tiling = ("<p class='limits'>† These sites share atoms" + (f" ({pairs})" if pairs else "")
                  + ". Site areas do not add up to the surface, and two overlapping "
                  "sites' means are not independent.</p>" if shared else
                  "<p class='limits'>These sites do not cover every atom, so their areas "
                  "do not add up to the surface.</p>")
    return (
        "<table><tr><th>site</th><th>mean</th><th>partition Δ</th><th>conformer range</th></tr>"
        + "".join(head) + "</table>" + tiling
        + f"<details><summary>Advanced: areas, tails, extrema, single atoms{first}</summary>"
        "<table><tr><th>site</th><th>area Å²</th><th>− / + area</th><th>p05 / p95</th>"
        "<th>V<sub>min</sub> / V<sub>max</sub></th></tr>" + "".join(adv) + "</table>"
        "<p class='limits'>Single atoms, nearest-nucleus / radius-scaled. These are "
        "<b>diagnostic</b>: assigning surface to one atom is a convention, and the two "
        "rules disagree by up to 11.9 kcal/mol on one atom. Quote the site.</p>"
        "<table><tr><th>site</th><th>atom</th><th>mean</th><th>area Å²</th></tr>"
        + "".join(atoms) + "</table></details>")


def page(title, subtitle, blocks):
    """blocks: list of HTML strings; viewers are ('viewer', div_id, js) tuples."""
    body, scripts = [], []
    for b in blocks:
        if isinstance(b, tuple):
            _, div_id, js = b
            body.append(f'<div class="panel"><div id="{div_id}" class="viewer"></div></div>')
            scripts.append(f"(function(){{\n{js}\n}})();")
        else:
            body.append(b)
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title><style>{CSS}</style>
<script src="{THREEDMOL}"></script></head><body>
<h1>{html.escape(title)}</h1><p class="sub">{html.escape(subtitle)}</p>
{''.join(body)}
<script>{''.join(scripts)}</script></body></html>"""
