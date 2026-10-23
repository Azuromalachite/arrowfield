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

from model import ComparisonStatus, compare

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
    return (f"\n<tr><td>… extended conformers only (n = {e['n'][0]} / {e['n'][1]}; "
            f"folded ones removed)</td><td>{_range(e['range_a'])} / {_range(e['range_b'])}</td></tr>")


def _extended_verdict(card):
    e = card.extended
    if e is None or card.status is not ComparisonStatus.COMPARABLE:
        return ""
    v = ("cannot judge (fewer than two extended conformers on one side)" if e["verdict"] is None
         else e["verdict"].value)
    return (f'<div class="verdict-ext">extended conformers only (n = {e["n"][0]} / {e["n"][1]}): '
            f"{html.escape(v)}</div>")


def card_html(card):
    title = (f'<div class="card-title">{html.escape(card.name_a)} → {html.escape(card.name_b)} '
             f'· site {html.escape(card.site)}<span class="unit">site mean, kcal/mol</span></div>')
    banner = (f'<div class="status status-{card.status.name.lower()}">'
              f'{card.status.value.upper()}</div>')
    if card.status is ComparisonStatus.NOT_COMPARABLE:
        return (f'<div class="card">{banner}{title}'
                f'<div class="verdict">{html.escape(card.verdict.value)}</div></div>')
    return f"""
<div class="card">{banner}{title}
<table>
<tr><td>difference of conformer means</td><td>{card.diff:+.2f}</td></tr>
<tr><td>conformer range ({html.escape(card.name_a)} / {html.escape(card.name_b)}, n = {card.n[0]} / {card.n[1]})</td>
<td>{_range(card.range_a)} / {_range(card.range_b)}</td></tr>{_extended_row(card)}
<tr><td>partition sensitivity (max)</td><td>{card.partition:.2f}</td></tr>
</table>
<div class="verdict">{html.escape(card.verdict.value)}</div>{_extended_verdict(card)}
<div class="scope">{_scope(card)}</div></div>"""


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
