"""Render the prototype pages from the records in records/.

    ESPMAP_PATH=/path/to/espmap python build_demo.py      # writes demo/*.html
"""
from pathlib import Path

import yaml

from record import Record, load_all, same_scale
from views import card_html, difference_card, page, site_table, svg_2d, viewer_js

HERE = Path(__file__).resolve().parent
REC, OUT = HERE / "records", HERE / "demo"
MAN = {m["id"]: m for m in yaml.safe_load((HERE / "molecules.yaml").read_text(encoding="utf-8"))["molecules"]}
LIMITS = ("<p class='limits'><b>What this shows:</b> local polarity and the electrostatic "
          "environment, on one stated scale. <b>Not:</b> which atom reacts, how fast, or by "
          "which mechanism.</p>")


def lesson_page(mol_id, site, out_name, caption):
    """site: a caller label from the manifest, or an atom map (an int) for a
    default site, whose label espmap derives from the map (#47)."""
    m = MAN[mol_id]
    rec = Record(REC / f"{mol_id}__s{12648430}")
    s = rec.site_by_map(site) if isinstance(site, int) else rec.site(site)
    site_label = s["label"]
    maps = [rec.atom_maps.get(i) for i in s["atoms"] if rec.atom_maps.get(i)]
    sc = rec.scale
    sub = (f"{m['display_name']} · site {site_label} highlighted · {sc['kind']} scale "
           f"{sc['vmin']:+.0f} … {sc['vmax']:+.0f} kcal/mol, {caption}")
    blocks = [
        '<div class="row"><div class="panel">' + svg_2d(m["mapped_smiles"], maps) + "</div>",
        ("viewer", f"{mol_id}_sticks", viewer_js(f"{mol_id}_sticks", rec, m["mapped_smiles"],
                                                 s["atoms"], show_surface=False)),
        ("viewer", f"{mol_id}_esp", viewer_js(f"{mol_id}_esp", rec, m["mapped_smiles"], s["atoms"])),
        "</div>", "<h2 style='font-size:16px'>Sites</h2>", site_table(rec, rec.site_labels()), LIMITS]
    (OUT / out_name).write_text(page(m["display_name"], sub, blocks), encoding="utf-8")
    return OUT / out_name


def card_page():
    e, z = load_all(REC, "e2_butene"), load_all(REC, "z2_butene")
    # the alkene carbon, chosen by its atom map (2), never by label: when
    # espmap #47 moved default labels from indices to maps, "C1" silently
    # changed from this carbon to the methyl and the card's answer flipped sign
    site = e[0].site_by_map(2)
    card = difference_card(e, z, site["label"], "(E)-but-2-ene", "(Z)-but-2-ene")
    pic = [("viewer", "e_esp", viewer_js("e_esp", e[0], MAN["e2_butene"]["mapped_smiles"],
                                         site["atoms"])),
           ("viewer", "z_esp", viewer_js("z_esp", z[0], MAN["z2_butene"]["mapped_smiles"],
                                         z[0].site_by_map(2)["atoms"]))]
    note = ("" if same_scale(e[0], z[0]) else
            "<p class='limits'>The two surfaces are on their own auto-scaled ranges (different "
            "render keys): compare the numbers in the card, not the colours.</p>")
    blocks = ['<div class="row">'] + pic + ["</div>", note, card_html(card), LIMITS]
    (OUT / "difference_card.html").write_text(
        page("E vs Z but-2-ene", f"site {site['label']} (the alkene carbon, map 2), 3 seeds each",
             blocks), encoding="utf-8")
    return card


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    print(lesson_page("acetaldehyde", "carbonyl O", "acetaldehyde.html", "warm = negative"))
    print(lesson_page("acetate", 3, "acetate.html", "darker = more negative"))  # an O, map 3
    print(card_page())
