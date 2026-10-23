"""Pages from espmap's P3 fixture records, read as Arrowfield will read them:
through record.py only, and with nothing from the manifest. If a page needs
anything the record does not carry, that is a finding about the contract
before it freezes (espmap #44, review 5).

    python build_p3.py /path/to/espmap-checkout-with-fixtures   # writes demo/p3_*.html
    python build_p3.py --records DIR                            # or a dir of records

Records from after espmap #47 carry geometry.bonds, map-based labels and the
site-overlap fields, all of which these pages use.
"""
import sys
from pathlib import Path

from record import load, same_scale
from views import page, site_table, svg_2d_view, viewer_js

HERE = Path(__file__).resolve().parent
FIX = (Path(sys.argv[2]).resolve() if sys.argv[1] == "--records"
       else Path(sys.argv[1]).resolve() / "fixtures" / "records")
OUT = HERE / "demo"
LIMITS = ("<p class='limits'><b>What this shows:</b> local polarity and the electrostatic "
          "environment, on one stated scale. <b>Not:</b> which atom reacts, how fast, or by "
          "which mechanism.</p>")


def caption(rec):
    sc = rec.scale
    warm = "warm = negative" if sc.kind == "diverging" else "darker = more negative"
    return (f"{rec.provenance.method} · {sc.kind} scale {sc.vmin:+.0f} … "
            f"{sc.vmax:+.0f} kcal/mol, {warm} · surface genus {rec.genus}")


def molecule_panels(name, rec, atom_map):
    """The site is chosen by atom map, never by label: labels follow the maps
    on a mapped input since espmap #47, and were indices before it."""
    s = rec.site_by_map(atom_map)
    smi = rec.mapped_smiles
    return s.label, [
        '<div class="row"><div class="panel">' + svg_2d_view(rec, s.atoms) + "</div>",
        ("viewer", f"{name}_sticks", viewer_js(f"{name}_sticks", rec, smi, s.atoms,
                                               show_surface=False)),
        ("viewer", f"{name}_esp", viewer_js(f"{name}_esp", rec, smi, s.atoms)), "</div>"]


def small(name, atom_map):
    rec = load(FIX / name)
    label, panels = molecule_panels(name, rec, atom_map)
    blocks = panels + ["<h2 style='font-size:16px'>Sites</h2>",
                       site_table(rec, rec.site_labels()), LIMITS]
    out = OUT / f"p3_{name}.html"
    out.write_text(page(f"{name} (P3 fixture)", f"site {label} highlighted · {caption(rec)}",
                        blocks), encoding="utf-8")
    return out


def c18(atom_map=3):
    ext, fold = load(FIX / "Z3Z13-18OAc_extended"), load(FIX / "Z3Z13-18OAc_folded")
    label, ext_panels = molecule_panels("ext", ext, atom_map)
    _, fold_panels = molecule_panels("fold", fold, atom_map)
    site = label
    labels = ext.site_labels()           # overlapping sites are marked, not hidden
    blocks = (["<h2 style='font-size:16px'>Extended (genus 0)</h2>"]
              + ext_panels
              + ["<h2 style='font-size:16px'>Folded (genus 3)</h2>"]
              + fold_panels
              + (["<p class='limits'><b>Different colour scales</b> (different render keys: "
                  f"{caption(ext).split(' · ')[1]} vs {caption(fold).split(' · ')[1]}). "
                  "Compare the numbers in the table, not the colours.</p>"]
                 if not same_scale(ext, fold) else [])
              + ["<p class='limits'>Two conformers of one molecule. The site table's "
                 "<b>conformer range</b> is the spread between these two shapes alone: a "
                 "lower bound on the real spread, not a measurement of it.</p>",
                 "<h2 style='font-size:16px'>Sites over both conformers</h2>",
                 site_table([ext, fold], labels), LIMITS])
    out = OUT / "p3_c18.html"
    out.write_text(page("(3Z,13Z)-octadeca-3,13-dien-1-yl acetate (P3 fixtures)",
                        f"site {site} highlighted · {caption(ext)}", blocks), encoding="utf-8")
    return out


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    # by atom map: the enolate's central C, a nitrite O, the thiocyanate S, a benzene C
    for name, atom_map in (("enolate", 2), ("nitrite", 2), ("thiocyanate", 1), ("benzene", 1)):
        print(small(name, atom_map))
    print(c18())
