"""The homework page: the plan (sesiidae-homework.md as data) and its sources
assembled into one page from the views already built.

Specified by tests/pedagogy/homework_page_set.json (#19). The page settles
what no single view can: one role per molecule, from the plan and checked
against the data; acts A-D in order, each prediction before its views and
the expected answers only for the instructor; one bundle, one results file,
one analysis key, named in a footer; the study's outcome words only inside
study cards, and a difference card at a tested site beside that site's study
card. Plan text is authored, shown verbatim; the page's own words join no
ESP difference to a biological result.

    page_html(plan, sources, audience="student") -> str, or model.Refused
"""
import html
import re
import unicodedata

from model import KNOWN_VERDICT_RULES, STUDY_FOLD_RULES, STUDY_KEYS, STUDY_LIMITS, Refused, compare, protocol, protocol_problems, protocol_words

ACTS = ("A", "B", "C", "D")
ROLE_WORDS = {"atlas_observed": "Atlas compound",
              "teaching_control": "teaching control: a constructed molecule, not an Atlas compound"}
VIEWS = ("study_card", "difference_card", "pair_view", "evidence_card", "contrast_card", "class_panel")
AUDIENCES = ("student", "instructor")
# identifier grammars (the round-8 probe on #21): ASCII, so nothing looks like something else
SLUG = re.compile(r"[a-z][a-z0-9_]*")
VERSION = re.compile(r"[0-9][0-9A-Za-z.+_-]*")
ANALYSIS_KEY = re.compile(r"[0-9a-f]{8,64}")
# the solvents espmap writes (esp_map.SOLVENTS, v0.2.1 and v0.3.0 alike; None is the gas phase),
# and "gas", which the views also read as gas: a closed vocabulary, so no synonym is guessed at
GAS_SOLVENTS = ("vacuum", "gas")
SOLVENTS = GAS_SOLVENTS + ("water", "dmso", "methanol", "acetonitrile", "chloroform", "toluene", "hexane")


def _e(x):
    return html.escape(str(x), quote=False)


def _v(x):
    """Authored or record text, verbatim and marked: the page's bans read only
    the words outside these spans."""
    return f'<span class="verbatim">{_e(x)}</span>'


def _attr(x):
    return html.escape(str(x), quote=True)


# Default_Ignorable_Code_Point (Unicode DerivedCoreProperties), and the symbols drawn as nothing
# by design (U+2800 BRAILLE PATTERN BLANK, U+1D159 MUSICAL SYMBOL NULL NOTEHEAD): characters a
# reader cannot see
_BLANK = ((0x00AD, 0x00AD), (0x034F, 0x034F), (0x061C, 0x061C), (0x115F, 0x1160), (0x17B4, 0x17B5),
          (0x180B, 0x180F), (0x200B, 0x200F), (0x202A, 0x202E), (0x2060, 0x206F), (0x2800, 0x2800),
          (0x3164, 0x3164), (0xFE00, 0xFE0F), (0xFEFF, 0xFEFF), (0xFFA0, 0xFFA0), (0xFFF0, 0xFFF8),
          (0x1BCA0, 0x1BCA3), (0x1D159, 0x1D159), (0x1D173, 0x1D17A), (0xE0000, 0xE0FFF))


def _visible(x):
    """Text with something to read: a letter, number, punctuation mark or symbol
    a reader can see. Whitespace, fillers, variation selectors, a lone combining
    mark or private-use character read like "" (the reruns on #21, 318da98 and
    2a67bdd)."""
    return isinstance(x, str) and any(
        unicodedata.category(ch)[0] in "LNPS" and not any(lo <= ord(ch) <= hi for lo, hi in _BLANK) for ch in x)


def _writable(x):
    """Text a page can carry: no lone surrogate (the page could not be written
    as UTF-8), no control character but tab, newline and carriage return, and no
    noncharacter (U+FDD0-FDEF, U+xFFFE/xFFFF: an HTML parse error), and no bidi
    embedding, override or isolate (U+202A-202E, U+2066-2069), which reorders the
    page's own words after it (the round-8 probe on #21)."""
    return not _UNCARRIABLE.search(x)


# the controls (Cc) but tab, newline and CR; the surrogates (Cs); the bidi embeddings, overrides
# and isolates; the noncharacters
_UNCARRIABLE = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f\x7f-\x9f\ud800-\udfff‪-‮⁦-⁩﷐-﷯"
                          + "".join(f"{chr(p * 0x10000 + 0xFFFE)}{chr(p * 0x10000 + 0xFFFF)}" for p in range(17))
                          + "]")


def _need(cond, why):
    if not cond:
        raise Refused(why)


def _check_plan(plan, sources, audience):
    _need(audience in AUDIENCES, f"audience must be one of {AUDIENCES}")
    _need(isinstance(plan, dict) and isinstance(plan.get("molecules"), dict) and isinstance(plan.get("acts"), list),
          "the plan must carry its molecules and acts")
    _need(isinstance(sources, dict), "the sources must be a mapping")
    bundle, study, ens = sources.get("bundle"), sources.get("study"), sources.get("ensembles")
    from model import StudyResults, TeachingBundle, bundle_problems
    _need(isinstance(bundle, TeachingBundle) and not bundle_problems(bundle), "the sources need one valid bundle")
    _need(isinstance(study, StudyResults), "the sources need one StudyResults")
    _need(study.verdict_rule in KNOWN_VERDICT_RULES,
          f"the study's verdict rule {study.verdict_rule!r} is not one Arrowfield knows: no study card can be shown")
    _need(isinstance(ens, dict) and isinstance(sources.get("conformers", {}), dict), "ensembles and conformers must be mappings")
    # source text the page shows (the round-7 probe on #21): each is text a page can carry, and
    # a name or reference has something to read
    for slug in plan["molecules"]:
        _need(isinstance(slug, str) and _writable(slug), f"molecule key {slug!r} is not text a page can carry")
        _need(_visible(slug) and "," not in slug,
              f"molecule key {slug!r} names no molecule a reader can read, or names two in data-molecules")
        # a slug is an identifier: lowercase ASCII, so no two molecules look alike on the page
        _need(SLUG.fullmatch(slug) is not None and slug != "none",
              f"molecule key {slug!r} is not a slug ({SLUG.pattern}): it could read as another molecule's")
    _need(isinstance(study.plan, str) and _writable(study.plan) and _visible(study.plan) and study.plan.strip() != "None",
          f"the study's plan reference {study.plan!r} is not one a page can show")
    for t in study.sites:
        _need(all(isinstance(f, str) and _writable(f) for f in (t.contrast, t.site, t.x, t.y)),
              f"a study row names {t.contrast!r} / {t.site!r} in text a page cannot carry")
    # roles: one per molecule, from the plan, checked against the data
    for slug, m in plan["molecules"].items():
        _need(isinstance(m, dict) and m.get("role") in ROLE_WORDS,
              f"{slug}: role {m.get('role') if isinstance(m, dict) else m!r} is not one of {tuple(ROLE_WORDS)}")
        src = ens.get(slug)
        _need(isinstance(src, tuple) and len(src) == 2, f"{slug}: no ensemble in the sources")
        source_id = src[0].atlas_id
        if m["role"] == "atlas_observed":
            _need(m.get("atlas_id") in bundle.compounds,
                  f"{slug}: called an Atlas compound, but {m.get('atlas_id')!r} is not a compound of the bundle")
            _need(source_id == m["atlas_id"],
                  f"{slug}: the plan says {m['atlas_id']}, its ensemble says {source_id}")
        else:
            _need(m.get("atlas_id") not in bundle.compounds,
                  f"{slug}: a teaching control the bundle carries as an Atlas compound is not a control")
            _need(isinstance(source_id, str) and source_id.startswith("arrowfield.control."),
                  f"{slug}: called a teaching control, but its ensemble's source is {source_id}")
    # acts: exactly A-D, in order; a prediction in each; act D holds no view
    acts = plan["acts"]
    _need(all(isinstance(a, dict) for a in acts) and [a.get("id") for a in acts] == list(ACTS),
          f"the acts must be {', '.join(ACTS)}, in that order and once each")
    # plan text is authored text: a field that is not text is refused, never shown as str(value)
    _need(plan.get("title") is None or isinstance(plan["title"], str), "plan.title is not text")
    _need(plan.get("title") is None or _writable(plan["title"]), "plan.title holds a character no page can carry")
    _need(plan.get("case_id") is None or isinstance(plan["case_id"], str), "plan.case_id is not text")
    _need(plan.get("case_id") is None or _writable(plan["case_id"]), "plan.case_id holds a character no page can carry")
    for a in acts:
        for field in ("title", "expected"):
            _need(a.get(field) is None or isinstance(a[field], str), f"act {a['id']}.{field} is not text")
            _need(a.get(field) is None or _writable(a[field]), f"act {a['id']}.{field} holds a character no page can carry")
        _need(isinstance(a.get("tasks", []), list) and all(isinstance(t, str) for t in a.get("tasks", [])),
              f"act {a['id']}.tasks must be a list of text")
        for i, t in enumerate(a.get("tasks", [])):
            _need(_writable(t), f"act {a['id']}.tasks[{i}] holds a character no page can carry")
            _need(_visible(t), f"act {a['id']}.tasks[{i}] has nothing to read")
        _need(not isinstance(a.get("predict"), str) or _writable(a["predict"]),
              f"act {a['id']}.predict holds a character no page can carry")
        _need(_visible(a.get("predict")),
              f"act {a['id']} has no prediction: the reveal would come first")
        _need(isinstance(a.get("blocks"), list) and isinstance(a.get("contrasts"), list), f"act {a['id']} is malformed")
        _need(all(isinstance(c, str) and _writable(c) and _visible(c) for c in a["contrasts"]),
              f"act {a['id']} names a contrast in text a page cannot show")
        _need(a["id"] != "D" or not a["blocks"], "act D is about what free-molecule ESP cannot show: it holds no view")
    # one analysis key across the page
    views = [v for _, vs in ens.values() for v in vs] + list((sources.get("conformers") or {}).values())
    for v in views:
        ver = v.provenance.espmap_version
        # identifiers have a grammar, so none holds the footer's separators or its "not recorded"
        _need(ver is None or (isinstance(ver, str) and VERSION.fullmatch(ver) is not None),
              f"{v.provenance.source}: espmap version {ver!r} is not a version ({VERSION.pattern})")
        _need(isinstance(v.analysis_key, str) and ANALYSIS_KEY.fullmatch(v.analysis_key) is not None,
              f"{v.provenance.source}: analysis key {v.analysis_key!r} is not one ({ANALYSIS_KEY.pattern})")
        pv = v.provenance
        # raw settings are provenance text (the user's decision after round 11): both stated, as
        # text with something to read, or neither; never checked as method names
        _need(all(isinstance(m, str) and _writable(m) and _visible(m) and m.strip() != "None" for m in (pv.xc, pv.basis))
              or (pv.xc, pv.basis) == (None, None),
              f"{pv.source}: settings {pv.xc!r}/{pv.basis!r} are not provenance a page can show (both, or neither)")
        _need(pv.solvent is None or pv.solvent in SOLVENTS,
              f"{pv.source}: solvent {pv.solvent!r} is not one espmap writes ({', '.join(SOLVENTS)}, or None for gas)")
        # the protocol is the canonical one its key digests, and the raw phase agrees with it
        problems = protocol_problems(pv)
        _need(not problems, problems[0] if problems else "")
        _need(protocol(pv).get("solvent") is None or protocol(pv)["solvent"] in SOLVENTS[len(GAS_SOLVENTS):],
              f"{pv.source}: its protocol's solvent {protocol(pv).get('solvent')!r} is not one espmap writes")
    for key, v in (sources.get("conformers") or {}).items():
        # the pair view names its records' method: a record it shows must state one
        _need(v.provenance.xc is not None, f"{key}: the record states no method, and the pair view names it")
    keys = {v.analysis_key for v in views}
    _need(len(keys) == 1, f"the page's molecules carry {len(keys)} analysis keys; the homework uses one throughout")
    # the records are the study's: its study cards scope its verdicts under its protocol, never
    # another's (the fresh probe on #21, 565e303)
    _need(keys == {STUDY_KEYS.get(study.verdict_rule)},
          f"the page's records are under analysis key {next(iter(keys))}, and {study.verdict_rule} was decided under "
          f"{STUDY_KEYS.get(study.verdict_rule)}: they are not that study's records")
    # and each of the study's Atlas molecules is shown on its member count: a count other than the
    # study's n_ok proves the members are not the study's (the user's decision after the fresh probe
    # on e17af42). An equal count proves no more than that; the results carry no member ids. The
    # teaching control is not one of the study's Atlas molecules
    for slug, m in plan["molecules"].items():
        if m["role"] != "atlas_observed":
            continue
        n_ok = (study.molecules.get(slug) or {}).get("n_ok")
        _need(isinstance(n_ok, int) and not isinstance(n_ok, bool),
              f"{slug}: {study.verdict_rule} reports no member count for this molecule ({n_ok!r})")
        _need(len(ens[slug][1]) == n_ok,
              f"{slug}: the page's ensemble has {len(ens[slug][1])} members; this study reports {n_ok} for this "
              "molecule, so its cards would not be the evidence the study decided on")
        # and on the study's folded / extended partition, which its extended-only rows refer to (the
        # user's decision after the fresh probe on 9e2a397). Equal counts prove no more: the results
        # carry no per-member flags
        n_folded = study.molecules[slug].get("n_folded")
        _need(isinstance(n_folded, int) and not isinstance(n_folded, bool),
              f"{slug}: {study.verdict_rule} reports no folded count for this molecule ({n_folded!r})")
        flags = [v.folded for v in ens[slug][1]]
        _need(all(isinstance(f, bool) for f in flags),
              f"{slug}: a member's folded state is not stated, so the page cannot count the study's partition")
        _need(sum(flags) == n_folded,
              f"{slug}: page has {sum(flags)} folded / {n_ok - sum(flags)} extended members; study reports "
              f"{n_folded} folded / {n_ok - n_folded} extended")
    # each molecule's ensemble is that molecule's: every record states its structure
    # (sites.identity.mapped_smiles) and an ensemble holds one; two molecules never share one (the
    # fresh probe on 9e2a397). Chemical identity, which the records state; not member identity
    structures = {}
    for slug, (_, vs) in ens.items():
        s = {v.mapped_smiles for v in vs}
        _need(len(s) == 1 and None not in s,
              f"{slug}: its ensemble's records state {len(s - {None})} structures"
              f"{' (and some none)' if None in s else ''}; one molecule's ensemble states one")
        structures.setdefault(next(iter(s)), []).append(slug)
    shared = [slugs for slugs in structures.values() if len(slugs) > 1]
    _need(not shared, f"{' and '.join(shared[0]) if shared else ''}: two molecules of the page with one structure")
    # and an Atlas molecule's records are its compound: their structure's InChIKey is the one the
    # bundle states (the fresh probe on dfeffcc: two molecules' records exchanged passed every check)
    for slug, m in plan["molecules"].items():
        if m["role"] == "atlas_observed":
            stated = bundle.compounds[m["atlas_id"]].inchikey
            got = _inchikey(next(iter({v.mapped_smiles for v in ens[slug][1]})))
            _need(isinstance(stated, str) and got == stated,
                  f"{slug}: its records' structure has InChIKey {got}; the bundle states {stated} for {m['atlas_id']}")
    return bundle, study, ens, keys, views


def _inchikey(mapped_smiles):
    """The InChIKey of a record's structure (its atom maps dropped), or None if it does not parse."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(mapped_smiles) if isinstance(mapped_smiles, str) else None
    if mol is None:
        return None
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToInchiKey(mol) or None


def _study_context(study, plan, x, y, views):
    # the plan commit only when the reference names exactly one: two are ambiguous, and
    # are never resolved by position (the probe on #21); the card then says it is not recorded
    commits = re.findall(r"@ ([0-9a-f]{7,40})\b", study.plan or "")
    # every record's protocol is verified and they share the study's key, so any record names the
    # study's protocol; its solvent is stated, never left "not stated"
    method, solvent = protocol_words(views[0].provenance)
    lim = {"method": method, "phase": "gas" if solvent is None else "condensed",
           "solvent": ("none: the gas phase only (solvent effects not studied)" if solvent is None
                       else f"PCM {solvent} only (other solvents and the gas phase not studied)"),
           **STUDY_LIMITS.get(study.verdict_rule, {})}
    return {"provenance": {"results_sha256": study.sha256, "plan_commit": commits[0] if len(commits) == 1 else None},
            "molecules": {s: plan["molecules"][s]["role"] for s in (x, y)}, "limitations": lim}


def _block(act, k, plan, bundle, study, ens, sources, views):
    """(html of the view, its molecules' slugs)."""
    import views as V
    _need(isinstance(k, dict) and k.get("view") in VIEWS, f"act {act['id']}: a block of unknown kind")
    _need(not ({"role", "roles"} & set(k)), f"act {act['id']}: a block carries its own role; roles are the plan's")
    mols = plan["molecules"]
    view = k["view"]
    if view in ("study_card", "difference_card"):
        _need(k.get("contrast") in act["contrasts"],
              f"act {act['id']}: a {view} for {k.get('contrast')!r}, not one of this act's contrasts")
        rows = [s for s in study.sites if s.contrast == k["contrast"] and s.population == "all"]
        _need(rows, f"the study has no contrast {k['contrast']!r}")
        x, y = rows[0].x, rows[0].y
    if view == "study_card":
        try:
            row = study.site(k["contrast"], k.get("site"))
        except KeyError:
            raise Refused(f"the study has no row for {k['contrast']} / {k.get('site')}") from None
        _need(x in mols and y in mols, f"{x} or {y} is not a molecule of the plan")
        return V.study_card_html(row, _study_context(study, plan, x, y, views)), (x, y)
    if view == "difference_card":
        _need((k.get("x"), k.get("y")) == (x, y),
              f"act {act['id']}: a card for {k.get('x')} -> {k.get('y')} under {k['contrast']}, which compares {x} -> {y}")
        _need(x in mols and y in mols, f"{x} or {y} is not a molecule of the plan")
        site = k.get("site")
        row = next((s for s in rows if s.site == site), None)
        if row is not None and row.tier in ("primary", "secondary"):
            _need(any(b.get("view") == "study_card" and b.get("contrast") == k["contrast"] and b.get("site") == site
                      for b in act["blocks"] if isinstance(b, dict)),
                  f"act {act['id']}: a difference card at {site}, which {k['contrast']} tests, needs that site's "
                  "study card in the act")
        _need(isinstance(site, str), f"act {act['id']}: a difference card at site {site!r}, which is not a site name")
        try:
            card = compare(ens[x][1], ens[y][1], site, x, y)
        except KeyError:
            raise Refused(f"act {act['id']}: a difference card at {site!r}, which {x} and {y} do not carry") from None
        ctx = {"study": study, "contrast": k["contrast"], "x": x, "y": y,
               "roles": {s: mols[s]["role"] for s in (x, y)}}
        return V.card_html(card, ctx), (x, y)
    if view == "pair_view":
        refs = (k.get("left"), k.get("right"))
        _need(all(isinstance(r, dict) and r.get("slug") in mols for r in refs), "a pair names a molecule the plan does not list")
        conf = sources.get("conformers") or {}
        keys = [f"{r['slug']}/{r.get('conformer')}" for r in refs]
        _need(all(c in conf for c in keys), f"a pair's conformer is not in the sources ({keys})")
        # each record is the member of its molecule's ensemble that it names: the same record file,
        # folded flag and structure, so the ensemble checks hold for the pair too (the fresh probe on
        # 844d3a6: a flag flipped, or another molecule's records, passed through the pair)
        for r, c in zip(refs, keys):
            v = conf[c]
            m = next((u for u in ens.get(r["slug"], (None, ()))[1] if u.provenance.source == r.get("conformer")), None)
            _need(m is not None, f"{c}: not a member of {r['slug']}'s ensemble")
            _need(v.provenance.sha256 is not None and v.provenance.sha256 == m.provenance.sha256,
                  f"{c}: not the record {r['slug']}'s ensemble holds as {r.get('conformer')}")
            _need(v.folded is m.folded, f"{c}: folded {v.folded}, where {r['slug']}'s ensemble says {m.folded}")
            _need(v.mapped_smiles == m.mapped_smiles, f"{c}: another structure than {r['slug']}'s ensemble")
        html_ = V.pair_view_html(conf[keys[0]], conf[keys[1]], {})
        rule = STUDY_FOLD_RULES.get(study.verdict_rule)
        if rule and any(isinstance(conf[c].folded, bool) for c in keys):
            # the labels stay "folded" / "extended"; what they mean is the study's rule, said beside them
            # (the user's browser check of the 9e2a397 preview: the folded conformer looked the longer; the
            # user's wording after round 17, since folded members are more compact on average)
            html_ += (f'<p class="limits">Folded here is this study\'s rule: {_e(rule)}. It is defined by contact, '
                      "not by length or compactness, so a folded conformer can look the longer of the two.</p>")
        return html_, (refs[0]["slug"], refs[1]["slug"])
    if view == "evidence_card":
        slug = k.get("molecule")
        _need(slug in mols, f"{slug!r} is not a molecule of the plan")
        st = bundle.observations.get(k.get("observation_id"))
        _need(st is not None, f"{k.get('observation_id')!r} is not an observation of the bundle")
        aid = mols[slug].get("atlas_id")
        ctx = {"bundle_object": bundle, "roles": {aid: mols[slug]["role"]} if aid else {}}
        return V.evidence_card_html(st, aid, ctx), (slug,)
    if view == "contrast_card":
        k_ = next((c for c in bundle.contrasts if c.label == k.get("label")), None)
        _need(k_ is not None, f"{k.get('label')!r} is not a contrast of the bundle")
        return V.contrast_card_html(k_, [bundle.observations[o] for o in k_.observation_ids],
                                    {"bundle_object": bundle}), ()
    slug = k.get("molecule")
    _need(slug in mols, f"{slug!r} is not a molecule of the plan")
    m = mols[slug]
    cid = m["atlas_id"] if m["role"] == "atlas_observed" else ens[slug][0].atlas_id
    ctx = {"bundle_object": bundle, "roles": {cid: m["role"]}}
    return V.evidence_class_panel_html(cid, k.get("taxon_id"), ctx), (slug,)


def page_html(plan, sources, audience="student"):
    """The homework as one HTML page, or model.Refused."""
    import views as V
    bundle, study, ens, keys, views = _check_plan(plan, sources, audience)
    mols = plan["molecules"]
    # an optional text field that is None or has nothing to read is absent: never "None", never
    # an empty heading
    body = [f'<h1>{_v(plan["title"])}</h1>'] if _visible(plan.get("title")) else []
    for act in plan["acts"]:
        parts = [f'<section class="act" data-act="{act["id"]}"><h2>Act {act["id"]}{": " + _v(act["title"]) if _visible(act.get("title")) else ""}</h2>',
                 f'<div class="predict"><b>Predict first:</b> {_v(act["predict"])}</div>']
        if act.get("tasks"):
            parts.append("<ol class='tasks'>" + "".join(f"<li>{_v(t)}</li>" for t in act["tasks"]) + "</ol>")
        for k in act["blocks"]:
            inner, slugs = _block(act, k, plan, bundle, study, ens, sources, views)
            roles = [mols[s]["role"] for s in slugs]
            legend = "".join(f'<div class="role">{_v(s)}: {_e(ROLE_WORDS[r])}</div>' for s, r in zip(slugs, roles))
            parts.append(f'<section class="block" data-view="{_attr(k["view"])}" '
                         f'data-molecules="{_attr(",".join(slugs))}" data-roles="{_attr(",".join(roles))}">'
                         f"{legend}{inner}</section>")
        if audience == "instructor" and _visible(act.get("expected")):
            parts.append(f'<div class="expected"><b>For the instructor:</b> {_v(act["expected"])}</div>')
        parts.append("</section>")
        body.append("".join(parts))
    versions = sorted({v.provenance.espmap_version or "not recorded" for v in views})
    commit = bundle.atlas_version
    atlas = (f"Atlas commit {_e(commit[:7])}" if re.fullmatch(r"[0-9a-f]{7,40}", commit)
             else "Provenance incomplete: the Atlas commit is not recorded")
    body.append(f'<div class="page-provenance" data-atlas="{_attr(commit)}" data-results-sha256="{_attr(study.sha256)}" '
                f'data-espmap="{_attr(",".join(versions))}" data-analysis-keys="{_attr(",".join(sorted(keys)))}">'
                f"{atlas} · results sha256 {_e(study.sha256[:12])}… · espmap {_e(', '.join(versions))} · "
                f"analysis key {_e(', '.join(k[:8] for k in sorted(keys)))}</div>")
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f"<title>{_e(plan['case_id'] if _visible(plan.get('case_id')) else 'homework')}</title><style>{V.CSS}</style>"
            f'<script src="{V.THREEDMOL}"></script></head><body>{"".join(body)}</body></html>')
    # the last line: whatever source text no check above names, the page refuses rather than
    # carry a character it cannot (a crash at the write, or an HTML parse error)
    bad = _UNCARRIABLE.search(page)
    _need(bad is None, f"a source field holds {bad.group()!r} (near {page[max(0, bad.start() - 40):bad.start()]!r}): "
                       "not text a page can carry" if bad is not None else "")
    # and every piece of source or plan text shown has something to read and is not the literal
    # "None", which a reader cannot tell from a missing value (the round-8 probe on #21)
    for m in re.finditer(r'<span class="verbatim">([^<]*)</span>', page):
        shown = html.unescape(m.group(1))
        _need(_visible(shown) and shown.strip() != "None",
              f"source text {shown!r} is shown (near {html.unescape(page[max(0, m.start() - 60):m.start()])!r}), "
              "but a reader cannot tell it from a missing value")
    return page
