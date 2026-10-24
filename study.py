"""Read a frozen study: its members and its own verdicts.

    from study import load_ensemble, load_results
    source, views = load_ensemble("controls/espmap/z3_z13_18_al/<protocol>/ensemble.json")
    results = load_results("espmap-sesiidae-001.results.json")
    results.site("C_position", "C4-H").outcome          # "resolved", the study's word

Two things this adapter will not do, because they are the study's:
- decide a verdict. A study's result says which rule decided it
  (verdict_rule); Arrowfield shows that verdict and never replaces it with
  model.compare()'s teaching verdict or espmap's generic resolved_under;
- decide what "folded" means. A member's folded flag is the producer's call,
  recorded in the ensemble index; it is passed into MoleculeView.folded as
  given (frozen study: contact OR genus > 0; espmap v0.3: contact).

Members are read from an Atlas-format index (spa.espmap_ensemble/0): the
Atlas's computed/espmap/ or Arrowfield's controls/espmap/. Each member is
checked against the index's sha256 and analysis key before record.py reads it.
"""
import hashlib
import json
import math
from dataclasses import replace
from pathlib import Path

from model import EnsembleSource, StudyResults, StudySite
from model import study_row_inconsistency as _consistent, study_row_malformed
from record import load as load_record

ENSEMBLE_SCHEMA = "spa.espmap_ensemble/0"
from model import KNOWN_VERDICT_RULES as VERDICT_RULES  # noqa: E402
OUTCOMES = ("resolved", "different on average", "unresolved", "not evaluable", "no data")
SITE_KEYS = {"outcome", "delta", "delta_other_rule", "partition_sensitivity",
             "partition_sensitivity_max_member", "range_x", "range_y", "mean_x", "mean_y",
             "bootstrap_interval", "n", "conditions", "why"}
class InvalidStudy(ValueError):
    """A study file Arrowfield will not show, with the reason."""


def _sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _strict(path):
    def no_nan(name):
        raise InvalidStudy(f"{path}: {name} is not a number JSON allows")
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=no_nan)
    except json.JSONDecodeError as e:
        raise InvalidStudy(f"{path}: not JSON ({e})") from None


def load_ensemble(index_path):
    """(EnsembleSource, [MoleculeView]) for one molecule's ensemble, members
    in index order, each MoleculeView.folded set from the index."""
    index_path = Path(index_path)
    ens = _strict(index_path)
    if ens.get("schema") != ENSEMBLE_SCHEMA:
        raise InvalidStudy(f"{index_path}: not a {ENSEMBLE_SCHEMA} index")
    d = index_path.parent
    views = []
    for m in ens["members"]:
        rec = d / m["record"]
        if _sha256(rec) != m["record_sha256"]:
            raise InvalidStudy(f"{rec}: sha256 differs from the index; the record changed after indexing")
        v = load_record(rec.with_suffix(""), name=f"{d.parent.name} {m['conformer']}")
        if v.analysis_key != ens["analysis_key"]:
            raise InvalidStudy(f"{rec}: analysis_key {v.analysis_key} is not the ensemble's {ens['analysis_key']}")
        if not isinstance(m.get("folded"), bool):
            raise InvalidStudy(f"{rec}: the index states no folded flag for {m['conformer']}")
        views.append(replace(v, folded=m["folded"]))
    sel = ens.get("selection", {})
    source = EnsembleSource(ens["atlas_id"], ens["protocol"], ens["analysis_key"], ens["espmap"]["version"],
                            _sha256(index_path), sel.get("method", ""), sel.get("notes"))
    return source, views


def _finite(v, where):
    if v is not None and not (isinstance(v, (int, float)) and math.isfinite(v)):
        raise InvalidStudy(f"{where}: {v!r} is not a finite number")
    return v


def _pair(v, where):
    if v is None:
        return None
    return tuple(_finite(x, where) for x in v)


def load_results(path):
    """StudyResults from a frozen study's results file (the Atlas's
    analyses/<study>/results.json). A partial (in-progress) file, an unknown
    verdict rule or outcome, or a non-finite number is refused."""
    j = _strict(path)
    rule = j.get("verdict_rule")
    if rule not in VERDICT_RULES:
        raise InvalidStudy(f"{path}: verdict_rule {rule!r} is not one Arrowfield knows how to show "
                           f"({VERDICT_RULES}); never show a study's verdicts under a rule it did not use")
    if j.get("partial"):
        raise InvalidStudy(f"{path}: a partial results file (members pending: {j.get('pending')}) "
                           "is a test of the scripts, never a result")
    sites = []
    for c in j["contrasts"]:
        for key, r in c["sites"].items():
            where = f"{path} {c['id']} {c['population']} {key}"
            extra = set(r) - SITE_KEYS
            if extra:
                raise InvalidStudy(f"{where}: field(s) {sorted(extra)} are not modelled")
            if r["outcome"] not in OUTCOMES:
                raise InvalidStudy(f"{where}: outcome {r['outcome']!r} is not one of {OUTCOMES}")
            why = _consistent(rule, r["outcome"], r.get("conditions") or {},
                              _pair(r.get("bootstrap_interval"), where))
            if why:
                raise InvalidStudy(f"{where}: an internally inconsistent row, refused as corrupt: {why}")
            tier, site = key.split(":", 1)
            sites.append(StudySite(
                verdict_rule=rule, contrast=c["id"], x=c["X"], y=c["Y"], population=c["population"],
                tier=tier, site=site, outcome=r["outcome"], n=tuple(r.get("n") or (0, 0)),
                delta=_finite(r.get("delta"), where), delta_other_rule=_finite(r.get("delta_other_rule"), where),
                partition=_finite(r.get("partition_sensitivity"), where),
                partition_max_member=_finite(r.get("partition_sensitivity_max_member"), where),
                range_x=_pair(r.get("range_x"), where), range_y=_pair(r.get("range_y"), where),
                mean_x=_finite(r.get("mean_x"), where), mean_y=_finite(r.get("mean_y"), where),
                bootstrap_interval=_pair(r.get("bootstrap_interval"), where),
                conditions=dict(r.get("conditions") or {})))
    for site in sites:
        why = study_row_malformed(site)
        if why:
            raise InvalidStudy(f"{path} {site.contrast} {site.population} {site.site}: a malformed row, "
                               f"refused as corrupt: {why}")
    return StudyResults(j["plan"], rule, _sha256(path), dict(j["molecules"]), tuple(sites))
