"""Read an Atlas teaching bundle (spa.teaching/0) into Arrowfield's model.

    from bundle import load
    b = load("sesiidae-homework-1.bundle.json")
    b.observations["spa.obs....gland"].display_state("spa.cmpd.z3_z13_18_oac")

The second of the two adapters (record.py reads espmap). Nothing here
interprets biology: every field the Atlas exports is kept as recorded, and a
card derives its broad category from those fields (EvidenceStatement.
display_state). A field this adapter does not know is refused, not dropped:
a new Atlas field must be modelled before Arrowfield shows a bundle that
carries it.
"""
import hashlib
import json
import math
from pathlib import Path

from model import (BundleCompound, Detection, EvidenceClass, EvidenceStatement, ObservedContrast,
                   Response, TeachingBundle)

BUNDLE = "spa.teaching/0"
SCHEMA_VERSIONS = ("0.1",)
TOP = {"bundle", "schema_version", "atlas_version", "case_id", "case_title", "comparisons",
       "compounds", "espmap", "observations"}
OBS = {"observation_id", "observation_type", "taxon", "statement", "context", "reference",
       "source_location", "evidence", "assay", "experiment_id", "materials", "treatment",
       "comparator", "dose", "trap", "dispenser", "outcome", "components"}
EVIDENCE = {"atlas_notes", "author_statement", "author_statement_form", "unresolved"}
OUTCOME = {"effect", "response", "significance"}
COMPONENT = {"atlas_id", "detected", "method", "identification", "abundance_class",
             "detection_limit", "measurements"}
RESPONSE = {"measure", "unit", "treatment_value", "comparator_value", "qualifier"}
COMPOUND = {"atlas_id", "shorthand", "systematic", "inchikey", "occurrences", "evidence_classes"}
OCCURRENCE = {"observation_id", "roles"}
# The Atlas's controlled words (its observation.schema.json and teaching-bundle.schema.json).
# A value outside them is refused: a view phrases each of these in its own words, and an
# unknown one would reach the page as Arrowfield's claim (probe round 3, #13).
EFFECTS = ("increased", "decreased", "no_effect", "caught_without_control")
SIGNIFICANCES = ("significant", "not_significant", "not_tested", "not_reported")
ROLES = ("treatment", "comparator", "detected", "not_detected", "not_tested", "antennal_active",
         "antennal_inconsistent", "antennal_inactive", "antennal_not_tested")
OBSERVATION_TYPES = ("chemical_measurement", "behavioural_response", "electrophysiology")
EVIDENCE_CLASSES = ("female_detected_gland", "female_detected_emission", "female_detected_combined",
                    "female_not_detected", "antennal_active", "antennal_inconsistent", "antennal_inactive",
                    "attractant_alone", "attractive_lure_component", "addition_increased",
                    "addition_decreased", "addition_no_effect", "female_pheromone_component")


def _number(v):
    """A finite real number; a bool is not one (True is an int to Python)."""
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def _count(v, nullable):
    return (v is None and nullable) or (isinstance(v, int) and not isinstance(v, bool) and v >= 0)
EVIDENCE_CLASS = {"evidence_class", "taxon_id", "observation_ids", "populations", "n_observations",
                  "n_significant", "n_tentative", "population_match"}


class InvalidBundle(ValueError):
    """A bundle Arrowfield will not show, with the reason."""


def _strict(path):
    def no_nan(name):
        raise InvalidBundle(f"{path}: {name} is not a number JSON allows")
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=no_nan)
    except json.JSONDecodeError as e:
        raise InvalidBundle(f"{path}: not JSON ({e})") from None


def _known(d, allowed, where):
    extra = set(d) - allowed
    if extra:
        raise InvalidBundle(f"{where}: field(s) {sorted(extra)} are not modelled by Arrowfield; "
                            "model them in model.py and bundle.py before showing this bundle")


def _response(r, where):
    if r is None:
        return None
    _known(r, RESPONSE, where)
    for k in ("treatment_value", "comparator_value"):
        v = r.get(k)
        if v is not None and not _number(v):
            raise InvalidBundle(f"{where}: {k} {v!r} is not a finite number")
    return Response(r["measure"], r["unit"], r.get("treatment_value"), r.get("comparator_value"),
                    r.get("qualifier"))


def _statement(o, roles):
    oid = o.get("observation_id", "?")
    _known(o, OBS, oid)
    ev = o.get("evidence") or {}
    _known(ev, EVIDENCE, f"{oid} evidence")
    out = o.get("outcome") or {}
    _known(out, OUTCOME, f"{oid} outcome")
    for key, words in (("effect", EFFECTS), ("significance", SIGNIFICANCES)):
        if out.get(key) is not None and out[key] not in words:
            raise InvalidBundle(f"{oid} outcome: {key} {out[key]!r} is not one of the Atlas's words {words}")
    dets = []
    for c in o.get("components") or []:
        _known(c, COMPONENT, f"{oid} component")
        if not isinstance(c["detected"], bool):
            raise InvalidBundle(f"{oid} component {c['atlas_id']}: detected {c['detected']!r} is not true or false")
        lim = c.get("detection_limit")
        if lim is not None and not (isinstance(lim, dict) and _number(lim.get("value"))):
            raise InvalidBundle(f"{oid} component {c['atlas_id']}: detection_limit {lim!r} has no numeric value")
        dets.append(Detection(c["atlas_id"], c["detected"], c.get("method"),
                              c.get("identification"), c.get("abundance_class"),
                              c.get("detection_limit"), tuple(c.get("measurements") or ())))
    if o["observation_type"] not in OBSERVATION_TYPES:
        raise InvalidBundle(f"{oid}: observation_type {o['observation_type']!r} is not one of {OBSERVATION_TYPES}")
    tax = o["taxon"]
    return EvidenceStatement(
        observation_id=oid, observation_type=o["observation_type"],
        taxon=tax["accepted_name"], taxon_as_reported=tax.get("name_as_reported"),
        statement=o["statement"], context=dict(o.get("context") or {}),
        reference=dict(o.get("reference") or {}), source_location=o.get("source_location"),
        assay=o.get("assay"), experiment_id=o.get("experiment_id"),
        materials=tuple(o.get("materials") or ()), treatment=tuple(o.get("treatment") or ()),
        comparator=o.get("comparator"), dose=o.get("dose"), trap=o.get("trap"),
        dispenser=o.get("dispenser"), effect=out.get("effect"), significance=out.get("significance"),
        response=_response(out.get("response"), f"{oid} response"), detections=tuple(dets),
        roles=roles.get(oid, {}), author_statement=ev.get("author_statement"),
        author_statement_form=ev.get("author_statement_form"), atlas_notes=ev.get("atlas_notes"),
        unresolved=tuple(ev.get("unresolved") or ()))


def load(path):
    """TeachingBundle from a spa.teaching/0 file; InvalidBundle otherwise."""
    j = _strict(path)
    if not isinstance(j, dict) or j.get("bundle") != BUNDLE:
        raise InvalidBundle(f"{path}: not a {BUNDLE} bundle")
    if j.get("schema_version") not in SCHEMA_VERSIONS:
        raise InvalidBundle(f"{path}: schema_version {j.get('schema_version')!r}; this reads "
                            f"{SCHEMA_VERSIONS}. Update bundle.py before showing it.")
    _known(j, TOP, str(path))
    if not j.get("atlas_version"):
        raise InvalidBundle(f"{path}: no atlas_version: a bundle must name the Atlas commit it came from")
    ids = [o.get("observation_id") for o in j["observations"]]
    if len(ids) != len(set(ids)):
        raise InvalidBundle(f"{path}: an observation id repeats")
    roles = {}
    compounds = {}
    for c in j["compounds"]:
        _known(c, COMPOUND, f"compound {c.get('atlas_id')}")
        classes = []
        for e in c.get("evidence_classes") or []:
            _known(e, EVIDENCE_CLASS, f"compound {c['atlas_id']} evidence class")
            where = f"compound {c['atlas_id']} evidence class {e.get('evidence_class')!r}"
            if e["evidence_class"] not in EVIDENCE_CLASSES:
                raise InvalidBundle(f"{where}: not one of the Atlas's classes {EVIDENCE_CLASSES}")
            if not (_count(e["n_observations"], False) and _count(e.get("n_significant"), True)
                    and _count(e.get("n_tentative"), True)):
                raise InvalidBundle(f"{where}: counts must be non-negative integers")
            classes.append(EvidenceClass(e["evidence_class"], e["taxon_id"], tuple(e["observation_ids"]),
                                         tuple(e.get("populations") or ()), e["n_observations"],
                                         e.get("n_significant"), e.get("n_tentative"),
                                         e.get("population_match")))
        compounds[c["atlas_id"]] = BundleCompound(c["atlas_id"], c.get("shorthand"), c.get("systematic"),
                                                  c.get("inchikey"), tuple(classes))
        for occ in c.get("occurrences") or []:
            _known(occ, OCCURRENCE, f"compound {c['atlas_id']} occurrence")
            if not (isinstance(occ["roles"], list) and all(isinstance(r, str) for r in occ["roles"])):
                raise InvalidBundle(f"{path}: {c['atlas_id']} in {occ['observation_id']}: roles must be a list of "
                                    f"role words, not {occ['roles']!r}")
            if occ["observation_id"] not in ids:
                raise InvalidBundle(f"{path}: {c['atlas_id']} occurs in {occ['observation_id']}, "
                                    "which the bundle does not carry")
            bad = [r for r in occ["roles"] if r not in ROLES]
            if bad:
                raise InvalidBundle(f"{path}: {c['atlas_id']} in {occ['observation_id']}: role(s) {bad} are not "
                                    f"the Atlas's role words {ROLES}")
            roles.setdefault(occ["observation_id"], {})[c["atlas_id"]] = tuple(occ["roles"])
    contrasts = []
    for k in j["comparisons"]:
        missing = [i for i in k["observation_ids"] if i not in ids]
        if missing:
            raise InvalidBundle(f"{path}: comparison {k.get('case_label')!r} names {missing}, "
                                "which the bundle does not carry")
        oc = k.get("observed_contrast") or {}
        contrasts.append(ObservedContrast(k["case_label"], tuple(k["observation_ids"]),
                                          oc.get("a"), oc.get("b")))
    obs = {o["observation_id"]: _statement(o, roles) for o in j["observations"]}
    return TeachingBundle(j["case_id"], j["case_title"], j["atlas_version"],
                          hashlib.sha256(Path(path).read_bytes()).hexdigest(), obs, compounds,
                          tuple(contrasts), tuple(j.get("espmap") or ()))
