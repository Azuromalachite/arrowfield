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

import model
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
# The Atlas's controlled words and every value check live in model.py
# (bundle_problems), shared with bundles built by hand: this module only
# parses, without coercing a value into the shape it should have had.
EFFECTS, SIGNIFICANCES, ROLES = model.EFFECTS, model.SIGNIFICANCES, model.ROLES
OBSERVATION_TYPES, EVIDENCE_CLASSES = model.OBSERVATION_TYPES, model.EVIDENCE_CLASSES
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


def _seq(v):
    """A JSON list as a tuple; anything else as it is, for the validator to
    refuse (tuple() of a string would make letters of it)."""
    return () if v is None else tuple(v) if isinstance(v, list) else v


def _response(r, where):
    if r is None:
        return None
    _known(r, RESPONSE, where)
    return Response(r["measure"], r["unit"], r.get("treatment_value"), r.get("comparator_value"),
                    r.get("qualifier"))


def _statement(o, roles):
    oid = o.get("observation_id", "?")
    _known(o, OBS, oid)
    ev = o.get("evidence") or {}
    _known(ev, EVIDENCE, f"{oid} evidence")
    out = o.get("outcome") or {}
    _known(out, OUTCOME, f"{oid} outcome")
    dets = []
    for c in o.get("components") or []:
        _known(c, COMPONENT, f"{oid} component")
        dets.append(Detection(c["atlas_id"], c["detected"], c.get("method"),
                              c.get("identification"), c.get("abundance_class"),
                              c.get("detection_limit"), _seq(c.get("measurements"))))
    tax = o["taxon"]
    return EvidenceStatement(
        observation_id=oid, observation_type=o["observation_type"],
        taxon=tax["accepted_name"], taxon_as_reported=tax.get("name_as_reported"),
        statement=o["statement"], context={} if o.get("context") is None else o["context"],
        reference={} if o.get("reference") is None else o["reference"], source_location=o.get("source_location"),
        assay=o.get("assay"), experiment_id=o.get("experiment_id"),
        materials=_seq(o.get("materials")), treatment=_seq(o.get("treatment")),
        comparator=o.get("comparator"), dose=o.get("dose"), trap=o.get("trap"),
        dispenser=o.get("dispenser"), effect=out.get("effect"), significance=out.get("significance"),
        response=_response(out.get("response"), f"{oid} response"), detections=tuple(dets),
        roles=roles.get(oid, {}), author_statement=ev.get("author_statement"),
        author_statement_form=ev.get("author_statement_form"), atlas_notes=ev.get("atlas_notes"),
        unresolved=_seq(ev.get("unresolved")))


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
            classes.append(EvidenceClass(e["evidence_class"], e["taxon_id"], _seq(e["observation_ids"]),
                                         _seq(e.get("populations")), e["n_observations"],
                                         e.get("n_significant"), e.get("n_tentative"),
                                         e.get("population_match")))
        compounds[c["atlas_id"]] = BundleCompound(c["atlas_id"], c.get("shorthand"), c.get("systematic"),
                                                  c.get("inchikey"), tuple(classes))
        for occ in c.get("occurrences") or []:
            _known(occ, OCCURRENCE, f"compound {c['atlas_id']} occurrence")
            if occ["observation_id"] not in ids:
                raise InvalidBundle(f"{path}: {c['atlas_id']} occurs in {occ['observation_id']}, "
                                    "which the bundle does not carry")
            roles.setdefault(occ["observation_id"], {})[c["atlas_id"]] = _seq(occ["roles"])
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
    b = TeachingBundle(j["case_id"], j["case_title"], j["atlas_version"],
                       hashlib.sha256(Path(path).read_bytes()).hexdigest(), obs, compounds,
                       tuple(contrasts), tuple(j.get("espmap") or ()))
    problems = model.bundle_problems(b)
    if problems:
        raise InvalidBundle(f"{path}: {problems[0]}" + (f" (and {len(problems) - 1} more)" if len(problems) > 1 else ""))
    return b
