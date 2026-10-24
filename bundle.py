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
    if extra:     # repr: a key need not be text, and the error path must not fail
        raise InvalidBundle(f"{where}: field(s) {sorted(map(repr, extra))} are not modelled by Arrowfield; "
                            "model them in model.py and bundle.py before showing this bundle")


# The parse boundary (probe round 5, part 2): a node's container type and its
# required keys are checked before anything indexes it, so a malformed shape
# is an InvalidBundle naming where, never a crash on the way to
# model.bundle_problems. Values are not checked here: the validator does that.
def _map(x, where, need=(), allowed=None):
    if not isinstance(x, dict):
        raise InvalidBundle(f"{where}: expected mapping, got {type(x).__name__}")
    missing = [k for k in need if k not in x]
    if missing:
        raise InvalidBundle(f"{where}.{missing[0]}: required field missing")
    if allowed is not None:
        _known(x, allowed, where)
    return x


def _optional_map(x, where, need=(), allowed=None):
    return {} if x is None else _map(x, where, need, allowed)


def _list(x, where):
    if not isinstance(x, list):
        raise InvalidBundle(f"{where}: expected list, got {type(x).__name__}")
    return x


def _optional_list(x, where):
    return [] if x is None else _list(x, where)


def _key(x, where):
    """An id used as a mapping key here: it must be text before it is used."""
    if not isinstance(x, str):
        raise InvalidBundle(f"{where}: expected text id, got {type(x).__name__}")
    return x


def _seq(v):
    """A JSON list as a tuple; anything else as it is, for the validator to
    refuse (tuple() of a string would make letters of it)."""
    return () if v is None else tuple(v) if isinstance(v, list) else v


def _response(r, where):
    if r is None:
        return None
    _map(r, where, ("measure", "unit"), RESPONSE)
    return Response(r["measure"], r["unit"], r.get("treatment_value"), r.get("comparator_value"),
                    r.get("qualifier"))


def _statement(o, roles, where):
    oid = o["observation_id"]
    ev = _optional_map(o.get("evidence"), f"{where}.evidence", allowed=EVIDENCE)
    out = _optional_map(o.get("outcome"), f"{where}.outcome", allowed=OUTCOME)
    dets = []
    for n, c in enumerate(_optional_list(o.get("components"), f"{where}.components")):
        _map(c, f"{where}.components[{n}]", ("atlas_id", "detected"), COMPONENT)
        dets.append(Detection(c["atlas_id"], c["detected"], c.get("method"),
                              c.get("identification"), c.get("abundance_class"),
                              c.get("detection_limit"), _seq(c.get("measurements"))))
    tax = _map(o["taxon"], f"{where}.taxon", ("accepted_name",))
    return EvidenceStatement(
        observation_id=oid, observation_type=o["observation_type"],
        taxon=tax["accepted_name"], taxon_as_reported=tax.get("name_as_reported"),
        statement=o["statement"], context={} if o.get("context") is None else o["context"],
        reference={} if o.get("reference") is None else o["reference"], source_location=o.get("source_location"),
        assay=o.get("assay"), experiment_id=o.get("experiment_id"),
        materials=_seq(o.get("materials")), treatment=_seq(o.get("treatment")),
        comparator=o.get("comparator"), dose=o.get("dose"), trap=o.get("trap"),
        dispenser=o.get("dispenser"), effect=out.get("effect"), significance=out.get("significance"),
        response=_response(out.get("response"), f"{where}.outcome.response"), detections=tuple(dets),
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
    _map(j, "bundle", ("case_id", "case_title", "observations", "compounds", "comparisons"), TOP)
    if not j.get("atlas_version"):
        raise InvalidBundle(f"{path}: no atlas_version: a bundle must name the Atlas commit it came from")
    observations = _list(j["observations"], "observations")
    for n, o in enumerate(observations):
        _map(o, f"observations[{n}]", ("observation_id", "observation_type", "taxon", "statement"), OBS)
        _key(o["observation_id"], f"observations[{n}].observation_id")
    ids = [o["observation_id"] for o in observations]
    if len(ids) != len(set(ids)):
        raise InvalidBundle(f"{path}: an observation id repeats")
    roles = {}
    compounds = {}
    for n, c in enumerate(_list(j["compounds"], "compounds")):
        _map(c, f"compounds[{n}]", ("atlas_id",), COMPOUND)
        cid = _key(c["atlas_id"], f"compounds[{n}].atlas_id")
        classes = []
        for m, e in enumerate(_optional_list(c.get("evidence_classes"), f"compounds[{n}].evidence_classes")):
            _map(e, f"compounds[{n}].evidence_classes[{m}]",
                 ("evidence_class", "taxon_id", "observation_ids", "n_observations"), EVIDENCE_CLASS)
            classes.append(EvidenceClass(e["evidence_class"], e["taxon_id"], _seq(e["observation_ids"]),
                                         _seq(e.get("populations")), e["n_observations"],
                                         e.get("n_significant"), e.get("n_tentative"),
                                         e.get("population_match")))
        compounds[cid] = BundleCompound(cid, c.get("shorthand"), c.get("systematic"),
                                        c.get("inchikey"), tuple(classes))
        for m, occ in enumerate(_optional_list(c.get("occurrences"), f"compounds[{n}].occurrences")):
            _map(occ, f"compounds[{n}].occurrences[{m}]", ("observation_id", "roles"), OCCURRENCE)
            oid = _key(occ["observation_id"], f"compounds[{n}].occurrences[{m}].observation_id")
            if oid not in ids:
                raise InvalidBundle(f"{path}: {cid} occurs in {oid}, which the bundle does not carry")
            roles.setdefault(oid, {})[cid] = _seq(occ["roles"])
    contrasts = []
    for n, k in enumerate(_list(j["comparisons"], "comparisons")):
        _map(k, f"comparisons[{n}]", ("case_label", "observation_ids"))
        kids = _list(k["observation_ids"], f"comparisons[{n}].observation_ids")
        missing = [i for i in kids if not isinstance(i, str) or i not in ids]
        if missing:
            raise InvalidBundle(f"{path}: comparison {k['case_label']!r} names {missing}, "
                                "which the bundle does not carry")
        oc = _optional_map(k.get("observed_contrast"), f"comparisons[{n}].observed_contrast")
        contrasts.append(ObservedContrast(k["case_label"], tuple(kids), oc.get("a"), oc.get("b")))
    obs = {o["observation_id"]: _statement(o, roles, f"observations[{n}]") for n, o in enumerate(observations)}
    b = TeachingBundle(j["case_id"], j["case_title"], j["atlas_version"],
                       hashlib.sha256(Path(path).read_bytes()).hexdigest(), obs, compounds,
                       tuple(contrasts), tuple(_optional_list(j.get("espmap"), "espmap")))
    problems = model.bundle_problems(b)
    if problems:
        raise InvalidBundle(f"{path}: {problems[0]}" + (f" (and {len(problems) - 1} more)" if len(problems) > 1 else ""))
    return b
