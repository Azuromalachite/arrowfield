"""Arrowfield's own data model: the only objects views and lessons consume.

Two narrow adapters fill it, and nothing else knows an upstream format:
  record.py   espmap records (JSON + NPZ)   -> MoleculeView
  bundle.py   the Atlas teaching bundle     -> TeachingBundle (EvidenceStatement,
              (spa.teaching/0)                 ObservedContrast)
A view never reads an espmap or Atlas dict. If either upstream changes its
storage, only its adapter changes.

The model also owns the judgements a view must not improvise: whether two
molecules may be compared (`compare`), and the evidence state of an
observation (`EvidenceState`). They are explicit values here, so a view can
show "not comparable" or "not tested" instead of failing or guessing.
"""
import hashlib
import json
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Optional

import numpy as np

RULES = ("nearest_nucleus", "radius_scaled")


# ---------------------------------------------------------------- provenance
@dataclass(frozen=True)
class Provenance:
    """Where a measurement came from, enough to trace a picture back."""
    source: str                     # a record's file stem, or a fixture's name
    analysis_key: str               # numbers compare only within one
    render_key: str                 # pictures compare only within one
    record_schema: Optional[str] = None
    espmap_version: Optional[str] = None
    code: Optional[str] = None      # the espmap commit or tag that made it
    sha256: Optional[str] = None    # of the record file
    xc: Optional[str] = None        # the settings as the run's command line spelled them
    basis: Optional[str] = None
    solvent: Optional[str] = None   # None or "vacuum": gas phase
    # key_fields.analysis, as sorted JSON: the canonical protocol analysis_key digests
    analysis_fields: Optional[str] = None

    @property
    def method(self):
        gas = self.solvent in (None, "vacuum", "gas")
        return f"{self.xc}/{self.basis}, {'gas' if gas else 'PCM ' + self.solvent}"


# a canonical functional or basis name in the shapes espmap's xc_key and basis_key give (b3lyp,
# def2svp, 631g**, 6311++g(2d,2p), libxc:hyb(0,0,0)+263*1+267*1): ASCII, no whitespace or "/", a
# comma only inside parentheses (#21's round 13), and no phase or solvent word anywhere
_NAME_CHARS = r"[A-Za-z0-9+*.:_-]"
CANONICAL_NAME = re.compile(rf"(?:{_NAME_CHARS}|\((?:{_NAME_CHARS}|,)*\))+")
PHASE_WORD = re.compile(r"gas|vacuum|pcm|solv|phase|water|dmso|methanol|acetonitrile|chloroform|toluene|hexane", re.I)


def analysis_digest(fields):
    """espmap's record_key: the sha1 of the fields' sorted JSON (v0.2.1 and v0.3.0 alike)."""
    return hashlib.sha1(json.dumps(fields, sort_keys=True).encode()).hexdigest()


def protocol(prov):
    """The canonical protocol a record's analysis_key digests (key_fields.analysis), or None.
    espmap decides what method names mean; Arrowfield reads its canonical fields and never
    normalises a functional or basis name itself."""
    try:
        f = json.loads(prov.analysis_fields) if isinstance(prov.analysis_fields, str) else None
    except ValueError:
        return None
    return f if isinstance(f, dict) else None


def protocol_problems(prov):
    """Why a record's protocol cannot be named, or []: no canonical fields, fields its key does
    not digest, or raw settings that contradict them where that is checkable without PySCF (the
    phase; and a solvent stated with no method)."""
    f = protocol(prov)
    if f is None:
        return [f"{prov.source}: no canonical analysis fields (key_fields.analysis) to name its protocol"]
    if analysis_digest(f) != prov.analysis_key:
        return [f"{prov.source}: analysis_key {prov.analysis_key} is not the digest of its canonical fields"]
    out = []
    # a caption prints "xc/basis, phase": each name in espmap's key shapes, never "None", never
    # naming a phase
    if not all(isinstance(f.get(k), str) and CANONICAL_NAME.fullmatch(f[k]) and f[k].casefold() != "none"
               and not PHASE_WORD.search(f[k]) for k in ("xc", "basis")):
        out.append(f"{prov.source}: its canonical fields name no functional and basis a caption can show "
                   f"({f.get('xc')!r}/{f.get('basis')!r})")
    if "solvent" not in f:
        out.append(f"{prov.source}: its canonical fields name no solvent (None is the gas phase)")
    if f.get("solvent") is not None and not isinstance(f["solvent"], str):
        out.append(f"{prov.source}: its canonical solvent {f['solvent']!r} is not a solvent name")
    raw = None if prov.solvent in (None, "vacuum", "gas") else prov.solvent
    if raw != f.get("solvent"):
        out.append(f"{prov.source}: its settings say solvent {prov.solvent!r}, its protocol "
                   f"{f.get('solvent') or 'gas phase'}")
    if prov.solvent is not None and prov.xc is None and prov.basis is None:
        out.append(f"{prov.source}: its settings state a solvent with no method")
    return out


def protocol_words(prov):
    """(method, solvent) from the canonical protocol: ("b3lyp/def2svp", None) in the gas phase,
    else the PCM solvent's name."""
    f = protocol(prov)
    return f"{f['xc']}/{f['basis']}", f.get("solvent")


# -------------------------------------------------------------- measurements
@dataclass(frozen=True)
class Scale:
    kind: str                       # "diverging" or "sequential"
    vmin: float
    vmax: float
    cmap: str


@dataclass(frozen=True)
class SiteMeasurement:
    """One site under one partition rule, kcal/mol and A^2."""
    label: str
    atoms: tuple                    # record atom indices
    rule: str
    mean: Optional[float]           # None: the site is buried (no surface)
    area: float
    area_pos: float
    area_neg: float
    partition: Optional[float]      # |difference between the two rules' means|
    p05: Optional[float] = None
    p95: Optional[float] = None
    vmin: Optional[float] = None
    vmax: Optional[float] = None
    h_generated_polar: bool = False
    overlaps: tuple = ()            # labels of sites sharing atoms with this one


@dataclass(frozen=True)
class Atom:
    """One atom, with its single-atom statistics. Diagnostic only: the
    partition rule moves one atom by up to 11.9 kcal/mol (espmap #42)."""
    index: int
    element: str
    atom_map: Optional[int]
    parent: Optional[int] = None    # a hydrogen's heavy atom
    by_rule: dict = field(default_factory=dict)   # rule -> (mean, area) or None


@dataclass(frozen=True, eq=False)
class Surface:
    verts: np.ndarray               # A
    faces: np.ndarray
    V: np.ndarray                   # kcal/mol
    owner: np.ndarray               # the atom owning each triangle, under one rule


@dataclass(frozen=True, eq=False)
class MoleculeView:
    """One computed geometry of one molecule: what every view draws from."""
    name: str
    symbols: tuple
    xyz: np.ndarray                 # A, record atom order
    bonds: Optional[tuple]          # ((i, j, order), ...) or None if unconfirmed
    charge: int
    atoms: tuple                    # Atom, by index
    sites: dict                     # label -> {rule: SiteMeasurement}
    sites_tile_surface: bool
    genus: int                      # surface topology: a diagnostic, never "folded"
    scale: Scale
    provenance: Provenance
    mapped_smiles: Optional[str] = None
    surface_loader: Optional[Callable] = None   # rule -> Surface
    # The producer's folded call, used as given and never re-derived from
    # genus: a macrocycle's surface has a hole and folds nowhere, and a chain
    # folded by contact can keep genus 0. Three definitions were in use when
    # this was inferred here (frozen study: contact OR genus > 0; espmap v0.3:
    # contact; Arrowfield: genus > 0). None: not stated, so no extended view.
    folded: Optional[bool] = None

    @property
    def atom_maps(self):
        return {a.index: a.atom_map for a in self.atoms if a.atom_map}

    @property
    def analysis_key(self):
        return self.provenance.analysis_key

    @property
    def render_key(self):
        return self.provenance.render_key

    def site_labels(self):
        return list(self.sites)

    def site(self, label, rule="nearest_nucleus"):
        if label not in self.sites:
            raise KeyError(f"{self.name}: no site {label!r}; sites are {self.site_labels()}")
        return self.sites[label][rule]

    def site_by_map(self, atom_map, rule="nearest_nucleus"):
        """The site whose heavy atom carries this atom map: 'the site of C3',
        whatever its label."""
        idx = next((i for i, m in self.atom_maps.items() if m == atom_map), None)
        if idx is None:
            raise KeyError(f"{self.name}: no atom with map {atom_map}")
        hits = [lab for lab, by in self.sites.items() if idx in by[rule].atoms]
        if not hits:
            raise KeyError(f"{self.name}: map {atom_map} is in no site")
        return self.site(hits[0], rule)

    def site_atoms(self, label):
        members = set(self.site(label).atoms)
        return [a for a in self.atoms if a.index in members]

    def surface(self, rule="nearest_nucleus"):
        if self.surface_loader is None:
            raise ValueError(f"{self.name}: no surface (no NPZ beside the record)")
        return self.surface_loader(rule)


def comparable(a, b):
    """Numbers may be compared only between measurements with one analysis key."""
    return a.analysis_key == b.analysis_key


def same_scale(a, b):
    """Pictures may be compared only between measurements with one render key."""
    return a.render_key == b.render_key


# ---------------------------------------------------------------- comparison
class ComparisonStatus(Enum):
    COMPARABLE = "comparable"
    NOT_COMPARABLE = "not comparable"
    SINGLE_CONFORMER = "single conformer only"


class Verdict(Enum):
    SEPARATED = "separated: every conformer differs in the same direction"
    SEPARATED_WITHIN_PARTITION = ("separated: every conformer differs in the same "
                                  "direction; but within the partition sensitivity")
    OVERLAPPING = "not attributable to the change: the conformer ranges overlap"
    # the user's wording (2026-10-24): what cannot be judged, and why
    CANNOT_JUDGE = "cannot judge conformational robustness: one geometry on a side"
    NOT_COMPARED = "not compared: the measurements have different analysis keys"
    RULE_DEPENDENT = ("direction depends on the partition rule: every conformer differs one way "
                      "under the nearest-nucleus rule, but the radius-scaled rule gives the opposite sign")


class Refused(LookupError):
    """A comparison that cannot be made, with the reason a student can act on
    (a site buried in a named conformer). Raised, never returned as numbers."""


@dataclass(frozen=True)
class Comparison:
    """One site, two molecules, each a list of conformers. The verdict is
    deliberately simple: separated only when the conformer ranges do not
    overlap, and never quoted without the partition sensitivity."""
    site: str
    name_a: str
    name_b: str
    status: ComparisonStatus
    verdict: Verdict
    n: tuple = (0, 0)
    diff: Optional[float] = None
    range_a: Optional[float] = None
    range_b: Optional[float] = None
    partition: Optional[float] = None
    extended: Optional[dict] = None  # the same ranges, and a verdict, over unfolded conformers only
    diff_other_rule: Optional[float] = None
    method: Optional[str] = None     # e.g. "b3lyp/def2svp, gas", shared by every measurement


# A difference under the other partition rule smaller than this (kcal/mol)
# has no sign: the card prints two decimals, so a flip below half its last
# digit could not be shown, and floating-point rounding (2e-16) must never
# read as a direction. -0.01 is a flip.
SIGN_TOL = 0.005


def _verdict(ma, mb, part, diff, diff_other):
    """The verdict for two lists of site means (at least two each)."""
    if not (min(mb) > max(ma) or max(mb) < min(ma)):
        return Verdict.OVERLAPPING
    # rounded first: 2.0 - 1.995 is 0.004999999999999893 in floats, and a
    # rule that says 'SIGN_TOL or more' must count exactly SIGN_TOL
    if diff_other is not None and diff * diff_other < 0 and round(abs(diff_other), 9) >= SIGN_TOL:
        return Verdict.RULE_DEPENDENT
    return Verdict.SEPARATED_WITHIN_PARTITION if abs(diff) <= part else Verdict.SEPARATED


def compare(views_a, views_b, label, name_a="A", name_b="B", rule="nearest_nucleus"):
    """Compare a site between two molecules' conformers. Measurements with
    different analysis keys are not compared at all: the result says so. A
    site buried in any conformer, under either rule, is refused (Refused),
    naming the conformer: dropping it would quietly change n; so is a
    non-finite mean, and an empty side. A difference whose sign flips under
    the other partition rule (by SIGN_TOL or more) is RULE_DEPENDENT. The
    extended view leaves out exactly the conformers stated folded."""
    import math
    empty = [n for n, vs in ((name_a, views_a), (name_b, views_b)) if not vs]
    if empty:
        raise Refused(f"no conformers for {' and '.join(map(repr, empty))}: there is nothing to "
                      "compare. Load at least one conformer per molecule.")
    if any(not comparable(views_a[0], v) for v in views_a[1:] + views_b):
        return Comparison(label, name_a, name_b, ComparisonStatus.NOT_COMPARABLE,
                          Verdict.NOT_COMPARED)
    other = next(r for r in RULES if r != rule)
    a = [v.site(label, rule) for v in views_a]
    b = [v.site(label, rule) for v in views_b]
    ao = [v.site(label, other) for v in views_a]
    bo = [v.site(label, other) for v in views_b]
    buried = [f"{v.name} ({r})" for vs, xs, r in ((views_a, a, rule), (views_b, b, rule),
                                                  (views_a, ao, other), (views_b, bo, other))
              for v, x in zip(vs, xs) if x.mean is None]
    if buried:
        raise Refused(f"site {label!r} has no surface in {', '.join(buried)}: it is buried in that "
                      "conformer, so there is no site mean to compare. Compare a site that is exposed "
                      "in every conformer, or read the folded conformer on its own.")
    bad = [f"{v.name} ({r}: {x.mean})" for vs, xs, r in ((views_a, a, rule), (views_b, b, rule),
                                                       (views_a, ao, other), (views_b, bo, other))
           for v, x in zip(vs, xs) if not math.isfinite(x.mean)]
    if bad:
        raise Refused(f"site {label!r} has a mean that is not a number in {', '.join(bad)}: that "
                      "conformer was not measured, and min/max would silently skip it. Recompute "
                      "or remove that record.")
    corrupt = [f"{v.name} ({x.rule}: {x.partition})" for vs, xs in ((views_a, a), (views_b, b))
               for v, x in zip(vs, xs)
               if x.partition is not None and not (math.isfinite(x.partition) and x.partition >= 0)]
    if corrupt:
        raise Refused(f"site {label!r} has a partition sensitivity that is not a finite, "
                      f"non-negative number in {', '.join(corrupt)}: it is the absolute difference "
                      "between the two rules' means, so that record is corrupt. Recompute it.")
    ma, mb = [x.mean for x in a], [x.mean for x in b]
    # The partition is |difference between the rules' means|, and the means
    # are the measurement: a stated partition never counts for less than they
    # show (probe rounds 3-4: unstated, or stated as 0 against means 8 apart,
    # it hid the rule disagreement).
    part = max(max(x.partition or 0.0, abs(x.mean - xo.mean)) for x, xo in zip(a + b, ao + bo))
    diff = float(np.mean(mb) - np.mean(ma))
    diff_other = float(np.mean([x.mean for x in bo]) - np.mean([x.mean for x in ao]))
    extended = None
    unstated = sum(v.folded is None for v in views_a + views_b)
    ext_a = [m for m, v in zip(ma, views_a) if v.folded is False]
    ext_b = [m for m, v in zip(mb, views_b) if v.folded is False]
    if unstated:
        extended = dict(n=None, range_a=None, range_b=None, verdict=None, unstated=unstated)
    elif len(ext_a) < len(ma) or len(ext_b) < len(mb):
        if min(len(ext_a), len(ext_b)) < 2:
            ext_verdict = None                       # fewer than two extended conformers
        else:
            eo_a = [x.mean for x, v in zip(ao, views_a) if v.folded is False]
            eo_b = [x.mean for x, v in zip(bo, views_b) if v.folded is False]
            ext_verdict = _verdict(ext_a, ext_b, part, float(np.mean(ext_b) - np.mean(ext_a)),
                                   float(np.mean(eo_b) - np.mean(eo_a)))
        extended = dict(n=(len(ext_a), len(ext_b)),
                        range_a=float(max(ext_a) - min(ext_a)) if len(ext_a) > 1 else None,
                        range_b=float(max(ext_b) - min(ext_b)) if len(ext_b) > 1 else None,
                        verdict=ext_verdict)
    if min(len(a), len(b)) < 2:
        status, verdict = ComparisonStatus.SINGLE_CONFORMER, Verdict.CANNOT_JUDGE
    else:
        status = ComparisonStatus.COMPARABLE
        verdict = _verdict(ma, mb, part, diff, diff_other)
    return Comparison(label, name_a, name_b, status, verdict, n=(len(a), len(b)), diff=diff,
                      range_a=float(max(ma) - min(ma)) if len(ma) > 1 else None,
                      range_b=float(max(mb) - min(mb)) if len(mb) > 1 else None,
                      partition=float(part), extended=extended, diff_other_rule=diff_other,
                      method=_method_words(views_a[0].provenance))


def _method_words(prov):
    """The method a caption names: the canonical protocol's when the record carries it and its
    key digests it, else the raw settings (a fixture without key_fields), else None."""
    if protocol(prov) is not None and not protocol_problems(prov):
        how, solvent = protocol_words(prov)
        return f"{how}, {'gas' if solvent is None else 'PCM ' + solvent}"
    return prov.method if prov.xc else None


# ------------------------------------------------------------------ evidence
class EvidenceState(Enum):
    """How an observation reads on a card. not_tested is never NEGATIVE: an
    absence of evidence is shown as such, never as evidence of absence."""
    OBSERVED = "observed"
    NEGATIVE = "negative"               # tested, and the effect was absent
    NOT_TESTED = "not tested"
    CONTESTED = "contested"             # studies disagree
    TEACHING_CONTROL = "teaching control"
    COMPUTATIONAL_MODEL = "computational model"


@dataclass(frozen=True)
class Detection:
    """One compound looked for in a chemical measurement, as recorded."""
    compound: str                       # Atlas id
    detected: bool
    method: Optional[str] = None
    identification: Optional[str] = None
    abundance_class: Optional[str] = None
    detection_limit: Optional[dict] = None     # {"value", "unit"} or None: not stated
    measurements: tuple = ()            # dicts, as recorded (value, unit, qualifier, spread)


@dataclass(frozen=True)
class Response:
    """A behavioural readout, as recorded: never rescaled or rounded."""
    measure: str
    unit: str
    treatment_value: Optional[float]
    comparator_value: Optional[float]
    qualifier: Optional[str] = None


@dataclass(frozen=True, eq=False)
class EvidenceStatement:
    """One Atlas observation, its raw fields first. Filled by bundle.py only.

    A gland-extract negative, a lure with no behavioural effect, an absent
    antennal response and an effect whose significance was never tested are
    different scientific statements; the Atlas keeps them apart and so does
    this. `statement` is the Atlas's own sentence for the record, shown
    verbatim. A card's broad category (`display_state`) is derived from these
    fields, per compound, and never stored in place of them."""
    observation_id: str
    observation_type: str               # chemical_measurement | behavioural_response | ...
    taxon: str                          # accepted name
    taxon_as_reported: Optional[str]
    statement: str                      # the Atlas's sentence, verbatim
    context: dict                       # population, locality, sex, ... as recorded
    reference: dict                     # citation, doi, atlas_id
    source_location: Optional[str]
    assay: Optional[str] = None
    experiment_id: Optional[str] = None
    materials: tuple = ()
    treatment: tuple = ()               # blend dicts, as recorded
    comparator: Optional[dict] = None
    dose: Optional[dict] = None
    trap: Optional[dict] = None
    dispenser: Optional[dict] = None
    effect: Optional[str] = None        # increased | decreased | no_effect | caught_without_control
    significance: Optional[str] = None  # significant | not_significant | not_tested | not_reported
    response: Optional[Response] = None
    detections: tuple = ()              # Detection
    roles: dict = field(default_factory=dict)   # compound -> roles in this observation
    author_statement: Optional[str] = None      # a quotation: the source's terms, not the Atlas's
    author_statement_form: Optional[str] = None
    atlas_notes: Optional[str] = None
    unresolved: tuple = ()

    def display_state(self, compound):
        """(EvidenceState, label) for one compound in this observation: the
        card's broad category plus the words that keep it specific. Absence
        of a test is never NEGATIVE; NEGATIVE always says what was absent."""
        det = next((d for d in self.detections if d.compound == compound), None)
        if self.observation_type == "chemical_measurement":
            if det is None:
                return EvidenceState.NOT_TESTED, "not looked for in this measurement"
            if det.detected:
                return EvidenceState.OBSERVED, "detected" + (f" ({det.abundance_class})"
                                                             if det.abundance_class else "")
            lim = det.detection_limit
            return EvidenceState.NEGATIVE, ("not detected, " + (
                f"detection limit {lim['value']} {lim['unit']}" if lim else "detection limit not stated"))
        roles = self.roles.get(compound, ())
        if not roles:
            return EvidenceState.NOT_TESTED, "not part of this experiment"
        if self.effect is None:
            return EvidenceState.NOT_TESTED, "no effect recorded"
        sig = {"significant": "significant", "not_significant": "not significant",
               "not_tested": "significance not tested", "not_reported": "significance not reported"}
        words = f"{self.effect.replace('_', ' ')}, {sig.get(self.significance, self.significance)}"
        if self.effect == "no_effect" and self.significance == "not_significant":
            return EvidenceState.NEGATIVE, f"{self.assay or 'assay'}: {words}"
        return EvidenceState.OBSERVED, f"{self.assay or 'assay'}: {words}"


@dataclass(frozen=True)
class ObservedContrast:
    """Two observations a teaching case sets side by side, as the Atlas
    states each: no verdict about the pair is implied."""
    label: str
    observation_ids: tuple
    a: str
    b: str


@dataclass(frozen=True)
class EvidenceClass:
    """One of the Atlas's own evidence classes for a compound in a taxon, as
    the Atlas derives it from its records: Arrowfield reads it, never asserts
    it. Atlas-wide, so its observation ids may lie outside the bundle."""
    evidence_class: str                 # e.g. addition_decreased, female_not_detected
    taxon_id: str
    observation_ids: tuple
    populations: tuple
    n_observations: int
    n_significant: Optional[int] = None
    n_tentative: Optional[int] = None
    population_match: Optional[str] = None


@dataclass(frozen=True)
class BundleCompound:
    atlas_id: str
    shorthand: Optional[str]
    systematic: Optional[str]
    inchikey: Optional[str]
    evidence_classes: tuple = ()        # EvidenceClass


@dataclass(frozen=True, eq=False)
class TeachingBundle:
    """An Atlas teaching bundle (spa.teaching/0), read by bundle.py."""
    case_id: str
    case_title: str
    atlas_version: str
    sha256: str                         # of the bundle file
    observations: dict                  # observation id -> EvidenceStatement
    compounds: dict                     # Atlas id -> BundleCompound
    contrasts: tuple                    # ObservedContrast
    espmap: tuple = ()                  # computed ensembles the bundle carries, as given


# --------------------------------------------------------------- pair view
# The pair view's three judgements, each its own function so a failure is
# attributable to one layer (tests/pedagogy/pair_view_set.json, #14 + #17).
class Relation(Enum):
    SAME_GEOMETRY = "same_geometry"
    CONFORMERS_OF_ONE_COMPOUND = "conformers_of_one_compound"
    TWO_COMPOUNDS = "two_compounds"


# All-atom RMSD after optimal superposition, minimised over symmetry-equivalent
# atoms, at or below which two geometries of one compound are one geometry.
# Measured before freezing (#14): the pilot's equivalent seed pairs are at
# most 3.4e-6 A apart; the Atlas ensembles' distinct conformers at least
# 0.98 A (heavy atoms) and 1.40 A (all atoms, closest five).
SAME_GEOMETRY_RMSD = 0.01
_MIN_ATOM_SEPARATION = 0.5      # A: closer than any bond; a geometry with such a pair is no geometry
_HEAVY_MATCHES = 1000           # heavy-atom symmetry matches tried (a C18 chain has a handful)
_ASSIGN_ROUNDS = 20             # superpose / reassign hydrogens until the assignment is stable


def _kabsch(P, Q):
    """(rmsd, rotation, centre of P, centre of Q) superposing P onto Q."""
    p0, q0 = P.mean(0), Q.mean(0)
    u, _, vt = np.linalg.svd((P - p0).T @ (Q - q0))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    rmsd = float(np.sqrt((((P - p0) @ r.T - (Q - q0)) ** 2).sum(1).mean()))
    return rmsd, r, p0, q0


def _best_rmsd(ma, xa, mb, xb):
    """All-atom RMSD after optimal superposition, minimised over symmetry-
    equivalent atoms, independent of how either input orders its atoms or
    bonds (the probe on #18): heavy atoms through the graph's own symmetry
    matches; each heavy atom's hydrogens by an optimal assignment among
    themselves (hydrogens on one parent are interchangeable: relabelling them
    changes neither the geometry nor its surface), alternated with a Kabsch
    superposition on all atoms until the assignment is stable."""
    from itertools import permutations
    from rdkit import Chem
    heavy_a = [a.GetIdx() for a in ma.GetAtoms() if a.GetAtomicNum() > 1]
    heavy_b = [a.GetIdx() for a in mb.GetAtoms() if a.GetAtomicNum() > 1]
    hyd = lambda m, i: sorted(n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1)  # noqa: E731
    matches = Chem.RemoveHs(mb).GetSubstructMatches(Chem.RemoveHs(ma), uniquify=False, useChirality=False,
                                                    maxMatches=_HEAVY_MATCHES)
    best = np.inf
    for match in matches:
        pairs = [(heavy_a[k], heavy_b[j]) for k, j in enumerate(match)]
        groups = [(hyd(ma, i), hyd(mb, j)) for i, j in pairs]
        if any(len(ga) != len(gb) for ga, gb in groups):
            continue
        mapping = pairs + [(h, g) for ga, gb in groups for h, g in zip(ga, gb)]
        for _ in range(_ASSIGN_ROUNDS):
            ia, ib = zip(*mapping)
            rmsd, r, p0, q0 = _kabsch(xa[list(ia)], xb[list(ib)])
            moved = (xa - p0) @ r.T + q0
            new = list(pairs)
            for ga, gb in groups:
                new += min((list(zip(ga, perm)) for perm in permutations(gb)),
                           key=lambda m: sum(float(((moved[h] - xb[g]) ** 2).sum()) for h, g in m))
            if sorted(new) == sorted(mapping):
                break
            mapping = new
        best = min(best, rmsd)
    return best


def _identity(v):
    """(canonical isomeric SMILES with charges, formal charge, RDKit mol with
    Hs) for one MoleculeView, or Refused when its constitution or 3D stereo
    cannot be established."""
    from rdkit import Chem
    import views
    xyz = np.asarray(v.xyz, dtype=float) if isinstance(v.xyz, np.ndarray) else None
    if xyz is None or xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise Refused(f"{v.name}: its coordinates are not a finite (n, 3) array")
    if len(xyz) > 1:
        d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
        if d[np.triu_indices(len(xyz), 1)].min() < _MIN_ATOM_SEPARATION:
            raise Refused(f"{v.name}: two atoms closer than {_MIN_ATOM_SEPARATION} A: no geometry to assign "
                          "stereo from")
    n = len(xyz)
    if v.bonds is not None and not all(
            isinstance(bd, (tuple, list)) and len(bd) == 3 and all(isinstance(k, int) and not isinstance(k, bool)
                                                                for k in bd[:2]) and 0 <= bd[0] < n and 0 <= bd[1] < n
            for bd in v.bonds):
        raise Refused(f"{v.name}: a bond names an atom the record does not have; its graph is not established")
    try:
        m = Chem.Mol(views.mol_from_view(v, v.mapped_smiles))
    except Exception as e:  # noqa: BLE001 - any failure to build the graph is 'not established'
        raise Refused(f"{v.name}: its constitution cannot be established ({e})") from None
    Chem.AssignStereochemistryFrom3D(m)
    unset = [s for s in Chem.FindPotentialStereo(m) if s.specified != Chem.StereoSpecified.Specified]
    if unset:
        raise Refused(f"{v.name}: {len(unset)} stereo element(s) cannot be assigned from its 3D coordinates")
    heavy = Chem.RemoveHs(m)
    for a in heavy.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(heavy, isomericSmiles=True), v.charge, m


def relation(a, b):
    """What a pair of MoleculeViews is, from their structures, never their
    names: one compound (same constitution, formal charge and 3D-assigned
    stereo) at one geometry (RMSD <= SAME_GEOMETRY_RMSD) or at two, or two
    compounds. Refused when either identity cannot be established."""
    for v in (a, b):
        if not isinstance(v, MoleculeView):
            raise Refused(f"relation compares MoleculeViews, not an object of type {type(v).__name__}")
    sa, ca, ma = _identity(a)
    sb, cb, mb = _identity(b)
    if (sa, ca) != (sb, cb):
        return Relation.TWO_COMPOUNDS
    rmsd = _best_rmsd(ma, np.asarray(a.xyz, float), mb, np.asarray(b.xyz, float))
    return Relation.SAME_GEOMETRY if rmsd <= SAME_GEOMETRY_RMSD else Relation.CONFORMERS_OF_ONE_COMPOUND


# espmap's DEFAULT_CMAP, frozen with its tables in vendor/ (palette.py)
PAIR_CMAP = {"diverging": "vik", "sequential": "lajolla"}


def shared_scale(a, b):
    """The one colour scale for a pair, from the two records' own scales:
    diverging (symmetric) if either is diverging or the pair crosses zero,
    else sequential over both. Refused for a scale that is not finite or not
    of a known kind, or a pair with no range to colour on."""
    import math
    for v in (a, b):
        s = v.scale
        if not (isinstance(s, Scale) and s.kind in PAIR_CMAP and all(
                isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)
                for x in (s.vmin, s.vmax)) and s.vmin <= s.vmax):
            raise Refused(f"{v.name}: its scale {getattr(v, 'scale', None)!r} is not a finite scale of a known kind")
    lo, hi = min(a.scale.vmin, b.scale.vmin), max(a.scale.vmax, b.scale.vmax)
    if "diverging" in (a.scale.kind, b.scale.kind) or lo < 0 < hi:
        m = max(abs(lo), abs(hi))
        kind, vmin, vmax = "diverging", -m, m
    else:
        kind, vmin, vmax = "sequential", lo, hi
    if not vmin < vmax:
        raise Refused(f"the pair's scales give no range to colour on ({vmin} to {vmax})")
    return Scale(kind, float(vmin), float(vmax), PAIR_CMAP[kind])


def fold_state(v):
    """The producer's folded flag as the page states it: "true", "false" or
    "unstated", never inferred from genus. Refused for a flag that is not a
    bool."""
    if v.folded is None:
        return "unstated"
    if not isinstance(v.folded, bool):
        raise Refused(f"{v.name}: folded {v.folded!r} is not true, false or unstated")
    return "true" if v.folded else "false"


# ------------------------------------------------------- bundle validation
# The Atlas's controlled words (its observation.schema.json and
# teaching-bundle.schema.json). One validator, bundle_problems, serves
# bundle.load and a bundle built by hand alike, so a view never meets a
# value it would have to trust (probe rounds 3-4 on #13; the user's
# direction: one shared semantic validator for both paths).
OBSERVATION_TYPES = ("chemical_measurement", "behavioural_response", "electrophysiology")
EFFECTS = ("increased", "decreased", "no_effect", "caught_without_control")
SIGNIFICANCES = ("significant", "not_significant", "not_tested", "not_reported")
ROLES = ("treatment", "comparator", "detected", "not_detected", "not_tested", "antennal_active",
         "antennal_inconsistent", "antennal_inactive", "antennal_not_tested")
EVIDENCE_CLASSES = ("female_detected_gland", "female_detected_emission", "female_detected_combined",
                    "female_not_detected", "antennal_active", "antennal_inconsistent", "antennal_inactive",
                    "attractant_alone", "attractive_lure_component", "addition_increased",
                    "addition_decreased", "addition_no_effect", "female_pheromone_component")
POPULATION_MATCH = ("yes", "no", "unknown")


def _behaviour(cls):
    """Classes that carry n_significant, as the Atlas builds them
    (scripts/build_evidence.py, rows())."""
    return cls.startswith(("attract", "addition")) or cls == "female_pheromone_component"


def _chemistry(cls):
    """Classes that carry n_tentative (the same function)."""
    return cls.startswith("female_detected") or cls == "female_pheromone_component"


def _is_real(v):
    """A finite real number; a bool is not one (True is an int to Python)."""
    import math
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def _is_count(v):
    return isinstance(v, int) and not isinstance(v, bool) and v >= 0


def _is_text(v, optional=False):
    return (v is None and optional) or isinstance(v, str)


def _is_texts(v):
    """A tuple of strings: never a bare string, which iterates as letters."""
    return isinstance(v, tuple) and all(isinstance(x, str) for x in v)


def _detection_problems(d, where):
    if not isinstance(d, Detection):
        return [f"{where}: a detection is a {type(d).__name__}, not a Detection"]
    out = []
    if not _is_text(d.compound):
        out.append(f"{where}: detection compound {d.compound!r} is not an id")
    if not isinstance(d.detected, bool):
        out.append(f"{where} {d.compound}: detected {d.detected!r} is not true or false")
    for name in ("method", "identification", "abundance_class"):
        if not _is_text(getattr(d, name), optional=True):
            out.append(f"{where} {d.compound}: {name} {getattr(d, name)!r} is not text")
    lim = d.detection_limit
    if lim is not None and not (isinstance(lim, dict) and _is_real(lim.get("value")) and lim["value"] >= 0
                                and _is_text(lim.get("unit"))):
        out.append(f"{where} {d.compound}: detection limit {lim!r} is not a non-negative number with a unit")
    if not isinstance(d.measurements, tuple):
        out.append(f"{where} {d.compound}: measurements are not a list")
    return out


def statement_problems(st):
    """Every way one EvidenceStatement breaks the Atlas's vocabulary, types
    or container shapes; empty when it may be shown."""
    if not isinstance(st, EvidenceStatement):
        return [f"an observation is a {type(st).__name__}, not an EvidenceStatement"]
    w = repr(st.observation_id)
    out = []
    if not _is_text(st.observation_id):
        out.append(f"observation id {w} is not text")
    if st.observation_type not in OBSERVATION_TYPES:
        out.append(f"{w}: observation_type {st.observation_type!r} is not one of {OBSERVATION_TYPES}")
    if st.effect is not None and st.effect not in EFFECTS:
        out.append(f"{w}: effect {st.effect!r} is not one of {EFFECTS}")
    if st.significance is not None and st.significance not in SIGNIFICANCES:
        out.append(f"{w}: significance {st.significance!r} is not one of {SIGNIFICANCES}")
    for name in ("taxon", "statement"):
        if not _is_text(getattr(st, name)):
            out.append(f"{w}: {name} is not text")
    for name in ("taxon_as_reported", "source_location", "assay", "author_statement"):
        if not _is_text(getattr(st, name), optional=True):
            out.append(f"{w}: {name} is not text")
    for name in ("context", "reference"):
        if not isinstance(getattr(st, name), dict):
            out.append(f"{w}: {name} is not a mapping")
    if isinstance(st.context, dict):
        for k in ("population", "locality"):
            if not _is_text(st.context.get(k), optional=True):
                out.append(f"{w}: context {k} is not text")
    if isinstance(st.reference, dict):
        for k in ("citation", "doi"):
            if not _is_text(st.reference.get(k), optional=True):
                out.append(f"{w}: reference {k} is not text")
    if not _is_texts(st.unresolved):
        out.append(f"{w}: unresolved items are not a list of text")
    for name in ("materials", "treatment"):
        if not isinstance(getattr(st, name), tuple):
            out.append(f"{w}: {name} is not a list")
    if not isinstance(st.detections, tuple):
        out.append(f"{w}: detections are a {type(st.detections).__name__}, not a list")
    else:
        for d in st.detections:
            out += _detection_problems(d, w)
    if not isinstance(st.roles, dict):
        out.append(f"{w}: roles are not a mapping")
    else:
        for comp, rs in st.roles.items():
            if not _is_texts(rs) or not rs:
                out.append(f"{w}: roles of {comp} must be a non-empty list of role words, not {rs!r}")
            else:
                bad = [r for r in rs if r not in ROLES]
                if bad:
                    out.append(f"{w}: role(s) {bad} of {comp} are not the Atlas's role words {ROLES}")
    r = st.response
    if r is not None:
        if not isinstance(r, Response):
            out.append(f"{w}: response is a {type(r).__name__}, not a Response")
        else:
            for name in ("measure", "unit"):
                if not _is_text(getattr(r, name)):
                    out.append(f"{w}: response {name} is not text")
            if not _is_text(r.qualifier, optional=True):
                out.append(f"{w}: response qualifier is not text")
            for name in ("treatment_value", "comparator_value"):
                v = getattr(r, name)
                if v is not None and not _is_real(v):
                    out.append(f"{w}: response {name} {v!r} is not a finite number")
    return out


def evidence_class_problems(e, where=""):
    """Every way one EvidenceClass breaks the Atlas's vocabulary, types or
    arithmetic: n_observations counts its distinct observation ids, and the
    significant and tentative counts are at most that."""
    if not isinstance(e, EvidenceClass):
        return [f"{where}: an evidence class is a {type(e).__name__}, not an EvidenceClass"]
    w = f"{where} class {e.evidence_class!r}"
    out = []
    if e.evidence_class not in EVIDENCE_CLASSES:
        out.append(f"{w}: not one of the Atlas's classes {EVIDENCE_CLASSES}")
    if not _is_text(e.taxon_id):
        out.append(f"{w}: taxon id is not text")
    if not _is_texts(e.observation_ids):
        out.append(f"{w}: observation ids are not a list of ids")
    elif len(set(e.observation_ids)) != len(e.observation_ids):
        out.append(f"{w}: an observation id repeats")
    if not _is_texts(e.populations):
        out.append(f"{w}: populations are not a list of text")
    fpc = e.evidence_class == "female_pheromone_component"
    if fpc and e.population_match not in POPULATION_MATCH:
        out.append(f"{w}: population_match {e.population_match!r} is not one of {POPULATION_MATCH}")
    if not fpc and e.population_match is not None:
        out.append(f"{w}: population_match is set, but only female_pheromone_component carries one")
    if not _is_count(e.n_observations):
        out.append(f"{w}: n_observations {e.n_observations!r} is not a non-negative integer")
    elif _is_texts(e.observation_ids) and e.n_observations != len(set(e.observation_ids)):
        out.append(f"{w}: n_observations {e.n_observations} but {len(set(e.observation_ids))} observation ids")
    for name, applies in (("n_significant", _behaviour), ("n_tentative", _chemistry)):
        v = getattr(e, name)
        if not applies(e.evidence_class):
            if v is not None:      # the Atlas: a count that does not apply is null
                out.append(f"{w}: {name} {v!r} given, but it does not apply to this class")
            continue
        if not _is_count(v):
            out.append(f"{w}: {name} {v!r} is not a non-negative integer")
        elif _is_count(e.n_observations) and v > e.n_observations:
            out.append(f"{w}: {name} {v} exceeds n_observations {e.n_observations}")
    return out


def bundle_problems(b):
    """Every problem in a TeachingBundle, however it was built: empty when
    the views may show it. bundle.load raises on these; the views refuse."""
    if not isinstance(b, TeachingBundle):
        return [f"the bundle is a {type(b).__name__}, not a TeachingBundle"]
    out = []
    if not _is_text(b.atlas_version):
        out.append("atlas_version is not text")
    if not isinstance(b.observations, dict):
        return out + ["observations are not a mapping"]
    for oid, st in b.observations.items():
        out += statement_problems(st)
        if isinstance(st, EvidenceStatement) and st.observation_id != oid:
            out.append(f"observation {oid!r} carries the id {st.observation_id!r}")
    if not isinstance(b.compounds, dict):
        return out + ["compounds are not a mapping"]
    for cid, c in b.compounds.items():
        if not isinstance(c, BundleCompound) or c.atlas_id != cid:
            out.append(f"compound {cid!r} is not a BundleCompound with that id")
            continue
        if not _is_text(c.shorthand, optional=True):
            out.append(f"compound {cid!r}: shorthand is not text")
        if not isinstance(c.evidence_classes, tuple):
            out.append(f"compound {cid!r}: evidence classes are not a list")
            continue
        for e in c.evidence_classes:
            out += evidence_class_problems(e, f"compound {cid!r}")
    if not isinstance(b.contrasts, tuple):
        out.append("contrasts are not a list")
    else:
        for k in b.contrasts:
            if not (isinstance(k, ObservedContrast) and _is_text(k.label) and _is_texts(k.observation_ids)
                    and all(i in b.observations for i in k.observation_ids)):
                out.append(f"contrast {getattr(k, 'label', k)!r} is malformed or names observations not carried")
    return out


@dataclass(frozen=True)
class CaseMolecule:
    id: str                             # an Atlas id, or arrowfield.control.*
    name: str
    evidence_status: str                # atlas_observed | teaching_control | computational_counterfactual


@dataclass(frozen=True)
class CaseContrast:
    id: str
    question: str
    members: tuple
    motivating: tuple = ()              # observation ids


@dataclass(frozen=True)
class TeachingCase:
    """A frozen study's case: its molecules with their evidence status, its
    contrasts, and the observations that motivate them."""
    case_id: str
    atlas_pin: str
    molecules: tuple
    contrasts: tuple
    evidence: tuple = ()                # EvidenceStatement


# ------------------------------------------------------------- frozen studies
# The verdict rules Arrowfield knows how to show; a study row under any other
# rule is refused, never shown as that study's.
KNOWN_VERDICT_RULES = ("espmap-sesiidae-001/range_partition_sign_v1",)
# The analysis key each frozen study was decided under. Its results file names none, so it is pinned
# here, as its rule is: espmap-sesiidae-001's plan (section 3 at Atlas 7e34658) fixes B3LYP/def2-SVP
# + RI, gas phase, one analysis_key, and its ensembles' indexes carry this one. A page's records under
# another key are not that study's, and its verdicts are never scoped under their protocol.
STUDY_KEYS = {"espmap-sesiidae-001/range_partition_sign_v1": "d0f5ed136cda12ada8c6324a79ec9f7e79134d2f"}
# What each frozen study says it did not run, in its own words, for its study cards' limits: the
# def2-TZVP check was optional (plan section 6) and its results.md reports it "not run".
STUDY_LIMITS = {"espmap-sesiidae-001/range_partition_sign_v1": {"basis_sensitivity": "not run (plan §6)"}}
# Each frozen study's rule for "folded", in its own terms (espmap-sesiidae-001's plan section 3, and
# each Atlas ensemble index's notes): a page that shows its folded / extended labels says what they mean
STUDY_FOLD_RULES = {"espmap-sesiidae-001/range_partition_sign_v1":
                    "an intramolecular contact or a surface of genus above zero (plan §3)"}
CONDITIONS = ("1_ranges_disjoint", "2_beyond_partition", "3_same_sign")


def study_row_inconsistency(rule, outcome, cond, interval):
    """None if a row's stated outcome follows from its own stated conditions
    under its rule, else why not. This checks a record against itself, as a
    charge that does not match its atoms would be: it never derives a verdict
    from the measurements, and a consistent row is shown as the study wrote it.
    range_partition_sign_v1 (the frozen plan, section 4): resolved iff all
    three conditions hold; different on average iff 2 and 3 hold, 1 fails and
    the bootstrap interval excludes 0; unresolved otherwise."""
    if outcome in ("not evaluable", "no data"):
        return None
    if rule == "espmap-sesiidae-001/range_partition_sign_v1":
        if set(cond) != set(CONDITIONS) or not all(isinstance(cond[k], bool) for k in CONDITIONS):
            return f"conditions {sorted(cond)} are not exactly {CONDITIONS}, each true or false"
        bad = _bad_pair(interval)
        if bad:
            return f"bootstrap interval {interval} is {bad}"
        c1, c2, c3 = (cond[k] for k in CONDITIONS)
        excl0 = interval is not None and (interval[0] > 0 or interval[1] < 0)
        follows = ("resolved" if c1 and c2 and c3 else
                   "different on average" if c2 and c3 and excl0 else "unresolved")
        if outcome != follows:
            return f"the row says {outcome!r}, but its own conditions {cond} and interval {interval} give {follows!r}"
    return None


def real_number(x):
    """A finite real measurement: an int, float or numpy float, never a bool
    (Python counts True as 1; probe round 2 of #11) and never NaN or inf."""
    import math
    import numbers
    return isinstance(x, numbers.Real) and not isinstance(x, bool) and math.isfinite(float(x))


def _bad_pair(p):
    """None for a finite (low, high) pair with low <= high, or for None;
    else what is wrong with it (probe of #11: a reversed interval read as
    'excluding 0')."""
    import math
    if p is None:
        return None
    if len(p) != 2 or not all(real_number(v) for v in p):
        return "not a pair of finite numbers"
    return "reversed (low above high)" if p[0] > p[1] else None


def study_row_malformed(row):
    """None if a verdict-bearing study row (resolved, different on average,
    unresolved) carries only finite numbers and ordered ranges, else what is
    wrong: such a row is a corrupt record, never rendered as measured (probe
    of #11: a NaN delta showed as '+nan'). Not evaluable / no data rows may
    carry anything; they show no difference."""
    import math
    if row.outcome in ("not evaluable", "no data"):
        return None
    for name in ("delta", "delta_other_rule", "partition", "partition_max_member", "mean_x", "mean_y"):
        v = getattr(row, name)
        if v is not None and not real_number(v):
            return f"{name} {v!r} is not a finite number"
    import numbers                                     # numpy integer counts are counts too
    n = row.n
    if not (isinstance(n, (tuple, list)) and len(n) == 2
            and all(isinstance(k, numbers.Integral) and not isinstance(k, (bool, np.bool_)) and k >= 0
                    for k in n)):
        return f"n {n!r} is not two non-negative integers"
    if row.partition is not None and row.partition < 0:
        return f"partition {row.partition} is negative"
    for name in ("range_x", "range_y", "bootstrap_interval"):
        bad = _bad_pair(getattr(row, name))
        if bad:
            return f"{name} {getattr(row, name)} is {bad}"
    return None


@dataclass(frozen=True)
class StudySite:
    """One site of one contrast in a frozen study, as the study decided it.
    Arrowfield shows this verdict; it never re-decides it, and never
    substitutes compare()'s teaching verdict or espmap's generic one for it.
    Every number is the study's, kcal/mol, delta = mean(y) - mean(x)."""
    verdict_rule: str                   # e.g. espmap-sesiidae-001/range_partition_sign_v1
    contrast: str
    x: str                              # slug of the first molecule
    y: str
    population: str                     # "all" | "extended"
    tier: str                           # "primary" | "secondary" | "exploratory"
    site: str
    outcome: str                        # resolved | different on average | unresolved | not evaluable | no data
    n: tuple = (0, 0)
    delta: Optional[float] = None
    delta_other_rule: Optional[float] = None
    partition: Optional[float] = None   # the study's definition (its plan, section 4)
    partition_max_member: Optional[float] = None   # reported beside, never decides
    range_x: Optional[tuple] = None
    range_y: Optional[tuple] = None
    mean_x: Optional[float] = None
    mean_y: Optional[float] = None
    # resamples of the SELECTED conformers: describes that set, not a
    # population, so it is not a confidence interval (Atlas 30c4a11)
    bootstrap_interval: Optional[tuple] = None
    conditions: dict = field(default_factory=dict)


@dataclass(frozen=True, eq=False)
class StudyResults:
    """A frozen study's results file, read by study.py."""
    plan: str
    verdict_rule: str
    sha256: str                         # of the results file
    molecules: dict                     # slug -> {n_ok, n_folded, excluded, evaluable}
    sites: tuple                        # StudySite

    def site(self, contrast, site, population="all"):
        hits = [s for s in self.sites if s.contrast == contrast and s.site == site
                and s.population == population]
        if len(hits) != 1:
            raise KeyError(f"{contrast} / {site} / {population}: {len(hits)} entries in the study")
        return hits[0]


@dataclass(frozen=True)
class EnsembleSource:
    """Where one molecule's ensemble came from: an Atlas-format index."""
    atlas_id: str                       # spa.cmpd.* or arrowfield.control.*
    protocol: str
    analysis_key: str
    espmap_version: str
    index_sha256: str
    selection: str
    folded_rule: Optional[str] = None   # the producer's statement of what "folded" means
