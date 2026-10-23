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
    xc: Optional[str] = None        # the method, for captions
    basis: Optional[str] = None
    solvent: Optional[str] = None   # None or "vacuum": gas phase

    @property
    def method(self):
        gas = self.solvent in (None, "vacuum", "gas")
        return f"{self.xc}/{self.basis}, {'gas' if gas else 'PCM ' + self.solvent}"


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
                      method=views_a[0].provenance.method if views_a[0].provenance.xc else None)


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


def _bad_pair(p):
    """None for a finite (low, high) pair with low <= high, or for None;
    else what is wrong with it (probe of #11: a reversed interval read as
    'excluding 0')."""
    import math
    if p is None:
        return None
    if len(p) != 2 or not all(isinstance(v, (int, float)) and math.isfinite(v) for v in p):
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
        if v is not None and not (isinstance(v, (int, float)) and math.isfinite(v)):
            return f"{name} {v!r} is not a finite number"
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
