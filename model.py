"""Arrowfield's own data model: the only objects views and lessons consume.

Two narrow adapters fill it, and nothing else knows an upstream format:
  record.py   espmap records (JSON + NPZ)   -> MoleculeView
  bundle.py   the Atlas teaching bundle     -> EvidenceStatement, TeachingCase
              (to come, once spa.teaching/0 is frozen)
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
    genus: int                      # > 0: the molecule folds onto itself
    scale: Scale
    provenance: Provenance
    mapped_smiles: Optional[str] = None
    surface_loader: Optional[Callable] = None   # rule -> Surface

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
    CANNOT_JUDGE = "single conformer: cannot judge"
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


def _verdict(ma, mb, part, diff, diff_other):
    """The verdict for two lists of site means (at least two each)."""
    if not (min(mb) > max(ma) or max(mb) < min(ma)):
        return Verdict.OVERLAPPING
    if diff_other is not None and diff * diff_other < 0:
        return Verdict.RULE_DEPENDENT
    return Verdict.SEPARATED_WITHIN_PARTITION if abs(diff) <= part else Verdict.SEPARATED


def compare(views_a, views_b, label, name_a="A", name_b="B", rule="nearest_nucleus"):
    """Compare a site between two molecules' conformers. Measurements with
    different analysis keys are not compared at all: the result says so. A
    site buried in any conformer, under either rule, is refused (Refused),
    naming the conformer: dropping it would quietly change n. A difference
    whose sign flips under the other partition rule is RULE_DEPENDENT."""
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
    ma, mb = [x.mean for x in a], [x.mean for x in b]
    part = max(x.partition or 0.0 for x in a + b)
    diff = float(np.mean(mb) - np.mean(ma))
    diff_other = float(np.mean([x.mean for x in bo]) - np.mean([x.mean for x in ao]))
    extended = None
    ext_a = [m for m, v in zip(ma, views_a) if v.genus == 0]
    ext_b = [m for m, v in zip(mb, views_b) if v.genus == 0]
    if len(ext_a) < len(ma) or len(ext_b) < len(mb):
        if min(len(ext_a), len(ext_b)) < 2:
            ext_verdict = None                       # fewer than two extended conformers
        else:
            eo_a = [x.mean for x, v in zip(ao, views_a) if v.genus == 0]
            eo_b = [x.mean for x, v in zip(bo, views_b) if v.genus == 0]
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
class EvidenceStatement:
    """One Atlas observation, as a card shows it. Filled by bundle.py only."""
    observation_id: str
    taxon: str
    compounds: tuple                    # Atlas compound ids involved
    assay: Optional[str]
    effect: Optional[str]               # as recorded, e.g. increased, no_effect
    significance: Optional[str]         # as recorded, e.g. significant, not_tested
    state: EvidenceState
    context: Optional[str] = None       # locality, population, sex, as recorded
    citation: Optional[str] = None


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
