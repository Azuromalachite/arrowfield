"""Build arrowfield-drafts/molecules.yaml: the curated teaching manifest.

Every entry is checked with RDKit before it is written: the SMILES parses,
every heavy atom carries a unique atom map (espmap's P1 gate rules), each
double bond that should have E/Z has it, and the formal charge matches. The
InChIKey and unmapped canonical SMILES are computed here, never typed.
"""
import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

OUT = Path(sys.argv[1])


def chain_ch2(a, b):
    return "".join(f"[CH2:{n}]" for n in range(a, b + 1))


# Sesiidae C18 skeleton: C1 bears O; maps 1-18 = chain, 31 = O on C1,
# 32 = carbonyl / aldehyde O, 21-22 = acetyl CH3 / C=O carbon.
TAIL = f"/[CH:13]=[CH:14]\\{chain_ch2(15, 17)}[CH3:18]"           # Z13
ACE = "[CH3:21][C:22](=[O:32])[O:31]"
CHAIN_SITES = [{"label": f"C{n}-H", "atoms": [f"{n}+H"]} for n in (1, 2, 3, 4, 5, 12, 13, 14, 15)]


def c18_sites(headgroup):
    """Chain sites plus head-group sites. "head O" has ONE definition in every
    molecule: the oxygen bonded to C1, the atom alone (map 31 in the alcohol
    and acetate, 32 in the aldehyde). The pilot compared "O31" across OH and
    OAc and got O+H against a bare ester O; the alcohol's O-H is its own site
    here, "hydroxyl", so the comparison is like-for-like (sesiidae-pilot.md)."""
    head = {"alcohol": [{"label": "head O", "atoms": [31]},
                        {"label": "hydroxyl", "atoms": ["31+H"]}],
            "aldehyde": [{"label": "head O", "atoms": [32]}],
            "acetate": [{"label": "head O", "atoms": [31]},
                        {"label": "carbonyl O", "atoms": [32]},
                        {"label": "acetyl", "atoms": ["21+H", 22, 32]}]}[headgroup]
    return CHAIN_SITES + head

ENTRIES = [
    # --- McMurry companion: small, cheap to compute live --------------------
    dict(id="ethane", name="ethane", smiles="[CH3:1][CH3:2]", lesson="polar-bonds",
         note="the non-polar control"),
    dict(id="ethanol", name="ethanol", smiles="[CH3:1][CH2:2][OH:3]", lesson="polar-bonds",
         sites=[{"label": "hydroxyl", "atoms": ["3+H"]}, {"label": "C-O carbon", "atoms": ["2+H"]}]),
    dict(id="dimethyl_ether", name="dimethyl ether", smiles="[CH3:1][O:2][CH3:3]",
         lesson="polar-bonds", note="same formula as ethanol, no O-H donor"),
    dict(id="acetaldehyde", name="acetaldehyde (ethanal)", smiles="[CH3:1][CH:2]=[O:3]",
         lesson="carbonyl",
         sites=[{"label": "carbonyl O", "atoms": [3]}, {"label": "carbonyl C-H", "atoms": ["2+H"]}]),
    dict(id="acetone", name="acetone (propanone)", smiles="[CH3:1][C:2](=[O:3])[CH3:4]",
         lesson="carbonyl",
         sites=[{"label": "carbonyl O", "atoms": [3]}, {"label": "carbonyl C", "atoms": [2]}]),
    dict(id="acetic_acid", name="acetic acid", smiles="[CH3:1][C:2](=[O:3])[OH:4]",
         lesson="acids-bases", solvent="water",
         note="compare with acetate in the same solvent setting"),
    dict(id="acetate", name="acetate", smiles="[CH3:1][C:2](=[O:3])[O-:4]", lesson="acids-bases",
         solvent="water", charge=-1,
         geometry="geometries/acetate.sdf (B3LYP/def2-TZVP, PCM water, make_geometries.py); "
                  "never SMILES + MMFF for small anions (up to 6.5 kcal/mol, H8)",
         note="an anion: PCM water is required (gas-phase small anions come out unbound)"),
    dict(id="ethyl_acetate", name="ethyl acetate", smiles="[CH3:1][CH2:2][O:3][C:4](=[O:5])[CH3:6]",
         lesson="esters",
         sites=[{"label": "ester O", "atoms": [3]}, {"label": "carbonyl O", "atoms": [5]}]),
    dict(id="ethene", name="ethene", smiles="[CH2:1]=[CH2:2]", lesson="alkenes"),
    dict(id="e2_butene", name="(E)-but-2-ene", smiles="[CH3:1]/[CH:2]=[CH:3]/[CH3:4]",
         lesson="alkenes", stereo={"C2=C3": "E"}),
    dict(id="z2_butene", name="(Z)-but-2-ene", smiles="[CH3:1]/[CH:2]=[CH:3]\\[CH3:4]",
         lesson="alkenes", stereo={"C2=C3": "Z"}),
    dict(id="enolate", name="acetaldehyde enolate", smiles="[CH2:1]=[CH:2][O-:3]",
         lesson="resonance", solvent="water", charge=-1,
         geometry="espmap bench/fixtures/enolate.sdf (B3LYP/def2-TZVP, PCM water); "
                  "never SMILES + MMFF for small anions",
         note="ESP points to O, the charge-controlled site; the alkylation outcome needs a "
              "citation with its conditions before a lesson states it"),
    # --- Sesiidae homework: C18. The Sesiidae Pheromone Atlas owns each
    # compound's identity and biology (atlas_compound); its ESP ensembles are
    # to live in the Atlas's computed/espmap and reach Arrowfield as a teaching
    # bundle. None is computed yet: espmap's only C18 surfaces are its two
    # Z3,Z13-18:OAc fixtures. role says whether the Atlas observes the compound
    # or Arrowfield constructed it as a control. ------------------------------
    dict(id="z3z13_18_oh", name="(3Z,13Z)-octadeca-3,13-dien-1-ol",
         smiles=f"[OH:31][CH2:1][CH2:2]/[CH:3]=[CH:4]\\{chain_ch2(5, 12)}{TAIL}",
         lesson="sesiidae-A", headgroup="alcohol", stereo={"C3=C4": "Z", "C13=C14": "Z"},
         sites=c18_sites("alcohol"), precompute=True,
         role="atlas_observed", atlas_compound="spa.cmpd.z3_z13_18_oh"),
    dict(id="z3z13_18_al", name="(3Z,13Z)-octadeca-3,13-dienal",
         smiles=f"[O:32]=[CH:1][CH2:2]/[CH:3]=[CH:4]\\{chain_ch2(5, 12)}{TAIL}",
         lesson="sesiidae-A", headgroup="aldehyde", stereo={"C3=C4": "Z", "C13=C14": "Z"},
         sites=c18_sites("aldehyde"), precompute=True,
         role="teaching_control", atlas_compound=None,
         note="a teaching control, not an Atlas-observed compound: built to keep the "
              "Z3,Z13 backbone fixed while only the head group changes (act A). The "
              "aldehyde the Atlas records is (2E,13Z)-octadeca-2,13-dienal, "
              "spa.cmpd.e2_z13_18_al"),
    dict(id="z3z13_18_ac", name="(3Z,13Z)-octadeca-3,13-dien-1-yl acetate",
         smiles=f"{ACE}[CH2:1][CH2:2]/[CH:3]=[CH:4]\\{chain_ch2(5, 12)}{TAIL}",
         lesson="sesiidae-A,B", headgroup="acetate", stereo={"C3=C4": "Z", "C13=C14": "Z"},
         sites=c18_sites("acetate"), precompute=True,
         role="atlas_observed", atlas_compound="spa.cmpd.z3_z13_18_oac"),
    dict(id="e3z13_18_ac", name="(3E,13Z)-octadeca-3,13-dien-1-yl acetate",
         smiles=f"{ACE}[CH2:1][CH2:2]/[CH:3]=[CH:4]/{chain_ch2(5, 12)}{TAIL}",
         lesson="sesiidae-B,C", headgroup="acetate", stereo={"C3=C4": "E", "C13=C14": "Z"},
         sites=c18_sites("acetate"), precompute=True,
         role="atlas_observed", atlas_compound="spa.cmpd.e3_z13_18_oac"),
    dict(id="e2z13_18_ac", name="(2E,13Z)-octadeca-2,13-dien-1-yl acetate",
         smiles=f"{ACE}[CH2:1]/[CH:2]=[CH:3]/{chain_ch2(4, 12)}{TAIL}",
         lesson="sesiidae-C", headgroup="acetate", stereo={"C2=C3": "E", "C13=C14": "Z"},
         sites=c18_sites("acetate"), precompute=True,
         role="atlas_observed", atlas_compound="spa.cmpd.e2_z13_18_oac"),
]

# The Atlas's InChIKey for each atlas_compound, checked against the structure
# here when the manifest is built (2026-10-23, Atlas 6775875). A mismatch means
# the two repositories disagree on what the compound is.
ATLAS_INCHIKEY = {
    "spa.cmpd.z3_z13_18_oh": "QBNCGBJHGBGHLS-NCCHCWSYSA-N",
    "spa.cmpd.z3_z13_18_oac": "VVJPJXKHBZNADP-DNNFRFAMSA-N",
    "spa.cmpd.e3_z13_18_oac": "VVJPJXKHBZNADP-ZSTZHTMOSA-N",
    "spa.cmpd.e2_z13_18_oac": "NWRSOYAGXTXEMK-NHRIVICHSA-N",
}


def check(e):
    m = Chem.MolFromSmiles(e["smiles"])
    assert m is not None, e["id"]
    heavy = [a for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    maps = [a.GetAtomMapNum() for a in heavy]
    assert all(maps) and len(set(maps)) == len(maps), f"{e['id']}: maps {maps}"
    assert Chem.GetFormalCharge(m) == e.get("charge", 0), e["id"]
    # every stereo double bond named in the entry has E/Z, and nothing else is open
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    want = e.get("stereo", {})
    got = {}
    by_map = {a.GetAtomMapNum(): a for a in heavy}
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() != Chem.BondStereo.STEREONONE:
            i, j = sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum()))
            got[f"C{i}=C{j}"] = {"STEREOE": "E", "STEREOTRANS": "E",
                                 "STEREOZ": "Z", "STEREOCIS": "Z"}[str(b.GetStereo())]
    assert got == want, f"{e['id']}: stereo {got} != {want}"
    from rdkit.Chem import FindPotentialStereo
    open_ = [s for s in FindPotentialStereo(m) if s.specified == Chem.StereoSpecified.Unspecified]
    assert not open_, f"{e['id']}: unspecified stereo"
    for g in e.get("sites", []):                                # site atoms exist
        for t in g["atoms"]:
            n = int(str(t).split("+")[0])
            assert n in by_map, f"{e['id']}: site {g['label']} names missing map {n}"
    plain = Chem.Mol(m)
    for a in plain.GetAtoms():
        a.SetAtomMapNum(0)
    key = Chem.MolToInchiKey(plain)
    if e["lesson"].startswith("sesiidae"):                      # who says what it is
        assert e.get("role") in ("atlas_observed", "teaching_control"), e["id"]
        if e["role"] == "atlas_observed":
            assert ATLAS_INCHIKEY.get(e["atlas_compound"]) == key, \
                f"{e['id']}: not the Atlas's {e['atlas_compound']}"
        else:
            assert e["atlas_compound"] is None and "teaching control" in e["note"], e["id"]
    return dict(canonical_smiles=Chem.MolToSmiles(plain), inchikey=key,
                formula=rdMolDescriptors.CalcMolFormula(plain))


rows = []
for e in ENTRIES:
    derived = check(e)
    row = dict(id=e["id"], display_name=e["name"], lesson=e["lesson"],
               mapped_smiles=e["smiles"], **derived, charge=e.get("charge", 0),
               solvent=e.get("solvent", "vacuum"))
    for k in ("role", "atlas_compound", "geometry", "headgroup", "stereo", "sites", "note",
              "precompute"):
        if k in e:
            row[k] = e[k]
    if e["lesson"].startswith("sesiidae"):
        row["biological_context"] = ("from the Atlas teaching bundle (not yet exported); "
                                     "never written here by hand")
    rows.append(row)

try:
    import yaml
    text = yaml.safe_dump({"schema": "arrowfield.molecules/0-draft", "molecules": rows},
                          sort_keys=False, allow_unicode=True, width=100)
    OUT.write_text(text, encoding="utf-8")
except ImportError:
    OUT = OUT.with_suffix(".json")
    OUT.write_text(json.dumps({"schema": "arrowfield.molecules/0-draft", "molecules": rows},
                              indent=1, ensure_ascii=False), encoding="utf-8")
print(f"wrote {OUT} ({len(rows)} molecules, all checked)")
