"""Analyse sesiidae_pilot.jsonl: for each comparison and site, is the change
between molecules larger than the spread across one molecule's conformers?

A site difference counts as SEPARATED only if the two molecules' conformer
ranges do not overlap: every conformer of one is above every conformer of the
other. Also reported: the difference of conformer means, each molecule's
range, and the partition-rule sensitivity at that site (spec section 5).
"""
import json
from pathlib import Path

SP = Path("/mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/"
          "669a3716-c2cd-43c5-9905-c4ad082b6f87/scratchpad")
rows = [json.loads(l) for l in (SP / "sesiidae_k10.jsonl").read_text().splitlines()]
by = {}
for r in rows:
    by.setdefault(r["molecule"], []).append(r)
COMPARE = [("stereo E3 vs Z3", "Z3Z13-OAc", "E3Z13-OAc"),
           ("position E2 vs E3", "E3Z13-OAc", "E2Z13-OAc"),
           ("head group OAc vs OH", "Z3Z13-OH", "Z3Z13-OAc")]
Q = "V_area_mean_kcal_mol"


def vals(mol, site, rule="nearest_nucleus", q=Q):
    return [r["sites"][site][rule][q] for r in by[mol] if r["sites"][site][rule][q] is not None]


def part(mol, site):
    return max(abs(r["sites"][site]["radius_scaled"][Q] - r["sites"][site]["nearest_nucleus"][Q])
               for r in by[mol] if r["sites"][site]["nearest_nucleus"][Q] is not None)


print("conformers:", {m: len(v) for m, v in by.items()},
      "| folded (genus>0):", {m: sum(r["genus"] > 0 for r in v) for m, v in by.items()})
print("global V_min / V_max range per molecule:")
for m, v in by.items():
    print(f"  {m:10s} V_min {min(r['V_min'] for r in v):+6.1f}..{max(r['V_min'] for r in v):+6.1f}"
          f"  V_max {min(r['V_max'] for r in v):+6.1f}..{max(r['V_max'] for r in v):+6.1f}")
out = {}
for title, a, b in COMPARE:
    print(f"\n=== {title}  ({b} minus {a}), site mean, kcal/mol")
    print(f"{'site':7s} {'diff of means':>13s} {'range A':>8s} {'range B':>8s} {'partition':>9s}  separated?")
    common = [s for s in by[a][0]["sites"] if s in by[b][0]["sites"]]
    sep_sites = []
    for s in common:
        A, B = vals(a, s), vals(b, s)
        if not A or not B:
            continue
        d = sum(B) / len(B) - sum(A) / len(A)
        ra, rb = max(A) - min(A), max(B) - min(B)
        sep = min(B) > max(A) or max(B) < min(A)
        if sep:
            sep_sites.append((s, d))
        p = max(part(a, s), part(b, s))
        print(f"{s:7s} {d:+13.2f} {ra:8.2f} {rb:8.2f} {p:9.2f}  {'YES' if sep else ''}")
    gmin = [sum(r[k] for r in by[b]) / len(by[b]) - sum(r[k] for r in by[a]) / len(by[a])
            for k in ("V_min", "V_max")]
    print(f"global: diff of mean V_min {gmin[0]:+.2f}, V_max {gmin[1]:+.2f}")
    out[title] = dict(separated=sep_sites, n_sites=len(common))
print("\nSUMMARY:", json.dumps(out))
