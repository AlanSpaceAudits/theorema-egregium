"""Checks on the transcribed triangles of Delambre & Méchain (1806-1810).

1. Sum: recompute each observed-angle sum and compare with the printed sum line (vol. 1).
2. Sine rule: side / sin(angle - E/3) must agree for the printed sides of a triangle (Legendre).
   Done with the vol. 2 spherical angles (the angles the sides were computed from) and,
   for information, with the vol. 1 observed angles.
3. Shared sides: a side printed in two triangles must have the same length.

Run: python3 check.py
"""
import json
import math
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent
R = json.loads((HERE / "review.json").read_text())
T1 = {}
for line in (HERE / "work" / "t1_transcription.txt").read_text().splitlines():
    if line.startswith("T "):
        n = int(line[2:].split("|")[0])
    elif line.startswith("SUM"):
        T1[n] = line[3:].split(";")[0].strip()


def sec(a):
    d, m, s = a.split()
    return int(d) * 3600 + int(m) * 60 + float(s)


print("1. Observed angle sums against the printed sum line (vol. 1)")
bad = 0
for t in R["triangles"]:
    n = int(t["id"].split(".")[1])
    calc = sum(a * 3600 for a in t["angles_observed_deg"]) - 648000
    pr = sec(T1[n]) - 648000
    if abs(calc - pr) > 0.006:
        bad += 1
        print(f"   {t['id']}: calc {calc:+.2f} printed {pr:+.2f}  {t.get('note', '')}")
print(f"   {len(R['triangles'])} triangles, {bad} mismatches\n")


def sine_check(t, angles_sec, E):
    idx = {v: i for i, v in enumerate(t["vertices"])}
    q = []
    for s in t["sides"]:
        opp = [v for v in t["vertices"] if v not in (s["a"], s["b"])][0]
        A = (angles_sec[idx[opp]] - E / 3) / 206264.806
        q.append(s["toises"] / math.sin(A))
    return (max(q) - min(q)) / (sum(q) / len(q)) if len(q) > 1 else None


print("2. Sine rule: side / sin(angle - E/3), relative spread")
worst2 = []
worst_obs = []
for t in R["triangles"]:
    if len(t["sides"]) < 2:
        continue
    if t.get("angles_spherical_vol2"):
        sp = [sec(a) for a in t["angles_spherical_vol2"]]
        E = sum(sp) - 648000
        worst2.append((sine_check(t, sp, E), t["id"]))
    if t["sum_minus_180_arcsec"] is not None:
        ob = [a * 3600 for a in t["angles_observed_deg"]]
        worst_obs.append((sine_check(t, ob, t["excess_printed_arcsec"]), t["id"]))
worst2.sort(reverse=True)
worst_obs.sort(reverse=True)
print(f"   vol. 2 spherical angles: {len(worst2)} triangles, max spread {worst2[0][0]:.1e} ({worst2[0][1]})")
for v, i in worst2:
    if v > 1e-5:
        print(f"   above 1e-5: {i} {v:.1e}")
print(f"   vol. 1 observed angles: {len(worst_obs)} triangles, max spread {worst_obs[0][0]:.1e} ({worst_obs[0][1]}), "
      f"median {sorted(v for v, _ in worst_obs)[len(worst_obs) // 2]:.1e}\n")

print("3. Sides printed in more than one triangle")
seen = defaultdict(list)
for t in R["triangles"]:
    for s in t["sides"]:
        seen[frozenset((s["a"], s["b"]))].append((s["toises"], t["id"]))
multi = {k: v for k, v in seen.items() if len(v) > 1}
bad3 = 0
for k, v in multi.items():
    vals = [x for x, _ in v]
    d = max(vals) - min(vals)
    if d > 0.0002:
        bad3 += 1
        print(f"   {' - '.join(sorted(k))}: {v}  difference {d:.4f} toise")
print(f"   {len(multi)} shared sides, {bad3} differ by more than 0.0002 toise")
