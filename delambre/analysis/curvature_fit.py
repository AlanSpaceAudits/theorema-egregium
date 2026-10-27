"""Gauss's measure of curvature from the observed angles of the Dunkerque-Barcelona chain
(Delambre & Méchain, 1806-1810).

The observed angle sum of each triangle, minus 180 degrees, is fitted against the triangle's
area by least squares through zero: excess = K x area. No radius, no latitude and no star enters.
Angles: the observed column of the Tableau général des triangles (vol. 1, pp. 513-542).
Areas: Heron's formula on the three printed sides (arcs in toises, vol. 2, Tableau complet,
pp. 800-839). A side not printed in its own triangle is taken from another triangle that prints it.
Triangles with a concluded (unobserved) angle have no observed sum and are left out.
1 toise = 10^7 / 5130740 m = 1.949036591 m (quarter meridian adopted by the commission, vol. 3, p. 107).

Run: python3 curvature_fit.py
"""
import json
import math
from pathlib import Path

TOISE = 1e7 / 5130740
RAD = math.pi / 180 / 3600
R_ = json.loads((Path(__file__).parent / "review.json").read_text())["triangles"]

side = {}
for t in R_:
    for s in t["sides"]:
        side.setdefault(frozenset((s["a"], s["b"])), s["toises"])

rows, skipped = [], []
for t in R_:
    if t["sum_minus_180_arcsec"] is None:
        continue
    v = t["vertices"]
    ab = [side.get(frozenset(p)) for p in ((v[0], v[1]), (v[1], v[2]), (v[0], v[2]))]
    if None in ab:
        skipped.append(t["id"])
        continue
    a, b, c = ab
    s = (a + b + c) / 2
    A = math.sqrt(s * (s - a) * (s - b) * (s - c))      # square toises
    rows.append((A, t["sum_minus_180_arcsec"] * RAD, t["excess_printed_arcsec"] * RAD, t["id"]))


def fit(r, i=1):
    K = sum(q[0] * q[i] for q in r) / sum(q[0] ** 2 for q in r)
    sd = math.sqrt(sum((q[i] - K * q[0]) ** 2 for q in r) / (len(r) - 1))
    return K, sd / math.sqrt(sum(q[0] ** 2 for q in r)), sd


def show(label, K, sK=None):
    Rt = 1 / math.sqrt(K)
    line = f"{label}: K = {K:.4e} per square toise = {K / TOISE ** 2:.4e} per square metre; R = {Rt:,.0f} toises = {Rt * TOISE / 1e3:,.1f} km"
    if sK:
        lo, hi = 1 / math.sqrt(K + sK), (1 / math.sqrt(K - sK) if K > sK else float("inf"))
        line += f" (1-sigma {lo:,.0f} to {hi:,.0f} toises = {lo * TOISE / 1e3:,.0f} to {hi * TOISE / 1e3:,.0f} km)"
    print(line)


n = len(rows)
print(f"{n} triangles with observed sums and three printed sides; left out for a missing printed side: {', '.join(skipped) or 'none'}")
print(f"mean area {sum(r[0] for r in rows) / n:,.0f} square toises = {sum(r[0] for r in rows) / n * TOISE ** 2 / 1e6:,.0f} km2")
K, sK, sd = fit(rows)
show("observed angle sums", K, sK)
print(f"  K / sigma = {K / sK:.1f} (distance from K = 0, a plane); scatter of sum about the fit {sd / RAD:.2f} arcsec")
Sx = sum(r[0] for r in rows); Sy = sum(r[1] for r in rows); Sxx = sum(r[0] ** 2 for r in rows); Sxy = sum(r[0] * r[1] for r in rows)
d = n * Sxx - Sx * Sx; k = (n * Sxy - Sx * Sy) / d; b = (Sy - k * Sx) / n
s2 = math.sqrt(sum((r[1] - b - k * r[0]) ** 2 for r in rows) / (n - 2)); sk = s2 * math.sqrt(n / d)
print(f"with a constant bias per triangle: bias {b / RAD:+.2f} arcsec")
show("  slope", k, sk)
show("the survey's own printed excess (for comparison only)", fit(rows, 2)[0])
