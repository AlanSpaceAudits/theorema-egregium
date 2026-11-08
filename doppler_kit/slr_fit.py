"""Speed and closest distance of a satellite from laser ranging alone.

Usage:
  python3 slr_fit.py FILE.np2 [--pass N] [--out DIR]
  python3 slr_fit.py https://edc.dgfi.tum.de/en/data/npt2/1166240/

Input: a CRD normal-point file (.np2 / .npt / .crd, from the EUROLAS Data Center or CDDIS; a daily file holds many passes),
or the EDC web page of one pass. Each "11" record is one measurement: time of day (s) and round-trip light time (s).

Method (no orbit and no station position is used):
  distance of each shot      r = round-trip time x c / 2
  straight path past station r(t)^2 = d^2 + (v (t - t0))^2      (Pythagoras)
  so for every shot          v = sqrt(r^2 - d^2) / |t - t0|
  d, v and t0 are fitted to all shots within +-3 minutes of the closest one.

Writes <name>_slr.csv (one row per pass) and <name>_pass<N>.png (range against time with the fit,
the triangle for one shot, and the speed every single shot gives).
"""
import argparse, os, re, urllib.request
from datetime import datetime, timedelta, timezone

import numpy as np
from scipy.optimize import least_squares

C_KM_S = 299792.458
WINDOW_S = 180.0


def read_crd(src):
    """Passes from a CRD file or an EDC page: list of dict(station, sat, day, t (s of day), r (km))."""
    if re.match(r"https?://", src):
        txt = urllib.request.urlopen(urllib.request.Request(src, headers={"User-Agent": "doppler-kit"}), timeout=60).read().decode("utf-8", "replace")
    else:
        txt = open(src, encoding="utf-8", errors="replace").read()
    txt = re.sub(r"<br\s*/?>", "\n", txt).replace("&nbsp;", " ")
    txt = re.sub(r"<[^>]+>", "\n", txt)
    passes, cur, st, sat = [], None, "", ""
    for line in txt.splitlines():
        w = line.split()
        if not w:
            continue
        if w[0] in ("h2", "H2"):
            st = w[1]
        elif w[0] in ("h3", "H3"):
            sat = w[1]
        elif w[0] in ("h4", "H4") and len(w) > 4:
            cur = dict(station=st, sat=sat, day=datetime(int(w[2]), int(w[3]), int(w[4]), tzinfo=timezone.utc), t=[], r=[])
            passes.append(cur)
        elif w[0] == "11" and cur is not None and len(w) > 2:
            cur["t"].append(float(w[1])); cur["r"].append(float(w[2]) * C_KM_S / 2)
    out = []
    for p in passes:
        if len(p["t"]) < 5:
            continue
        t = np.array(p["t"]); t[t < t[0] - 43200] += 86400       # pass across midnight
        o = np.argsort(t); p["t"], p["r"] = t[o], np.array(p["r"])[o]
        out.append(p)
    return out


def fit(t, r):
    m = np.abs(t - t[np.argmin(r)]) <= WINDOW_S
    T, R = t[m], r[m]
    q = least_squares(lambda p: np.hypot(p[0], p[1] * (T - p[2])) - R, [R.min(), 7.0, T[np.argmin(R)]], x_scale="jac")
    d, v, t0 = abs(q.x[0]), abs(q.x[1]), q.x[2]
    return dict(d_km=d, v_km_s=v, t0=t0, rms_m=float(np.sqrt(np.mean(q.fun ** 2)) * 1e3), n=int(m.sum()), used=m)


def figure(p, F, png):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    t, r = p["t"] - F["t0"], p["r"]; d, v = F["d_km"], F["v_km_s"]
    fig, ax = plt.subplots(2, 1, figsize=(10, 8.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 2]))
    tt = np.linspace(t.min(), t.max(), 600)
    ax[0].plot(t, r, "o", ms=4, color="#2b6cb0", label="measured: round-trip time × c/2")
    ax[0].plot(tt, np.hypot(d, v * tt), color="#c53030", lw=1.5, label="straight path: r² = d² + (v·t)²")
    ax[0].axhline(d, color="#555", ls="--", lw=1); ax[0].text(t.min(), d, f" d = {d:.1f} km (closest)", va="bottom", fontsize=9)
    k = int(np.argmax(np.abs(t) * F["used"]))
    tk, rk = t[k], r[k]; along = np.sqrt(max(rk * rk - d * d, 0))
    ax[0].plot([tk, tk], [d, rk], color="#975a16", lw=2)
    ax[0].annotate(f"shot at t = {tk:+.1f} s\nr = {rk:.1f} km\n√(r² − d²) = {along:.1f} km\nv = {along:.1f} / {abs(tk):.1f} = {along / abs(tk):.4f} km/s",
                   (tk, (d + rk) / 2), (-15 if tk > 0 else 15, 0), textcoords="offset points", ha="right" if tk > 0 else "left", va="center",
                   fontsize=9, color="#975a16", bbox=dict(boxstyle="round", fc="white", ec="#d69e2e"))
    ax[0].set_ylabel("distance to the satellite (km)"); ax[0].legend(loc="upper center", fontsize=9)
    when = (p["day"] + timedelta(seconds=float(F["t0"]))).strftime("%Y-%m-%d %H:%M:%S.%f")[:-5]
    ax[0].set_title(f"{p['sat']} past {p['station']}: v {v:.4f} km/s, d {d:.1f} km, closest {when} UTC, {F['n']} shots, fit miss {F['rms_m']:.0f} m rms", fontsize=10)
    far = np.abs(t) > 20
    ax[1].plot(t[far], np.sqrt(np.clip(r[far] ** 2 - d * d, 0, None)) / np.abs(t[far]), "o", ms=4, color="#2b6cb0", label="v = √(r² − d²) / |t|, each shot")
    ax[1].axhline(v, color="#c53030", lw=1.5, label=f"fit {v:.4f} km/s")
    ax[1].set_ylabel("speed (km/s)"); ax[1].set_xlabel("time from closest approach (s)"); ax[1].legend(fontsize=9)
    fig.tight_layout(); fig.savefig(png, dpi=120); plt.close(fig)


def run(src, which=0, out=None, log=print):
    passes = read_crd(src)
    assert passes, "no passes with 5 or more '11' records found"
    name = "edc_" + re.sub(r"\D", "", src.rstrip("/").split("/")[-1]) if re.match(r"https?://", src) else os.path.splitext(os.path.basename(src))[0]
    out = out or (os.getcwd() if re.match(r"https?://", src) else os.path.dirname(os.path.abspath(src))); os.makedirs(out, exist_ok=True)
    rows = []
    for i, p in enumerate(passes):
        F = fit(p["t"], p["r"]); p["F"] = F
        when = (p["day"] + timedelta(seconds=float(F["t0"]))).isoformat(timespec="milliseconds")
        rows.append(f"{i},{p['sat']},{p['station']},{when},{F['v_km_s']:.4f},{F['d_km']:.1f},{F['n']},{F['rms_m']:.0f}")
        log(f"pass {i}: {p['sat']} past {p['station']}, closest {when}, v {F['v_km_s']:.4f} km/s, d {F['d_km']:.1f} km, {F['n']} shots, miss {F['rms_m']:.0f} m rms")
    csvp = os.path.join(out, f"{name}_slr.csv")
    open(csvp, "w").write("pass,satellite,station,closest_utc,v_km_s,d_km,shots,fit_miss_m_rms\n" + "\n".join(rows) + "\n")
    png = os.path.join(out, f"{name}_pass{which}.png")
    figure(passes[which], passes[which]["F"], png)
    return dict(passes=[dict(sat=p["sat"], station=p["station"], **{k: v for k, v in p["F"].items() if k != "used"}) for p in passes],
                files=dict(csv=csvp, plot=png))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("src", help="CRD normal-point file or EDC dataset page URL")
    ap.add_argument("--pass", dest="which", type=int, default=0, help="which pass to plot (default 0)")
    ap.add_argument("--out", help="output folder")
    a = ap.parse_args()
    run(a.src, a.which, a.out)


if __name__ == "__main__":
    main()
