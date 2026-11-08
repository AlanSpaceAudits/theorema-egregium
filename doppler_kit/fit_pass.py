"""Speed, closest distance, time of closest approach and transmitted frequency from a Doppler trace alone.

Usage:
  python3 fit_pass.py TRACE.csv [--c 299792458] [--f0 HZ] [--out DIR]

TRACE.csv needs columns t_s and f_Hz (as written by iq_to_trace.py). Any source works: a trace read
off any waterfall, an SDR log, a spreadsheet, as long as it has time in seconds and received frequency in Hz.

Model: a source moving in a straight line at constant speed v, closest distance d at time t0,
transmitting f0. Received frequency, with the light (or sound) travel time included:
  f(t) = f0 / (1 + v^2 tau / (c sqrt(d^2 + v^2 tau^2)))
where tau is the emission time measured from closest approach. No station position, no orbit, no map.

Two estimates:
  A  plateaus and half-plateau width: v = c (f+ - f-)/(f+ + f-), f0 = 2 f+ f- / (f+ + f-),
     d = dt * v * sqrt(3)/2 with dt the time between the two points where the range rate is +-v/2.
     (Needs the trace to reach both plateaus; a pass cut off by the horizon gives a low v here.)
  D  least squares of the full model above, free t0, v, d, f0, started from A.
Path test: a straight path leaves only scatter around the fit; a bending path leaves a smooth pattern.
  A degree-6 polynomial through the leftover, F-test with the independent point count: p < 0.001 = bending,
  otherwise no bending detected. detectable_rms_Hz: the smallest smooth-pattern rms this test calls bending.
Writes <name>_fit.json, <name>_scurve.png (the S-curve with plateaus, middle, half-plateau points and the width marked,
  plus what the fit leaves over) and, when iq_to_trace.py saved a waterfall, <name>_waterfall_marked.png.
"""
import argparse, csv, json, os

import numpy as np
from scipy.optimize import least_squares


def model(t, t0, v, d, f0, c):
    tp = t - t0
    te = tp.copy()
    for _ in range(6):
        te = tp - (np.sqrt(d * d + v * v * te * te) - d) / c
    return f0 / (1 + v * v * te / np.sqrt(d * d + v * v * te * te) / c)


def method_a(t, f, c, f0=None):
    n = len(t); k = max(3, round(n * 0.02))
    fp, fm = f[:k].mean(), f[-k:].mean()
    v = c * (fp - fm) / (fp + fm)
    f0 = f0 or 2 * fp * fm / (fp + fm)
    i = np.argmin(np.abs(f - f0)); t0 = t[i]
    def cross(F):
        s = np.where(np.diff(np.sign(f - F)) != 0)[0]
        if not len(s):
            return np.nan
        j = s[np.argmin(np.abs(t[s] - t0))]
        return t[j] + (t[j + 1] - t[j]) * (F - f[j]) / (f[j + 1] - f[j])
    dt = cross(f0 / (1 + v / (2 * c))) - cross(f0 / (1 - v / (2 * c)))
    return dict(v=v, f0=f0, t0=t0, d=abs(dt) * v * np.sqrt(3) / 2, dt=abs(dt))


def straightness(t, left, indep, k=6):
    """Straight path or bending path? A straight path leaves only scatter; a bending one leaves a smooth pattern.
    F-test of a degree-k polynomial through what the fit leaves over, with the effective number of independent points."""
    from scipy.stats import f as fdist
    x = (t - t.mean()) / np.ptp(t)
    smooth = np.polyval(np.polyfit(x, left, k), x)
    n = len(t) / indep
    df1, df2 = k + 1, max(1.0, n - k - 1)
    rss0, rss1 = np.sum(left ** 2), np.sum((left - smooth) ** 2)
    F = ((rss0 - rss1) / df1) / (rss1 / df2)
    p = float(fdist.sf(F, df1, df2))
    amp = float(np.sqrt(np.mean(smooth ** 2)))
    scatter = float(np.sqrt(rss1 / len(t)))
    # Detection limit. N = len(t) points, A = pattern rms, s = scatter rms.
    #   the polynomial is a least-squares projection, so rss0 - rss1 = sum(smooth^2) = N A^2
    #   rss1 = N s^2
    #   F = (N A^2 / df1) / (N s^2 / df2)
    #   F = (A^2 / s^2) (df2 / df1)
    #   bending when F >= F_crit = F.ppf(0.999; df1, df2)
    #   A^2 = s^2 F_crit df1 / df2
    #   A_min = s sqrt(F_crit df1 / df2)
    detect = scatter * float(np.sqrt(fdist.ppf(0.999, df1, df2) * df1 / df2))
    return dict(path="bending" if p < 0.001 else "no bending detected", p_value=p, pattern_rms_Hz=amp,
                scatter_rms_Hz=scatter, detectable_rms_Hz=detect)


def run(trace, c=299792458.0, f0=None, out=None, log=print):
    """Fit one trace. Returns the result dict, with the output file paths under "files"."""
    rows = list(csv.DictReader(open(trace)))
    t = np.array([float(r["t_s"]) for r in rows]); f = np.array([float(r["f_Hz"]) for r in rows])
    o = np.argsort(t); t, f = t[o], f[o]
    A = method_a(t, f, c, f0)
    p0 = [A["t0"], A["v"], A["d"] if np.isfinite(A["d"]) and A["d"] > 0 else A["v"] * 30, A["f0"]]
    free = [0, 1, 2] if f0 else [0, 1, 2, 3]
    def res(q):
        p = list(p0); [p.__setitem__(k, q[i]) for i, k in enumerate(free)]
        return model(t, p[0], abs(p[1]), abs(p[2]), p[3], c) - f
    r = least_squares(res, [p0[k] for k in free], x_scale="jac")
    p = list(p0); [p.__setitem__(k, r.x[i]) for i, k in enumerate(free)]
    p[1], p[2] = abs(p[1]), abs(p[2])
    left = f - model(t, *p, c)
    rms = float(np.sqrt(np.mean(left ** 2)))
    step = float(np.median(np.diff(t)))
    indep = max(1.0, 1.0 / step) if step < 1 else 1.0          # iq_to_trace: 1 s windows every 0.25 s, so 4 points share one window
    J = r.jac
    cov = np.linalg.pinv(J.T @ J) * rms ** 2 * indep
    err = dict(zip([["t0", "v", "d", "f0"][k] for k in free], np.sqrt(np.diag(cov))))
    S = straightness(t, left, indep)
    rng = [float(np.hypot(p[2], p[1] * (T - p[0]))) for T in (t[0], t[-1])]
    D = dict(t0_s=p[0], v_m_s=p[1], d_m=p[2], f0_Hz=p[3], rms_Hz=rms, n=len(t),
             err={k + ("_s" if k == "t0" else "_m_s" if k == "v" else "_m" if k == "d" else "_Hz"): float(v) for k, v in err.items()},
             range_first_m=rng[0], range_last_m=rng[1])
    if rows and rows[0].get("utc"):
        from datetime import datetime, timedelta
        D["t0_utc"] = (datetime.fromisoformat(rows[0]["utc"]) + timedelta(seconds=p[0] - float(rows[0]["t_s"]))).isoformat(timespec="milliseconds")
    out = out or os.path.dirname(os.path.abspath(trace)); name = os.path.splitext(os.path.basename(trace))[0].replace("_trace", "")
    M = marks(p, c)
    res_all = dict(c_m_s=c, method_A={k: float(v) for k, v in A.items()}, method_D=D, path=S,
                   marks=dict(f_plus_Hz=M["fp"], f_minus_Hz=M["fm"], half_plateau_Hz=M["fh"], half_plateau_t_s=M["th"], dt_s=M["dt"],
                              slope_Hz_s=M["slope"], v_plateaus_m_s=M["v_plat"], d_width_m=M["d_width"], d_slope_m=M["d_slope"]))
    files = dict(fit=os.path.join(out, f"{name}_fit.json"), scurve=os.path.join(out, f"{name}_scurve.png"))
    json.dump(res_all, open(files["fit"], "w"), indent=1)
    log(f"A (plateaus): v {A['v'] / 1e3:.3f} km/s, d {A['d'] / 1e3:.1f} km, f0 {A['f0']:.1f} Hz")
    log(f"D (full fit): v {p[1] / 1e3:.3f} ± {err['v'] / 1e3:.3f} km/s, d {p[2] / 1e3:.1f} ± {err['d'] / 1e3:.1f} km, "
        f"t0 {p[0]:.2f} s, f0 {p[3]:.1f} Hz, left over {rms:.1f} Hz rms over {len(t)} points")
    log(f"Path: {S['path']} (smooth pattern {S['pattern_rms_Hz']:.1f} Hz rms against scatter {S['scatter_rms_Hz']:.1f} Hz, p = {S['p_value']:.2g}; "
        f"detectable from {S['detectable_rms_Hz']:.1f} Hz rms)")
    centre = None
    if rows and rows[0].get("f_offset_Hz"):
        centre = float(rows[0]["f_Hz"]) - float(rows[0]["f_offset_Hz"])
    if figures(t, f, p, c, left, out, name, centre, S):
        files["waterfall_marked"] = os.path.join(out, f"{name}_waterfall_marked.png")
    res_all["files"] = files
    return res_all


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("trace"); ap.add_argument("--c", type=float, default=299792458.0, help="wave speed, m/s (radio/light default; 343 for sound)")
    ap.add_argument("--f0", type=float, help="fix the transmitted frequency, Hz (otherwise fitted)")
    ap.add_argument("--out", help="output folder (default: next to the trace)")
    a = ap.parse_args()
    run(a.trace, a.c, a.f0, a.out)


def marks(p, c):
    """Plateaus, middle, half-plateau points and the width, all from the fitted curve."""
    t0, v, d, f0 = p
    fp, fm = f0 / (1 - v / c), f0 / (1 + v / c)
    fh = [f0 / (1 - v / (2 * c)), f0 / (1 + v / (2 * c))]
    tt = np.linspace(t0 - 50 * d / v, t0 + 50 * d / v, 200001)
    ff = model(tt, *p, c)
    th = [float(np.interp(-F, -ff, tt)) for F in fh]            # ff falls with time
    slope = float(np.gradient(model(np.array([t0 - 0.5, t0, t0 + 0.5]), *p, c), 0.5)[1])
    return dict(fp=fp, fm=fm, f0=f0, fh=fh, th=th, dt=th[1] - th[0], slope=slope,
                d_width=(th[1] - th[0]) * v * np.sqrt(3) / 2, d_slope=f0 * v * v / (c * abs(slope)), v_plat=c * (fp - fm) / (fp + fm))


def figures(t, f, p, c, left, out, name, centre, S=None):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    M = marks(p, c); t0, v, d, f0 = p
    halo = [pe.withStroke(linewidth=3, foreground="white")]
    # 1: the S-curve, with what the fit leaves over
    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
    lo, hi = min(t[0], M["th"][0]) - 20, max(t[-1], M["th"][1]) + 20
    tt = np.linspace(lo, hi, 1500)
    k = lambda F: (F - f0) / 1e3
    ax[0].plot(t - t0, k(f), ".", ms=2.5, color="#2b6cb0", label="measured")
    ax[0].plot(tt - t0, k(model(tt, *p, c)), color="#c53030", lw=1.5, label="straight-pass fit")
    for F, lab in ((M["fp"], "plateau f₊ (coming)"), (M["fm"], "plateau f₋ (going)"), (f0, "middle f₀")):
        ax[0].axhline(k(F), color="#555", lw=1, ls="--")
        ax[0].text(hi - t0, k(F), " " + lab, va="center", ha="left", fontsize=9)
    ax[0].axvline(0, color="#555", lw=0.8, ls=":")
    for i, (T, F) in enumerate(zip(M["th"], M["fh"])):
        ax[0].plot(T - t0, k(F), "o", ms=13, mfc="none", mec="#d69e2e", mew=2.5)
        ax[0].annotate(f"half-plateau {i + 1}\nrange rate {'−' if i == 0 else '+'}v/2", (T - t0, k(F)), (20 if i else -20, 18 if i else -30),
                       textcoords="offset points", ha="left" if i else "right", fontsize=9, color="#975a16", path_effects=halo)
    yb = k(M["fp"]) * 1.12 + 0.1
    ax[0].annotate("", (M["th"][0] - t0, yb), (M["th"][1] - t0, yb), arrowprops=dict(arrowstyle="<->", color="#975a16", lw=2))
    ax[0].text((M["th"][0] + M["th"][1]) / 2 - t0, yb, f"Δt = {M['dt']:.2f} s", ha="center", va="bottom", fontsize=11, color="#975a16")
    for T in M["th"]:
        ax[0].plot([T - t0] * 2, [k(f0 / (1 - v / (2 * c))) if T == M["th"][0] else k(f0 / (1 + v / (2 * c))), yb], color="#975a16", lw=1, ls=":")
    txt = (f"v = c (f₊ − f₋)/(f₊ + f₋) = {M['v_plat'] / 1e3:.3f} km/s\n"
           f"d = Δt · v · √3/2 = {M['d_width'] / 1e3:.1f} km\n"
           f"d = f₀ v² / (c · |slope at middle|) = {M['d_slope'] / 1e3:.1f} km   (slope {M['slope']:.1f} Hz/s)\n"
           f"full fit: v {v / 1e3:.3f} km/s, d {d / 1e3:.1f} km, f₀ {f0:,.1f} Hz")
    ax[0].text(0.02, 0.03, txt, transform=ax[0].transAxes, fontsize=9.5, family="monospace", va="bottom",
               bbox=dict(boxstyle="round", fc="white", ec="#aaa"))
    ax[0].set_ylabel("received − f₀ (kHz)"); ax[0].legend(loc="upper right"); ax[0].set_xlim(lo - t0, hi - t0 + (hi - lo) * 0.18)
    ax[0].set_ylim(k(M["fm"]) * 1.25 - 0.1, yb + abs(k(M["fp"])) * 0.15 + 0.3)
    ax[0].set_title("Doppler S-curve: plateaus, middle, half-plateau points and width")
    ax[1].plot(t - t0, left, ".", ms=2.5, color="#2b6cb0"); ax[1].axhline(0, color="k", lw=0.8)
    if S:
        ax[1].set_title(f"left over: {S['path']} (smooth pattern {S['pattern_rms_Hz']:.1f} Hz rms, scatter {S['scatter_rms_Hz']:.1f} Hz rms, "
                        f"detectable from {S['detectable_rms_Hz']:.1f} Hz rms)", fontsize=10)
    ax[1].set_ylabel("left over (Hz)"); ax[1].set_xlabel("time from closest approach (s)")
    fig.tight_layout(); fig.savefig(os.path.join(out, f"{name}_scurve.png"), dpi=120); plt.close(fig)
    # 2: the same marks on the waterfall, if iq_to_trace.py saved one
    wf = os.path.join(out, f"{name}_waterfall.npz")
    if centre is None or not os.path.exists(wf):
        return False
    W = np.load(wf)
    fig, ax = plt.subplots(figsize=(9, 11))
    ax.imshow(W["db"].astype(np.float32), aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=15,
              extent=[W["f"][0] / 1e3, W["f"][-1] / 1e3, W["t"][0], W["t"][-1]])
    x = lambda F: (F - centre) / 1e3
    ty = np.linspace(W["t"][0], W["t"][-1], 1500)
    ax.plot(x(model(ty, *p, c)), ty, color="#ff5a5a", lw=1, alpha=0.8, label="straight-pass fit")
    for F, lab in ((M["fp"], "f₊"), (M["fm"], "f₋"), (f0, "f₀")):
        if W["f"][0] < F - centre < W["f"][-1]:
            ax.axvline(x(F), color="white", lw=1, ls="--", alpha=0.8)
            ax.text(x(F) + 0.15, W["t"][0] + 0.02 * (W["t"][-1] - W["t"][0]), lab, color="white", ha="left", va="bottom", fontsize=13,
                    path_effects=[pe.withStroke(linewidth=3, foreground="black")])
    for i, (T, F) in enumerate(zip(M["th"], M["fh"])):
        ax.plot(x(F), T, "o", ms=16, mfc="none", mec="#ffd23f", mew=2.5)
        ax.text(x(F) + (0.6 if i else -0.6), T, f"half-plateau {i + 1}", color="#ffd23f", ha="left" if i else "right", va="center",
                fontsize=10, path_effects=[pe.withStroke(linewidth=3, foreground="black")])
    xb = W["f"][-1] / 1e3 * 0.85
    for T, F in zip(M["th"], M["fh"]):
        ax.plot([x(F), xb], [T, T], color="#ffd23f", lw=1, ls=":")
    ax.annotate("", (xb, M["th"][0]), (xb, M["th"][1]), arrowprops=dict(arrowstyle="<->", color="#ffd23f", lw=2))
    ax.text(xb - 0.3, sum(M["th"]) / 2, f"Δt = {M['dt']:.2f} s", color="#ffd23f", ha="right", va="center", fontsize=12,
            path_effects=[pe.withStroke(linewidth=3, foreground="black")])
    ax.set_xlabel("offset from centre (kHz)"); ax.set_ylabel("time from start of recording (s)")
    ax.set_title(f"Waterfall with the half-plateau points: v {v / 1e3:.3f} km/s, d {d / 1e3:.1f} km")
    fig.tight_layout(); fig.savefig(os.path.join(out, f"{name}_waterfall_marked.png"), dpi=110); plt.close(fig)
    return True

if __name__ == "__main__":
    main()
