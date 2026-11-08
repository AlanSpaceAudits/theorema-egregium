"""Doppler kit window: pick a recording, get the waterfall, the S-curve, the marked waterfall and every number.

Usage: python3 doppler_gui.py   (Python 3 with tkinter, numpy, scipy, matplotlib; ffmpeg for compressed audio)

Takes IQ recordings (.wav 2-channel, .raw/.cs16/.cf32/.cu8/.sigmf-data), audio (.wav 1-channel, .ogg/.mp3/.flac/.m4a/.opus)
and laser-ranging normal points (.np2/.npt/.crd). Runs iq_to_trace.py and fit_pass.py, or slr_fit.py, and shows the results.
"""
import contextlib, io, os, queue, re, threading, traceback
import tkinter as tk
from tkinter import filedialog, ttk

import matplotlib.image as mpimg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import fit_pass, iq_to_trace, slr_fit

LASER = (".np2", ".npt", ".crd")
WAVES = {"radio / light (299,792,458 m/s)": 299792458.0, "sound in air (343 m/s)": 343.0}


class Log(io.TextIOBase):
    def __init__(self, q): self.q = q
    def write(self, s): self.q.put(("log", s)); return len(s)


class App:
    def __init__(self, root):
        self.root, self.q = root, queue.Queue()
        root.title("Doppler kit"); root.geometry("1400x900")
        left = ttk.Frame(root, padding=8); left.pack(side="left", fill="y")
        self.v = {k: tk.StringVar(value=d) for k, d in dict(file="", centre="", start="", f0="", span="", rate="", fmt="auto",
                                                            out="", wave=list(WAVES)[0], c="", which="0", lo="", hi="").items()}
        self.swap, self.burst = tk.BooleanVar(), tk.BooleanVar()
        row = 0
        self.auto = {}
        def add(label, w, key=None):
            nonlocal row
            ttk.Label(left, text=label).grid(row=row, column=0, sticky="w", pady=2); w.grid(row=row, column=1, sticky="we", pady=2)
            if key:
                self.auto[key] = ttk.Label(left, text="", foreground="#777", wraplength=260)
                self.auto[key].grid(row=row, column=2, sticky="w", padx=4)
            row += 1
        f = ttk.Frame(left); ttk.Entry(f, textvariable=self.v["file"], width=34).pack(side="left"); ttk.Button(f, text="Browse…", command=self.pick).pack(side="left")
        add("Recording", f)
        add("Centre / dial frequency (Hz)", ttk.Entry(left, textvariable=self.v["centre"]), "centre")
        add("Start time UTC (optional)", ttk.Entry(left, textvariable=self.v["start"]), "start")
        add("Wave", ttk.Combobox(left, textvariable=self.v["wave"], values=list(WAVES) + ["other (m/s below)"], state="readonly", width=30))
        add("Wave speed m/s (other)", ttk.Entry(left, textvariable=self.v["c"]))
        add("Fix transmitted f₀ (Hz, optional)", ttk.Entry(left, textvariable=self.v["f0"]))
        add("Search ± Hz (blank: scan)", ttk.Entry(left, textvariable=self.v["span"]), "span")
        add("Tone band low (Hz)", ttk.Entry(left, textvariable=self.v["lo"]), "band")
        add("Tone band high (Hz)", ttk.Entry(left, textvariable=self.v["hi"]))
        add("Burst signal (FSK satellites)", ttk.Checkbutton(left, variable=self.burst), "burst")
        add("Raw IQ sample rate (S/s)", ttk.Entry(left, textvariable=self.v["rate"]), "rate")
        add("Raw IQ format", ttk.Combobox(left, textvariable=self.v["fmt"], values=["auto", "cs16", "cf32", "cu8"], state="readonly"), "format")
        add("Q stored before I", ttk.Checkbutton(left, variable=self.swap))
        add("Laser file: pass to plot", ttk.Entry(left, textvariable=self.v["which"]))
        f = ttk.Frame(left); ttk.Entry(f, textvariable=self.v["out"], width=34).pack(side="left"); ttk.Button(f, text="…", width=3, command=self.pick_out).pack(side="left")
        add("Output folder (blank: next to file)", f)
        self.go = ttk.Button(left, text="Run", command=self.start); self.go.grid(row=row, column=0, columnspan=3, sticky="we", pady=8); row += 1
        self.log = tk.Text(left, width=60, height=18, wrap="word", font=("monospace", 9)); self.log.grid(row=row, column=0, columnspan=3, sticky="nswe")
        left.rowconfigure(row, weight=1)
        self.nb = ttk.Notebook(root); self.nb.pack(side="left", fill="both", expand=True)
        self.res = tk.Text(self.nb, wrap="word", font=("monospace", 11)); self.nb.add(self.res, text="Results")
        self.tabs = {}
        root.after(100, self.poll)

    def pick(self):
        p = filedialog.askopenfilename(filetypes=[("Recordings", "*.wav *.raw *.cs16 *.cf32 *.cu8 *.sigmf-data *.ogg *.oga *.mp3 *.flac *.m4a *.opus *.aac"),
                                                  ("Laser ranging", "*.np2 *.npt *.crd"), ("All files", "*")])
        if p:
            self.v["file"].set(p)

    def pick_out(self):
        p = filedialog.askdirectory()
        if p:
            self.v["out"].set(p)

    def num(self, k):
        s = self.v[k].get().strip().replace(",", "")
        return float(s) if s else None

    def start(self):
        self.go.state(["disabled"]); self.log.delete("1.0", "end"); self.res.delete("1.0", "end")
        for w in self.auto.values():
            w.config(text="")
        threading.Thread(target=self.work, daemon=True).start()

    def work(self):
        try:
            path = self.v["file"].get().strip()
            assert path, "pick a recording first"
            out = self.v["out"].get().strip() or None
            with contextlib.redirect_stdout(Log(self.q)):
                if path.lower().endswith(LASER) or re.match(r"https?://", path):
                    r = slr_fit.run(path, int(self.num("which") or 0), out)
                    self.q.put(("laser", r))
                    return
                c = WAVES.get(self.v["wave"].get()) or self.num("c")
                assert c, "give the wave speed"
                lo, hi = self.num("lo"), self.num("hi")
                assert (lo is None) == (hi is None), "give both tone band limits, or neither"
                fmt = self.v["fmt"].get()
                S = iq_to_trace.run(path, self.num("centre"), self.v["start"].get().strip() or None, self.num("span"),
                                    self.swap.get(), out, self.num("rate"), None if fmt == "auto" else fmt,
                                    band=(lo, hi) if lo is not None else None, burst=self.burst.get())
                self.q.put(("auto", S["auto"]))
                self.q.put(("image", ("Waterfall", S["files"]["waterfall"])))
                if S["points"] < 20:
                    raise RuntimeError("too few points on the trace to fit")
                F = fit_pass.run(S["files"]["trace"], c, self.num("f0"), out)
                self.q.put(("done", (S, F)))
        except AssertionError as e:
            self.q.put(("log", f"Stopped: {e}\n"))
        except Exception:
            self.q.put(("log", traceback.format_exc()))
        finally:
            self.q.put(("end", None))

    def poll(self):
        while not self.q.empty():
            kind, x = self.q.get()
            if kind == "log":
                self.log.insert("end", x); self.log.see("end")
            elif kind == "auto":
                for k, text in x.items():
                    if k in self.auto:
                        self.auto[k].config(text="auto: " + text)
            elif kind == "image":
                self.show(*x)
            elif kind == "done":
                S, F = x
                self.show("S-curve", F["files"]["scurve"])
                if "waterfall_marked" in F["files"]:
                    self.show("Marked waterfall", F["files"]["waterfall_marked"])
                self.res.insert("end", report(S, F)); self.nb.select(0)
            elif kind == "laser":
                self.show("Laser ranging", x["files"]["plot"])
                self.res.insert("end", laser_report(x)); self.nb.select(0)
            elif kind == "end":
                self.go.state(["!disabled"])
        self.root.after(100, self.poll)

    def show(self, name, png):
        if name not in self.tabs:
            fr = ttk.Frame(self.nb); self.nb.add(fr, text=name)
            fig = Figure(figsize=(8, 8)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
            cv = FigureCanvasTkAgg(fig, fr); NavigationToolbar2Tk(cv, fr).update(); cv.get_tk_widget().pack(fill="both", expand=True)
            self.tabs[name] = (fr, ax, cv)
        fr, ax, cv = self.tabs[name]
        ax.clear(); ax.axis("off"); ax.imshow(mpimg.imread(png)); cv.draw()
        self.nb.select(fr)


def report(S, F):
    D, A, M, P = F["method_D"], F["method_A"], F["marks"], F["path"]
    e = D["err"]; c = F["c_m_s"]
    km = c > 1e5
    L, V = (lambda m: f"{m / 1e3:,.2f} km") if km else (lambda m: f"{m:,.2f} m"), (lambda m: f"{m / 1e3:.4f} km/s") if km else (lambda m: f"{m:.3f} m/s")
    lines = [
        f"File            {S['file']}",
        f"Points          {S['points']} (point noise {S['point_noise_Hz']} Hz, swing {S['swing_Hz']} Hz)",
        "",
        "FULL FIT (straight pass, constant speed, signal travel time included)",
        f"  speed v                 {V(D['v_m_s'])}  ± {V(e.get('v_m_s', 0))}",
        f"  closest distance d      {L(D['d_m'])}  ± {L(e.get('d_m', 0))}",
        f"  closest approach t₀     {D['t0_s']:.2f} s from start" + (f"  ({D['t0_utc']} UTC)" if D.get("t0_utc") else ""),
        f"  transmitted f₀          {D['f0_Hz']:,.1f} Hz",
        f"  distance at first point {L(D['range_first_m'])}",
        f"  distance at last point  {L(D['range_last_m'])}",
        f"  fit leaves over         {D['rms_Hz']:.1f} Hz rms",
        "",
        "READ OFF THE CURVE",
        f"  plateau f₊ (coming)     {M['f_plus_Hz']:,.1f} Hz",
        f"  plateau f₋ (going)      {M['f_minus_Hz']:,.1f} Hz",
        f"  half-plateau 1          {M['half_plateau_Hz'][0]:,.1f} Hz at {M['half_plateau_t_s'][0]:.2f} s",
        f"  half-plateau 2          {M['half_plateau_Hz'][1]:,.1f} Hz at {M['half_plateau_t_s'][1]:.2f} s",
        f"  width Δt                {M['dt_s']:.2f} s",
        f"  slope at the middle     {M['slope_Hz_s']:.2f} Hz/s",
        f"  v = c (f₊ − f₋)/(f₊ + f₋)      {V(M['v_plateaus_m_s'])}",
        f"  d = Δt · v · √3/2              {L(M['d_width_m'])}",
        f"  d = f₀ v² / (c · |slope|)       {L(M['d_slope_m'])}",
        f"  plateaus from the raw ends:    v {V(A['v'])}, d {L(A['d'])}   (low when the pass is cut short)",
        "",
        "PATH",
        f"  {P['path'].upper()}: smooth pattern left by the fit {P['pattern_rms_Hz']:.1f} Hz rms against scatter {P['scatter_rms_Hz']:.1f} Hz rms (p = {P['p_value']:.2g})",
        f"  This recording can detect a smooth pattern of {P['detectable_rms_Hz']:.1f} Hz rms or more (p below 0.001).",
        "  A straight path leaves only scatter. A bending path leaves a smooth wave.",
    ]
    if S.get("auto"):
        lines += ["", "AUTO VALUES (from the file; a value you type always wins)"] + [f"  {v}" for v in S["auto"].values()]
    for k in ("note",):
        if S.get(k):
            lines += ["", "NOTE: " + S[k]]
    lines += ["", "Files:"] + [f"  {p}" for p in list(S["files"].values()) + list(F["files"].values())]
    return "\n".join(lines) + "\n"


def laser_report(r):
    lines = ["LASER RANGING: distance of each shot = round-trip time × c/2; straight path r² = d² + (v·t)²", ""]
    for i, p in enumerate(r["passes"]):
        lines.append(f"pass {i}: {p['sat']} past {p['station']}: v {p['v_km_s']:.4f} km/s, d {p['d_km']:.1f} km, {p['n']} shots, fit miss {p['rms_m']:.0f} m rms")
    lines += ["", "Files:"] + [f"  {p}" for p in r["files"].values()]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    root = tk.Tk(); App(root); root.mainloop()
