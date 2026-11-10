"""Self-check: write a synthetic IQ recording of a straight pass, run both scripts on it, compare with the known answer.

Usage: python3 make_test_iq.py [folder]   (about 70 MB, deleted afterwards). Also checks slr_fit.py on a synthetic laser pass.
Pass: v 7,660 m/s, d 420 km, f0 437.800 MHz + 1,200 Hz, 300 s at 96 kS/s, 16-bit, with noise and one fixed receiver line.
  Run with the file alone (centre from the name, span from the scan), then with --span 15000: the same v and d.
Burst pass (--burst): v 7,400 m/s, d 600 km, f0 436.650 MHz, 300 s at 24 kS/s cf32, a 1 s 2-FSK burst every 10 s
(+-500 Hz, 1,200 baud) under a steady line 5 times stronger. Run again without --burst: burst mode comes on by itself.
"""
import json, os, subprocess, sys

import numpy as np

C, FC, FS, T = 299792458.0, 437.8e6, 96000, 300.0
V, D, T0, OFF = 7660.0, 420e3, 150.0, 1200.0


def doppler(t, t0, v, d, f0):
    tp = t - t0; te = tp.copy()
    for _ in range(6):
        te = tp - (np.sqrt(d * d + v * v * te * te) - d) / C
    return f0 / (1 + v * v * te / np.sqrt(d * d + v * v * te * te) / C)


def burst_check(d, here):
    fs, fc, v, dist, t0 = 24000, 436.65e6, 7400.0, 600e3, 150.0
    n = int(T * fs); t = np.arange(n) / fs
    rng = np.random.default_rng(2)
    bits = rng.integers(0, 2, int(T * 1200) + 1)[(t * 1200).astype(int)] * 2 - 1
    on = (t % 10.0) < 1.0
    f = doppler(t, t0, v, dist, fc) - fc + 500.0 * bits
    x = 0.1 * on * np.exp(2j * np.pi * np.cumsum(f) / fs) + 0.5 * np.exp(2j * np.pi * 3000.0 * t)
    x += 0.1 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    raw = os.path.join(d, "burst_test.cf32")
    x.astype(np.complex64).tofile(raw)
    args = [sys.executable, os.path.join(here, "iq_to_trace.py"), raw, "--rate", str(fs), "--center", str(fc), "--format", "cf32", "--span", "11500"]
    subprocess.run(args + ["--burst", "--out", d], check=True)
    auto = os.path.join(d, "auto_burst"); os.makedirs(auto, exist_ok=True)
    subprocess.run(args + ["--out", auto], check=True)
    os.remove(raw)
    subprocess.run([sys.executable, os.path.join(here, "fit_pass.py"), os.path.join(d, "burst_test_trace.csv")], check=True)
    r = json.load(open(os.path.join(d, "burst_test_fit.json")))["method_D"]
    assert json.load(open(os.path.join(auto, "burst_test_summary.json")))["burst"], "burst mode did not come on by itself"
    assert open(os.path.join(auto, "burst_test_trace.csv")).read() == open(os.path.join(d, "burst_test_trace.csv")).read(), "auto burst trace differs"
    print(f"known: v {v / 1e3:.3f} km/s, d {dist / 1e3:.1f} km, t0 {t0:.2f} s, f0 {fc:.1f} Hz", flush=True)
    assert abs(r["v_m_s"] - v) / v < 0.01 and abs(r["d_m"] - dist) / dist < 0.03 and abs(r["t0_s"] - t0) < 1.0 and abs(r["f0_Hz"] - fc) < 30, r


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "test_out"
    os.makedirs(d, exist_ok=True)
    wav = os.path.join(d, "test_437800000Hz.wav")
    n = int(T * FS)
    t = np.arange(n) / FS
    f = doppler(t, T0, V, D, FC + OFF) - FC
    rng = np.random.default_rng(1)
    x = 0.05 * np.exp(2j * np.pi * np.cumsum(f) / FS) + 0.02 * np.exp(2j * np.pi * 5306.0 * t)
    x += 0.2 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    iq = np.empty((n, 2), "<i2"); iq[:, 0] = np.clip(x.real * 8000, -32767, 32767); iq[:, 1] = np.clip(x.imag * 8000, -32767, 32767)
    hdr = b"RIFF" + (36 + iq.nbytes).to_bytes(4, "little") + b"WAVEfmt " + (16).to_bytes(4, "little") + \
          (1).to_bytes(2, "little") + (2).to_bytes(2, "little") + FS.to_bytes(4, "little") + (FS * 4).to_bytes(4, "little") + \
          (4).to_bytes(2, "little") + (16).to_bytes(2, "little") + b"data" + iq.nbytes.to_bytes(4, "little")
    with open(wav, "wb") as fh:
        fh.write(hdr); fh.write(iq.tobytes())
    here = os.path.dirname(os.path.abspath(__file__))
    subprocess.run([sys.executable, os.path.join(here, "iq_to_trace.py"), wav, "--out", d], check=True)
    fixed = os.path.join(d, "span15000"); os.makedirs(fixed, exist_ok=True)
    subprocess.run([sys.executable, os.path.join(here, "iq_to_trace.py"), wav, "--span", "15000", "--out", fixed], check=True)
    os.remove(wav)
    for o in (d, fixed):
        subprocess.run([sys.executable, os.path.join(here, "fit_pass.py"), os.path.join(o, "test_437800000Hz_trace.csv")], check=True)
    fit = json.load(open(os.path.join(d, "test_437800000Hz_fit.json")))
    r, r2 = fit["method_D"], json.load(open(os.path.join(fixed, "test_437800000Hz_fit.json")))["method_D"]
    print(f"known: v {V / 1e3:.3f} km/s, d {D / 1e3:.1f} km, t0 {T0:.2f} s, f0 {FC + OFF:.1f} Hz", flush=True)
    assert abs(r["v_m_s"] - V) / V < 0.01 and abs(r["d_m"] - D) / D < 0.02 and abs(r["t0_s"] - T0) < 0.5 and abs(r["f0_Hz"] - FC - OFF) < 20, r
    assert (r["v_m_s"], r["d_m"]) == (r2["v_m_s"], r2["d_m"]), ("file alone and --span 15000 differ", r, r2)
    assert fit["path"]["path"] == "no bending detected", fit["path"]
    assert fit["pass_found"], fit["no_pass_reasons"]
    flat = os.path.join(d, "flat_trace.csv")                     # a fixed tone is no pass
    with open(flat, "w") as fh:
        fh.write("t_s,f_Hz\n" + "".join(f"{i * 0.25},{1.4e7 + (3, -2, 5, -4)[i % 4]}\n" for i in range(400)))
    import fit_pass
    assert not fit_pass.run(flat, log=lambda *_: None)["pass_found"]
    burst_check(d, here)
    # laser ranging: a synthetic CRD pass, 5 s shots, d 500 km, v 7.4 km/s
    import slr_fit
    crd = os.path.join(d, "test.np2")
    with open(crd, "w") as fh:
        fh.write("H2 TEST 0000\nH3 testsat 0\nH4 1 2026 1 1 0 0 0\n")
        for tt in np.arange(-150.0, 150.1, 5.0):
            fh.write(f"11 {43200 + tt:.6f} {2 * np.hypot(500.0, 7.4 * tt) / slr_fit.C_KM_S:.12f} PS 2\n")
    q = slr_fit.run(crd, 0, d)["passes"][0]
    assert abs(q["v_km_s"] - 7.4) < 1e-3 and abs(q["d_km"] - 500.0) < 0.1, q
    print("self-check passed")


if __name__ == "__main__":
    main()
