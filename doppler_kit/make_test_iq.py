"""Self-check: write a synthetic IQ recording of a straight pass, run both scripts on it, compare with the known answer.

Usage: python3 make_test_iq.py [folder]   (about 70 MB, deleted afterwards). Also checks slr_fit.py on a synthetic laser pass.
Pass: v 7,660 m/s, d 420 km, f0 437.800 MHz + 1,200 Hz, 300 s at 96 kS/s, 16-bit, with noise and one fixed receiver line.
"""
import json, os, subprocess, sys

import numpy as np

C, FC, FS, T = 299792458.0, 437.8e6, 96000, 300.0
V, D, T0, OFF = 7660.0, 420e3, 150.0, 1200.0


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "test_out"
    os.makedirs(d, exist_ok=True)
    wav = os.path.join(d, "test_437800000Hz.wav")
    n = int(T * FS)
    t = np.arange(n) / FS
    tp = t - T0; te = tp.copy()
    for _ in range(6):
        te = tp - (np.sqrt(D * D + V * V * te * te) - D) / C
    f = (FC + OFF) / (1 + V * V * te / np.sqrt(D * D + V * V * te * te) / C) - FC
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
    subprocess.run([sys.executable, os.path.join(here, "fit_pass.py"), os.path.join(d, "test_437800000Hz_trace.csv")], check=True)
    os.remove(wav)
    r = json.load(open(os.path.join(d, "test_437800000Hz_fit.json")))["method_D"]
    print(f"known: v {V / 1e3:.3f} km/s, d {D / 1e3:.1f} km, t0 {T0:.2f} s, f0 {FC + OFF:.1f} Hz")
    assert abs(r["v_m_s"] - V) / V < 0.01 and abs(r["d_m"] - D) / D < 0.02 and abs(r["t0_s"] - T0) < 0.5 and abs(r["f0_Hz"] - FC - OFF) < 20, r
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
