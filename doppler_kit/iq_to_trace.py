"""Received frequency against time from an IQ or audio recording.

Accepts a WAV file (2 channels I/Q, or 1 channel audio; 8/16-bit int or 32-bit float), a headerless raw IQ file
(.raw/.cs16/.cf32/.cu8/.sigmf-data: give --rate, and --format if not 16-bit signed), or any audio file ffmpeg reads
(.ogg/.mp3/.flac/.m4a/.opus). For audio, --center is the receiver dial frequency (0 = the audio frequency itself,
for sound recordings), and the strongest tone is tracked.

Usage:
  python3 iq_to_trace.py RECORDING.wav [--center HZ] [--start 2026-06-14T18:10:09Z] [--span 15000] [--swap-iq] [--out DIR]
  python3 iq_to_trace.py RECORDING.raw --rate 48000 --center 437800000 [--format cs16|cf32|cu8]
  python3 iq_to_trace.py RECORDING.wav --band LO HI     (audio: tone between LO and HI Hz; IQ: offsets LO..HI from the centre)
  python3 iq_to_trace.py RECORDING.cf32 --rate 48000 --center 436650000 --format cf32 --burst   (short FSK/packet bursts)

Writes to DIR (default: next to the recording):
  <name>_trace.csv      t_s, utc, f_offset_Hz, f_Hz, snr_dB   (accepted points only)
  <name>_waterfall.png  0.1 s FFT rows, +-span around the centre, with the trace drawn on top
  <name>_summary.json   settings, point noise, number of points

Method (no orbit and no station position is used anywhere):
  1. Decimate to about 48 kS/s (polyphase FIR), remove the DC line (10 s running mean).
  2. Find the signal: 0.25 s spectra, every bin and frame median-normalised, fixed receiver lines masked,
     300 Hz boxcar peak per frame, a smooth S-shaped curve through it (used only to steer the search).
  3. Measure: multiply by the conjugate phase of that curve (dechirp) so the carrier sits near 0 Hz,
     then 1.0 s Hann windows every 0.25 s on a 1 Hz grid, search +-400, +-100, then +-150 Hz,
     each pass with a better reference. Final value = noise-subtracted power centroid in the band.
  4. Keep a point if its peak is more than 14 dB above the frame median, it is not on a fixed
     receiver line, and it is within 4 sigma of a smoothing spline (sigma from second differences).

--band LO HI: one FFT of the whole signal, bins outside LO..HI set to 0 (audio: also the negative bins, the rest x2),
  inverse FFT, band centre (LO+HI)/2 moved to 0 Hz, search +-(HI-LO)/2. For IQ the DC blank moves to -(LO+HI)/2.
  Audio has no DC line and no receiver lines: no DC removal, no DC blank, no fixed-line mask.
  In the dechirp passes, bins whose source frequency (bin + reference) lies outside the band are set to 0; the frame
  floor is the median of the other bins, with no per-bin floor.
--burst (for short bursts under steady lines): 0.5 s Hann frames, 100 Hz boxcar across frequency, each bin divided by
  its running median over +-20 s (steady lines go to 1), excess E = ratio - 1. Strongest 2000 Hz band of E per frame;
  keep the frame if 10 log10(mean E in the band + 1) >= 6 dB; value = centroid of max(E, 0) over that band (for 2-FSK
  the midpoint of the two tones). Outliers: soft-L1 fit of the pass shape from a grid of starts, drop points off by > 800 Hz.
"""
import argparse, json, os, re, struct
from datetime import datetime, timedelta, timezone

import numpy as np
from scipy.interpolate import make_smoothing_spline
from scipy.ndimage import uniform_filter1d, median_filter
from scipy.optimize import least_squares
from scipy.signal import resample_poly

TARGET = 48000
EDGE_S = 1.0
COARSE_THR_DB, A_THR_DB, C_THR_DB = 3.0, 15.0, 14.0
SEARCH_B, SEARCH_C, SEARCH_D = 400.0, 100.0, 150.0
SPUR_HALF = 40.0
CLIP_K = 4.0
BURST_FRAME_S, BURST_SMOOTH_HZ, BURST_MEDIAN_S, BURST_BAND_HZ, BURST_THR_DB = 0.5, 100.0, 20.0, 2000.0, 6.0
BURST_FSCALE_HZ, BURST_GATE_HZ = 300.0, 800.0


# ---------- reading ----------
def wav_info(path):
    """Channels, rate, sample format and data offset. Length comes from the file size (SDR++ size fields overflow past 4 GB)."""
    with open(path, "rb") as fh:
        head = fh.read(1 << 16)
    assert head[:4] == b"RIFF" and head[8:12] == b"WAVE", "not a WAV file"
    i, fmt, data = 12, None, None
    while i + 8 <= len(head):
        cid, size = head[i:i + 4], struct.unpack("<I", head[i + 4:i + 8])[0]
        if cid == b"fmt ":
            tag, ch, rate, _, _, bits = struct.unpack("<HHIIHH", head[i + 8:i + 24])
            fmt = (tag, ch, rate, bits)
        if cid == b"data":
            data = i + 8
            break
        i += 8 + size + (size & 1)
    assert fmt and data, "no fmt/data chunk"
    tag, ch, rate, bits = fmt
    assert ch in (1, 2), "need 2 channels (I and Q) or 1 channel (audio)"
    dt = {8: np.uint8, 16: np.dtype("<i2"), 32: np.dtype("<f4") if tag == 3 else np.dtype("<i4")}[bits]
    n = (os.path.getsize(path) - data) // (ch * np.dtype(dt).itemsize)
    return dict(rate=rate, dtype=dt, offset=data, n=n, bits=bits, float=(tag == 3), ch=ch)


def raw_info(path, rate, fmt):
    dt, bits, flt = {"cs16": (np.dtype("<i2"), 16, False), "cf32": (np.dtype("<f4"), 32, True), "cu8": (np.uint8, 8, False)}[fmt]
    return dict(rate=int(rate), dtype=dt, offset=0, n=os.path.getsize(path) // (2 * np.dtype(dt).itemsize), bits=bits, float=flt, ch=2)


def to_complex(raw, info, swap):
    x = raw.astype(np.float32)
    if info["bits"] == 8:
        x = (x - 127.5) / 128.0
    elif not info["float"]:
        x /= float(2 ** (info["bits"] - 1))
    if info["ch"] == 1:
        return x[:, 0]
    return (x[:, 1] + 1j * x[:, 0]) if swap else (x[:, 0] + 1j * x[:, 1])


def decimate(path, info, swap):
    fs, n = info["rate"], info["n"]
    dec = max(1, int(round(fs / TARGET)))
    x = np.memmap(path, dtype=info["dtype"], mode="r", offset=info["offset"], shape=(n, info["ch"]))
    ch = dec * TARGET            # about 1 s per chunk, a multiple of dec
    marg = dec * 200
    out = []
    for s in range(0, n, ch):
        a, b = max(0, s - marg), min(n, s + ch + marg)
        y = resample_poly(to_complex(np.asarray(x[a:b]), info, swap), 1, dec) if dec > 1 else to_complex(np.asarray(x[a:b]), info, swap)
        i0 = (s - a) // dec
        out.append(y[i0:i0 + (min(s + ch, n) - s) // dec])
    return np.concatenate(out), fs / dec


def to_analytic(y, fsd):
    """Audio: keep the positive frequencies only and move the strongest tone to 0 Hz. Returns the signal and the tone (Hz)."""
    from scipy.signal import hilbert, welch
    f, P = welch(y, fsd, nperseg=8192)
    P[f < 50] = 0
    tone = float(f[np.argmax(P)])
    return (hilbert(y) * np.exp(-2j * np.pi * tone * np.arange(len(y)) / fsd)).astype(np.complex64), tone


def band_shift(y, fsd, band, audio):
    """Zero everything outside LO..HI in one FFT of the whole signal (for audio this also drops the negative
    frequencies, x2 keeps the amplitude) and move the band centre to 0 Hz. Returns the signal and the centre (Hz).
    The cut must come first: the dechirp passes shift the spectrum by up to the swing and would pull outside
    sound (engine rumble) into the search window."""
    f = np.fft.fftfreq(len(y), 1 / fsd)
    z = np.fft.ifft(np.fft.fft(y) * (((f >= band[0]) & (f <= band[1])) * (2.0 if audio else 1.0)))
    c = (band[0] + band[1]) / 2
    return (z * np.exp(-2j * np.pi * c * np.arange(len(y)) / fsd)).astype(np.complex64), c


def audio_to_wav(path, out):
    """Any audio format ffmpeg reads -> a 48 kHz mono 16-bit WAV next to the outputs."""
    import subprocess
    wav = os.path.join(out, os.path.splitext(os.path.basename(path))[0] + "_audio48k.wav")
    subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-i", path, "-ac", "1", "-ar", "48000", "-c:a", "pcm_s16le", wav], check=True)
    return wav


# ---------- helpers ----------
def stft(y, fsd, win, hop, span, nfft):
    nfr = (len(y) - win) // hop + 1
    w = np.hanning(win).astype(np.float32)
    f = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fsd))
    keep = np.abs(f) <= span
    P = np.empty((nfr, keep.sum()), np.float32)
    for a in range(0, nfr, 128):
        idx = np.arange(a, min(nfr, a + 128))
        seg = np.stack([y[i * hop:i * hop + win] for i in idx]) * w
        X = np.fft.fftshift(np.fft.fft(seg, nfft, axis=1), axes=1)[:, keep]
        P[idx] = X.real ** 2 + X.imag ** 2
    return (np.arange(nfr) * hop + win / 2) / fsd, f[keep], P


def remove_dc(y, h):
    c = np.concatenate([[0], np.cumsum(y.astype(np.complex128))])
    i = np.arange(len(y)); a = np.clip(i - h, 0, len(y)); b = np.clip(i + h + 1, 0, len(y))
    return (y - (c[b] - c[a]) / (b - a)).astype(np.complex64)


def fixed_lines(f, P):
    """Receiver spurs: bins whose time-median power stands > 6 dB over the local floor."""
    m = 10 * np.log10(np.median(P, axis=0) + 1e-30)
    base = median_filter(m, size=max(5, int(2000 / (f[1] - f[0]))), mode="nearest")
    hot = np.where(m - base > 6)[0]
    lines, grp = [], []
    for j in hot:
        if grp and j - grp[-1] > 2:
            lines.append(float(np.mean(f[grp]))); grp = []
        grp.append(j)
    if grp:
        lines.append(float(np.mean(f[grp])))
    return np.array(lines)


def near(fr, lines):
    if len(lines) == 0:
        return np.zeros(len(fr), bool)
    return np.min(np.abs(np.asarray(fr)[:, None] - lines[None, :]), axis=1) < SPUR_HALF


def second_diff_sigma(t, f, dt=1.0):
    lut = {round(x * 8): v for x, v in zip(t, f)}
    d2 = [v - 2 * lut[round((x + dt) * 8)] + lut[round((x + 2 * dt) * 8)] for x, v in zip(t, f)
          if round((x + dt) * 8) in lut and round((x + 2 * dt) * 8) in lut]
    d2 = np.array(d2)
    return (1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6) if len(d2) > 5 else np.nan), len(d2)


def robust_spline(t, f, k=CLIP_K, sigma=None, iters=10):
    def fit(a, b):
        if sigma is None:
            return make_smoothing_spline(a, b)
        best = make_smoothing_spline(a, b, lam=1e-8)
        for lam in np.logspace(-8, 6, 57):
            sp = make_smoothing_spline(a, b, lam=lam); r = b - sp(a)
            if 1.4826 * np.median(np.abs(r - np.median(r))) > sigma:
                break
            best = sp
        return best
    ok = np.ones(len(t), bool)
    for _ in range(iters):
        sp = fit(t[ok], f[ok]); r = f - sp(t)
        mad = 1.4826 * np.median(np.abs(r[ok] - np.median(r[ok])))
        new = np.abs(r) < k * max(mad if sigma is None else sigma, 2.0)
        if np.array_equal(new, ok):
            break
        ok = new
    return fit(t[ok], f[ok]), ok


def pass_model(t, a, b, t0, d):
    u = t - t0
    return a - b * u / np.sqrt(d * d + u * u)


def fit_model(t, f, span):
    tm = t[np.argmin(np.abs(f - np.median(f)))]
    p0 = [float(np.clip(np.median(f), -span * 0.9, span * 0.9)), max(abs(np.ptp(f)) / 2, 100.0), tm, 60.0]
    r = least_squares(lambda p: pass_model(t, *p) - f, p0, loss="soft_l1", f_scale=200.0,
                      bounds=([-span, 10, t.min() - 300, 2], [span, 2 * span, t.max() + 300, 2000]))
    return r.x


def peak_interp(row, j):
    a, b, c = np.log(row[j - 1:j + 2] + 1e-30)
    den = a - 2 * b + c
    return 0.5 * (a - c) / den if den != 0 else 0.0


# ---------- tracking ----------
def pass_a(y, fsd, T, span, dc):
    win = int(round(fsd / 4))
    t, f, P = stft(y, fsd, win, win // 2, span, win)                  # 4 Hz grid
    if dc is None:
        lines = np.array([])            # audio with --band: no receiver lines; a source on one plateau is not a line
    else:
        lines = np.append(fixed_lines(f, P), dc)
    R = P / np.median(P, axis=0); R /= np.median(R, axis=1, keepdims=True)
    R[:, near(f, lines)] = 1.0
    df = f[1] - f[0]
    inside = (t > EDGE_S) & (t < T - EDGE_S)
    Rs = uniform_filter1d(R, max(3, int(min(300.0, span / 2) / df)), axis=1)   # narrower in a narrow --band
    ci = np.argmax(Rs, axis=1)
    csnr = 10 * np.log10(Rs[np.arange(len(t)), ci] / np.median(Rs, axis=1))
    cdet = (csnr > COARSE_THR_DB) & inside
    sp, _ = robust_spline(t[cdet], f[ci][cdet], k=3.0)
    ref = sp(t)
    pk = np.full(len(t), np.nan); snr = np.full(len(t), np.nan)
    for i in np.where(inside)[0]:
        j0, j1 = np.searchsorted(f, ref[i] - 600), np.searchsorted(f, ref[i] + 600)
        if j1 - j0 < 3:
            continue
        j = int(np.argmax(R[i, j0:j1])) + j0
        if 0 < j < len(f) - 1:
            pk[i] = f[j] + peak_interp(R[i], j) * df; snr[i] = 10 * np.log10(R[i, j])
    det = snr > A_THR_DB
    assert det.sum() > 20, "signal not found: check --center, --span and --swap-iq"
    return fit_model(t[det], pk[det], span), lines


def dechirp_pass(y, fsd, T, ref, search, span, cut=False):
    ts = np.arange(len(y)) / fsd
    z = (y * np.exp(-2j * np.pi * np.cumsum(ref(ts)) / fsd)).astype(np.complex64)
    win = int(round(fsd))
    t, f, P = stft(z, fsd, win, win // 4, span, win)                 # 1 s windows, 0.25 s hop, 1 Hz grid
    df = f[1] - f[0]
    j0, j1 = np.searchsorted(f, -search), np.searchsorted(f, search)
    fref = ref(t)
    if cut:
        # --band: bins that came from outside the cut band are empty, so the floor is the median of the bins
        # inside it, and they cannot hold the peak or the centroid. No per-bin floor: in a narrow band the
        # tone sits in the same dechirped bins all the time and that floor would subtract the tone itself.
        out = np.abs(f[None, :] + fref[:, None]) > span
        med = np.array([np.median(P[i, ~out[i]]) if (~out[i]).any() else np.inf for i in range(len(t))])
        nbin = np.ones(len(f))
        P[out] = 0.0
    else:
        med = np.median(P, axis=1)
        nbin = np.median(P / med[:, None], axis=0)
    pk = np.full(len(t), np.nan); snr = np.full(len(t), np.nan); lc = np.full(len(t), np.nan)
    for i in range(len(t)):
        j = int(np.argmax(P[i, j0:j1])) + j0
        pk[i] = fref[i] + f[j] + peak_interp(P[i], j) * df
        snr[i] = 10 * np.log10(P[i, j] / med[i])
        e = np.clip(P[i, j0:j1] - med[i] * nbin[j0:j1], 0, None)
        if e.sum() > 0:
            lc[i] = fref[i] + (e * f[j0:j1]).sum() / e.sum()
    inside = (t - 0.5 > EDGE_S) & (t + 0.5 < T - EDGE_S)
    return dict(t=t, pk=pk, snr=snr, lc=lc, inside=inside)


def track(y, fsd, span, dc, cut=False):
    """dc: offset of the receiver DC line (blanked), or None for audio. cut: the signal was cut to +-span (--band)."""
    T = len(y) / fsd
    pA, lines = pass_a(y, fsd, T, span, dc)
    D = dechirp_pass(y, fsd, T, lambda tt: pass_model(np.asarray(tt, float), *pA), SEARCH_B, span, cut)
    det = (D["snr"] > C_THR_DB) & D["inside"] & ~near(D["pk"], lines)
    _, ok = robust_spline(D["t"][det], D["pk"][det])
    tb, fb = D["t"][det][ok], D["pk"][det][ok]
    pB = fit_model(tb, fb, span)
    resB = fb - pass_model(tb, *pB)
    cx, cy = [], []
    for e0 in np.arange(0, T + 10, 10.0):
        m = (tb >= e0) & (tb < e0 + 10)
        if m.sum() >= 3:
            cx.append(tb[m].mean()); cy.append(np.median(resB[m]))
    cx, cy = np.array(cx), np.array(cy)
    refC = lambda tt: pass_model(np.asarray(tt, float), *pB) + np.interp(tt, cx, cy)
    D = dechirp_pass(y, fsd, T, refC, SEARCH_C, span, cut)
    det = (D["snr"] > C_THR_DB) & D["inside"] & ~near(D["pk"], lines)
    _, ok0 = robust_spline(D["t"][det], D["pk"][det])
    sig, _ = second_diff_sigma(D["t"][det][ok0], D["pk"][det][ok0])
    sp, ok = robust_spline(D["t"][det], D["pk"][det], sigma=sig)
    tc = D["t"][det][ok]; c0, c1 = tc[0], tc[-1]
    refD = lambda tt: refC(tt) + sp(np.clip(np.asarray(tt, float), c0, c1)) - refC(np.clip(np.asarray(tt, float), c0, c1))
    D = dechirp_pass(y, fsd, T, refD, SEARCH_D, span, cut)
    det = (D["snr"] > C_THR_DB) & D["inside"] & ~near(D["pk"], lines) & np.isfinite(D["lc"])
    _, ok0 = robust_spline(D["t"][det], D["lc"][det])
    sig, nsig = second_diff_sigma(D["t"][det][ok0], D["lc"][det][ok0])
    _, ok = robust_spline(D["t"][det], D["lc"][det], sigma=sig)
    acc = np.zeros(len(D["t"]), bool); acc[np.where(det)[0][ok]] = True
    return D, acc, sig, nsig, lines


def robust_pass(t, f, T, span):
    """Soft-L1 fit of the pass shape from 15 x 3 starts (t0 across the trace, width T/30, T/10, T/3); the lowest cost wins."""
    best = None
    for t0 in np.linspace(t.min(), t.max(), 15):
        for w in (T / 30, T / 10, T / 3):
            r = least_squares(lambda p: pass_model(t, *p) - f, [float(np.median(f)), max(np.ptp(f) / 2, 100.0), t0, max(w, 1.0)],
                              loss="soft_l1", f_scale=BURST_FSCALE_HZ,
                              bounds=([-span, 10, t.min() - T, 0.5], [span, 2 * span, t.max() + T, 10 * T]))
            if best is None or r.cost < best.cost:
                best = r
    return best.x


def burst_track(y, fsd, span):
    T = len(y) / fsd
    win = int(BURST_FRAME_S * fsd)
    t, f, S = stft(y, fsd, win, win, span, win)
    df = f[1] - f[0]
    S = uniform_filter1d(S, max(1, int(BURST_SMOOTH_HZ / df)), axis=1)
    base = median_filter(S, size=(2 * int(BURST_MEDIAN_S / BURST_FRAME_S) + 1, 1), mode="nearest")
    E = S / (base + 1e-30) - 1.0
    nb = max(1, int(BURST_BAND_HZ / df))
    B = uniform_filter1d(E, nb, axis=1)
    lc = np.full(len(t), np.nan); snr = np.full(len(t), np.nan)
    for i in range(len(t)):
        j = int(np.argmax(B[i]))
        snr[i] = 10 * np.log10(max(B[i, j], 1e-6) + 1.0)
        lo, hi = max(0, j - nb // 2), min(len(f), j + nb // 2 + 1)
        e = np.clip(E[i, lo:hi], 0, None)
        if snr[i] >= BURST_THR_DB and e.sum() > 0:
            lc[i] = (e * f[lo:hi]).sum() / e.sum()
    det = np.isfinite(lc)
    assert det.sum() > 10, "no bursts found: check --center, --span and --band"
    p = robust_pass(t[det], lc[det], T, span)
    acc = det.copy(); acc[det] = np.abs(lc[det] - pass_model(t[det], *p)) < BURST_GATE_HZ
    r = lc[acc] - pass_model(t[acc], *p)
    sig = 1.4826 * float(np.median(np.abs(r - np.median(r))))
    return dict(t=t, lc=lc, snr=snr), acc, sig, int(acc.sum()), np.array([])


def waterfall(y, fsd, span, t, f, acc, png, title):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    win = int(round(fsd / 10))
    tt, ff, P = stft(y, fsd, win, win * 5, span, win)                # 0.1 s FFT, 10 Hz bins, one row every 0.5 s
    R = P / np.median(P, axis=0); R /= np.median(R, axis=1, keepdims=True)
    fig, ax = plt.subplots(figsize=(8, 10))
    np.savez_compressed(png.replace(".png", ".npz"), t=tt.astype(np.float32), f=ff.astype(np.float32), db=(10 * np.log10(R)).astype(np.float16))
    ax.imshow(10 * np.log10(R), aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=15,
              extent=[ff[0] / 1e3, ff[-1] / 1e3, tt[0], tt[-1]])
    ax.plot(f[acc] / 1e3, t[acc], ",", color="red", alpha=0.6)
    ax.set_xlabel("offset from centre (kHz)"); ax.set_ylabel("time from start of recording (s)"); ax.set_title(title, fontsize=10)
    fig.tight_layout(); fig.savefig(png, dpi=110); plt.close(fig)


def run(path, center=None, start=None, span=15000.0, swap_iq=False, out=None, rate=None, fmt="cs16", log=print, band=None, burst=False):
    """Trace one recording. band = (lo, hi) Hz or None; burst = short-burst mode.
    Returns the summary dict, with the output file paths under "files"."""
    name, orig = os.path.splitext(os.path.basename(path))[0], os.path.basename(path)
    out = out or os.path.dirname(os.path.abspath(path)); os.makedirs(out, exist_ok=True)
    fc = center
    if fc is None:
        m = re.search(r"(\d{6,})Hz", name)
        fc = float(m.group(1)) if m else None
    t_start = datetime.fromisoformat(start.replace("Z", "+00:00")) if start else None
    ext = os.path.splitext(path)[1].lower()
    tmp = None
    if ext in (".ogg", ".oga", ".mp3", ".flac", ".m4a", ".opus", ".aac"):
        path = tmp = audio_to_wav(path, out); ext = ".wav"
    if ext == ".wav":
        info = wav_info(path)
    else:
        if not rate and re.fullmatch(r"iq_\d+", name):
            rate = 48000.0                                          # SatNOGS IQ files: 16-bit IQ at 48,000 S/s
            log("SatNOGS IQ file: 48,000 S/s, 16-bit assumed")
        assert rate, "give the sample rate for a raw IQ file"
        info = raw_info(path, rate, fmt)
    audio = info["ch"] == 1
    if fc is None:
        assert audio, "give the centre frequency (Hz)"
        fc = 0.0
    log(f"{name}: {info['rate']} S/s, {info['bits']}-bit, {'audio' if audio else 'IQ'}, {info['n'] / info['rate']:.1f} s, centre {fc / 1e6:.6f} MHz")
    y, fsd = decimate(path, info, swap_iq)
    if tmp:
        os.remove(tmp)
    tone, dc = 0.0, None
    if band:
        band = (float(min(band)), float(max(band)))
        assert band[1] - band[0] >= 50, "the band must be at least 50 Hz wide"
        assert (0 <= band[0] if audio else -fsd / 2 < band[0]) and band[1] < fsd / 2, f"the band must lie inside the recording (+-{fsd / 2:.0f} Hz)"
    if audio and band:
        y, tone = band_shift(y, fsd, band, True)
        span = (band[1] - band[0]) / 2
        log(f"  audio band {band[0]:.0f} to {band[1]:.0f} Hz")
    elif audio:
        y, tone = to_analytic(y, fsd)
        y, dc, span = remove_dc(y, int(5 * fsd)), 0.0, min(span, fsd / 2 - 100)
        log(f"  strongest audio tone {tone:.0f} Hz")
    else:
        y, dc = remove_dc(y.astype(np.complex64), int(5 * fsd)), 0.0
        if band:
            y, tone = band_shift(y, fsd, band, False)
            span, dc = (band[1] - band[0]) / 2, -tone
            log(f"  offsets {band[0]:.0f} to {band[1]:.0f} Hz from the centre")
    fc += tone
    if burst:
        log("  burst mode")
        D, acc, sig, nsig, lines = burst_track(y, fsd, span)
    else:
        D, acc, sig, nsig, lines = track(y, fsd, span, dc, cut=bool(band))
    t, f = D["t"], D["lc"]
    files = dict(trace=os.path.join(out, f"{name}_trace.csv"), waterfall=os.path.join(out, f"{name}_waterfall.png"),
                 summary=os.path.join(out, f"{name}_summary.json"))
    with open(files["trace"], "w") as fh:
        fh.write("t_s,utc,f_offset_Hz,f_Hz,snr_dB\n")
        for i in np.where(acc)[0]:
            utc = (t_start + timedelta(seconds=float(t[i]))).isoformat(timespec="milliseconds") if t_start else ""
            fh.write(f"{t[i]:.3f},{utc},{f[i]:.2f},{fc + f[i]:.2f},{D['snr'][i]:.1f}\n")
    waterfall(y, fsd, span, t, f, acc, files["waterfall"], f"{name}\ncentre {fc / 1e6:.6f} MHz")
    swing = float(np.ptp(f[acc])) if acc.any() else 0.0
    summ = dict(file=orig, centre_Hz=fc, audio_tone_Hz=tone if audio else None, band_Hz=list(band) if band else None,
                burst=bool(burst), sample_rate=info["rate"],
                decimated_rate=fsd, duration_s=info["n"] / info["rate"], start_utc=start, points=int(acc.sum()),
                point_noise_Hz=round(float(sig), 2), point_noise_pairs=nsig, swing_Hz=round(swing, 1),
                fixed_receiver_lines_Hz=[round(float(v), 1) for v in lines if v != dc],
                window="0.5 s Hann frames, no overlap" if burst else "1.0 s Hann, 0.25 s hop (only every 4th point is independent)")
    json.dump(summ, open(files["summary"], "w"), indent=1)
    log(f"  {acc.sum()} points, point noise {sig:.1f} Hz, swing {swing:.0f} Hz -> {out}")
    if not audio and swing < 0.2 * span:
        summ["note"] = ("The signal moves only %.0f Hz over the recording. A satellite pass at a fixed tuning moves several kHz. "
                        "This recording looks Doppler-corrected (the receiver followed a predicted orbit), so it holds only the leftover. "
                        "The fit needs a recording made at a fixed frequency with Doppler correction off." % swing)
        log("  NOTE: " + summ["note"])
    if not audio and not burst:
        cover = acc.sum() / max(1, D["inside"].sum())
        if swing < 0.2 * span or cover < 0.3:
            summ["burst_hint"] = ("The trace covers %.0f%% of the recording and moves %.0f Hz. If the satellite sends short bursts "
                                  "(FSK or packets) under steady lines, run again with --burst." % (100 * cover, swing))
            log("  NOTE: " + summ["burst_hint"])
    summ["files"] = files
    return summ


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("wav"); ap.add_argument("--center", type=float, help="centre (IQ) or dial (audio) frequency, Hz (default: read from an SDR++ file name; 0 for audio)")
    ap.add_argument("--start", help="UTC time of the first sample, e.g. 2026-06-14T18:10:09Z (optional)")
    ap.add_argument("--span", type=float, default=15000, help="search +- this many Hz around the centre (default 15000)")
    ap.add_argument("--swap-iq", action="store_true", help="the recording stores Q before I")
    ap.add_argument("--out", help="output folder (default: next to the recording)")
    ap.add_argument("--rate", type=float, help="sample rate, for a raw IQ file")
    ap.add_argument("--format", default="cs16", choices=["cs16", "cf32", "cu8"], help="raw IQ sample type (default cs16)")
    ap.add_argument("--band", type=float, nargs=2, metavar=("LO", "HI"),
                    help="audio: pick the tone between LO and HI Hz; IQ: search offsets LO to HI Hz from the centre")
    ap.add_argument("--burst", action="store_true", help="the signal comes in short bursts (FSK/packet satellites)")
    a = ap.parse_args()
    run(a.wav, a.center, a.start, a.span, a.swap_iq, a.out, a.rate, a.format, band=a.band, burst=a.burst)


if __name__ == "__main__":
    main()
