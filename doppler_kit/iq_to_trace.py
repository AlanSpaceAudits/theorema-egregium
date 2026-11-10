"""Received frequency against time from an IQ or audio recording.

Accepts a WAV file (2 channels I/Q, or 1 channel audio; 8/16-bit int or 32-bit float), a headerless raw IQ file
(.raw/.cs16/.cf32/.cu8/.sigmf-data: give --rate, and --format if not 16-bit signed), or any audio file ffmpeg reads
(.ogg/.mp3/.flac/.m4a/.opus). For audio, --center is the receiver dial frequency (0 = the audio frequency itself,
for sound recordings), and the strongest tone is tracked.

Usage:
  python3 iq_to_trace.py RECORDING.wav [--center HZ] [--start 2026-06-01T12:00:00Z] [--span HZ] [--swap-iq] [--out DIR]
  python3 iq_to_trace.py RECORDING.raw --rate 48000 --center 437800000 [--format cs16|cf32|cu8]
  python3 iq_to_trace.py RECORDING.wav --band LO HI     (audio: tone between LO and HI Hz; IQ: offsets LO..HI from the centre)
  python3 iq_to_trace.py RECORDING.cf32 --rate 48000 --center 436650000 --format cf32 --burst   (short FSK/packet bursts)

Blank values are worked out from the file; a value given always wins, and every auto value is logged with its source:
  rate, format, centre, start: SigMF meta file (FILE.sigmf-meta beside FILE.sigmf-data), WAV header, an SDR++ file name
  (...437800000Hz...), a SatNOGS iq_<id>.raw (48,000 S/s cs16; centre and start from network.satnogs.org), and for the
  start the file's modified time minus the recording length (check it).
  span (IQ): a scan of the whole recorded band finds the falling line (steady lines removed) and sets the span and,
  when needed, moves the centre. band (audio): the scan picks the strongest tone that moves.
  burst: when the trace looks like bursts (covers < 30%, moves < span/5, or does not fit a pass cleanly), burst mode
  runs too and the clean pass that covers more of the recording is kept.

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
--fm (FM carriers such as ISS SSTV or voice): 1.0 s Hann frames every 0.5 s; frame noise = median bin power; the
  strongest 2000 Hz band (boxcar) must be >= 15 dB over the noise; value = centroid of max(P - 5 noise, 0) within
  +-10000 Hz, the window started on that band and moved to the centroid 3 times (the mean of an FM signal's frequency
  is its carrier; the strongest band is often a sideband). Outliers as in --burst.
"""
import argparse, json, math, os, re, struct
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
FM_FRAME_S, FM_HOP_S, FM_BAND_HZ, FM_HALF_HZ, FM_NOISE_K, FM_THR_DB, FM_ITER = 1.0, 0.5, 2000.0, 10000.0, 5.0, 15.0, 3


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


def raw_info(path, rate, fmt="cs16"):
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


def decimate(path, info, swap, shift=0, target=TARGET):
    """Mix the recording down by `shift` Hz (an integer), then decimate to about `target` S/s."""
    fs, n = info["rate"], info["n"]
    dec = max(1, int(round(fs / target)))
    x = np.memmap(path, dtype=info["dtype"], mode="r", offset=info["offset"], shape=(n, info["ch"]))
    ch = dec * TARGET            # about 1 s per chunk, a multiple of dec
    marg = dec * 200
    out = []
    for s in range(0, n, ch):
        a, b = max(0, s - marg), min(n, s + ch + marg)
        z = to_complex(np.asarray(x[a:b]), info, swap)
        if shift:
            ph = (np.arange(a, b, dtype=np.int64) * int(shift)) % int(fs)        # exact phase index, no float drift
            z = (z * np.exp(-2j * np.pi * ph / fs)).astype(np.complex64)
        y = resample_poly(z, 1, dec) if dec > 1 else z
        i0 = (s - a) // dec
        out.append(y[i0:i0 + (min(s + ch, n) - s) // dec])
    return np.concatenate(out), fs / dec


def strongest_tone(y, fsd):
    """Audio: the strongest tone above 50 Hz in the average spectrum (Welch, 8192 samples)."""
    from scipy.signal import welch
    f, P = welch(y, fsd, nperseg=8192)
    P[f < 50] = 0
    return float(f[np.argmax(P)])


def to_analytic(y, fsd, tone=None):
    """Audio: keep the positive frequencies only and move the strongest tone to 0 Hz. Returns the signal and the tone (Hz)."""
    from scipy.signal import hilbert
    tone = strongest_tone(y, fsd) if tone is None else tone
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


# ---------- auto values: metadata ----------
SIGMF_FORMATS = {"ci16_le": "cs16", "ci16": "cs16", "cf32_le": "cf32", "cf32": "cf32", "cu8": "cu8", "cu8_le": "cu8"}


def utc_text(ts):
    """Unix time (s) -> 2025-01-02T03:04:05Z, rounded to the nearest second."""
    return datetime.fromtimestamp(math.floor(ts + 0.5), timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def iso_utc(s):
    """An ISO time from a metadata file, as a start time the kit reads: fraction cut to 6 digits, no zone = UTC."""
    s = s.strip()
    m = re.fullmatch(r"([0-9-]+T[0-9:]+)(?:[.,]([0-9]+))?(Z|z|[+-][0-9:]+)?", s)
    if not m:
        return None
    zone = m.group(3) or "Z"
    return m.group(1) + ("." + m.group(2)[:6] if m.group(2) else "") + ("Z" if zone in "Zz" else zone)


def sigmf_files(path):
    """(data file, meta file or None) for a SigMF pair; any other file: (path, None)."""
    if path.endswith(".sigmf-meta"):
        return path[:-len("meta")] + "data", path
    if path.endswith(".sigmf-data") and os.path.exists(path[:-len("data")] + "meta"):
        return path, path[:-len("data")] + "meta"
    return path, None


def sigmf_values(meta):
    """rate, fmt, centre and start from a SigMF meta file (only the fields present)."""
    d = json.load(open(meta, encoding="utf-8"))
    g, c = d.get("global", {}), (d.get("captures") or [{}])[0]
    out = {}
    if isinstance(g.get("core:sample_rate"), (int, float)):
        out["rate"] = float(g["core:sample_rate"])
    dt = g.get("core:datatype")
    if isinstance(dt, str):
        assert dt.lower() in SIGMF_FORMATS, f"SigMF datatype {dt} is not supported (ci16_le, cf32_le or cu8)"
        out["fmt"] = SIGMF_FORMATS[dt.lower()]
    if isinstance(c.get("core:frequency"), (int, float)):
        out["centre"] = float(c["core:frequency"])
    if isinstance(c.get("core:datetime"), str) and iso_utc(c["core:datetime"]):
        out["start"] = iso_utc(c["core:datetime"])
    return out


def satnogs_values(obs, log):
    """Start and centre of SatNOGS observation `obs` from network.satnogs.org (None when offline or on failure)."""
    import urllib.request
    url = f"https://network.satnogs.org/api/observations/{obs}/"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "doppler-kit"})
        d = json.load(urllib.request.urlopen(req, timeout=10))
    except Exception as e:
        log(f"auto: SatNOGS observation {obs} not read ({type(e).__name__}): give the centre and start time")
        return {}
    out = {}
    fc, src = d.get("center_frequency"), "centre frequency"
    if fc is None:
        fc, src = d.get("transmitter_downlink_low"), "transmitter downlink"
    if isinstance(fc, (int, float)):
        out["centre"] = (float(fc), src)
    if isinstance(d.get("start"), str) and iso_utc(d["start"]):
        out["start"] = iso_utc(d["start"])
    return out


# ---------- auto values: signal scan ----------
SCAN_FRAMES, SCAN_SAMPLES, SCAN_SLEW_HZ_S = 400, 1 << 26, 500.0
SCAN_ON_DB, SCAN_CAP_DB = 3.0, 20.0
DEFAULT_SPAN = 15000.0


def scan_spectra(path, info, swap, audio):
    """Coarse spectra of the full recorded band, streamed: frames of nblk Hann blocks (power averaged) at M times
    spread evenly over the file. IQ: bins about fs/64 wide, frames up to 1 s; audio: about fs/16, up to 0.25 s."""
    fs, n = info["rate"], info["n"]
    nfft = 1 << int(round(np.log2(fs / (16.0 if audio else 64.0))))
    nblk = max(1, min(16, int(fs * (0.25 if audio else 1.0)) // nfft))
    L = nblk * nfft
    if nfft < 16 or n < 10 * L:
        return None
    M = int(min(n // L, max(SCAN_FRAMES, SCAN_SAMPLES // L)))
    x = np.memmap(path, dtype=info["dtype"], mode="r", offset=info["offset"], shape=(n, info["ch"]))
    w = np.hanning(nfft)
    P = np.empty((M, nfft)); t = np.empty(M)
    for k in range(M):
        s = k * (n - L) // (M - 1)
        z = to_complex(np.asarray(x[s:s + L]), info, swap).astype(np.complex128).reshape(nblk, nfft) * w
        X = np.fft.fft(z, axis=1)
        P[k] = np.fft.fftshift((X.real ** 2 + X.imag ** 2).mean(axis=0))
        t[k] = (s + L / 2) / fs
    return t, np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fs)), P


def strongest_path(v, D):
    """Column index per row of the path with the largest sum of v that never moves up and moves down at most
    D columns from one row to the next (the received frequency of a passing source only falls)."""
    M, N = v.shape
    S = v[0].copy(); back = np.zeros((M, N), np.int32); idx = np.arange(N)
    for k in range(1, M):
        best = np.full(N, -np.inf); arg = idx.copy()
        for d in range(D + 1):                  # predecessor j + d; the smallest d wins a tie
            cand = np.full(N, -np.inf); cand[:N - d] = S[d:]
            better = cand > best
            best[better] = cand[better]; arg[better] = idx[better] + d
        S = v[k] + best; back[k] = arg
    j = np.empty(M, np.int64); j[-1] = int(np.argmax(S))
    for k in range(M - 1, 0, -1):
        j[k - 1] = back[k, j[k]]
    return j


def moving_line(t, f, v, keep, min_move, mask_half):
    """The strongest falling path over the bins in `keep`; a path that moves less than min_move (5th to 95th percentile
    of its strong frames) is masked (+-mask_half around its median) and the search runs again, 3 times at most.
    Returns (lo, hi, strong frames) or None."""
    cols = np.where(keep)[0]
    if len(cols) < 3:
        return None
    c0, c1 = cols[0], cols[-1] + 1
    vv = np.where(keep[None, c0:c1], v[:, c0:c1], 0.0); ff = f[c0:c1]
    D = max(1, int(np.ceil(SCAN_SLEW_HZ_S * (t[1] - t[0]) / (f[1] - f[0]))))
    for _ in range(3):
        j = strongest_path(vv, D)
        on = vv[np.arange(len(t)), j] > SCAN_ON_DB
        if on.sum() < 8:
            return None
        fo = ff[j][on]
        lo, hi = np.percentile(fo, [5, 95])
        if hi - lo >= min_move:
            return float(lo), float(hi), int(on.sum())
        vv[:, np.abs(ff - np.median(fo)) < mask_half] = 0.0
    return None


def scan_iq(t, f, P, fs):
    """IQ: every bin over its median in time (steady lines go to 1), every frame over its median, 300 Hz boxcar.
    First the default window (+-21.6 kHz), then the full band. Returns (lo, hi, frames) in Hz from the centre, or None."""
    df = f[1] - f[0]
    R = P / np.median(P, axis=0); R /= np.median(R, axis=1, keepdims=True)
    Rs = uniform_filter1d(R, max(1, int(round(300.0 / df))), axis=1)
    v = np.clip(10 * np.log10(Rs) - SCAN_ON_DB, 0.0, SCAN_CAP_DB)
    near_lim = min(0.45 * fs, 0.45 * TARGET)
    for lim in [near_lim] + ([0.45 * fs] if 0.45 * fs > near_lim else []):
        r = moving_line(t, f, v, np.abs(f) <= lim, max(10 * df, 1000.0), max(10 * df, 1000.0))
        if r:
            return r
    return None


def scan_audio(t, f, P, fs):
    """Audio: positive frequencies 50 Hz to 0.45 fs, each spectrum over its running median across 1 kHz (tones stand
    out, broadband sound does not), every bin over its median in time (steady tones go to 1), 50 Hz boxcar."""
    from scipy.ndimage import median_filter
    df = f[1] - f[0]
    keep = (f >= 50.0) & (f <= 0.45 * fs)
    f, P = f[keep], P[:, keep]
    if len(f) < 3:
        return None
    C = P / median_filter(P, size=(1, int(round(1000.0 / df)) // 2 * 2 + 1), mode="nearest")
    R = C / np.median(C, axis=0)
    Rs = uniform_filter1d(R, max(1, int(round(50.0 / df))), axis=1)
    v = np.clip(10 * np.log10(Rs) - SCAN_ON_DB, 0.0, SCAN_CAP_DB)
    return moving_line(t, f, v, np.ones(len(f), bool), max(5 * df, 100.0), max(5 * df, 200.0))


def span_from_scan(lo, hi):
    """(centre shift, search span) covering lo..hi (Hz from the centre) plus a margin; the default span when it fits."""
    m = max(2000.0, 0.1 * (hi - lo))
    a, b = lo - m, hi + m
    if a >= -DEFAULT_SPAN and b <= DEFAULT_SPAN:
        return 0, DEFAULT_SPAN
    shift = int(1000 * math.floor((a + b) / 2000 + 0.5))
    return shift, max(DEFAULT_SPAN, 500.0 * math.ceil(max(b - shift, shift - a) / 500))


def band_from_scan(lo, hi, fs):
    m = max(50.0, 0.25 * (hi - lo))
    return max(0.0, 10.0 * math.floor((lo - m) / 10)), min(10.0 * math.ceil((hi + m) / 10), float(math.floor(0.45 * fs)))


def pass_quality(t, f, acc, T, span, sig):
    """Accepted points of a trace against the pass shape: clean when >= 20 points, the swing is at least a fifth of
    the span (not a steady line) and the soft-L1 pass fit leaves at most 3 x the point noise (rms)."""
    if acc.sum() < 20:
        return dict(clean=False, why="fewer than 20 points", cover=0.0)
    tt, ff = t[acc], f[acc]
    swing = float(np.ptp(ff)); cover = float((tt.max() - tt.min()) / T)
    if swing < 0.2 * span:
        return dict(clean=False, why=f"moves only {swing:.0f} Hz", cover=cover)
    p = robust_pass(tt, ff, T, span)
    rms = float(np.sqrt(np.mean((ff - pass_model(tt, *p)) ** 2)))
    noise = max(float(sig), 1.0) if np.isfinite(sig) else 1.0
    if rms > 3 * noise:
        return dict(clean=False, why=f"pass fit leaves {rms:.0f} Hz rms against point noise {noise:.0f} Hz", cover=cover)
    return dict(clean=True, why=f"covers {100 * cover:.0f}% of the recording", cover=cover)


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


def fm_track(y, fsd, span):
    T = len(y) / fsd
    win = int(FM_FRAME_S * fsd)
    t, f, S = stft(y, fsd, win, int(FM_HOP_S * fsd), span, win)
    df = f[1] - f[0]
    noise = np.median(S, axis=1)
    B = uniform_filter1d(S, max(1, int(FM_BAND_HZ / df)), axis=1)
    h = int(FM_HALF_HZ / df)
    lc = np.full(len(t), np.nan); snr = np.full(len(t), np.nan)
    for i in range(len(t)):
        j = int(np.argmax(B[i]))
        snr[i] = 10 * np.log10(B[i, j] / (noise[i] + 1e-30))
        if snr[i] < FM_THR_DB:
            continue
        w = np.clip(S[i] - FM_NOISE_K * noise[i], 0, None).astype(np.float64)
        for _ in range(FM_ITER):
            lo, hi = max(0, j - h), min(len(f), j + h + 1)
            if w[lo:hi].sum() <= 0:
                break
            c = (w[lo:hi] * f[lo:hi]).sum() / w[lo:hi].sum()
            j = int(np.argmin(np.abs(f - c)))
        else:
            lc[i] = c
    det = np.isfinite(lc)
    assert det.sum() > 10, "no FM carrier found: check --center, --span and --band"
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


def run(path, center=None, start=None, span=None, swap_iq=False, out=None, rate=None, fmt=None, log=print, band=None, burst=False, fm=False):
    """Trace one recording. band = (lo, hi) Hz or None; burst = short-burst mode (on), else tried when the trace
    looks like bursts. Blank values (None) are worked out from the file; a value given here always wins.
    Returns the summary dict, with the output file paths under "files" and the auto values under "auto"."""
    auto = {}

    def note(key, text, indent=""):
        auto[key] = text
        log(f"{indent}auto: {text}")

    path, meta = sigmf_files(path)
    name, orig = os.path.splitext(os.path.basename(path))[0], os.path.basename(path)
    out = out or os.path.dirname(os.path.abspath(path)); os.makedirs(out, exist_ok=True)
    fc = center
    sm = sigmf_values(meta) if meta else {}
    if meta:
        mname = os.path.basename(meta)
        if rate is None and "rate" in sm:
            rate = sm["rate"]; note("rate", f"sample rate {rate:,.0f} S/s (SigMF {mname})")
        if fmt is None and "fmt" in sm:
            fmt = sm["fmt"]; note("format", f"format {fmt} (SigMF {mname})")
        if fc is None and "centre" in sm:
            fc = sm["centre"]; note("centre", f"centre {fc:,.0f} Hz (SigMF {mname})")
        if start is None and "start" in sm:
            start = sm["start"]; note("start", f"start {start} (SigMF {mname})")
    if fc is None:
        m = re.search(r"(\d{6,})Hz", name)
        if m:
            fc = float(m.group(1)); note("centre", f"centre {fc:,.0f} Hz (from file name)")
    ext = os.path.splitext(path)[1].lower()
    satnogs = re.fullmatch(r"iq_(\d+)", name)
    if ext in (".wav", ".ogg", ".oga", ".mp3", ".flac", ".m4a", ".opus", ".aac"):
        satnogs = None
    if satnogs and (fc is None or start is None):
        sv = satnogs_values(satnogs.group(1), log)
        if fc is None and "centre" in sv:
            fc = sv["centre"][0]; note("centre", f"centre {fc:,.0f} Hz (SatNOGS observation {satnogs.group(1)}, {sv['centre'][1]})")
        if start is None and "start" in sv:
            start = sv["start"]; note("start", f"start {start} (SatNOGS observation {satnogs.group(1)})")
    t_start = datetime.fromisoformat(start.replace("Z", "+00:00")) if start else None
    mtime = os.path.getmtime(path)
    tmp = None
    if ext in (".ogg", ".oga", ".mp3", ".flac", ".m4a", ".opus", ".aac"):
        path = tmp = audio_to_wav(path, out); ext = ".wav"
    if ext == ".wav":
        info = wav_info(path)
    else:
        if not rate and satnogs:
            rate = 48000.0                                          # SatNOGS IQ files: 16-bit IQ at 48,000 S/s
            fmt = fmt or "cs16"
            note("rate", "sample rate 48,000 S/s, cs16 (SatNOGS IQ file)")
        assert rate, "give the sample rate for a raw IQ file"
        info = raw_info(path, rate, fmt or "cs16")
    audio = info["ch"] == 1
    T = info["n"] / info["rate"]
    if start is None:
        start = utc_text(mtime - T); t_start = datetime.fromisoformat(start.replace("Z", "+00:00"))
        note("start", f"start {start} (file modified time minus length: check it)")
    if fc is None:
        assert audio, "give the centre frequency (Hz): it is not in the file name, a SigMF file or a SatNOGS observation"
        fc = 0.0
    log(f"{name}: {info['rate']} S/s, {info['bits']}-bit, {'audio' if audio else 'IQ'}, {T:.1f} s, centre {fc / 1e6:.6f} MHz")
    if band:
        band = (float(min(band)), float(max(band)))
    shift, target, auto_band, scanned = 0, TARGET, None, False
    if not audio and band is None and span is None:
        sc = scan_spectra(path, info, swap_iq, False)
        r = scan_iq(*sc, info["rate"]) if sc else None
        if r:
            shift, span = span_from_scan(r[0], r[1])
            if shift:
                note("span", f"signal found at {r[0] / 1e3:+.1f} to {r[1] / 1e3:+.1f} kHz: centre moved by {shift:+,} Hz, "
                             f"search \u00b1{span:,.0f} Hz (scan)", "  ")
            else:
                note("span", f"signal found at {r[0] / 1e3:+.1f} to {r[1] / 1e3:+.1f} kHz, search \u00b1{span:,.0f} Hz (scan)", "  ")
            target, scanned = max(TARGET, 2.5 * span), True
        else:
            span = DEFAULT_SPAN
            note("span", f"no moving signal found: search \u00b1{span:,.0f} Hz (scan)", "  ")
    elif audio and band is None:
        sc = scan_spectra(path, info, False, True)
        r = scan_audio(*sc, info["rate"]) if sc else None
        if r:
            auto_band = band_from_scan(r[0], r[1], info["rate"])
        else:
            note("band", "no moving tone found: the strongest tone is tracked (scan)", "  ")
    span = DEFAULT_SPAN if span is None else span
    y, fsd = decimate(path, info, swap_iq, shift, target)
    if tmp:
        os.remove(tmp)
    if scanned:
        span = min(span, 0.45 * fsd)
    tone, dc = 0.0, None
    if auto_band:
        strong = strongest_tone(y, fsd)
        if auto_band[0] <= strong <= auto_band[1]:
            note("band", f"the strongest tone {strong:.0f} Hz moves: no band needed (scan)", "  ")
        else:
            band = auto_band
            note("band", f"band {band[0]:.0f} to {band[1]:.0f} Hz around the moving tone; the strongest tone {strong:.0f} Hz is steady (scan)", "  ")
    if band:
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
        y, dc = remove_dc(y.astype(np.complex64), int(5 * fsd)), (-float(shift) if shift else 0.0)
        if band:
            y, tone = band_shift(y, fsd, band, False)
            span, dc = (band[1] - band[0]) / 2, -tone
            log(f"  offsets {band[0]:.0f} to {band[1]:.0f} Hz from the centre")
    fc += tone + shift
    hint = None
    if fm:
        log("  FM mode")
        D, acc, sig, nsig, lines = fm_track(y, fsd, span)
        burst = False
    elif burst:
        log("  burst mode")
        D, acc, sig, nsig, lines = burst_track(y, fsd, span)
    else:
        try:
            D, acc, sig, nsig, lines = track(y, fsd, span, dc, cut=bool(band))
            normal_err = None
        except Exception as e:                  # the normal tracer can fail in many ways on a burst signal
            if audio:
                raise
            normal_err = e
        qn = None
        if normal_err is None and not audio:
            swing0 = float(np.ptp(D["lc"][acc])) if acc.any() else 0.0
            cover0 = acc.sum() / max(1, D["inside"].sum())
            hint = swing0 < 0.2 * span or cover0 < 0.3
            qn = pass_quality(D["t"], D["lc"], acc, T, span, sig)
        if normal_err is not None or hint or (qn and not qn["clean"]):
            try:
                B = burst_track(y, fsd, span)
            except Exception:
                B = None
            qb = pass_quality(B[0]["t"], B[0]["lc"], B[1], T, span, B[2]) if B else dict(clean=False, why="no bursts found")
            # FM is tried when the normal trace is not a clean pass. Against a clean burst trace, FM needs a lower point
            # noise and at least as many points.
            normal_clean = bool(qn and qn["clean"])
            F = None
            if not normal_clean:
                try:
                    F = fm_track(y, fsd, span)
                except Exception:
                    F = None
            qf = pass_quality(F[0]["t"], F[0]["lc"], F[1], T, span, F[2]) if F else dict(clean=False, why="no FM carrier found")
            fm_wins = qf["clean"] and (not (B and qb["clean"]) or (F[2] < B[2] and F[3] >= B[3]))
            if not normal_clean and not fm_wins:
                why = "burst trace has lower point noise or more points" if qf["clean"] else f"FM trace {qf['why']}"
                note("fm", f"FM mode tried and not used ({why})", "  ")
            if fm_wins:
                D, acc, sig, nsig, lines = F
                fm = True
                why = f"normal trace {qn['why']}" if qn else "the normal trace found no signal"
                vs = f"; point noise {F[2]:.0f} Hz over {F[3]} points, burst mode {B[2]:.0f} Hz over {B[3]}" if B and qb["clean"] else ""
                note("fm", f"FM mode on ({why}{vs})", "  ")
            elif normal_err is not None:
                if not qb["clean"]:
                    raise normal_err
                D, acc, sig, nsig, lines = B
                burst = True
                note("burst", "burst mode on (the normal trace found no signal)", "  ")
            elif qb["clean"] and (not qn["clean"] or qb["cover"] > qn["cover"]):
                D, acc, sig, nsig, lines = B
                burst = True
                why = f"normal trace {qn['why']}" if not qn["clean"] else f"burst trace {qb['why']}, normal trace {100 * qn['cover']:.0f}%"
                note("burst", f"burst mode on (steady lines stronger than the signal: {why})", "  ")
            else:
                note("burst", f"burst mode tried and not used (burst trace {qb['why']})", "  ")
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
                burst=bool(burst), fm=bool(fm), sample_rate=info["rate"],
                decimated_rate=fsd, duration_s=T, start_utc=start, points=int(acc.sum()),
                point_noise_Hz=round(float(sig), 2), point_noise_pairs=nsig, swing_Hz=round(swing, 1),
                fixed_receiver_lines_Hz=[round(float(v), 1) for v in lines if v != dc],
                window="1.0 s Hann, 0.5 s hop (only every 2nd point is independent)" if fm else "0.5 s Hann frames, no overlap" if burst else "1.0 s Hann, 0.25 s hop (only every 4th point is independent)",
                auto=auto)
    json.dump(summ, open(files["summary"], "w"), indent=1)
    log(f"  {acc.sum()} points, point noise {sig:.1f} Hz, swing {swing:.0f} Hz -> {out}")
    if not audio and swing < 0.2 * span:
        summ["note"] = ("The signal moves only %.0f Hz over the recording. A satellite pass at a fixed tuning moves several kHz. "
                        "This recording looks Doppler-corrected (the receiver followed a predicted orbit), so it holds only the leftover. "
                        "The fit needs a recording made at a fixed frequency with Doppler correction off." % swing)
        log("  NOTE: " + summ["note"])
    summ["files"] = files
    return summ


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("wav"); ap.add_argument("--center", type=float, help="centre (IQ) or dial (audio) frequency, Hz (default: read from an SDR++ file name; 0 for audio)")
    ap.add_argument("--start", help="UTC time of the first sample, e.g. 2026-06-01T12:00:00Z (optional)")
    ap.add_argument("--span", type=float, help="search +- this many Hz around the centre (default: found by a scan of the whole band)")
    ap.add_argument("--swap-iq", action="store_true", help="the recording stores Q before I")
    ap.add_argument("--out", help="output folder (default: next to the recording)")
    ap.add_argument("--rate", type=float, help="sample rate, for a raw IQ file")
    ap.add_argument("--format", choices=["cs16", "cf32", "cu8"], help="raw IQ sample type (default: from SigMF, else cs16)")
    ap.add_argument("--band", type=float, nargs=2, metavar=("LO", "HI"),
                    help="audio: pick the tone between LO and HI Hz; IQ: search offsets LO to HI Hz from the centre")
    ap.add_argument("--burst", action="store_true", help="the signal comes in short bursts (FSK/packet satellites); without it, tried when the trace looks like bursts")
    ap.add_argument("--fm", action="store_true", help="the signal is an FM carrier (ISS SSTV or voice): trace the centre of its energy; without it, tried when the normal trace is not a clean pass")
    a = ap.parse_args()
    run(a.wav, a.center, a.start, a.span, a.swap_iq, a.out, a.rate, a.format, band=a.band, burst=a.burst, fm=a.fm)


if __name__ == "__main__":
    main()
