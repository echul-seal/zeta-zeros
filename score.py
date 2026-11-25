"""zeta-zeros: score code.

The piece sonifies the explicit formula. For the nontrivial zeros rho = 1/2 + i*gamma,
    sum_gamma cos(gamma * u)  =  (1/2) e^{u/2}  -  (1/2) e^{-u/2} sum_n Lambda(n) delta(u - log n)
(in the sense of distributions), so a sum of cosines at the zero heights, read along
u = log x, piles up into spikes at every prime power.

Structure:
- The piece is a series of passes ("bars"). Each pass scans u = log x from log x0 up to
  log x1, so prime powers arrive as an accelerating run.
- Zeros enter one at a time, on a log clock in their unfolded height. Zero 1 is alone
  for several seconds; by the end thousands pour in per bar.
- x1 for each pass grows with the resolution the present zeros allow (p ~ gamma_max/pi),
  so every bar reaches as far into the primes as the zeros can see.
- The complex sum Z(u) = sum a_k e^{i gamma_k u} is put onto a carrier whose pitch is
  tied to u (f = f0 * x^alpha), single sideband. Every prime power therefore sounds at its
  own fixed pitch, bar after bar. Early on Z is incoherent and each bar is a rising wash;
  as zeros accumulate, the wash collapses into struck notes at the primes.
- A tick voice strikes once per zero entry. The entry clock is a smooth warp of the
  unfolded heights, so the ticks keep the zeros' local spacing: GUE repulsion, by ear.
- Coda: once every zero is in, the last passes slow down, so each prime rings on its own.

Usage: python3 score.py params.json out.wav
"""
import json, sys
import numpy as np
from multiprocessing import Pool
from scipy.signal import resample_poly
from scipy.io import wavfile

P = json.load(open(sys.argv[1]))
OUT = sys.argv[2]
SR = 48000
ZSR = P.get("zsr", 4000)
N = P["N"]
DUR = P["dur"]

g = np.loadtxt("zeros1")[:N]


def unfold(t):
    return t / (2 * np.pi) * np.log(t / (2 * np.pi * np.e)) + 7 / 8


def smooth(x):
    x = np.clip(x, 0, 1)
    return x * x * (3 - 2 * x)


nu = unfold(g)
nu = nu - nu[0]
Tk = P["t0"] + P["A"] * np.log1p(nu / P["nu0"])
RAMP = P["ramp"]


def amps(T):
    return smooth((T[:, None] - Tk[None, :]) / RAMP)


# ---- pass schedule -------------------------------------------------------
u0 = np.log(P["x0"])
passes = []          # (t_start, t_end, u_end, rate)
ts_ = 0.0
coda = P.get("coda_rates", [])
ci = 0
while True:
    live = g[Tk <= ts_ + 0.5 * RAMP]
    gmax = live.max() if len(live) else g[0]
    x1 = float(np.clip(P["reach"] * gmax / np.pi, P["x1_min"], P["x1_max"]))
    u1 = np.log(x1)
    if ts_ >= Tk[-1] + RAMP:
        if ci >= len(coda):
            break
        r = coda[ci]; ci += 1
    else:
        r = P["r"]
    d = (u1 - u0) / r + P["gap"]
    passes.append((ts_, ts_ + d, u1, r))
    ts_ += d
DUR = passes[-1][1] + P.get("tail", 4.0)

t = np.arange(int(DUR * ZSR)) / ZSR
starts = np.array([p[0] for p in passes])
pidx = np.clip(np.searchsorted(starts, t, side="right") - 1, 0, len(passes) - 1)
pstart = starts[pidx]
pu1 = np.array([p[2] for p in passes])[pidx]
prate = np.array([p[3] for p in passes])[pidx]
ph = t - pstart
u = np.minimum(u0 + prate * ph, pu1 + 0.5)
inside = (u0 + prate * ph) <= pu1

B = 2000
blocks = [(i, min(i + B, len(t))) for i in range(0, len(t), B)]


def work(bl):
    i, j = bl
    a = amps(t[i:j]).astype(np.float32)
    act = np.nonzero(a.max(0) > 0)[0]
    if len(act) == 0:
        return i, np.zeros(j - i, np.complex64), np.zeros(j - i, np.float32)
    a = a[:, act]
    E = np.exp(1j * np.outer(u[i:j], g[act]).astype(np.float64)).astype(np.complex64)
    Z = (E * a).sum(1)
    W2 = (a * a).sum(1)
    return i, Z, W2


if __name__ == "__main__":
    Z = np.zeros(len(t), np.complex128)
    W2 = np.zeros(len(t))
    with Pool(P.get("procs", 8)) as pool:
        for i, z, w in pool.imap_unordered(work, blocks, chunksize=2):
            Z[i:i + len(z)] = z
            W2[i:i + len(w)] = w

    neff = np.maximum(W2, 1e-9)
    Zn = Z / np.sqrt(neff)
    Zn *= neff ** (-P["tilt"] / 2)
    edge = P["edge"]
    w = smooth(ph / edge) * smooth((pu1 - (u0 + prate * ph)) / (prate * edge))
    w[~inside] = 0
    Zn *= w

    # per-pass gain: hold the loudest struck primes at a constant level, so the
    # floor recedes under them instead of everything fading together
    m = np.abs(Zn)
    gp = np.ones(len(passes))
    for k in range(len(passes)):
        sel = (pidx == k) & inside
        if sel.sum() > 10 and m[sel].max() > 0:
            gp[k] = 1 / max(np.percentile(m[sel], P["ref_pct"]), 1e-6)
    gp = np.minimum(gp, P["max_gain"] * np.median(gp))
    gs = np.convolve(np.pad(gp, 1, mode="edge"), [0.25, 0.5, 0.25], mode="valid")
    Zn *= gs[pidx]
    m = np.abs(Zn)
    k = P["knee"]
    Zn = Zn / np.sqrt(1 + (m / k) ** 2)
    np.save(OUT + ".Z.npy", np.stack([t, u, Zn.real, Zn.imag, pidx]).astype(np.float32))

    up = SR // ZSR
    Zr = resample_poly(Zn.real, up, 1) + 1j * resample_poly(Zn.imag, up, 1)
    ts = np.arange(len(Zr)) / SR
    tu = np.interp(ts, t, u)
    tr = np.interp(ts, t, prate)

    al = P["alpha"]
    L = np.zeros(len(Zr))
    Rch = np.zeros(len(Zr))
    for h, amp, pan in P["partials"]:
        f0 = P["f0"] * h
        inst = f0 * np.exp(al * tu)
        phi = 2 * np.pi * np.cumsum(inst) / SR
        s = (Zr * np.exp(1j * phi)).real * amp
        s *= 1 / (1 + (inst / P["roll_hz"]) ** 2)
        L += s * np.cos(pan * np.pi / 2)
        Rch += s * np.sin(pan * np.pi / 2)

    # tick voice
    tick = np.zeros(len(ts))
    rate = 1 / np.maximum(np.gradient(Tk), 1e-6)
    tl = int(P["tick_len"] * SR)
    tt = np.arange(tl) / SR
    env = np.exp(-tt / P["tick_decay"]) * smooth(tt / 0.004)
    tw = (np.sin(2 * np.pi * P["tick_f"] * tt) + 0.25 * np.sin(2 * np.pi * P["tick_f"] * 2 * tt) * np.exp(-tt / 0.02)) * env
    for kk in range(N):
        i = int(Tk[kk] * SR)
        if i + tl >= len(tick):
            continue
        a = P["tick_amp"] * float(np.clip(1 - rate[kk] / P["tick_rate_max"], 0, 1)) ** 1.5
        if a <= 1e-4:
            continue
        tick[i:i + tl] += a * tw
    tp = P["tick_pan"]
    L += tick * np.cos(tp * np.pi / 2)
    Rch += tick * np.sin(tp * np.pi / 2)

    st = np.stack([L, Rch], 1)
    # tame the top: a gentle low-pass, then even out loudness across the piece
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, P["lp_hz"], fs=SR, output="sos")
    st = sosfiltfilt(sos, st, axis=0)
    win = int(P["lev_win"] * SR)
    pw = (st ** 2).mean(1)
    cs = np.concatenate([[0], np.cumsum(pw)])
    lo = np.clip(np.arange(len(pw)) - win // 2, 0, len(pw)); hi = np.clip(np.arange(len(pw)) + win // 2, 0, len(pw))
    e = np.sqrt((cs[hi] - cs[lo]) / np.maximum(hi - lo, 1) + 1e-12)
    ref = np.percentile(e, 90)
    mk = np.clip(ref / e, 1, 10 ** (P["lev_max_db"] / 20)) ** P["lev_amt"]
    st *= mk[:, None]
    fi, fo = int(P["fade_in"] * SR), int(P["fade_out"] * SR)
    st[:fi] *= smooth(np.arange(fi) / fi)[:, None]
    st[-fo:] *= smooth(np.arange(fo)[::-1] / fo)[:, None]
    st -= st.mean(0)
    st /= np.abs(st).max()
    dr = P.get("drive", 1.0)
    if dr > 1:
        st = np.tanh(dr * st) / np.tanh(dr)
        st = sosfiltfilt(sos, st, axis=0)   # filter again: the drive makes new harmonics
    st -= st.mean(0)
    st *= 10 ** (-1 / 20) / np.abs(st).max()
    wavfile.write(OUT, SR, (st * 32767).astype(np.int16))
    info = {"dur": DUR, "n_passes": len(passes),
            "passes": [[round(a, 2), round(b, 2), round(float(np.exp(c)), 1), r] for a, b, c, r in passes],
            "gamma_max": float(g[-1]), "Tk": [float(Tk[i]) for i in [0, 1, 2, 9, 99, 999, N - 1]]}
    json.dump(info, open(OUT + ".info.json", "w"), indent=1)
    np.save(OUT + ".Tk.npy", Tk)
    np.save(OUT + ".passes.npy", np.array([p[:3] for p in passes]))
    print(json.dumps(info))
