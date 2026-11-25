"""Listen by numbers: levels, clipping, DC, spectral balance per 10 s, and how clearly
the primes stand out of the sum in each pass (prime-to-floor contrast)."""
import sys, json
import numpy as np
from scipy.io import wavfile
from scipy.signal import stft
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

wav = sys.argv[1]
sr, x = wavfile.read(wav)
x = x.astype(np.float64) / 32768
m = x.mean(1)
out = {}
out["dur"] = len(x) / sr
out["sr"]=int(sr)
out["peak_dbfs"] = float(20 * np.log10(np.abs(x).max()))
out["clipped_samples"] = int((np.abs(x) > 0.999).sum())
out["dc"] = x.mean(0).tolist()
out["rms_dbfs"] = float(20 * np.log10(np.sqrt((m ** 2).mean())))
seg = int(10 * sr)
rows = []
f, tt, S = stft(m, sr, nperseg=4096, noverlap=2048)
P = np.abs(S) ** 2
for i in range(0, len(m), seg):
    s = m[i:i + seg]
    if len(s) < sr: break
    j0, j1 = np.searchsorted(tt, [i / sr, (i + seg) / sr])
    p = P[:, j0:j1].mean(1)
    cen = float((f * p).sum() / (p.sum()+1e-30))
    hi = float(10 * np.log10(p[f > 4000].sum() / (p.sum()+1e-30) + 1e-12))
    rms = float(20 * np.log10(np.sqrt((s ** 2).mean()) + 1e-12))
    crest = float(20 * np.log10(np.abs(s).max() / (np.sqrt((s ** 2).mean()) + 1e-12)))
    rows.append([int(i // sr), round(rms, 1), round(crest, 1), round(cen), round(hi, 1)])
out["per10s [t, rms dB, crest dB, centroid Hz, >4k dB]"] = rows

# prime contrast from the envelope
try:
    t, u, zr, zi, pidx = np.load(wav + ".Z.npy")
    mag = np.hypot(zr, zi)
    pp = []
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    for p_ in np.unique(pidx).astype(int):
        sel = pidx == p_
        uu, mm = u[sel], mag[sel]
        near = np.zeros_like(uu, bool)
        for q in primes + [4, 8, 9, 16, 25, 27, 32]:
            near |= np.abs(uu - np.log(q)) < 0.01
        if near.sum() == 0 or (~near).sum() == 0: continue
        on = np.percentile(mm[near], 95); off = np.median(mm[~near]) + 1e-9
        pp.append([int(p_), round(float(t[sel][0]), 1), round(float(20 * np.log10((on+1e-9) / off)), 1)])
    out["pass contrast [pass, t, primes vs floor dB]"] = pp
except Exception as e:
    out["contrast_err"] = str(e)

plt.figure(figsize=(14, 5))
plt.imshow(10 * np.log10(P[f < 6000] + 1e-12), aspect="auto", origin="lower",
           extent=[0, tt[-1], 0, 6000], vmin=-110, vmax=-30, cmap="magma")
plt.colorbar(); plt.title(wav)
plt.savefig(wav + ".spec.png", dpi=80)
print(json.dumps(out))
