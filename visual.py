"""zeta-zeros: the picture.

Black ground. Left: the critical line, each zero a point at its height (log scale),
lighting as it enters. Right: the current bar, drawn as it is played: the envelope
|Z(u)| along log x, a thin trace that is noise early and a row of spikes at the primes
late. Below the trace, faint marks at the prime powers the trace has hit. A playhead.

Usage: python3 visual.py params.json out.wav out_video.mp4
"""
import sys, json, subprocess
import numpy as np
from PIL import Image, ImageDraw, ImageFont

P = json.load(open(sys.argv[1]))
WAV = sys.argv[2]
OUT = sys.argv[3]
W, H, FPS = 1920, 1080, 30
t, u, zr, zi, pidx = np.load(WAV + ".Z.npy")
mag = np.hypot(zr, zi)
Tk = np.load(WAV + ".Tk.npy")
passes = np.load(WAV + ".passes.npy")
g = np.loadtxt("zeros1")[:len(Tk)]
DUR = float(json.load(open(WAV + ".info.json"))["dur"])
u0 = np.log(P["x0"])
UMAX = np.log(P["x1_max"])

# layout
LX = 170                     # critical line x
PX0, PX1 = 360, 1840         # plot area for log x
PY = 640                     # baseline of trace
PH = 380                     # trace height


def ux(uu):
    return PX0 + (uu - u0) / (UMAX - u0) * (PX1 - PX0)


gl = np.log(g)
gy = lambda gg: 1000 - (np.log(gg) - np.log(10)) / (np.log(g[-1] * 1.05) - np.log(10)) * 920

pp = []
n = 2
while n <= P["x1_max"]:
    m, p = n, None
    for q in range(2, int(n ** 0.5) + 1):
        if n % q == 0:
            p = q
            break
    if p is None:
        pp.append((n, n, 1.0))
    else:
        k = n
        while k % p == 0:
            k //= p
        if k == 1:
            pp.append((n, p, 0.5))
    n += 1

try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 22)
except Exception:
    font = ImageFont.load_default()

ff = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
                       "-s", f"{W}x{H}", "-r", str(FPS), "-i", "-", "-i", WAV,
                       "-c:v", "libx264", "-preset", "medium", "-crf", "18", "-pix_fmt", "yuv420p",
                       "-c:a", "aac", "-b:a", "256k", "-shortest", "-movflags", "+faststart", OUT],
                      stdin=subprocess.PIPE)

nframes = int(DUR * FPS)
ZSR = 1 / (t[1] - t[0])
for fr in range(nframes):
    T = fr / FPS
    im = Image.new("RGB", (W, H), (0, 0, 0))
    d = ImageDraw.Draw(im)
    # critical line
    d.line([(LX, 70), (LX, 1010)], fill=(38, 38, 38), width=1)
    live = Tk <= T
    nl = int(live.sum())
    if nl:
        age = T - Tk[:nl]
        ys = gy(g[:nl])
        br = np.clip(80 + 175 * np.exp(-age / 1.5), 80, 255)
        for y, b in zip(ys[-1500:], br[-1500:]):
            d.point((LX, int(y)), fill=(int(b),) * 3)
            d.point((LX + 1, int(y)), fill=(int(b * .6),) * 3)
        for k in range(max(0, nl - 40), nl):
            if age[k] < 1.2:
                r = 2 + 4 * (1 - age[k] / 1.2)
                y = gy(g[k])
                c = int(255 * (1 - age[k] / 1.2))
                d.ellipse([LX - r, y - r, LX + r, y + r], outline=(c, c, c))
    # counter
    d.text((LX - 60, 1030), f"{nl:>5d}", fill=(110, 110, 110), font=font)

    # current pass trace
    i = int(T * ZSR)
    i = min(i, len(t) - 1)
    k = int(pidx[i])
    sel = np.nonzero(pidx == k)[0]
    sel = sel[sel <= i]
    if len(sel) > 2:
        uu = u[sel]
        mm = mag[sel]
        ok = uu <= passes[k][2]
        uu, mm = uu[ok], mm[ok]
        if len(uu) > 2:
            # one value per pixel column, max-pooled, so no spike is skipped
            cols = np.floor(ux(uu)).astype(int)
            c0, c1 = cols[0], cols[-1]
            mx = np.full(c1 - c0 + 1, 0.0)
            np.maximum.at(mx, cols - c0, mm)
            xs = np.arange(c0, c1 + 1)
            ys = PY - np.clip(mx, 0, 1.15) * PH
            d.line(list(zip(xs.tolist(), ys.tolist())), fill=(225, 225, 225), width=2)
            # ghost of the previous pass
        cur_u = uu[-1] if len(uu) else u0
        hx = ux(min(cur_u, passes[k][2]))
        d.line([(hx, PY + 6), (hx, PY + 20)], fill=(90, 90, 90), width=1)
    # baseline and prime marks: a prime is lit once the zeros are fine enough to resolve it
    d.line([(PX0, PY + 1), (PX1, PY + 1)], fill=(30, 30, 30), width=1)
    gmax = g[nl - 1] if nl else 0
    for n_, p_, w_ in pp:
        x = ux(np.log(n_))
        res = gmax * np.log(n_ + 1) - gmax * np.log(n_) if gmax else 0
        on = np.clip((res - 1.5) / 3, 0, 1)
        if on <= 0:
            continue
        c = int(30 + 150 * on * w_)
        hgt = 10 if w_ == 1 else 5
        d.line([(x, PY + 30), (x, PY + 30 + hgt)], fill=(c, c, c), width=1)
    for n_ in [2, 3, 5, 7, 11, 13, 101, 1009]:
        if n_ <= P["x1_max"]:
            x = ux(np.log(n_))
            on = np.clip(((gmax * (np.log(n_ + 1) - np.log(n_)) if gmax else 0) - 1.5) / 3, 0, 1)
            if on > 0:
                c = int(90 * on)
                d.text((x - 6, PY + 52), str(n_), fill=(c, c, c), font=font)

    # gentle fades
    fade = min(1, T / 3, (DUR - T) / 4)
    if fade < 1:
        a = np.asarray(im, dtype=np.float32) * max(fade, 0)
        im = Image.fromarray(a.astype(np.uint8))
    ff.stdin.write(im.tobytes())
ff.stdin.close()
ff.wait()
print("ok", nframes)
