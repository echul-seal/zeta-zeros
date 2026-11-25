# Version notes

Each version was rendered on a rented box, then "listened to" by numbers: per-10-second RMS, crest factor, spectral centroid, energy above 4 kHz, DC, clipping, and a per-bar measure of how far the prime powers stand above the floor of the sum (`analyze.py`). Contact sheets of the video were checked by eye.

## v1: the idea, badly tuned
Sum of `e^{iγu}` over the first 4,000 zeros, scanned along `u = log x` in repeating bars, shifted onto a fixed 110 Hz carrier. Zeros entered on a straight time warp.
- The fixed carrier made every bar the same smear: the primes had no pitch of their own, so you could not hear *which* prime lit up, only that something did.
- The straight warp spent the first minute with a handful of zeros and then dumped thousands at once. No arc.

## v2: primes get their own pitch
- Carrier pitch now follows x (`f = f0·x^α`), so each prime power rings at its own fixed note, the same note in every bar. The bar becomes a rising scale whose notes are the primes.
- Zero entry moved to a log clock in the zeros' unfolded height: zero 1 is alone for 1.5 s, zero 10 arrives at 16 s, zero 100 at 56 s, zero 1,000 at 116 s, zero 10,000 at 180 s.
- Each bar reaches as far into the primes as the present zeros can resolve (x up to about γ_max/π), so the scale grows longer as the piece goes on, from 12 to 2,400.
- Went to 10,000 zeros. 4,000 barely got past x ≈ 1,000.
- Problem: prime contrast climbed nicely (3 dB in bar 2, 34 dB by bar 22) but overall level sagged 13 dB from the start to the end, and above 160 s the carrier for high primes pushed real energy above 4 kHz (from −72 dB to −25 dB). Harsh, and quieter exactly where it should be most present. Also ran 282 s with a long, slow coda.

## v3: tame the top
- Lowered α from 0.33 to 0.27, so x = 2,400 lands near 900 Hz instead of 1.4 kHz.
- 1/(1+(f/1600)²) rolloff per partial, then a 4th-order 3.2 kHz low-pass on the mix.
- Slow loudness rider (6 s window, up to +9 dB) to keep the late bars from sinking.
- Coda shortened to two bars (half and quarter speed). 251 s.
- Result: >4 kHz energy stayed under −71 dB everywhere. Level drift fell from 13 dB to 6 dB.

## v4: level
- Loudness rider pushed harder (90%, up to +10 dB). Drift across the body is now about 4 dB; the remaining dips are the quiet bar lines of the slow coda, which I want.
- First full video. Contact sheet showed the trace was being decimated to 2,500 points per bar, which silently skipped narrow prime spikes in the late bars: the most important picture in the piece was undersampled.

## v5: drive
- Trace is max-pooled per pixel column, so every spike that sounds is also drawn.
- A gentle tanh drive (2.0) on the master, so the struck primes sit forward of the wash. Loudness came up about 6 dB (RMS -27 to -22 dBFS) and the body now holds within 4 dB.
- Problem: the drive made its own harmonics after the low-pass. Energy above 4 kHz climbed back from -84 dB to -42 dB in the densest bars (150-190 s). The exact harshness v3 removed, reintroduced by the master.

## v6: final
- Low-pass applied again after the drive.
- Same score otherwise. Peak -1 dBFS, zero clipped samples, DC around 3e-6.
