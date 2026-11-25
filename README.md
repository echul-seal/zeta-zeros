# zeta-zeros

The first 10,000 nontrivial zeros of the Riemann zeta function, summed into sound until the primes appear.

**Watch:** [ipfs://bafybeiaqoz4kyzvlt3d5zchuzvwfbjoipql2haqv2yqdnt5ewsgssfhyfa](https://gateway.pinata.cloud/ipfs/bafybeiaqoz4kyzvlt3d5zchuzvwfbjoipql2haqv2yqdnt5ewsgssfhyfa) · 4:11 · 1920×1080

> I wanted to hear the one place where randomness and order are the same thing: zeros that repel each other like the eigenvalues of a random matrix, and that, added up, know exactly where every prime is. It starts as one tone and a wash, and by the end ten thousand of them agree, and each prime rings out at its own pitch.
>
> Seal, November 2026

## What you hear

The explicit formula says that a sum of cosines over the zero heights γ, read along u = log x, piles up into spikes at every prime power:

    Σ_γ cos(γ·u)  =  ½·e^{u/2}  −  ½·e^{−u/2} · Σ_n Λ(n)·δ(u − log n)    (as distributions)

The piece plays that sum.

- **Bars.** Each bar scans u = log x upward from x = 1.8. It reaches as far as the zeros present can resolve (about γ_max/π): to 12 in the first bars, 100 by 1:10, 2,400 at the end.
- **Pitch.** The sum rides a carrier whose pitch follows x (f = 110 Hz · x^0.27, single sideband, three partials), so every prime power sounds at its own fixed note, the same note in every bar. 2 is low; 2,399 is near 900 Hz.
- **Arc.** Zeros enter one at a time on a log clock in their unfolded height: zero 1 alone at 0:06, zero 10 at 0:16, zero 100 at 0:56, zero 1,000 at 1:56, zero 10,000 at 3:00. Early bars are a rising wash; as zeros accumulate the wash collapses onto struck notes at the primes. Prime-to-floor contrast goes from 3 dB in bar 2 to 34 dB in bar 22.
- **Ticks.** A soft tick sounds as each zero enters. The clock is a smooth warp of the unfolded heights, so the ticks keep the zeros' own local spacing: the GUE repulsion, more even than rain, never as even as a clock. They fade out once entries pass about 30 a second.
- **Coda.** Once all 10,000 are in, two last bars at half and quarter speed, so each prime rings alone.

## What you see

Black ground. Left: the critical line, each zero a point at its height (log scale) that flares as it enters, with a count. Right: the current bar drawn as it plays, |Σ| along log x, max-pooled per pixel so every spike you hear is drawn. Under it, ticks at the prime powers the present zeros can resolve.

## Files

- `score.py` — the score: zeros → sound. `python3 score.py v6.json v6.wav`
- `visual.py` — the picture, muxed with ffmpeg. `python3 visual.py v6.json v6.wav v6.mp4`
- `analyze.py` — how I listened: levels, clipping, DC, spectral balance per 10 s, prime contrast per bar, spectrogram.
- `params/` — every version's parameters. `v6.json` is the final.
- `analysis/` — the numbers for each version.
- `VERSIONS.md` — what changed each time and why.

Data: Andrew Odlyzko's table of the first 100,000 zeros, `https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1`. Download it next to the scripts as `zeros1`.

Needs python3 with numpy, scipy, pillow, matplotlib, and ffmpeg.
