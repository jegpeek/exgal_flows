# How the radial reversal was found

*A record of the reasoning path that led to the "CO–divergence relation
reverses sign with galactocentric radius" result. Written because the
finding was not something anyone set out to look for — it fell out of
chasing down an unexpected sign flip during a routine cross-check.*

**Date:** 2 July 2026
**Context:** refactoring the exploratory M33 notebook into the
`exgal_flows` package, then re-running the headline results on the new
VLA A+B+C+D+GBT HI data (3.05 pc/pixel) as part of the requested
new-vs-old epoch comparison.

---

## The one-line summary

While confirming that "CO-bright gas sits at converging (colliding) flows"
still held on the new high-resolution data, the number came out with the
**wrong sign**. Running that discrepancy to ground — through two wrong
guesses about its cause — revealed that the CO–divergence relation is not a
single number at all: it **changes sign with galactocentric radius**.
Diverging flows dominate the inner ~1 kpc; colliding flows dominate the
outer disk, where 90% of the CO-bright gas actually lives. The notebook's
original negative result was the whole-disk average; a fixed *pixel*
aperture had silently been re-measuring only the inner kpc on the finer new
grid.

## The trigger: a sign that shouldn't have flipped

The plan was mundane — run `result_divergence_of_co` on the new data and
compare to the old epoch. The expectation, from the original notebook, was
that CO-bright pixels sit at **negative** divergence (converging), more so
for brighter CO. The new data returned the opposite:

```
CO 2x-3x: +2.82   CO 3x-4x: +2.17   CO 4x-7x: +4.37   CO >7x: +6.94  km/s/kpc
```

Positive, and *steepening* positive with CO brightness. That is a clean,
confident, wrong-signed result — exactly the kind of thing that is
dangerous to wave away. There were two honest possibilities:

- **(a)** the higher-resolution A-array data genuinely changes the science
  (which is precisely what a new-vs-old comparison is *for*), or
- **(b)** the refactor introduced an arithmetic difference from the
  notebook somewhere.

The rule followed here: **you cannot report (a) until you have ruled out
(b).** A refactor that changes a sign has to be presumed guilty.

## Wrong guess #1: "my package flipped a sign"

To isolate the package, the notebook's divergence computation was
reproduced *inline* — straight from the cached velocity field, not through
any `exgal_flows` code — on the same pixels. That gave:

```
notebook-style, from cached velocity:  -2.28 -> -2.54 -> -2.73 -> -5.74  (steepening negative)
```

Negative and steepening — the original result. And the cached `div` field
(produced by the package pipeline) *agreed* with it. So the raw divergence
field was fine, yet `result_divergence_of_co` reported positive. The
natural conclusion at that point was: **the sign flip lives inside the
statistics function, not the data.** That felt like a located bug.

## Wrong guess #2: "the cached arrays are internally inconsistent"

Inspecting the stat function deepened the confusion instead of resolving
it. Using the *same cached `div` array*, the statistic returned **+8.14**
for the `>7x` bin, while the manual recompute on what looked like the
identical `co > 4.2` mask gave **−5.74**. Same field, same threshold,
opposite sign. The only way that seemed possible was if the cached `div`
and cached `velocity` arrays disagreed in sign somewhere — a corrupted
cache.

So that was tested directly, pixel-by-pixel:

```
corr(cached div, gradient-from-cached-velocity) = 0.9999999999999
median ratio div/gx                              = 1.0000000
```

Perfectly consistent. The cache was not the problem. Both "bug" hypotheses
were now dead — and that is what forced the actual explanation into view.

## The realization: it was the aperture, and the aperture is a radius

The two manual tests had **not** been measuring the same pixels. One
selected `rad < 916` and the other `rad < 300`, and — the crucial
oversight — `rad` is stored in **pixels**, not parsecs. The `< 916` test
had quietly integrated out to ~2795 pc; the stat's `< 300 px` covered only
~915 pc. They disagreed in sign because **they sampled different parts of
the galaxy.**

That reframed everything. The question was never "which number is right"
— both were right, for their respective regions:

> **The sign is radius-dependent.** Inner ~915 pc → positive divergence
> (diverging); larger radii → negative (colliding).

At that point the sign flip stopped being a bug to fix and became a result
to map.

## Turning it into a measurement

The radial dependence was measured properly with annuli from 0 to 4000 pc,
for three CO-brightness thresholds. The relation is a smooth, monotonic
gradient, not a threshold artifact:

| CO bin | inner disk (r < 916 pc) | whole disk (r < 3600 pc) |
|--------|------------------------:|-------------------------:|
| 2x–3x  | +3.25 | −1.53 |
| > 7x   | +8.03 | **−7.26** [−12.0, −1.4] |

- **r < ~1 kpc:** CO-bright gas sits at mildly *positive* divergence. This
  inner region is dominated by M33's central molecular structures and
  strong rotation-curve shear.
- **r > ~1.5 kpc:** CO-bright gas sits at *negative* divergence, and the
  brightest CO reaches ≈ −9.5 km/s/kpc — the "CO at colliding flows"
  result, cleanly recovered, brightest bin significant.
- The zero-crossing (turnover) sits near ~1.0–1.2 kpc.

A radial-extent check settled the default aperture: half of all CO>2×
pixels lie beyond 2000 pc, and 90% within 4000 pc. A fixed 300-pixel
aperture captures only the inner ~10% of the CO disk on the new data.

## The reconciliation (why the notebook was never "wrong")

The last piece was showing the old and new results are the *same physics*
seen through different apertures:

| Data | pc / pixel | `rmax = 300 px` covers | Which regime |
|------|-----------:|-----------------------:|--------------|
| Old (~12 pc/pix) | ~12 | ~3600 pc | ~whole disk → **negative** (colliding) |
| New (3.05 pc/pix) | 3.05 | ~916 pc | inner kpc → **positive** (diverging) |

The notebook's `rmax = 300` was written in pixels against ~12 pc/pixel
data, where it happened to encompass the whole molecular disk. Ported to
3 pc/pixel data, the *same 300 pixels* silently shrank 4× in physical size
and started reporting the inner-kpc regime instead. The result didn't
break; the aperture changed underneath it. This is exactly the failure mode
the pixels→parsecs refactor was meant to prevent — it just took a
sign flip to notice the trap was real.

## What made this findable

A few habits did the actual work here, and they are worth naming because
none of them are specific to this dataset:

1. **Treating a confident wrong-signed number as a stop, not a nuisance.**
   The positive divergence could have been rounded away as "close to zero."
   It wasn't.
2. **Ruling out the tool before believing the science.** The inline,
   package-free recompute is what proved the data were sound and pointed
   the search inward.
3. **Distrusting a located "bug" that doesn't fully add up.** The
   "it's in the stat function" conclusion was abandoned the moment the
   same cached field gave two signs — because a real explanation has to
   account for *all* the evidence, not most of it.
4. **Checking units on the thing you're slicing by.** The whole reversal
   hid behind `rad` being in pixels. One `corr = 1.0` sanity check on the
   cache is what exposed it.

The reversal is a real astrophysical gradient in M33 — inner rotation/shear
regime vs. outer colliding-flow regime — and it was surfaced by refusing to
explain away a sign that shouldn't have flipped.

---

## Where to look

- `outputs/radial_divergence.png` — the radius-resolved divergence profile
  (the figure that came directly out of this investigation).
- `outputs/radial_blockboot.png` — a later spatial block-bootstrap check
  confirming the reversal survives correlated-noise-aware errors
  (iid errors understate the scatter ~10×; the inner-vs-outer difference is
  still +8.65 km/s/kpc, ~2.9σ).
- `FINDINGS.md` — the aperture caveat and all three headline results in
  their final, reconciled form.
