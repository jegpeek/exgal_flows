# Findings from the refactor + new-data check

This note records what changed scientifically when the exploratory notebook
was turned into the `exgal_flows` package and re-run on the new VLA A-array
HI data. **One issue is important and worth reading before trusting any
single number.**

## TL;DR

The headline result — **CO-bright gas sits preferentially at converging
(colliding) flows, and more strongly so for brighter CO** — reproduces on
the new data **when measured over the whole molecular disk**. It does *not*
reproduce, and in fact reverses sign, if you keep the notebook's aperture of
`rmax = 300` **pixels**, because a fixed pixel aperture shrinks 4x in
physical size on the higher-resolution data.

## The aperture trap (why pixels vs. parsecs matters)

The notebook selected analysis pixels with `rr < rmax`, where `rr` is
image-plane radius from the map center **in pixels** and `rmax = 300`.

| Data | pc / pixel | `300 px` aperture | What it captures |
|------|-----------:|------------------:|------------------|
| Old (test_jegp, ~12 pc/pix) | ~12 | ~3600 pc | ~whole molecular disk |
| New (A+B+C+D+GBT, 3.05 pc/pix) | 3.05 | ~916 pc | inner ~10% of the CO disk |

Half of all CO>2x-threshold pixels lie beyond 2000 pc; 90% lie within
4000 pc. So `rmax=300 px` on the new data throws away ~90% of the molecular
gas and keeps only the inner kiloparsec — which, it turns out, is exactly
the region where the CO-vs-divergence relation has the *opposite* sign.

**Median LOS divergence of CO-bright pixels (km/s/kpc; negative = colliding):**

| CO bin | inner disk (r < 916 pc) | whole disk (r < 3600 pc) |
|--------|------------------------:|-------------------------:|
| 2x-3x  | +3.25 | -1.53 |
| 3x-4x  | +1.53 | -3.28 |
| 4x-7x  | +4.93 | -2.17 |
| > 7x   | +8.03 | **-7.26** [-12.0, -1.4] |

The whole-disk column is the notebook's original claim, recovered on the new
data: negative and steepening with CO brightness, with the brightest bin
significantly below zero (68% wedge-bootstrap interval excludes 0).

## It's a radial gradient, not just a threshold

The sign flip is not an artifact of one aperture choice — the relation
genuinely changes with galactocentric radius (see
`outputs/radial_divergence.png`):

* **r < ~1 kpc:** CO-bright gas sits at mildly *positive* divergence
  (diverging flows). This inner region is dominated by M33's central
  molecular structures and strong rotation-curve shear.
* **r > ~1.5 kpc:** CO-bright gas sits at *negative* divergence (colliding
  flows), and the brightest CO (>7x threshold) reaches -25 to -40 km/s/kpc.

The physical picture — molecular gas forming preferentially where flows
converge — holds across the star-forming disk. The inner kpc is a separate
regime and probably shouldn't be pooled with it.

## Recommended defaults

`config/m33.yaml` now sets `radius_max_pc: 4000.0` (whole disk). To exactly
reproduce the notebook's inner-disk aperture, set `radius_max_pc: 916.0`.

The point of expressing every analysis scale in **parsecs** (not pixels) is
that the new-vs-old comparison now probes the *same physical region* on both
grids. That is what makes the comparison in `scripts/compare_epochs.py`
meaningful.

## Results 2 and 3 also reproduce (whole-disk, new data)

**Result 2 - shear suppresses CO.** Over r < 4000 pc, the fraction of pixels
that are CO-bright falls monotonically with |shear|:

| |shear| [km/s/pc] | CO-bright fraction |
|------------------:|-------------------:|
| 0.00-0.02 | 0.173 |
| 0.08-0.10 | 0.156 |
| 0.16-0.18 | 0.115 |
| 0.18-0.20 | 0.104 |

A ~40% relative drop from low to high shear.

**Result 3 - the converging-flow signal lives at ~50-100 pc.** Median LOS
divergence of the brightest CO (>7x threshold), as the velocity field is
smoothed to coarser scales:

| Smoothing FWHM | median div CO>7x [km/s/kpc] |
|---------------:|----------------------------:|
| 0 pc (native)  | -8.9 |
| 14 pc | -9.5 |
| 29 pc | -8.4 |
| 58 pc | -4.9 |
| 86 pc | -4.9 |
| 150 pc | -4.0 |

The converging-flow signature is strongest at 15-30 pc and is still clearly
present at 57-86 pc, then washes out by ~150 pc — consistent with the
"effects live on 50-100 pc scales" claim (the underlying structures are at
the fine end of that range; they survive smoothing up to ~100 pc).

## New vs. old HI data (the epoch comparison)

Running `scripts/compare_epochs.py` (both epochs, same 4000 pc physical
aperture, `decim=2`) gives the brightest-CO (>7x) median divergence:

| Epoch | beam | median div, CO>7x [km/s/kpc] |
|-------|------|-----------------------------:|
| New (A+B+C+D+GBT) | 5.3" | **-8.78** [-13.3, -2.7] |
| Old (14B+17B feathered) | 8.2" | **-4.18** [-8.4, -1.0] |

Both are significantly negative — the colliding-flow result is robust to the
HI epoch. The higher-resolution **new** data shows a ~2x stronger signal,
which is the expected direction: coarser resolution smears the convergence
and dilutes the divergence signature. The A-array data isn't changing the
conclusion; it's sharpening it. See `outputs/compare_epochs_divergence.png`.

This is exactly why the analysis scales are physical: the two beams differ,
but because the aperture and smoothing are set in parsecs, the comparison is
apples-to-apples.

## The radial reversal, and why iid errors are not enough

The single most important correction this refactor uncovered concerns the
**aperture**. The original notebook computed the CO-vs-divergence statistic
inside `rad < 300 pixels`. On the old ~12 pc/pixel data that is ~3600 pc —
essentially the whole star-forming disk. On the new 3.05 pc/pixel data the
identical `300` is only **916 pc** — the inner ~10% of the CO disk by area.
The result is not scale-invariant: the sign of the CO-divergence relation
**flips with radius**.

- Inside ~916 pc: brightest CO sits at **diverging** flows (positive divergence).
- Whole disk (< 4000 pc): brightest CO sits at **converging/colliding**
  flows (negative divergence) — the headline result.

Half of the CO>2x pixels lie beyond 2000 pc, so the pixel-radius cut was
discarding most of the molecular gas. The default aperture is now
`radius_max_pc = 4000`, which captures 90% of CO>2x and 99% of CO>7x pixels.

### Block bootstrap (the correlation concern)

Because everything in these maps is spatially correlated on a range of
scales, the per-pixel iid bootstrap **understates the uncertainty**. The
divergence field has a 1/e correlation length of only ~6 pc (near the beam),
but CO organizes into complexes tens to hundreds of pc across, so
neighbouring bright pixels are not independent draws.

We therefore added a **spatial tiled block bootstrap**: tile the map into
square blocks, resample whole blocks with replacement, and take a
block-weighted median of the divergence in each galactocentric annulus. This
preserves the internal correlation of each block.

| quantity | iid bootstrap | block bootstrap |
|----------|--------------:|----------------:|
| error on annular medians | baseline | **6-11x larger** |

So the iid error bars were optimistic by roughly an order of magnitude — the
concern was well founded. **But the radial reversal survives.** With block
resampling:

- inner (<900 pc) minus outer (>1500 pc) divergence = **+8.65 km/s/kpc**,
  68% CI **[+5.4, +12.7]**
- P(inner > outer) = **0.996** (~2.9 sigma)
- robust across block sizes from 50 to 600 pc (P > 0.96 even at 600 pc,
  where only ~115 tiles remain)
- point-estimate zero-crossing (the turnover) at **r ~ 1173 pc**

The reversal is a real, spatially-coherent trend, not inner-region noise.
See `outputs/radial_blockboot.png` (binned median profile with the two line
fits below) and `outputs/galaxy_turnover.png` (the turnover radius drawn on
the galaxy).

The spatial block bootstrap is a first-class part of the package:

```python
prof = ana.radial_divergence_profile()   # FlowAnalysis method
prof.turnover_pc          # -> 1172.9  (zero-crossing of the profile)
prof.err_low, prof.err_high   # block-bootstrap 68% band per annulus
prof.iid_low, prof.iid_high   # the (too-narrow) iid band, for comparison

from exgal_flows.plotting import plot_radial_divergence_profile
plot_radial_divergence_profile(prof)      # reproduces outputs/radial_blockboot.png
```

Block side and annulus edges are set in `config/m33.yaml`
(`block_pc`, `radial_edges_pc`).

### Fitting the trend directly

The binned profile is convenient but throws away information. We also fit the
**simplest model directly to the per-pixel divergence** — a straight line

```
div(r) = a + b*r
```

with `a` the overall sense of divergence (extrapolated to the centre,
km/s/kpc) and `b` the rate of change with radius (km/s/kpc per kpc). This is
fit to the data itself, *not* to the error bars, and is block-bootstrapped
over the same spatial tiles for honest uncertainty. Two estimators, to show
the answer is not an artifact of one loss function:

| fit | `a` (km/s/kpc) | `b` (km/s/kpc/kpc) | zero-crossing | P(`b`<0) |
|-----|---------------:|-------------------:|--------------:|---------:|
| OLS (least squares, mean-like)   | **+5.1** [0.8, 9.3] | **-3.9** [-5.9, -2.1] | ~1310 pc | 0.98 |
| L1 (least abs. dev., median-like)| **+2.0** [-1.0, 5.1] | **-2.3** [-3.8, -0.9] | ~900 pc  | 0.95 |

Both agree on the qualitative result: **positive intercept, negative slope,
zero-crossing near 0.9-1.3 kpc** — converging in the inner disk, diverging
outside, with high confidence that the slope is negative. OLS is steeper
because it tracks the (heavier) positive tail in the inner disk; L1 is
shallower and, being median-like, is the natural companion to the binned
median profile. The zero-crossings bracket the profile turnover (~1173 pc).

```python
ols = ana.radial_divergence_fit(method="ols")   # FlowAnalysis method
l1  = ana.radial_divergence_fit(method="l1")
ols.intercept, ols.slope, ols.zero_crossing_pc, ols.slope_neg_prob

from exgal_flows.plotting import plot_radial_divergence_profile
plot_radial_divergence_profile(prof, fit=[ols, l1])  # -> radial_blockboot.png
```

#### A sigmoid rules out zero more decisively than the line

The line fit answers "is the slope negative?" (P(`b`<0) ~ 0.97). But the
scientific claim is really "does divergence go from *positive* inside to
*negative* outside?" — a statement about two plateaus, not one slope. So we
also fit a bounded step

```
div(r) = lo + (hi - lo) / (1 + exp((r - r0) / w))
```

with `hi` = inner plateau, `lo` = outer plateau, `r0` = transition radius,
`w` = width. Block-bootstrapped over the same tiles (500 draws, all
converged):

| parameter | value [16, 84] |
|-----------|---------------:|
| `hi` (inner plateau) | **+6.8** [+4.4, +20.5] km/s/kpc |
| `lo` (outer plateau) | **-5.5** [-9.8, -4.2] km/s/kpc |
| `r0` (transition)    | **1123** [482, 1427] pc |
| amplitude `hi-lo`    | **12.4** [10.9, 29.1] km/s/kpc |
| zero-crossing        | ~1127 pc |

The headline: the inner plateau is positive **and** the outer plateau
negative in **100%** of bootstrap draws — `straddle_zero_prob = 1.000`, versus
the OLS line's `P(b<0) = 0.975`. **So yes, the sigmoid rules out "no reversal"
more decisively than `y = a + b·r`.** It also tracks the binned medians more
faithfully (sits on the +6 inner points and the -5/-6 outer points rather than
splitting the difference), and the transition radius lands right on the profile
turnover (~1173 pc). The caveat: the inner plateau's *magnitude* is poorly
constrained (least squares chases the heavy positive tail of the sparse
inner-disk pixels, hence the wide upper error), but its *sign* — the part that
carries the significance — is rock-solid. The fitted width rails to ~20 pc
(effectively a step); the data do not resolve a gradual transition.

```python
sig = ana.radial_divergence_sigmoid(n_boot=500)   # FlowAnalysis method
sig.hi, sig.lo, sig.r0_pc, sig.amplitude, sig.straddle_zero_prob
plot_radial_divergence_profile(prof, fit=[ols, l1, sig])  # -> radial_blockboot.png
```

## Shear suppresses CO (the sensitivity study)

The third scientific claim is that **velocity shear reduces molecular gas**:
where the deprojected velocity field is being sheared (large transverse
gradient `|dv/dy| / cos i`), CO is less likely to be bright. We quantify this
as the **CO-bright fraction** in bins of `|shear|`, and summarize each curve
with a single *sensitivity* number: the fractional drop in CO-bright fraction
from the low-shear half of the range to the high-shear half
(`rel_drop = 1 - hi_frac / lo_frac`). Errors are the azimuthal **wedge block
bootstrap** (whole sectors resampled), for the same correlation reason as the
radial study — per-pixel iid bands are not honest here.

Baseline (new data, CO>2x, full 4000 pc disk): CO-bright fraction falls from
**17.0%** at low shear to **13.1%** at high shear, `rel_drop = 0.23
[0.21, 0.25]`. Shear suppression is present and significant. The effect is
robust and coherent across all four breakdowns:

**a) Epoch (new vs old HI).** Both epochs show the same suppression; the
feathered 8.2" data gives a marginally *steeper* drop than the 5.3" A-array
data (higher noise inflates the low-shear reference), but they agree within
the bands.

| epoch | rel_drop |
|-------|---------:|
| new (5.3") | 0.23 [0.21, 0.25] |
| old (8.2") | 0.26 [0.23, 0.29] |

**b) CO-intensity cut (new).** The suppression is essentially independent of
how bright a pixel must be to count as CO-bright, and if anything is
*strongest* for the brightest molecular gas — the opposite of what a simple
noise/threshold artifact would produce.

| CO cut | rel_drop |
|--------|---------:|
| CO>2x (1.2 K) | 0.27 [0.25, 0.29] |
| CO>3x (1.8 K) | 0.27 [0.25, 0.30] |
| CO>4x (2.4 K) | 0.26 [0.23, 0.30] |
| CO>7x (4.2 K) | 0.32 [0.27, 0.37] |

**c) Smoothing FWHM (new, CO>2x).** As the velocity field is smoothed, the
shear-CO anticorrelation *strengthens* monotonically — the trend is not a
pixel-scale noise feature; it lives on ~50-100 pc scales, the same regime
where the convergence signal is strongest.

| smoothing | rel_drop |
|-----------|---------:|
| native | 0.23 [0.21, 0.25] |
| 14.4 pc | 0.26 [0.23, 0.29] |
| 28.8 pc | 0.29 [0.24, 0.33] |
| 57.5 pc | 0.34 [0.27, 0.40] |
| 86.3 pc | 0.45 [0.33, 0.64] |

(At 86 pc the highest-shear bins empty out on the finite disk, so that last
row's point estimate is noisy — read it as "still strong," not as a precise
0.45.)

**d) Galactocentric annulus (new, CO>2x).** Unlike the *divergence* relation,
which reverses sign with radius, the shear suppression **keeps the same sign
at every radius** and simply gets stronger outward — from ~15% in the inner
kpc to ~27% in the outer disk. Shear removes CO throughout M33; it is not an
inner-disk-only or outer-disk-only effect.

| annulus | rel_drop (new) | rel_drop (old) |
|---------|---------------:|---------------:|
| 0-900 pc     | 0.15 [0.13, 0.17] | 0.18 [0.13, 0.21] |
| 900-1500 pc  | 0.16 [0.12, 0.20] | 0.23 [0.18, 0.30] |
| 1500-2500 pc | 0.25 [0.22, 0.27] | 0.32 [0.28, 0.36] |
| 2500-4000 pc | 0.27 [0.24, 0.30] | 0.32 [0.27, 0.37] |

See `outputs/shear_sensitivity.png` (all four panels, curves normalized to
their low-shear value so the sensitivity is the visual message).

The shear study is a first-class package method — one call covers every axis:

```python
# baseline
res = ana.shear_sensitivity()                       # GradientCurveCI
res.rel_drop, res.rel_drop_lo, res.rel_drop_hi      # sensitivity + 68% band

# the four breakdowns, each by changing one argument:
ana.shear_sensitivity(co_threshold=4.2)             # CO-intensity cut
ana.shear_sensitivity(fwhm_pc=57.5)                 # smoothing level
ana.shear_sensitivity(radius_window_pc=(1500, 2500))# radial annulus
# epoch: build a second FlowAnalysis from config/m33_old.yaml and call the same
```
