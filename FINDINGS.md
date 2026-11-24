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

## Is the |gradient| result just beam smearing of patchy HI? (null test)

**Concern.** The CO-vs-|shear| and CO-vs-|divergence| statistics use absolute
values. If HI is patchy and the beam smears it, velocity gradients are largest
*between* HI clumps (low S/N, beam-mixed spectra), where CO is rare anyway.
That alone would produce "less CO at high |gradient|", with no physics at all.

**What we found, in brief**
1. *The concern is valid.* A forward simulation with **no** link between flows
   and CO (patchy lognormal HI with holes, an independent velocity field,
   observed through the real deprojected 5.3" beam, mom0/mom1 noise matched to
   the data, CO placed by true HI only) reproduces raw drops of 0.26-0.41, as large as
   or larger than observed (0.22 |shear|, 0.39 |div|).
2. *The Laplacian split does not diagnose it.* In those null simulations the
   drop is the same in HI peaks (Lap<0) and troughs (Lap>0), just as in the
   data (fig. panel a). Persistence across the split is therefore not evidence.
3. *An HI-only null does.* A gradient-boosted classifier predicts CO-brightness
   from HI structure alone (column + Laplacian at 0-98 pc, + radius; spatially
   cross-fitted on 300 pc tiles, out-of-fold AUC 0.85). The drop it implies is
   what HI structure, and any smearing coupled to it, can explain. In the null
   simulations the observed-minus-null excess is ~0 at all smoothing scales
   (10 realisations). In a positive control with CO really suppressed by 57 pc
   shear, the excess is large and grows with smoothing (0.08 -> 0.49).
4. **|divergence|: fully explained by HI structure.** Excess 0.02 at native
   resolution, ~0 at 14-86 pc; 94% of the drop is reproduced by the HI-only
   null. The |div| trend cannot be distinguished from the smearing artefact.
   (Signed convergence of CO-bright gas is a different statistic. Observed
   mean div -1.31 vs HI-only null -0.08 km/s/kpc, difference -1.23 [-2.07,
   -0.33], P(>=0)=0.09: weak, but it is not touched by this |.| concern.)
5. **|shear|: modest residual signal.** Excess over the HI-only null is
   positive at every scale and grows with smoothing, like the positive control:
   0.023 (native), 0.042 (14 pc), 0.076 (29 pc), 0.130 (58 pc), 0.166 (86 pc).
   Averaged over 14-58 pc: 0.083 vs null 0.00 +/- 0.057 (z ~ 1.5; 0/10 null
   realisations reach it). Suggestive, not decisive; the null scatter at
   >= 58 pc is large because few independent cells remain.

**Caveats.** The HI-only null is conservative: any real pathway through HI
structure (flows compress HI -> CO) is absorbed into it. The simulation's
"worst-case" noise assigns *all* sub-beam velocity roughness in faint HI to
noise; with radiometer-only noise the artefact nearly vanishes (raw |shear|
drop ~0). Roughly 1/3 of the raw |shear| drop disappears if pixels within
61 pc of blanked velocity are excluded (0.22 -> 0.17), showing mask edges matter.

**Code.** `exgal_flows.nulls`: `hi_structure_features`, `hi_only_co_null`,
`relative_drop_vs_null`, `SmearingSimParams`, `simulate_observation`,
`simulated_gradients` (6 tests in `tests/test_nulls.py`). Figure:
`outputs/beam_smearing_null.png`; tables `outputs/null_test_*.csv`.

### Can a better |gradient| estimator avoid the artefact? (benchmark)

We tested five ways of measuring the gradients at 14-86 pc, using the same null
simulations (no flow-CO link) and positive controls (CO really suppressed by
57 pc shear or divergence):

| estimator | what changes |
|---|---|
| current | `G*(mom0 v)/G*mom0`, then finite difference; NaNs spread through the filter |
| NaN-safe smooth, uniform / mom0² weights | normalised convolution, blanked pixels get zero weight (`flows.normalized_smooth`) |
| plane fit, uniform / mom0² weights | weighted local plane fit in a Gaussian aperture (`flows.local_plane_gradient`); first-order unbiased when weights vary across the aperture, mom0² is inverse-variance for radiometer noise |

Results (`outputs/gradient_estimator_benchmark.png`, `.csv`):
1. **No estimator removes the artefact in the raw statistic.** Null-simulation
   raw |shear| drops are 0.48-0.59 at 14 pc and 0.07-0.31 at 58-86 pc for all
   five estimators. The plane fit with mom0² weights is the least affected at
   small scales (0.48 vs 0.59), but the difference is minor. The artefact is
   driven by noise that scales with S/N, not by the weighting bias.
2. **The fix is the statistic: compare at fixed HI column.** Computing the
   drop within mom0 deciles and recombining (standardised to the high-gradient
   population) gives ~0 in the null simulations (|mean| <= 0.09, scatter
   0.04-0.16) for every estimator, while keeping most of the positive-control
   signal (0.2-0.8). It is a transparent, ML-free version of the HI-only null.
3. **M33 under this statistic** (plane fit, mom0² weights; tile-bootstrap 68%):
   |shear| 0.057 [0.021, 0.091] at 14 pc, 0.106 [0.054, 0.161] at 29 pc,
   0.132 [0.050, 0.212] at 58 pc, 0.075 [-0.043, 0.198] at 86 pc. Positive, but
   about 1 sigma of the null-simulation scatter at each scale. |div|: consistent
   with zero at all scales.
4. **NaN handling matters for what is measured.** The current estimator drops
   pixels near blanked velocity as smoothing grows (1.21M -> 0.31M pixels at
   86 pc), preferentially at HI holes. NaN-safe estimators keep ~1.20M. With
   them, the |shear| excess over the HI-only classifier null is ~0.03-0.07 and
   flat with scale, rather than rising to 0.17. Part of the earlier rise with
   smoothing was this selection.
5. The worst-case-noise simulations give a *larger* raw |shear| drop at
   14-29 pc (0.5-0.6) than M33 shows (0.26-0.30). So the worst-case noise
   model overstates the artefact in the data; treat it as an upper bound.

Not tested: per-pixel noise debiasing of |grad|^2, cube-based (peak or
Gaussian-fit) velocities, and HI-mass weighting of pixels.

### Top third vs bottom third (supersedes the fixed-threshold drops above)

The fixed windows (low < 100, high 100-190 km/s/kpc) select different
populations at each scale. At native resolution they exclude the top 20% of
|shear| entirely; at 58-86 pc the "high" window is only the top 4-8% of
pixels. We therefore redefine the drop by rank, at every scale:

    drop = 1 - f(top third of |gradient|) / f(bottom third)

The "fixed HI column" version stratifies by mom0 decile (functions
`nulls.tercile_windows`, `nulls.quantile_strata`, `nulls.window_drop_boot`).
The null simulations have 6 realisations and the 57 pc positive controls 2
each. Errors are 68% tile-bootstrap intervals. Figure:
`outputs/tercile_drop.png`; table: `outputs/tercile_drop_summary.csv`.

Plane-fit estimator (mom0² weights); native = finite difference:

| scale | \|shear\| raw | null raw | \|shear\| fixed-HI | null fixed-HI (mean ± sd) | positive control |
|---|---|---|---|---|---|
| native | 0.51 | 0.69 ± 0.06 | 0.100 [0.071, 0.131] | 0.02 ± 0.07 | 0.20 |
| 14 pc | 0.37 | 0.55 ± 0.06 | 0.055 [0.014, 0.094] | 0.03 ± 0.07 | 0.21 |
| 29 pc | 0.29 | 0.41 ± 0.06 | 0.025 [-0.028, 0.074] | 0.05 ± 0.08 | 0.27 |
| 58 pc | 0.22 | 0.21 ± 0.07 | 0.012 [-0.049, 0.071] | 0.06 ± 0.10 | 0.36 |
| 86 pc | 0.17 | 0.12 ± 0.09 | 0.032 [-0.031, 0.093] | 0.06 ± 0.10 | 0.38 |

| scale | \|div\| raw | null raw | \|div\| fixed-HI | null fixed-HI | positive control |
|---|---|---|---|---|---|
| native | 0.55 | 0.53 ± 0.08 | 0.121 [0.092, 0.152] | 0.00 ± 0.06 | 0.24 |
| 14 pc | 0.41 | 0.40 ± 0.08 | 0.062 [0.022, 0.103] | 0.01 ± 0.07 | 0.25 |
| 29 pc | 0.30 | 0.29 ± 0.10 | 0.005 [-0.044, 0.055] | 0.01 ± 0.07 | 0.29 |
| 58 pc | 0.18 | 0.11 ± 0.11 | -0.062 [-0.123, 0.005] | 0.02 ± 0.07 | 0.36 |
| 86 pc | 0.14 | 0.03 ± 0.13 | -0.062 [-0.131, 0.010] | 0.04 ± 0.06 | 0.37 |

Reading:
1. Raw tercile drops are large (0.14-0.55) but no larger than in the no-link
   null simulations, so on their own they carry no information.
2. At fixed HI column, neither |shear| nor |div| shows a suppression of CO
   beyond the null scatter at 14-86 pc. The positive controls show what a
   real effect of the injected size looks like (0.2-0.38), and the data are
   well below it. Only at native resolution (finite differences on 3 pc
   pixels, sub-beam and noise-dominated) are both ~0.1, i.e. 1-2 null sigma.
3. At 58-86 pc the |div| drop at fixed HI column goes *negative*: high-|div|
   pixels have more CO than expected (-0.06 with the plane fit, -0.15 to
   -0.24 with the current estimator). This is consistent with the
   signed result: the most converging gas is CO-rich, and |div| mixes it in
   with the most diverging gas.
4. The fixed-threshold shear signal at 29-58 pc (0.11-0.13) does not survive
   the tercile definition, so it came from the high window being a thin tail.
   Conclusion: once HI structure is controlled for, the evidence that
   |shear| suppresses CO on 14-86 pc scales is not significant in these data.

### Signed divergence: converging third vs diverging third

The colliding-flows claim is about the sign of the divergence, so we compare
the most converging third (most negative divergence) with the most diverging
third: contrast = 1 - f(diverging) / f(converging), positive when converging
gas is more CO-rich. The setup is in `scripts/null_test_setup.py`, and the
simulation parameters are in `nulls.m33_calibrated_params`. The positive
control uses CO multiplied by exp(-div_57pc / median|div|) (`inject_signed`).
Figure: `outputs/signed_divergence_tercile.png`; table:
`outputs/signed_divergence_tercile.csv`.

Plane-fit estimator (mom0² weights); native = finite difference:

| scale | f_conv / f_div | contrast | at fixed HI column | null fixed-HI (worst-case, mean ± sd, n=6) | excess over HI-only classifier null |
|---|---|---|---|---|---|
| native | 0.122 / 0.118 | 0.039 [0.006, 0.075] | 0.039 [0.010, 0.069] | 0.01 ± 0.07 | 0.040 [0.014, 0.065] |
| 14 pc | 0.131 / 0.123 | 0.062 [0.013, 0.115] | 0.058 [0.013, 0.101] | 0.01 ± 0.07 | 0.058 [0.021, 0.096] |
| 29 pc | 0.136 / 0.126 | 0.078 [0.020, 0.139] | 0.069 [0.014, 0.120] | 0.00 ± 0.08 | 0.068 [0.022, 0.113] |
| 58 pc | 0.142 / 0.128 | 0.098 [0.029, 0.167] | 0.088 [0.027, 0.148] | 0.00 ± 0.09 | 0.095 [0.046, 0.147] |
| 86 pc | 0.146 / 0.125 | 0.141 [0.068, 0.214] | 0.121 [0.056, 0.184] | 0.01 ± 0.11 | 0.139 [0.090, 0.194] |

Reading:
1. Converging gas is more CO-rich than diverging gas at every scale, and the
   contrast grows with scale, from 0.04 (native) to 0.12-0.14 (86 pc). The
   current estimator gives 0.18-0.20 at 86 pc. This is the opposite scale
   dependence from the |gradient| artefact, which weakens with smoothing.
2. The worst-case-noise null is centred on zero at fixed HI column. Smearing
   does not create a preferred sign, so the signed statistic does not have the
   rectification problem that |div| has.
3. **The significance depends on the error model.** Against the spatial tile
   bootstrap, the excess over the HI-only null excludes zero at every scale
   (about 2.5-3σ at 58-86 pc). Against the realisation-to-realisation scatter
   of the null simulations it is only about 0.7-1.1σ per scale (1.6σ with the
   current estimator at 86 pc). The two radiometer-noise nulls sit at
   +0.04 to +0.08, so a small positive bias cannot yet be excluded (n=2).
4. The positive control (0.62-0.77) is far stronger than the data, so it shows
   the test has power but does not calibrate the expected size of a real effect.

Conclusion: of the three statistics, signed divergence is the only one with
the right sign and the right scale dependence once HI structure is controlled
for, but it is a ~10% contrast whose significance is ~1σ to ~3σ depending on
how spatial correlation is modelled.

#### Calibrating the signed-divergence significance (40 nulls, realistic controls, combined test)

This supersedes the significance statements in the previous subsection.
Same statistic (converging vs diverging third, at fixed HI column, plane-fit
estimator); scales 14, 29, 58 and 86 pc. Figure:
`outputs/signed_divergence_significance.png`; tables:
`outputs/signed_divergence_ensemble.csv`, `signed_divergence_combined.csv`,
`signed_divergence_joint_fit.csv`. The joint bootstrap across scales is
`nulls.window_drop_joint`, and the GLS helper is `nulls.gls_amplitude`.

1. **Null ensemble (30 worst-case + 10 radiometer noise).** The means are
   -0.015 to -0.03 at every scale, so there is no positive bias. The +0.04 to
   +0.08 seen in the first two radiometer runs was a fluctuation. The scatter
   between realisations is 0.08-0.11, and it is **1.28x** the tile-bootstrap
   error measured *within* a simulation. The 300 pc tile bootstrap therefore
   underestimates the error by ~28%, because it misses correlations on scales
   larger than a tile. The earlier "2.5-3 sigma" (the excess over the
   classifier null with the raw tile bootstrap) was over-stated for this
   reason.
2. **Per scale.** The fraction of null realisations at or above M33 is 0.20,
   0.17, 0.17 and 0.07 at 14, 29, 58 and 86 pc.
3. **Combined over scales** (empirical p from the 40 nulls):
   - mean contrast over 14-86 pc: 0.084 vs null -0.018 +/- 0.094, z = 1.1,
     p = 0.15;
   - rise with scale (linear slope over 14-86 pc, per 86 pc): 0.074 vs null
     -0.007 +/- 0.054, z = 1.5, p = 0.05 (1 of 40 nulls exceeds it).
   - A GLS amplitude fit with the null covariance gives a *negative* constant
     amplitude. The scales are so strongly correlated (r ~ 0.9) that GLS
     weights the *differences* between scales, so for this data it is not a
     usable summary. The model-free mean and slope are reported instead.
4. **Power.** In realistic positive controls the CO emissivity is scaled by
   exp(-div_57pc / (24 x median|div|)), giving a contrast of 0.05-0.07 with
   20 realisations, slightly *weaker* than M33's 0.06-0.12. These are
   detected at p < 0.05 only **23%** (mean) and **27%** (slope) of the time.
   The realisation-to-realisation scatter (0.07-0.10) is as large as the
   effect, so one galaxy of M33's size cannot decide it.

Conclusion: M33 shows the expected sign (converging gas more CO-rich) and a
contrast that rises toward 86 pc, but once the error is calibrated on 40 null
galaxies the evidence is ~1-1.5 sigma (p = 0.15 for the level, p = 0.05 for
the rise). This is consistent with a real effect of roughly the observed
size, which the test would detect only about a quarter of the time. A firm
result needs more independent area: the LMC, NGC 6822, or stacking galaxies.

#### Simulations vs block bootstrap: which error to trust?

Maps: `outputs/simulation_maps.png` (3.1 kpc cutouts: M33, a null
realisation, and a control with a strong injected effect).

They measure different things. The block bootstrap estimates the scatter of
M33's own statistic by resampling its tiles. It assumes the tiles are
independent and interchangeable, and it says nothing about bias. The
simulations give the full distribution under "no flow-CO link", including
any bias from smearing (large for |div|, ~0 for the signed statistic), but
only for galaxies that look like the model.

Checks:
1. **Block size, in the data.** The bootstrap error at 86 pc grows with
   block size: 0.054 (150 pc), 0.066 (300), 0.073 (600), 0.082 (1000), 0.077
   (1500). So M33 has correlations beyond 300 pc; at >= 1000 pc blocks the
   error is ~1.2x the 300 pc value. This is independent of any simulation and
   agrees with the 1.28x inflation.
2. **Bootstrap vs truth, in the simulations.** Within a simulation the
   bootstrap error plateaus at ~600 pc (0.089 at 86 pc), but the true scatter
   across realisations is 0.116. The bootstrap still misses variance from
   structure on the scale of the whole box.
3. **Area mismatch.** The simulated galaxies have ~70% of M33's unblanked area
   (median 420 vs 612 300-pc tiles), so their scatter is ~1.2x too large for
   M33. Rescaled to M33's area, the null sd at 14-86 pc is 0.064-0.096, and the
   combined tests become: mean z = 1.3, p = 0.10; slope z = 1.8, p = 0.024.
4. **Morphology mismatch (from the maps).** The simulated HI is smoother and
   dominated by a few large holes, whereas M33's HI is filamentary with many
   ~200-500 pc holes. The simulated divergence noise is more anisotropic
   (streaky along y), because worst-case mom1 noise is correlated along the
   stretched deprojected beam. More large-scale power means more scatter
   between realisations, so if anything the simulations overstate the error.

Bottom line: the 300 pc bootstrap is too small (shown by the data), and the
simulations as run are somewhat too large (area and morphology). The best
estimate of the scatter lies between, close to the >= 1000 pc block bootstrap
and the area-rescaled simulations, which agree (0.077-0.082 vs 0.096 at 86 pc).
Signed-divergence significance: ~1.3-1.8 sigma.

#### Wedge blocks and excluding the centre

Errors here come from fixed-angle azimuthal wedges around the M33 centre
(`nulls.wedge_index`), using a delete-one-wedge jackknife
(`window_drop_joint(..., resample="jackknife")`). The jackknife is more
stable than a bootstrap with only 8-24 groups. Figure:
`outputs/signed_divergence_wedges.png`; tables: `signed_divergence_wedges.csv`,
`signed_divergence_wedge_nulls.csv`.

**The signal is in the outer disk.** Converging-vs-diverging CO contrast at
fixed HI column, plane-fit estimator, 12-wedge jackknife sd in brackets:

| sample | 14 pc | 29 pc | 58 pc | 86 pc |
|---|---|---|---|---|
| full disk | 0.058 (0.050) | 0.069 (0.057) | 0.088 (0.062) | 0.121 (0.065) |
| R > 0.9 kpc | 0.111 | 0.126 | 0.158 | 0.201 |
| **R > 1.2 kpc** | **0.134 (0.037)** | **0.153 (0.040)** | **0.194 (0.037)** | **0.248 (0.034)** |
| R > 1.5 kpc | 0.112 | 0.128 | 0.165 | 0.215 |
| R < 1.2 kpc | -0.056 (0.08) | -0.058 (0.10) | -0.047 (0.11) | -0.044 (0.11) |

- The inner 1.2 kpc holds 9% of the pixels but 24% of the CO-bright ones
  (36% CO-bright vs 12% outside). Its contrast is slightly negative but
  consistent with zero. Because the drop is a ratio of CO counts, this CO-rich
  region carries a large weight and pulls the full-disk value down.
- Neither the tercile boundaries nor an extra radius stratification explains
  the difference. Outer pixels give 0.25 at 86 pc with full-disk or
  outer-only boundaries, and stratifying the full disk by mom0 x radius still
  gives 0.12. The inner and outer disk genuinely differ, as in the earlier
  radial reversal of mean divergence.
- The result is not fine-tuned to 1.2 kpc: cuts at 0.9 and 1.5 kpc give 0.20
  and 0.215 at 86 pc.

**Significance for R > 1.2 kpc**
- *Null simulations (same annulus and wedges, 40 realisations).* None reaches
  M33 at any scale. Mean over 14-86 pc: z = 2.1 (as run), 2.6 (rescaled to
  M33 area). Rise with scale: z = 2.7 / 3.2. p < 0.025, limited by the number
  of realisations.
- *Wedge jackknife on the data.* Mean over scales 0.182 +/- 0.035-0.054
  (8-24 wedges). In the simulations the wedge jackknife under-estimates the
  true scatter by 1.12-1.29x. After that correction z = 2.8-4.3 depending on
  wedge count (4.0 for 12 wedges).
- M33's outer disk is about half as variable as the simulated galaxies
  (wedge jackknife 0.034-0.040 vs 0.066-0.104 in the simulations), consistent
  with the morphology mismatch. The simulation-based z is the conservative
  bound.

**Caveat: the cut was chosen after looking.** Excluding the centre was
motivated by the radial reversal found earlier in these same data, and
it doubles the signal. The full-disk statistic, which needs no choice, remains
~1.1-1.8 sigma. A clean confirmation needs the R > 1.2 kpc, signed-tercile,
fixed-HI test specified in advance and applied to independent data: the old
8.2" M33 epoch is not independent; NGC 6822 or the LMC would be.

#### |shear| and |div| over the same outer-disk area

Method: the same test (top vs bottom third, fixed HI column, R > 1.2 kpc,
12-wedge jackknife, 40 nulls under the same geometry). Figure:
`outputs/same_area_gradients.png`; tables: `same_area_gradients.csv`,
`same_area_div_decomposition.csv`, `same_area_abs_nulls.csv`.

| R > 1.2 kpc | 14 pc | 29 pc | 58 pc | 86 pc | mean, z vs nulls |
|---|---|---|---|---|---|
| \|shear\| (top vs bottom third) | 0.031 | -0.014 | -0.028 | -0.001 | -0.003, z = 0.0 |
| \|div\| (top vs bottom third) | 0.050 | 0.002 | -0.054 | -0.032 | -0.009, z = -0.2 |
| signed div (converging vs diverging third) | 0.134 | 0.153 | 0.194 | 0.248 | 0.182, z = 2.6 |

**|shear|.** Nothing, in either region: the full disk gives 0.03 (z = 0.6) and
the outer disk 0.00. At fixed HI column we find no sign of shear suppressing
CO on 14-86 pc scales.

**|div| is ~0 because the signed effect cancels.** Splitting the top-|div| third
by sign (outer disk, relative to the low-|div| third at fixed HI column; error
bars are wedge jackknife values):
- strongly *converging* gas: +0.06, +0.14, +0.23, +0.24 (+/- 0.10-0.15) more
  often CO-bright;
- strongly *diverging* gas: -0.16, -0.15, -0.14, -0.20 (+/- 0.05-0.07) less
  often CO-bright.

So CO fraction falls monotonically from converging, through quiescent, to
diverging gas. Taking |div| adds the two tails together, and they cancel. The
deficit in diverging gas is the more tightly measured half (~2.5-3 sigma per
scale by jackknife). Signed divergence is the right statistic, and |div| is
not.

Note: on the data, the wedge jackknife for the signed statistic in the outer
disk (0.034-0.040) is about half the jackknife for |div| and |shear|
(0.064-0.115). In the simulations the two are similar (~0.07-0.10). We have not
explained why the signed result is so consistent from wedge to wedge.

#### Does the |div| artefact dilute the signed asymmetry? (paired test)

**Hypothesis.** Whatever spuriously suppresses CO at high |div| (beam smearing,
low S/N in faint HI) would also drive down the measured converging-vs-diverging
asymmetry. Two routes:
- *Symmetric suppression.* Both tails lose CO equally. This does not reduce
  the asymmetry: a multiplicative loss cancels in f_div/f_conv, and an
  additive one *inflates* the contrast.
- *Sign scrambling.* Velocity noise flips the sign of the divergence in faint
  gas, which dilutes any real asymmetry toward zero.

**Test.** Positive-control simulations (CO depends on the true 57 pc signed
divergence; strengths s = 24 and 6; 6 seeds each). Each seed was re-observed
with *no* velocity noise, radiometer noise and worst-case noise, so the true
galaxy and true CO are identical across noise models. The recovered contrast
was compared with the contrast computed from the true velocity field at the
same scale (29, 58, 86 pc). Figure: `outputs/asymmetry_attenuation_test.png`;
table: `attenuation_paired_sims.csv`.

**Result: no meaningful dilution at 29-86 pc.** The slope of recovered vs true
contrast is 1.04 (no noise), 1.04 (radiometer) and 1.03 (worst case). Worst-case
noise lowers the weak-effect contrast by ~0.003-0.006 (about 20% of a 0.02
effect), and the strong effect by <3%. With the plane-fit estimator at
>= 29 pc, velocity noise is averaged down far enough that the sign of the
divergence is rarely flipped. The |div| drop in these controls stays small
(0.01-0.04), as expected when the effect cancels. The scrambling route could
matter at native and 14 pc resolution, which we did not test here.

**Matched-|div| statistics on the data** (R > 1.2 kpc, fixed HI column;
these cancel any symmetric |div| effect by construction):

| statistic | 29 pc | 58 pc | 86 pc | z vs nulls (area-rescaled, n=30) |
|---|---|---|---|---|
| converging vs diverging third | 0.153 | 0.194 | 0.248 | 2.9 / 3.1 / 3.3 |
| sign, within top |div| third | 0.245 | 0.297 | 0.347 | 2.9 / 2.9 / 2.9 |
| sign, at matched |div| (all thirds) | 0.106 | 0.139 | 0.171 | 2.7 / 3.0 / 3.1 |

Comparing converging with diverging gas at the *same* |div| gives the same
answer, and no null realisation reaches the data for any of these
statistics. Table: `matched_div_sign_test.csv`. (These z values are per scale,
against 30 nulls; the R > 1.2 kpc caveat above still applies.)

#### Asymmetry vs |div| (quintiles): is a high-|div| suppression hiding signal?

**Concern.** The asymmetry is larger in the top |div| third (0.25-0.35) than
averaged over all |div| (0.11-0.17). If CO is suppressed at high |div|, might
the true signed signal be even larger?

**Two reasons it is not.**
1. The asymmetry is a *ratio within a |div| bin*, 1 - f(div>0)/f(div<0). A
   suppression that acts equally on both signs cancels if it is multiplicative,
   and enlarges the contrast if it is additive. Only a raw difference
   f_conv - f_div would be shrunk by it.
2. In M33 there is no common-mode suppression left to act. At fixed HI
   column (R > 1.2 kpc), the CO-bright fraction in each |div| quintile relative
   to the lowest quintile is -0.02 to +0.07, consistent with zero (wedge
   jackknife errors 0.02-0.16).

**Dose-response.** The asymmetry rises steadily with |div|. At 86 pc it goes
0.02, 0.03, 0.17, 0.27, 0.39 (+/- 0.05-0.07) from the weakest to the
strongest quintile; at 58 pc 0.03 -> 0.37; at 29 pc 0.02 -> 0.31. The
all-thirds average is lower because the weakest |div| bins, where the sign of
the divergence carries little information, dilute it. This is what a real
flow-CO link predicts.

**Simulation check.** With a strong injected effect at 57 pc and 6 paired
galaxies, the per-quintile asymmetry curve is the same for the true velocity
field, for observation with no velocity noise and for worst-case noise
(top quintile 0.22 / 0.24 / 0.24). The artefact does not flatten the
dose-response. Figure: `outputs/div_quintile_asymmetry.png`; table:
`div_quintile_asymmetry.csv`.

#### Toward the true amplitude of convergence onto molecular gas (R > 1.2 kpc)

**What "divergence" is in this code.** `geometry.deproject` multiplies the
line-of-sight residual velocity by cos i, and `flows.velocity_gradients`
divides by cos i. So `divergence` = d v_los / d x, with x the in-plane
minor-axis coordinate (km/s/pc). For in-plane motion v_los = sin i u_x, so
the in-plane gradient component is **du_x/dx = divergence / sin i
(= 1.22 x divergence)**. (The cos i factors cancel and do not deproject the
velocity.)

**Orientation.** Only du_x/dx is observable, because u_y (along the major axis)
does not project onto the line of sight. At deprojected azimuth theta (from
the minor axis), x = cos theta e_R - sin theta e_phi, so

    du_x/dx(theta) = A + B cos 2theta + C sin 2theta,
    A = (div u)/2,  B = (G_RR - G_pp)/2,  C = -(G_Rp + G_pR)/2,

where G_RR = du_R/dR, G_pp is the azimuthal term and G_Rp, G_pR are the
cross terms. Fitting this across 24 wedges gives the full in-plane
convergence **div u = 2A** without assuming the flows are isotropic. It also
separates radial (A + B) from azimuthal (A - B) convergence. The only
assumption is that flow statistics do not depend on azimuth.

**Estimator.** For each wedge we take the mean du_x/dx at CO-bright pixels
minus the mean over all pixels in the same wedge x ring (1.2-2, 2-3, 3-4 kpc)
x HI-decile stratum. Comparing within the same location removes
rotation-model errors (which themselves produce cos 2theta patterns) and
anything tied to HI column. The fit is CO-count weighted; errors are 24-wedge
jackknife values. In 20 null simulations the estimator has a small positive
bias, A = +0.33 to +0.45 +/- 0.9-1.25, which is subtracted below. Table:
`outputs/convergence_amplitude.csv`; figure: `convergence_amplitude.png`.

| scale | div u at CO-bright gas (km/s/kpc) | z | radial A+B | azimuthal A-B | 1/|div u| (Myr) | |div u| x L (km/s) |
|---|---|---|---|---|---|---|
| native (22 pc beam) | -4.5 +/- 2.6 | 1.7 | +0.7 | -5.2 | 220 | 0.10 |
| 14 pc | -4.5 +/- 2.6 | 1.8 | +0.4 | -5.0 | 215 | 0.07 |
| 29 pc | -4.7 +/- 2.5 | 1.9 | +0.1 | -4.8 | 205 | 0.14 |
| 58 pc | -5.3 +/- 2.5 | 2.1 | -0.5 | -4.8 | 185 | 0.30 |
| 86 pc | -5.5 +/- 2.4 | 2.3 | -1.0 | -4.5 | 180 | 0.47 |

(Error = max(1.2 x jackknife, null scatter); radial/azimuthal errors ~1.7-2.1.)

Reading:
- CO-bright gas sits in net converging HI flows of ~5 km/s/kpc. That is a
  compression time of ~200 Myr, or a ~0.1-0.5 km/s velocity difference across
  the smoothing scale. This is ~2 sigma, weaker than the CO-fraction tests,
  because mean velocity gradients have heavy tails while the tercile
  statistics are rank-based.
- The convergence appears to be mostly **azimuthal** (A - B ~ -5, ~2.5
  sigma), with little radial convergence (A + B ~ 0). This is tentative,
  since B is only ~2 sigma. It suggests gas piling up along orbits rather
  than radial inflow onto clouds.

**What still separates this from the "true" amplitude**
1. *Resolution and smoothing.* We measure the convergence averaged over the
   beam plus smoothing kernel. If inflow is concentrated on cloud scales
   (~10-30 pc), the true cloud-scale rate is higher. This is a lower limit.
2. *Dilution.* All CO-bright pixels are averaged, accreting or not. This is
   the population mean.
3. *Reference level.* Convergence is measured relative to HI of the same
   column at the same location. Any convergence shared by all dense HI is
   absorbed into the reference (conservative).
4. *Tracer.* mom1 is HI-intensity weighted, so this is the HI envelope's flow,
   not the molecular gas's.
5. *Vertical motions.* cos i x u_z also contributes to v_los. Assumed zero.
6. *Azimuthal symmetry.* M33's spiral arms are an m = 2 pattern, which could
   alias into the cos 2theta term. Needs checking.

The next step for (1)-(3) is a forward-model calibration with 2-D *vector*
velocity fields: generate u = (u_x, u_y) with a known compressive part, place
CO where the true div u is converging, observe v_los = sin i u_x through the
beam and noise, and apply this estimator. The recovered/true ratio at each
scale is the correction factor. The current simulations use a scalar
line-of-sight field, so they cannot test this.

#### Calibrating the convergence amplitude with vector-flow simulations

This supersedes the radial/azimuthal interpretation in the previous
subsection.

**Model.** `nulls.simulate_vector_observation` builds an in-plane vector velocity field
u = grad(phi) + curl(psi) with a chosen compressive fraction (0.5 by
default) on the calibrated patchy HI. It observes v_los = sin i u_x through the
real deprojected beam, with mom1 noise and blanking matched to M33, places CO
by HI and optionally boosts CO where the true div u (or only G_RR or G_phiphi)
is converging, on a chosen scale. The *same* estimator
(`nulls.azimuthal_excess` + `nulls.azimuthal_fit`) is run on it, and the
truth (HI-matched excess of the true div u at CO-bright gas) is computed
alongside. Driver: `scripts/convergence_calibration.py`, with 54 simulated
galaxies over 11 configurations. Tables: `outputs/convergence_calibration_sims.csv`,
`convergence_calibration_summary.csv`, `convergence_amplitude_calibrated.csv`;
figure: `outputs/convergence_calibration.png`. Tests: 3 new (40 pass).

**What the simulations show**
1. *Bias.* The estimator reads ~0.6-0.8 km/s/kpc too far toward divergence
   (mean over all sims, +/- 0.3). Its scatter per galaxy is 1.6-2.2, comparable
   to the data's calibrated jackknife error, so the error bars are realistic.
2. *Effect on the measured scale or larger* (57 or 150 pc, mixed flows):
   recovered/true = 0.9-1.1. The measurement is essentially unbiased.
3. *Effect on cloud scales (20 pc):* recovered/true = 0.20-0.30. The beam
   (22 pc on the sky, ~38 pc along x after deprojection) washes it out.
4. *Compressive fraction matters:* 1.0 at 0.25, ~1.0 at 0.5, 0.46-0.69 at
   0.75. With mostly compressive flows, d u_x/dx has its power along k_x,
   which is the direction the deprojected beam smooths most.
5. *Radial vs azimuthal split (B) is unbiased but noisy:* per-galaxy scatter
   1.6-2.0. M33's B = 2.0 is only ~1.1 sigma, so **we cannot tell whether the
   convergence is radial, azimuthal or isotropic.** (The earlier "mostly
   azimuthal" reading used the jackknife error alone.) Single-component
   injections are recovered at 0.6-0.8 in div u.
6. Large-scale m = 1, 2 residual-velocity patterns of 5 km/s (rotation-model
   errors, spiral streaming) do not bias A beyond the noise, because each
   wedge x ring stratum subtracts its own mean.
7. *Scale of the effect.* In M33 the recovered convergence *strengthens* from
   native to 86 pc (change -1.24 +/- 0.94). In the simulations it is flat
   for a 150 pc effect (+0.03 +/- 0.54), weakens for 57 pc (+0.4 to +0.6) and
   is noisy for 20 pc (+0.4 +/- 1.5). This weakly favours a CO-linked
   convergence on scales of >~ 100 pc, not cloud-scale, at ~1-2 sigma.

**Calibrated amplitude (R > 1.2 kpc, in-plane div u at CO-bright gas, relative to HI-matched gas)**

| measured at | debiased div u | z | if compressive 0.75 | if 20 pc effect | 1/|div u| |
|---|---|---|---|---|---|
| native | -4.4 +/- 2.6 | 1.7 | -5.2 | -12 | 220 Myr |
| 29 pc | -4.6 +/- 2.5 | 1.8 | -6.1 | -15 | 210 Myr |
| 58 pc | -5.2 +/- 2.5 | 2.1 | -8.3 | -20 | 190 Myr |
| 86 pc | -5.4 +/- 2.4 | 2.3 | -10.6 | -24 | 180 Myr |

Bottom line: the best estimate is **div u ~ -5 +/- 2.5 km/s/kpc on
~60-90 pc scales** (compression time ~200 Myr), assuming the effect lives on
those scales or larger and the flows are mixed (consistent with the scale
trend). Systematics push it more negative, never less: up to ~ -10 if the
flows are mostly compressive, ~ -15 to -24 if the convergence is really
concentrated on ~20 pc cloud scales (weakly disfavoured). The HI-matched
reference (convergence shared by all dense HI is subtracted), dilution by
non-accreting CO, and the HI-weighted tracer all also make this a lower
limit. Vertical motions are assumed zero.
