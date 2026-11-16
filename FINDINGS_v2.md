# Findings v2: rerun with the 8" CO(2-1) integrated-intensity map

This note summarises `notebooks/reproduce_results_check_co_mom0.ipynb`
(branch `check`) and lists the observational biases that still need testing.
It supplements `FINDINGS.md`; it does not replace it.

## What was run

| Input | File |
|-------|------|
| HI (new) | `HI/lglbs_ABCD/m33_A+B+C+D+GBT_hi21cm_strict_mom{0,1}.fits` (5.3", 3.05 pc/pix) |
| HI (old, shear epoch panel only) | `HI/17B/M33_14B_17B_..._GBT_feathered_K.mom{0,1}.fits` (8.2") |
| Rotation model | `rot_model/rad.mod.fits` (DiskFit, peak velocities, no warp, no radial flow) |
| CO | `CO/m33_all_co21_2p5kms_8p0as_aligned_strict_mom0.fits` (8", K km/s) |
| Config | `config/m33.yaml` |

Two things to know about the CO selection in this run:

- The notebook loads `config/m33.yaml`, so the threshold is the peak-map value
  of 0.6, now applied to a map in K km/s. The bins are therefore 1.2, 1.8, 2.4
  and 4.2 K km/s on the `cos i`-scaled map (2.1, 3.1, 4.2 and 7.3 K km/s on the
  sky). `config/m33_co_mom0.yaml` (threshold 2.0) exists but is commented out.
- The CO-intensity panel of the shear figure hard-codes `co_threshold=2.0 * f`
  (4, 6, 8, 14 K km/s), so that one panel uses different cuts from everything
  else. Its `CO>7x` curve has too few pixels to read.

## Results

"7m peak" is the same notebook run with the original
`m33_7m_co21_2p5kms_pbcorr_round.image.max.fits` map
(`notebooks/reproduce_results_check.ipynb`), which reproduces `FINDINGS.md`.

**Median divergence of CO-bright pixels** (km/s/kpc, r < 4000 pc, `decim=2`,
68% wedge-bootstrap interval; negative = converging):

| CO bin | 8" mom0 | 7m peak |
|--------|--------:|--------:|
| 2x-3x | -2.53 [-4.40, -0.55] | -1.78 [-3.77, +0.07] |
| 3x-4x | -2.25 [-4.65, -0.29] | -2.87 [-5.59, -0.69] |
| 4x-7x | -1.48 [-4.35, +1.53] | -2.14 [-4.78, +0.57] |
| > 7x  | -3.21 [-6.08, -0.39] | -8.78 [-13.29, -2.70] |

- The sign is negative in every bin with both CO maps.
- The steepening with CO brightness is not seen with the mom0 map: all four
  bins agree within their intervals at -1.5 to -3.2.

**Radial trend of CO-bright divergence** (CO > 2x, tiled block bootstrap):

| Fit | 8" mom0 | 7m peak |
|-----|--------:|--------:|
| OLS intercept `a` [km/s/kpc] | +5.50 [0.54, 10.45] | +5.07 [0.81, 9.34] |
| OLS slope `b` [km/s/kpc/kpc] | -4.32 [-6.60, -1.79] | -3.87 [-5.85, -2.08] |
| OLS zero crossing | 1274 pc | 1308 pc |
| L1 intercept `a` | +4.50 [0.51, 8.01] | +2.06 [-1.05, 5.13] |
| L1 slope `b` | -3.69 [-5.37, -1.69] | -2.28 [-3.81, -0.94] |
| L1 zero crossing | 1218 pc | 904 pc |
| Sigmoid inner plateau | +7.0 [2.4, 14.2] | +6.8 [4.4, 20.5] |
| Sigmoid outer plateau | -6.1 [-10.8, -4.4] | -5.5 [-9.8, -4.2] |
| Sigmoid transition `r0` | 1150 pc | 1123 pc |
| `straddle_zero_prob` | 1.000 | 1.000 |

- The sign reversal near 1.2 kpc is reproduced with the mom0 map, at almost
  the same radius and amplitude.
- The binned medians are positive in all four annuli inside 1200 pc, negative
  from 1200 to 3200 pc, and near zero in the 3200-4000 pc annulus.
- Sample: 11751 CO-bright pixels in 231 tiles at `decim=4` (7m peak: 16971
  pixels in 267 tiles).

**Shear** (CO > 1x threshold):

- The CO-bright fraction falls from about 0.10 at low |shear| to about 0.05 at
  0.2 km/s/pc.
- Sensitivity `rel_drop` = 0.24 [0.22, 0.26] (7m peak: 0.23 [0.21, 0.25]).
- The old 8.2" HI epoch falls slightly more steeply than the new one, as in
  `FINDINGS.md`. All four radial annuli show the decline.

**Resolution study** (CO-bright fraction against divergence):

- At native resolution the fraction peaks near zero divergence (about 0.105)
  and falls to about 0.035 at both -200 and +200 km/s/kpc.
- This two-sided fall is larger than the fall with |shear|. See bias 4 below.
- The curves at 58 and 86 pc smoothing are noisy and the outer bins are empty.

## Observational biases and the tests that would check them

Ordered by how much each could change the conclusions. "Control" means the
median divergence of all pixels with CO coverage in the same annulus.

| # | Possible bias | Test with the data |
|---|---------------|--------------------|
| 1 | **Sign convention.** Only the line-of-sight velocity is used, so whether negative `dv/dx` means converging depends on which side of M33 is the near side after rotation. The code does not check this. | Establish the near side independently (dust lanes, or the sense of rotation with trailing arms) and confirm the sign of `dv/dx` for a known inflow. |
| 2 | **Large-scale residuals.** The rotation model has one PA and inclination, no radial flow, and was fit to peak velocities but subtracted from moment-1. Each leaves a smooth gradient of 1-2 km/s/kpc that varies with radius, which could produce the 1.2 kpc reversal in all gas. | (a) Subtract the control from the CO medians in each annulus. (b) Repeat with the high-pass field (`set_velocity_mode("highpass")`) and with no subtraction. (c) Repeat with DiskFit models that include radial flow and a warp. |
| 3 | **Only one divergence term is measured.** `dv/dx` sees convergence along the minor-axis direction only, so sensitivity depends on azimuth and on arm orientation. | Split the CO-bright medians into wedges about the minor and major axes. Arm-driven convergence should be strongest near the minor axis; a `sin 2θ` pattern points to bias 2. |
| 4 | **Noisy gradients mimic suppression.** Large |gradient| pixels are likely low signal-to-noise HI, where CO is also rarer. The CO fraction falls with |divergence| on both sides at least as fast as with |shear|, which suggests the decline is not specific to shear. | (a) Compare the |shear| and |divergence| curves directly. (b) Repeat after cutting on the moment-1 error map (`strict_emom1`). (c) Repeat within narrow bins of HI moment-0. |
| 5 | **Chance alignment.** CO and the gradient field share large-scale structure (arms, radius). | Null test: mirror or rotate the CO map by 180° about the centre and rerun. Radial distributions are preserved, local correspondence is not, so the CO medians should fall to the control. |
| 6 | **CO threshold.** The cut is a fixed value tuned on another map, in different units, and scaled by `cos i`. Fixed cuts on a map with varying noise select by sensitivity. | Define CO-bright by signal-to-noise using `strict_emom0`, and repeat with the 7m peak map, the 8" `tpeak` map and the 8" mom0 map at matched pixel fractions. |
| 7 | **CO coverage in the denominator.** `result_co_vs_shear` and `resolution_study` mask on radius only, so pixels without CO coverage count as CO-dark. | Recompute both with a finite-CO mask and compare the fractions. |
| 8 | **Resolution mismatch.** CO is 8" (33 pc), new HI is 5.3" (22 pc), and deprojection stretches the beam 1.74x along x. Smoothing kernels of 14 and 29 pc are at or below the HI beam. | Use the 8.2" HI epoch as the matched-resolution case for all results, not only the shear panel. Quote smoothing scales as the quadrature sum with the beam. |
| 9 | **Mask-edge artefacts.** Blanked pixels are zero-filled before interpolation, which creates false gradients at footprint edges and moment-1 holes. | Rerun with `edge_erode_pix=3` and compare. |
| 10 | **Moment-1 as the velocity tracer.** Multiple components, self-absorption and feedback shells shift moment-1 where CO is bright, and vertical motions enter with weight `cos i`. | (a) Repeat with an HI peak-velocity or fitted-centroid map. (b) Mask known HII regions and shells and rerun. (c) Compare CO moment-1 with HI moment-1 at CO-bright pixels. |

Amplitudes are line-of-sight gradients. In-plane values are larger by
`1/sin i` (about 1.22); signs are unaffected.

## Status of the checks

- Test 2(b): run on all three HI epochs with the 8" mom0 map, see below. The
  rotation model is not needed for the divergence results. For shear, some
  large-scale subtraction is needed, but a high-pass works as well as the
  model. An earlier run with the 8" `broad_tpeak` map is superseded (see the
  edge-rim note under Test 2(b)).
- CO coverage-edge rim in `broad_tpeak`: found and quantified at 8.2" in
  rotation-model mode. `strict_mom0` and the 7m peak map are not affected.
  The other epochs were not rerun with `tpeak`.
- Test 9: `edge_erode_pix` is implemented but not yet compared.
- Test 5 (chance alignment): run, see below.
- All other tests are not started.

### Test 5: chance-alignment null (8" CO mom0, `decim=2`)

`FlowAnalysis.null_test_chance_alignment(mode)` reruns the CO-bright results
with the deprojected CO map displaced (`rot180`, `mirror_major`,
`mirror_minor`) and compares with the control (median divergence of all
pixels with CO coverage in the same annulus, point estimate). Driven from
`notebooks/reproduce_results_check_co_mom0.ipynb`. One realisation of each
displacement, so there is no p-value.

**Median divergence per CO bin** [km/s/kpc]. Control median over the
aperture is -1.7.

| CO bin | real | `rot180` | `mirror_major` | `mirror_minor` |
|--------|-----:|---------:|---------------:|---------------:|
| 2x-3x | -2.5 [-4.4, -0.6] | -2.4 [-3.9, -0.8] | -2.0 [-3.7, -0.1] | -0.7 [-3.4, +1.7] |
| 3x-4x | -2.3 [-4.7, -0.3] | -2.0 [-3.9, +0.3] | +1.3 [-0.9, +3.2] | -3.4 [-5.5, -1.1] |
| 4x-7x | -1.5 [-4.4, +1.5] | -1.2 [-3.4, +1.0] | -0.4 [-2.5, +2.2] | -1.0 [-4.2, +2.2] |
| > 7x  | -3.2 [-6.1, -0.4] | +2.5 [-2.2, +7.0] | +6.6 [+3.8, +10.8] | -3.7 [-7.1, 0.0] |

**Radial profile, median minus control** [km/s/kpc]:

| Annulus [pc] | control | real | `rot180` | `mirror_major` | `mirror_minor` |
|--------------|--------:|-----:|---------:|---------------:|---------------:|
| 0-300 | +3.3 | -0.5 | -15.6 | -2.3 | -6.0 |
| 300-600 | -4.5 | +9.6 | -1.6 | -2.9 | +4.6 |
| 600-900 | +5.4 | -0.4 | +8.1 | +6.7 | -0.2 |
| 900-1200 | +2.5 | +3.9 | +0.5 | +6.6 | +4.1 |
| 1200-1500 | -8.6 | +0.1 | +3.8 | +2.3 | -18.3 |
| 1500-2000 | -1.4 | -2.3 | +0.9 | +1.3 | +4.1 |
| 2000-2500 | -4.0 | +0.5 | -2.4 | 0.0 | +3.3 |
| 2500-3200 | -0.9 | -10.3 | -1.2 | -0.8 | +0.1 |
| 3200-4000 | -1.2 | +1.8 | +0.3 | -0.3 | -8.3 |

Assessment:

- The local CO-to-flow correspondence is not isolated by this test. The
  radial sign reversal near 1.2 kpc is in the control (+2.5 at 0.9-1.2 kpc,
  -8.6 at 1.2-1.5 kpc) and in the displaced maps. Real minus control at
  1.2-1.5 kpc is +0.1. Test 2(b) below shows the step survives with no
  rotation model, so it is not a model residual; it is a feature of the HI
  velocity field that CO-covered gas shares.
- Real minus control has an RMS of 5.0 km/s/kpc over the nine annuli. The
  three displaced maps give 6.1, 3.5 and 7.5.
- Only two annuli stand out against the nulls: 2.5-3.2 kpc (real -10.3,
  nulls -1.2 to +0.1) and 300-600 pc (real +9.6, nulls -2.9 to +4.6). Three
  nulls and no error on the control are not enough to call either
  significant.
- The negative sign in the 2x-4x bins is reproduced by `rot180`, which keeps
  the radial distribution of CO. The real bins sit within about 1 km/s/kpc of
  the aperture control. Only the > 7x bin separates, and only partly:
  `mirror_minor` reproduces the real value and the `rot180` interval overlaps
  it.
- The control medians swing by 4-9 km/s/kpc between annuli, larger than the
  1-2 km/s/kpc quoted for bias 2.
- The aperture control is not radially matched to CO-bright pixels. Use
  `rot180` as the matched comparison for the bins.
- Next: many displacements (arbitrary angles, random shifts) to build a null
  distribution, tests 2(a) and 2(b), and an error band on the control.

### Test 2(b): is the rotation model needed? (`notebooks/highpass_check.ipynb`)

Four large-scale subtractions on three HI epochs: rotation model (default),
none (moment-1 as observed), and a Gaussian high-pass of 500 pc and 1000 pc.
Same `edge_erode_pix=3` for every mode. CO is the 8" `strict_mom0` map
(K km/s) with a base threshold of 0.1 on the `cos i`-scaled map, so the bins
are 0.2, 0.3, 0.4 and 0.7 on that map and the numbers are not comparable with
the main run (threshold 0.6) or with Test 5. `decim` is 1, 3 and 4 for the
21.4", 8.2" and 5.3" epochs. Results are cached in `outputs/highpass_check/`.

| Epoch | Subtraction | All-pixel median | CO > 7x median | Shear `rel_drop` | Turnover [pc] |
|-------|-------------|-----------------:|---------------:|-----------------:|--------------:|
| 21.4" | rotation model | -0.12 | -1.0 [-2.0, +0.2] | 0.45 [0.37, 0.53] | 465 |
| 21.4" | none | +0.04 | -2.0 [-3.8, +0.1] | 0.15 [0.06, 0.25] | 404 |
| 21.4" | high-pass 500 pc | -0.24 | -1.0 [-1.7, -0.2] | 0.57 [0.51, 0.64] | 406 |
| 21.4" | high-pass 1000 pc | -0.25 | -1.5 [-2.3, -0.7] | 0.49 [0.43, 0.56] | 268 |
| 8.2" | rotation model | -1.16 | -2.3 [-3.8, -0.9] | 0.17 [0.14, 0.20] | 429 |
| 8.2" | none | -0.18 | -2.2 [-3.9, -0.6] | 0.09 [0.06, 0.13] | 326 |
| 8.2" | high-pass 500 pc | -0.94 | -1.4 [-2.8, -0.0] | 0.21 [0.18, 0.24] | 290 |
| 8.2" | high-pass 1000 pc | -1.08 | -2.0 [-3.2, -0.8] | 0.19 [0.16, 0.23] | 230 |
| 5.3" | rotation model | -1.63 | -3.0 [-4.7, -1.3] | 0.17 [0.15, 0.18] | 1137 |
| 5.3" | none | -0.89 | -2.6 [-4.3, -1.2] | 0.12 [0.10, 0.13] | 317 |
| 5.3" | high-pass 500 pc | -1.16 | -1.8 [-3.6, -0.2] | 0.18 [0.17, 0.20] | 424 |
| 5.3" | high-pass 1000 pc | -1.45 | -2.6 [-4.1, -1.3] | 0.17 [0.16, 0.19] | 266 |

Divergence in km/s/kpc, with the 68% wedge-bootstrap interval. The all-pixel
median is the control over the whole aperture.

**Step near 1.2 kpc** (5.3", median in the 900-1200 pc then 1200-1500 pc
annulus, km/s/kpc):

| Subtraction | CO-bright | control (all pixels) |
|-------------|----------:|---------------------:|
| rotation model | +5.4 to -13.0 | +3.2 to -9.2 |
| none | +4.6 to -13.2 | +3.2 to -8.6 |
| high-pass 500 pc | +6.1 to -6.1 | +2.7 to -3.7 |
| high-pass 1000 pc | +6.1 to -11.1 | +4.2 to -7.0 |

In all 12 epoch and mode combinations the CO-bright median goes from positive
to negative across this boundary (drop of 8 to 18 km/s/kpc), and the control
drops too (3 to 12).

Assessment:

- With the mom0 map the CO-bright divergence does not steepen with CO
  brightness. The four bins are within about 1 km/s/kpc of each other in every
  mode (5.3" model: -2.8, -3.5, -3.7, -3.0) and the > 7x bin is 0.5 to 2.1
  km/s/kpc below the all-pixel control in all 12 combinations. The CO-bright
  gas is converging a little more than the average pixel, with no dose
  response. The earlier steepening (-6 to -10 at > 7x) came from the `tpeak`
  edge rim, see below.
- Shear sensitivity agrees between the model, high-pass 500 pc and high-pass
  1000 pc at 5.3" and 8.2" (5.3": 0.17, 0.18, 0.17; 8.2": 0.17, 0.21, 0.19).
  At 21.4" the three span 0.45 to 0.57. With no subtraction it is lower than
  with the model by 29% (5.3"), 47% (8.2") and 67% (21.4"), so some
  large-scale removal is needed, but the DiskFit model is not required.
- The step near 1.2 kpc is in the control and survives without a model, down
  to the 500 pc scale. It is not a rotation-model residual, and it is not
  specific to CO. What causes it is not identified.
- The outer disk is negative from 1.5 to 3.2 kpc in all modes. At 2.5-3.2 kpc
  the CO-bright median is below the control in all 12 combinations (1.3 to 5.9
  km/s/kpc; 4.2 to 5.9 at 5.3").
- Model-dependent results: the turnover radius (first zero crossing) is 230 to
  465 pc in 11 of 12 cases and 1137 pc for the 5.3" model. The 0-600 pc medians
  are within a few km/s/kpc of zero, so this crossing is not well defined; the
  1.2 kpc step is the better-defined feature. The 600-900 pc annulus at 5.3" is
  +5.9 with the model, +4.7 with none, +3.7 for high-pass 1000 pc and +0.1 for
  high-pass 500 pc. The 900-1200 pc annulus is positive in all four (+4.6 to
  +6.1).
- The 0-600 pc medians depend on the CO map: at 8.2" with the model they are
  +5.2 and +5.4 with `tpeak` and +2.0 and -0.2 with mom0. A positive
  divergence inside 600 pc is not robust.
- The all-pixel median is negative in every mode at 8.2" and 5.3" (-0.2 to
  -1.6) and more negative at higher resolution. Even the high-pass gives
  -0.9 to -1.5 at 5.3" and 8.2", so this offset is not a model residual either.
- Shear sensitivity depends on resolution: with the model `rel_drop` is 0.45,
  0.17 and 0.17 at 21.4", 8.2" and 5.3". The drop is between 21.4" and 8.2";
  it does not continue to 5.3". This is bias 8. With `tpeak` the 8.2" value was
  0.26, and the 5.3" value was unchanged. I have not separated the effect of
  the rim from the change of threshold units.
- The high-pass is not valid within about one kernel width of the centre, so
  the inner 500-1000 pc cannot be tested without the model. It also removes
  real flows coherent over more than the kernel, so agreement shows the signal
  is on scales below 500 pc, not that the model is correct.
- The bootstrap intervals reflect sampling only (about +/-0.01 on `rel_drop`
  at 5.3"). The spread between modes is the estimate of the systematic.

### Edge rim in the `broad_tpeak` CO map

The first Test 2(b) run used the 8" `broad_tpeak` map (K) with a base
threshold of 0.1 K. It showed divergence steepening with CO brightness at 8.2"
and 5.3" (> 7x: -5.9 to -10.5 km/s/kpc in all four modes). Rerunning with
`strict_mom0` removed it. The cause is the footprint edge of the `tpeak` map.

- `tpeak` is a peak (maximum) statistic, so noise biases it upward, and the
  noise is largest at the coverage edge. The median `tpeak` is 1.05 K in the
  outermost pixel and 0.08 K in the interior; 98% of pixels within 2 pixels
  of the edge exceed the 2x cut and 33% of the outermost ring exceeds the 7x
  cut. `broad_etpeak` is about 8 times higher in the rim (0.24 K against
  0.03 K). `strict_mom0` is zero there, with under 0.4% of the outermost ring
  above either cut. The two maps have the same footprint.
- The rim supplies a large share of the CO-bright sample (8.2", r < 4 kpc),
  against 4.9% of all CO-covered pixels:

| Map | Selection | Pixels | Within 100 pc of the edge |
|-----|-----------|-------:|--------------------------:|
| `tpeak` | CO > 2x | 17516 | 33% |
| `tpeak` | CO > 7x | 747 | 27% |
| `strict_mom0` | CO > 2x | 47582 | 0.2% |
| `strict_mom0` | CO > 7x | 20738 | 0.3% |

- The rim pixels have a much more negative divergence. For `tpeak` > 7x the
  median is -11.5 km/s/kpc within 100 pc of the edge (n=194) and -3.5 beyond
  (n=546). This is consistent with gradients at the edge of the HI map, or
  the HI and CO footprints, but I did not test the cause.
- Excluding CO within 100 pc of the edge (8.2", rotation model, median
  divergence per bin, km/s/kpc):

| CO bin | `tpeak` | `tpeak`, edge excluded | `strict_mom0` |
|--------|--------:|-----------------------:|--------------:|
| 2x-3x | -1.0 | -2.1 | -1.9 |
| 3x-4x | -1.4 | -2.9 | -2.7 |
| 4x-7x | -3.4 | -1.5 | -2.6 |
| > 7x | -6.1 | -3.5 [-8.0, -0.3] | -2.3 |

  The `tpeak` steepening disappears and the bins match mom0. Excluding 200 pc
  gives the same, and mom0 is unchanged by the exclusion.
- The 7m peak map used in `FINDINGS.md` has no rim: 0.4% of its CO > 2x pixels
  and none of its CO > 7x pixels lie within 100 pc of the edge. The steepening
  in those results is not explained by this artefact, and it has not been
  retested.
- Limits: this was run only at 8.2" in rotation-model mode, as an ad hoc
  script that is not in the repo. Any peak-type CO map used with a fixed
  threshold needs its footprint edge masked. This also bears on bias 6
  (threshold) and bias 9 (mask edges).
