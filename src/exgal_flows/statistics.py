"""The statistical results: CO fraction vs. flow, and its significance.

Three functions, one per headline result:

* :func:`co_fraction_vs_gradient` -- fraction of pixels that are CO-bright
  as a function of |shear| or divergence (the resolution-study curves).
* :func:`divergence_of_co_bright` -- the median divergence of CO-bright
  pixels, with a wedge block-bootstrap confidence interval (the "CO sits
  at converging flows" result).
* :func:`wedge_bootstrap_median` -- the reusable block bootstrap used by
  the above, resampling azimuthal wedges to respect the spatial
  correlation of the maps.

Everything operates on flattened (raveled) pixel arrays so the caller
controls masking (radius cut, finite values) and decimation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def assign_wedges(azimuth_rav: np.ndarray, n_wedges: int) -> np.ndarray:
    """Assign each pixel to an azimuthal wedge.

    Wedges are defined over the full range of azimuth so that a given wedge
    index maps to the same physical sector for every subsequent subset.

    Parameters
    ----------
    azimuth_rav : numpy.ndarray
        Raveled azimuth array (radians).
    n_wedges : int
        Number of wedges.

    Returns
    -------
    numpy.ndarray
        Integer wedge index in ``[0, n_wedges)`` for each pixel.
    """
    a_min, a_max = np.nanmin(azimuth_rav), np.nanmax(azimuth_rav)
    edges = np.linspace(a_min, a_max, n_wedges + 1)
    return np.clip(np.digitize(azimuth_rav, edges) - 1, 0, n_wedges - 1)


def wedge_bootstrap_median(
    values: np.ndarray,
    wedge_indices: np.ndarray,
    n_wedges: int,
    n_boot: int = 1000,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Block-bootstrap distribution of the median by resampling wedges.

    Draws ``n_wedges`` wedges with replacement, pools their pixels, and
    takes the median; repeats ``n_boot`` times.  Resampling whole wedges
    (rather than individual pixels) preserves the spatial correlation
    within a sector, giving an honest uncertainty on the median.

    Parameters
    ----------
    values : numpy.ndarray
        The quantity whose median is bootstrapped (e.g. divergence in
        km/s/kpc), one entry per pixel in the subset.
    wedge_indices : numpy.ndarray
        Wedge index for each entry of ``values`` (same length).
    n_wedges : int
        Number of distinct wedges.
    n_boot : int
        Number of bootstrap resamples.
    rng : numpy.random.Generator, optional
        Random generator for reproducibility.  Defaults to a fresh default
        generator.

    Returns
    -------
    numpy.ndarray
        Length-``n_boot`` array of bootstrap medians (may contain NaN if a
        draw happens to pool no data).
    """
    if rng is None:
        rng = np.random.default_rng()

    # Pre-split the values by wedge once, so each resample is a cheap gather.
    by_wedge = [values[wedge_indices == w] for w in range(n_wedges)]

    boot = np.empty(n_boot)
    for b in range(n_boot):
        drawn = rng.integers(0, n_wedges, size=n_wedges)
        pooled = [by_wedge[w] for w in drawn if by_wedge[w].size > 0]
        boot[b] = np.nanmedian(np.concatenate(pooled)) if pooled else np.nan
    return boot


@dataclass
class GradientCurve:
    """CO-bright fraction as a function of a velocity-gradient axis.

    Attributes
    ----------
    centers : numpy.ndarray
        Bin centres (same units as the input gradient).
    fraction : numpy.ndarray
        Fraction of pixels in each bin that exceed the CO threshold.
    """

    centers: np.ndarray
    fraction: np.ndarray


def co_fraction_vs_gradient(
    gradient: np.ndarray,
    co: np.ndarray,
    co_threshold: float,
    bin_edges: np.ndarray,
    mask: np.ndarray | None = None,
) -> GradientCurve:
    """Fraction of CO-bright pixels as a function of a gradient axis.

    Parameters
    ----------
    gradient : numpy.ndarray
        Raveled gradient values (e.g. ``|shear|`` or divergence, km/s/pc).
    co : numpy.ndarray
        Raveled CO peak-intensity values (same length).
    co_threshold : float
        A pixel is "CO-bright" if ``co > co_threshold``.
    bin_edges : numpy.ndarray
        Left edges of the gradient bins.  Bin ``i`` spans
        ``[bin_edges[i], bin_edges[i] + d)`` where ``d`` is the (uniform)
        spacing.
    mask : numpy.ndarray, optional
        Boolean pre-selection (e.g. radius cut) applied to every bin.

    Returns
    -------
    GradientCurve
    """
    d = bin_edges[1] - bin_edges[0]
    if mask is None:
        mask = np.ones_like(gradient, dtype=bool)

    is_bright = co > co_threshold
    frac = np.full(len(bin_edges), np.nan)
    for i, lo in enumerate(bin_edges):
        sel = (gradient > lo) & (gradient < lo + d) & mask
        if np.any(sel):
            frac[i] = np.nanmean(is_bright[sel])
    return GradientCurve(centers=bin_edges + d / 2, fraction=frac)


@dataclass
class DivergenceResult:
    """Median divergence of CO-bright pixels, per CO-brightness bin.

    Attributes
    ----------
    labels : list[str]
        Human-readable CO-brightness bin labels (e.g. ``"3x - 4x"``).
    medians : numpy.ndarray
        Median divergence (km/s/kpc) in each CO-brightness bin.
    err_low, err_high : numpy.ndarray
        16th / 84th percentile of the wedge-bootstrap median (km/s/kpc).
    """

    labels: list[str]
    medians: np.ndarray
    err_low: np.ndarray
    err_high: np.ndarray


def divergence_of_co_bright(
    divergence: np.ndarray,
    co: np.ndarray,
    azimuth: np.ndarray,
    co_threshold: float,
    factors=(2.0, 3.0, 4.0, 7.0),
    n_wedges: int = 8,
    n_boot: int = 1000,
    seed: int = 0,
) -> DivergenceResult:
    """Median divergence in CO-brightness bins, with bootstrap errors.

    For each CO-brightness bin (defined as multiples of ``co_threshold``),
    compute the median line-of-sight divergence of those pixels and a
    wedge block-bootstrap 68% interval.  A negative median means CO
    preferentially sits at converging / colliding flows.

    Parameters
    ----------
    divergence : numpy.ndarray
        Raveled divergence values in km/s/pc (converted to km/s/kpc
        internally by the ``*1000`` factor, matching the notebook axes).
    co : numpy.ndarray
        Raveled CO peak-intensity values.
    azimuth : numpy.ndarray
        Raveled azimuth (radians) for wedge assignment.
    co_threshold : float
        Base CO threshold; bins are ``factors[i] * co_threshold``.
    factors : sequence of float
        Lower edges of the CO-brightness bins, in units of ``co_threshold``.
        The last bin is open-ended (``> factors[-1]``).
    n_wedges, n_boot : int
        Bootstrap configuration.
    seed : int
        Seed for reproducible resampling.

    Returns
    -------
    DivergenceResult
    """
    rng = np.random.default_rng(seed)
    factors = list(factors)

    wedge_all = assign_wedges(azimuth, n_wedges)

    labels: list[str] = []
    medians, err_low, err_high = [], [], []

    for i, f in enumerate(factors):
        lower = co_threshold * f
        if i < len(factors) - 1:
            upper = co_threshold * factors[i + 1]
            sel = (co > lower) & (co <= upper)
            label = f"{f:g}x - {factors[i + 1]:g}x"
        else:
            sel = co > lower
            label = f"> {f:g}x"

        dvs_sub = divergence[sel] * 1000.0  # km/s/kpc
        wedge_sub = wedge_all[sel]

        labels.append(label)
        medians.append(np.nanmedian(dvs_sub))

        boot = wedge_bootstrap_median(dvs_sub, wedge_sub, n_wedges, n_boot, rng=rng)
        err_low.append(np.nanpercentile(boot, 16))
        err_high.append(np.nanpercentile(boot, 84))

    return DivergenceResult(
        labels=labels,
        medians=np.array(medians),
        err_low=np.array(err_low),
        err_high=np.array(err_high),
    )


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """Weighted median of ``values`` with integer or float ``weights``.

    Used by the tiled block bootstrap, where each pixel's weight is the
    number of times its block was drawn in a given resample.
    """
    if values.size == 0:
        return np.nan
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cw = np.cumsum(w)
    total = cw[-1]
    if total <= 0:
        return np.nan
    # First value at which the cumulative weight reaches half the total.
    return float(v[np.searchsorted(cw, 0.5 * total)])


def assign_tiles(i_pix, j_pix, block_pix: int):
    """Assign pixels to square spatial tiles for the block bootstrap.

    Tiles the image plane into ``block_pix`` x ``block_pix`` squares and
    returns, for each pixel, the index of the (occupied) tile it belongs to.
    This is the exact keying used by
    :func:`tiled_block_bootstrap_radial_profile`, factored out so a figure
    can draw the same tiles the bootstrap resamples.

    Parameters
    ----------
    i_pix, j_pix : numpy.ndarray
        Integer pixel coordinates (row, col).
    block_pix : int
        Tile side length in pixels (``>= 1``).

    Returns
    -------
    tile_idx : numpy.ndarray
        Index into ``0 .. n_tiles-1`` of the tile each pixel falls in.
    n_tiles : int
        Number of distinct occupied tiles.
    block_pix : int
        The (clamped) tile side actually used, in pixels.
    """
    block_pix = max(1, int(block_pix))
    bi = (np.asarray(i_pix) // block_pix).astype(np.int64)
    bj = (np.asarray(j_pix) // block_pix).astype(np.int64)
    tile_key = bi * (bj.max() + 1) + bj
    uniq, tile_idx = np.unique(tile_key, return_inverse=True)
    return tile_idx, uniq.size, block_pix


@dataclass
class RadialProfile:
    """Radial profile of the CO-bright divergence with block-bootstrap errors.

    Attributes
    ----------
    r_centers : numpy.ndarray
        Annulus centres (pc).
    r_edges : numpy.ndarray
        Annulus edges (pc), length ``len(r_centers) + 1``.
    medians : numpy.ndarray
        Point-estimate median divergence (km/s/kpc) in each annulus.
    err_low, err_high : numpy.ndarray
        16th/84th percentile of the block-bootstrap median in each annulus.
    iid_low, iid_high : numpy.ndarray
        16th/84th percentile of the per-pixel (iid) bootstrap, for comparison.
    turnover_pc : float
        Radius (pc) where the point-estimate profile crosses zero divergence
        (linearly interpolated), or NaN if it never crosses.
    """

    r_centers: np.ndarray
    r_edges: np.ndarray
    medians: np.ndarray
    err_low: np.ndarray
    err_high: np.ndarray
    iid_low: np.ndarray
    iid_high: np.ndarray
    turnover_pc: float


def tiled_block_bootstrap_radial_profile(
    i_pix: np.ndarray,
    j_pix: np.ndarray,
    divergence: np.ndarray,
    radius_pc: np.ndarray,
    pixel_scale_pc: float,
    r_edges: np.ndarray,
    block_pc: float = 300.0,
    n_boot: int = 1000,
    seed: int | None = 0,
) -> RadialProfile:
    """Radial profile of CO-bright divergence with a spatial block bootstrap.

    The per-pixel (iid) bootstrap badly understates the uncertainty here
    because the maps are correlated on a range of scales and CO organises
    into complexes tens to hundreds of pc across.  This routine tiles the
    image plane into square blocks ``block_pc`` on a side, resamples whole
    blocks with replacement, and takes a block-weighted median of the
    divergence in each galactocentric annulus -- preserving the internal
    correlation of each block.

    Parameters
    ----------
    i_pix, j_pix : numpy.ndarray
        Integer pixel coordinates (row, col) of the CO-bright pixels.
    divergence : numpy.ndarray
        Line-of-sight divergence per pixel, in km/s/kpc (already scaled).
    radius_pc : numpy.ndarray
        Galactocentric (deprojected image-plane) radius per pixel, in pc.
    pixel_scale_pc : float
        Physical size of one pixel, in pc.
    r_edges : numpy.ndarray
        Annulus edges in pc (monotonic increasing).
    block_pc : float
        Block side length in pc.
    n_boot : int
        Number of bootstrap resamples.
    seed : int, optional
        Seed for reproducibility.

    Returns
    -------
    RadialProfile
    """
    r_edges = np.asarray(r_edges, float)
    rng = np.random.default_rng(seed)

    block_pix = max(1, int(round(block_pc / pixel_scale_pc)))
    # Tile id per pixel (same keying a figure can reuse via assign_tiles).
    tile_idx, n_tiles, block_pix = assign_tiles(i_pix, j_pix, block_pix)

    ann = np.digitize(radius_pc, r_edges) - 1  # annulus index per pixel
    n_ann = len(r_edges) - 1

    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    medians = np.full(n_ann, np.nan)
    err_low = np.full(n_ann, np.nan)
    err_high = np.full(n_ann, np.nan)
    iid_low = np.full(n_ann, np.nan)
    iid_high = np.full(n_ann, np.nan)

    for a in range(n_ann):
        in_ann = ann == a
        if not np.any(in_ann):
            continue
        dv = divergence[in_ann]
        medians[a] = np.nanmedian(dv)

        # iid bootstrap (per-pixel) for comparison.
        n = dv.size
        iid = np.empty(n_boot)
        for b in range(n_boot):
            iid[b] = np.nanmedian(dv[rng.integers(0, n, size=n)])
        iid_low[a], iid_high[a] = np.nanpercentile(iid, [16, 84])

        # Spatial block bootstrap: resample tiles, weight pixels by how many
        # times their tile was drawn.
        t_ann = tile_idx[in_ann]
        boot = np.empty(n_boot)
        for b in range(n_boot):
            counts = np.bincount(rng.integers(0, n_tiles, size=n_tiles),
                                 minlength=n_tiles)
            w = counts[t_ann].astype(float)
            boot[b] = _weighted_median(dv, w)
        err_low[a], err_high[a] = np.nanpercentile(boot, [16, 84])

    # Turnover: first zero-crossing of the point-estimate profile.
    turnover = np.nan
    good = np.isfinite(medians)
    rc = r_centers[good]
    mc = medians[good]
    for k in range(len(mc) - 1):
        if mc[k] == 0:
            turnover = rc[k]
            break
        if mc[k] * mc[k + 1] < 0:  # sign change between k and k+1
            frac = mc[k] / (mc[k] - mc[k + 1])
            turnover = rc[k] + frac * (rc[k + 1] - rc[k])
            break

    return RadialProfile(
        r_centers=r_centers,
        r_edges=r_edges,
        medians=medians,
        err_low=err_low,
        err_high=err_high,
        iid_low=iid_low,
        iid_high=iid_high,
        turnover_pc=float(turnover),
    )


@dataclass
class GradientCurveCI:
    """CO-bright fraction vs a gradient axis, with a block-bootstrap band.

    Like :class:`GradientCurve` but carries a wedge block-bootstrap
    confidence interval on the fraction in each bin, plus the pixel count
    per bin and a scalar *sensitivity* (the fractional drop from the low-
    gradient to the high-gradient end).

    Attributes
    ----------
    centers : numpy.ndarray
        Bin centres (same units as the input gradient).
    fraction : numpy.ndarray
        CO-bright fraction (point estimate) in each bin.
    err_low, err_high : numpy.ndarray
        16th / 84th percentile of the block-bootstrap fraction per bin.
    counts : numpy.ndarray
        Number of pixels contributing to each bin.
    lo_frac, hi_frac : float
        Point-estimate fraction in the low- and high-gradient reference
        windows (see :func:`co_fraction_vs_gradient_boot`).
    rel_drop : float
        ``1 - hi_frac / lo_frac`` -- the relative suppression of CO from
        low to high gradient (the "sensitivity").  NaN if ``lo_frac`` is 0.
    rel_drop_lo, rel_drop_hi : float
        16th / 84th percentile of ``rel_drop`` under the block bootstrap.
    """

    centers: np.ndarray
    fraction: np.ndarray
    err_low: np.ndarray
    err_high: np.ndarray
    counts: np.ndarray
    lo_frac: float
    hi_frac: float
    rel_drop: float
    rel_drop_lo: float
    rel_drop_hi: float


def co_fraction_vs_gradient_boot(
    gradient: np.ndarray,
    co: np.ndarray,
    azimuth: np.ndarray,
    co_threshold: float,
    bin_edges: np.ndarray,
    n_wedges: int = 8,
    n_boot: int = 1000,
    seed: int | None = 0,
    mask: np.ndarray | None = None,
    lo_hi_split: float | None = None,
) -> GradientCurveCI:
    """CO-bright fraction vs a gradient axis with a wedge block bootstrap.

    Identical binning to :func:`co_fraction_vs_gradient`, but resamples
    azimuthal wedges (whole sectors, with replacement) to put an honest
    error band on the fraction in each bin -- the per-pixel iid error is
    invalid here because the maps are spatially correlated.

    A scalar *sensitivity* is also returned: the relative drop in CO-bright
    fraction between a low-gradient and a high-gradient reference window,
    split at ``lo_hi_split`` (defaults to the midpoint of ``bin_edges``).

    Parameters
    ----------
    gradient : numpy.ndarray
        Raveled gradient values (e.g. ``|shear|``, km/s/pc).
    co : numpy.ndarray
        Raveled CO peak-intensity values (same length).
    azimuth : numpy.ndarray
        Raveled azimuth (radians), for wedge assignment.
    co_threshold : float
        A pixel is CO-bright if ``co > co_threshold``.
    bin_edges : numpy.ndarray
        Left edges of the (uniform) gradient bins.
    n_wedges : int
        Number of azimuthal wedges to resample.
    n_boot : int
        Number of bootstrap resamples.
    seed : int, optional
        Seed for reproducibility.
    mask : numpy.ndarray, optional
        Boolean pre-selection (e.g. a radius window) applied everywhere.
    lo_hi_split : float, optional
        Gradient value separating the low/high reference windows for the
        sensitivity metric.  Defaults to the midpoint of ``bin_edges``.

    Returns
    -------
    GradientCurveCI
    """
    d = bin_edges[1] - bin_edges[0]
    rng = np.random.default_rng(seed)
    if mask is None:
        mask = np.ones_like(gradient, dtype=bool)

    g = gradient[mask]
    bright = (co[mask] > co_threshold)
    wedge = assign_wedges(azimuth[mask], n_wedges)

    # bin index per pixel (-1 outside the range)
    bidx = np.floor((g - bin_edges[0]) / d).astype(int)
    in_range = (bidx >= 0) & (bidx < len(bin_edges))

    centers = bin_edges + d / 2
    counts = np.array([np.count_nonzero(in_range & (bidx == i))
                       for i in range(len(bin_edges))])

    def _profile(sel_pix):
        frac = np.full(len(bin_edges), np.nan)
        for i in range(len(bin_edges)):
            s = sel_pix & (bidx == i)
            if np.any(s):
                frac[i] = np.mean(bright[s])
        return frac

    fraction = _profile(in_range)

    if lo_hi_split is None:
        lo_hi_split = 0.5 * (bin_edges[0] + bin_edges[-1])
    lo_win = in_range & (g < lo_hi_split)
    hi_win = in_range & (g >= lo_hi_split)

    def _drop(sel_pix):
        lo = np.mean(bright[sel_pix & lo_win]) if np.any(sel_pix & lo_win) else np.nan
        hi = np.mean(bright[sel_pix & hi_win]) if np.any(sel_pix & hi_win) else np.nan
        rd = 1.0 - hi / lo if (lo and np.isfinite(lo) and lo > 0) else np.nan
        return lo, hi, rd

    lo_frac, hi_frac, rel_drop = _drop(in_range)

    # wedge block bootstrap
    by_wedge = [np.nonzero(wedge == w)[0] for w in range(n_wedges)]
    boot_frac = np.full((n_boot, len(bin_edges)), np.nan)
    boot_drop = np.full(n_boot, np.nan)
    for b in range(n_boot):
        drawn = rng.integers(0, n_wedges, size=n_wedges)
        idx = np.concatenate([by_wedge[w] for w in drawn if by_wedge[w].size])
        sel = np.zeros(g.shape[0], dtype=bool)
        # a pixel can be drawn multiple times; for a fraction the multiplicity
        # cancels, so membership is enough
        sel[idx] = True
        boot_frac[b] = _profile(in_range & sel)
        _, _, boot_drop[b] = _drop(in_range & sel)

    err_low = np.nanpercentile(boot_frac, 16, axis=0)
    err_high = np.nanpercentile(boot_frac, 84, axis=0)
    rel_drop_lo = float(np.nanpercentile(boot_drop, 16))
    rel_drop_hi = float(np.nanpercentile(boot_drop, 84))

    return GradientCurveCI(
        centers=centers,
        fraction=fraction,
        err_low=err_low,
        err_high=err_high,
        counts=counts,
        lo_frac=float(lo_frac),
        hi_frac=float(hi_frac),
        rel_drop=float(rel_drop),
        rel_drop_lo=rel_drop_lo,
        rel_drop_hi=rel_drop_hi,
    )


@dataclass
class LinearRadialFit:
    """A straight-line fit ``div(r) = a + b * r`` to CO-bright divergence.

    The simplest model of the radial reversal: a constant *sense of
    divergence* ``a`` (the value extrapolated to the galaxy centre) plus a
    constant *rate of change with radius* ``b``.  Fit directly to the
    per-pixel divergence by (weighted) least squares, with a spatial block
    bootstrap over tiles for honest errors.

    Attributes
    ----------
    intercept, slope : float
        Point-estimate ``a`` (km/s/kpc) and ``b`` (km/s/kpc per kpc).
    intercept_lo, intercept_hi, slope_lo, slope_hi : float
        16th / 84th percentile of ``a`` and ``b`` under the block bootstrap.
    zero_crossing_pc : float
        Radius where the fitted line crosses zero, ``-a / b`` (pc).
    zero_crossing_lo_pc, zero_crossing_hi_pc : float
        16th / 84th percentile of the zero-crossing (pc).
    r_grid_pc : numpy.ndarray
        Radii (pc) at which the fitted line and its band are tabulated.
    line, band_lo, band_hi : numpy.ndarray
        Point-estimate fitted line and its 16th / 84th percentile band
        (km/s/kpc) on ``r_grid_pc``.
    slope_neg_prob : float
        Fraction of bootstrap draws with ``b < 0`` (evidence the divergence
        decreases with radius, i.e. the reversal is real).
    method : str
        ``"ols"`` (least squares) or ``"l1"`` (least absolute deviations).
    label : str
        Legend label with the fitted parameters.
    """

    intercept: float
    slope: float
    intercept_lo: float
    intercept_hi: float
    slope_lo: float
    slope_hi: float
    zero_crossing_pc: float
    zero_crossing_lo_pc: float
    zero_crossing_hi_pc: float
    r_grid_pc: np.ndarray
    line: np.ndarray
    band_lo: np.ndarray
    band_hi: np.ndarray
    slope_neg_prob: float
    method: str = "ols"
    label: str = ""


def _wls_line(r, y, w):
    """Weighted least-squares slope/intercept of ``y = a + b r``."""
    sw = w.sum()
    if sw <= 0:
        return np.nan, np.nan
    mr = (w * r).sum() / sw
    my = (w * y).sum() / sw
    srr = (w * (r - mr) ** 2).sum()
    if srr <= 0:
        return np.nan, np.nan
    b = (w * (r - mr) * (y - my)).sum() / srr
    a = my - b * mr
    return a, b


def _l1_line(r, y, w, n_iter: int = 30, eps: float = 1e-3):
    """Weighted least-absolute-deviations line via IRLS.

    Minimises ``sum_i w_i |y_i - a - b r_i|`` by iteratively reweighted
    least squares: each iteration solves a weighted-LS problem with weights
    ``w_i / max(|residual_i|, eps)``, which converges to the L1 (median-like)
    solution.  This is the direct-data analogue of the per-annulus *median*
    profile and is robust to the heavy tails of per-pixel divergence.
    """
    a, b = _wls_line(r, y, w)  # OLS warm start
    if not (np.isfinite(a) and np.isfinite(b)):
        return a, b
    for _ in range(n_iter):
        resid = np.abs(y - a - b * r)
        ww = w / np.maximum(resid, eps)
        a_new, b_new = _wls_line(r, y, ww)
        if not (np.isfinite(a_new) and np.isfinite(b_new)):
            break
        if abs(a_new - a) < 1e-4 and abs(b_new - b) < 1e-4:
            a, b = a_new, b_new
            break
        a, b = a_new, b_new
    return a, b


_LINE_FITTERS = {"ols": _wls_line, "l1": _l1_line}


def block_bootstrap_line_fit(
    radius_pc: np.ndarray,
    divergence: np.ndarray,
    tile_idx: np.ndarray,
    n_tiles: int,
    n_boot: int = 1000,
    seed: int | None = 0,
    r_grid_pc: np.ndarray | None = None,
    method: str = "ols",
) -> LinearRadialFit:
    """Fit ``div(r) = a + b r`` to pixels, block-bootstrapped over tiles.

    The point estimate is a straight line through the per-pixel divergence
    (radius in kpc internally, so ``b`` is km/s/kpc *per kpc*), fit either by
    ordinary least squares (``method="ols"``) or by least absolute
    deviations (``method="l1"``, the robust median-like fit).  Errors come
    from resampling whole spatial tiles with replacement (see
    :func:`assign_tiles`) and refitting with per-pixel weights equal to the
    number of times each pixel's tile was drawn -- the same block scheme used
    for the radial profile, so the line's uncertainty respects the maps'
    spatial correlation rather than treating every pixel as independent.

    Parameters
    ----------
    radius_pc : numpy.ndarray
        Galactocentric radius per pixel (pc).
    divergence : numpy.ndarray
        Line-of-sight divergence per pixel (km/s/kpc).
    tile_idx : numpy.ndarray
        Tile index per pixel from :func:`assign_tiles`.
    n_tiles : int
        Number of occupied tiles.
    n_boot : int
        Number of bootstrap resamples.
    seed : int, optional
        Seed for reproducibility.
    r_grid_pc : numpy.ndarray, optional
        Radii (pc) at which to tabulate the fitted line/band.  Defaults to
        100 points spanning the data.

    Returns
    -------
    LinearRadialFit
    """
    rng = np.random.default_rng(seed)
    r_kpc = np.asarray(radius_pc, float) / 1000.0
    y = np.asarray(divergence, float)
    good = np.isfinite(r_kpc) & np.isfinite(y)
    r_kpc, y, tile_idx = r_kpc[good], y[good], np.asarray(tile_idx)[good]

    if r_grid_pc is None:
        r_grid_pc = np.linspace(radius_pc.min(), radius_pc.max(), 100)
    r_grid_pc = np.asarray(r_grid_pc, float)
    rg_kpc = r_grid_pc / 1000.0

    try:
        fitter = _LINE_FITTERS[method]
    except KeyError:
        raise ValueError(f"method must be 'ols' or 'l1', got {method!r}")

    a0, b0 = fitter(r_kpc, y, np.ones_like(y))
    line = a0 + b0 * rg_kpc

    a_s = np.empty(n_boot)
    b_s = np.empty(n_boot)
    band = np.empty((n_boot, rg_kpc.size))
    for k in range(n_boot):
        counts = np.bincount(rng.integers(0, n_tiles, size=n_tiles),
                             minlength=n_tiles)
        w = counts[tile_idx].astype(float)
        a_s[k], b_s[k] = fitter(r_kpc, y, w)
        band[k] = a_s[k] + b_s[k] * rg_kpc

    a_lo, a_hi = np.nanpercentile(a_s, [16, 84])
    b_lo, b_hi = np.nanpercentile(b_s, [16, 84])
    band_lo = np.nanpercentile(band, 16, axis=0)
    band_hi = np.nanpercentile(band, 84, axis=0)

    zc = -a0 / b0 * 1000.0 if (np.isfinite(b0) and b0 != 0) else np.nan
    zc_s = np.where(b_s != 0, -a_s / b_s * 1000.0, np.nan)
    zc_lo, zc_hi = np.nanpercentile(zc_s, [16, 84])
    slope_neg_prob = float(np.mean(b_s < 0))

    return LinearRadialFit(
        intercept=float(a0), slope=float(b0),
        intercept_lo=float(a_lo), intercept_hi=float(a_hi),
        slope_lo=float(b_lo), slope_hi=float(b_hi),
        zero_crossing_pc=float(zc),
        zero_crossing_lo_pc=float(zc_lo), zero_crossing_hi_pc=float(zc_hi),
        r_grid_pc=r_grid_pc, line=line, band_lo=band_lo, band_hi=band_hi,
        slope_neg_prob=slope_neg_prob, method=method,
        label=f"{method.upper()} fit: a={a0:.1f}, b={b0:.1f}",
    )


@dataclass
class SigmoidRadialFit:
    """A saturating (sigmoid) fit to CO-bright divergence vs radius.

    A bounded alternative to :class:`LinearRadialFit`: instead of a line that
    grows without limit, the divergence saturates to an inner plateau ``hi``
    and an outer plateau ``lo``, switching over a transition radius ``r0`` with
    a characteristic ``width``::

        div(r) = lo + (hi - lo) / (1 + exp((r - r0) / width))

    (a *decreasing* step, since divergence goes from positive/converging in the
    inner disk to negative/diverging outside).  Fit directly to the per-pixel
    data by least squares, block-bootstrapped over spatial tiles.

    Attributes
    ----------
    hi, lo : float
        Inner and outer plateau of the divergence (km/s/kpc).
    r0_pc : float
        Transition radius (pc).
    width_pc : float
        Transition width (pc); small width -> sharp step.
    hi_lo, hi_hi, lo_lo, lo_hi, r0_lo_pc, r0_hi_pc, width_lo_pc, width_hi_pc : float
        16th / 84th percentiles of each parameter under the block bootstrap.
    amplitude : float
        ``hi - lo`` -- the size of the reversal (km/s/kpc).
    amplitude_lo, amplitude_hi : float
        16th / 84th percentile of the amplitude.
    zero_crossing_pc : float
        Radius where the curve crosses zero (pc), or NaN if it never does.
    zero_crossing_lo_pc, zero_crossing_hi_pc : float
        16th / 84th percentile of the zero-crossing.
    straddle_zero_prob : float
        Fraction of bootstrap draws with ``hi > 0`` **and** ``lo < 0`` -- i.e.
        the probability that the reversal genuinely straddles zero.  This is
        the sigmoid analogue of :attr:`LinearRadialFit.slope_neg_prob`.
    r_grid_pc, line, band_lo, band_hi : numpy.ndarray
        Point-estimate curve and its 16th / 84th percentile band on
        ``r_grid_pc`` (``line`` named to match :class:`LinearRadialFit` so the
        two plug into the same plotting routine).
    n_converged : int
        Number of bootstrap draws whose fit converged.
    method : str
        Always ``"sigmoid"`` (lets the plot route on the same field).
    label : str
        Legend label with the fitted plateaus.
    """

    hi: float
    lo: float
    r0_pc: float
    width_pc: float
    hi_lo: float
    hi_hi: float
    lo_lo: float
    lo_hi: float
    r0_lo_pc: float
    r0_hi_pc: float
    width_lo_pc: float
    width_hi_pc: float
    amplitude: float
    amplitude_lo: float
    amplitude_hi: float
    zero_crossing_pc: float
    zero_crossing_lo_pc: float
    zero_crossing_hi_pc: float
    straddle_zero_prob: float
    r_grid_pc: np.ndarray
    line: np.ndarray
    band_lo: np.ndarray
    band_hi: np.ndarray
    n_converged: int
    method: str = "sigmoid"
    label: str = ""


def _sigmoid(r_kpc, hi, lo, r0, w):
    """Decreasing sigmoid in kpc: hi at small r, lo at large r."""
    return lo + (hi - lo) / (1.0 + np.exp((r_kpc - r0) / w))


def _sigmoid_zero_pc(hi, lo, r0, w):
    """Radius (pc) where the sigmoid crosses zero, or NaN."""
    if hi > 0 and lo < 0:
        return (r0 + w * np.log(-hi / lo)) * 1000.0
    return np.nan


def block_bootstrap_sigmoid_fit(
    radius_pc: np.ndarray,
    divergence: np.ndarray,
    tile_idx: np.ndarray,
    n_tiles: int,
    n_boot: int = 1000,
    seed: int | None = 0,
    r_grid_pc: np.ndarray | None = None,
    p0=(5.0, -5.0, 1.2, 0.3),
    bounds=((0.0, -30.0, 0.0, 0.02), (30.0, 0.0, 4.0, 3.0)),
    maxfev: int = 8000,
) -> SigmoidRadialFit:
    """Fit a saturating sigmoid ``div(r)`` to pixels, block-bootstrapped.

    Least-squares fit of :func:`_sigmoid` to the per-pixel divergence (radius
    in kpc internally), with errors from resampling whole spatial tiles (see
    :func:`assign_tiles`).  The headline significance is
    :attr:`SigmoidRadialFit.straddle_zero_prob` -- the fraction of draws whose
    inner plateau is positive *and* outer plateau negative, i.e. the reversal
    rules out zero.  This concentrates the sign information into two plateaus
    rather than the intercept+slope of the line fit, so it typically rejects
    "no reversal" more decisively.

    Parameters
    ----------
    radius_pc, divergence : numpy.ndarray
        Per-pixel radius (pc) and divergence (km/s/kpc).
    tile_idx, n_tiles : numpy.ndarray, int
        Tile assignment from :func:`assign_tiles`.
    n_boot, seed : int
        Bootstrap resamples and RNG seed.
    r_grid_pc : numpy.ndarray, optional
        Radii (pc) at which to tabulate the curve/band.
    p0, bounds, maxfev :
        Passed to :func:`scipy.optimize.curve_fit` (initial guess
        ``(hi, lo, r0_kpc, width_kpc)``, box bounds, and iteration cap).

    Returns
    -------
    SigmoidRadialFit
    """
    from scipy.optimize import curve_fit

    rng = np.random.default_rng(seed)
    r_kpc = np.asarray(radius_pc, float) / 1000.0
    y = np.asarray(divergence, float)
    tile_idx = np.asarray(tile_idx)
    good = np.isfinite(r_kpc) & np.isfinite(y)
    r_kpc, y, tile_idx = r_kpc[good], y[good], tile_idx[good]

    if r_grid_pc is None:
        r_grid_pc = np.linspace(radius_pc.min(), radius_pc.max(), 100)
    r_grid_pc = np.asarray(r_grid_pc, float)
    rg_kpc = r_grid_pc / 1000.0

    lb = np.asarray(bounds[0], float)
    ub = np.asarray(bounds[1], float)

    def _fit(rr, yy):
        po, _ = curve_fit(_sigmoid, rr, yy, p0=p0, bounds=(lb, ub), maxfev=maxfev)
        return po

    p_hat = _fit(r_kpc, y)
    hi0, lo0, r00, w0 = p_hat
    line = _sigmoid(rg_kpc, *p_hat)

    by_tile = [np.nonzero(tile_idx == k)[0] for k in range(n_tiles)]
    params = np.full((n_boot, 4), np.nan)
    band = np.full((n_boot, rg_kpc.size), np.nan)
    for b in range(n_boot):
        drawn = rng.integers(0, n_tiles, size=n_tiles)
        idx = np.concatenate([by_tile[k] for k in drawn if by_tile[k].size])
        try:
            po = _fit(r_kpc[idx], y[idx])
            params[b] = po
            band[b] = _sigmoid(rg_kpc, *po)
        except Exception:
            pass

    ok = np.isfinite(params[:, 0])
    P = params[ok]
    amp = P[:, 0] - P[:, 1]
    zc_s = np.array([_sigmoid_zero_pc(*row) for row in P])

    def pct(a, q):
        return float(np.nanpercentile(a, q)) if np.any(np.isfinite(a)) else np.nan

    return SigmoidRadialFit(
        hi=float(hi0), lo=float(lo0), r0_pc=float(r00 * 1000.0),
        width_pc=float(w0 * 1000.0),
        hi_lo=pct(P[:, 0], 16), hi_hi=pct(P[:, 0], 84),
        lo_lo=pct(P[:, 1], 16), lo_hi=pct(P[:, 1], 84),
        r0_lo_pc=pct(P[:, 2] * 1000.0, 16), r0_hi_pc=pct(P[:, 2] * 1000.0, 84),
        width_lo_pc=pct(P[:, 3] * 1000.0, 16), width_hi_pc=pct(P[:, 3] * 1000.0, 84),
        amplitude=float(hi0 - lo0), amplitude_lo=pct(amp, 16), amplitude_hi=pct(amp, 84),
        zero_crossing_pc=float(_sigmoid_zero_pc(*p_hat)),
        zero_crossing_lo_pc=pct(zc_s, 16), zero_crossing_hi_pc=pct(zc_s, 84),
        straddle_zero_prob=float(np.mean((P[:, 0] > 0) & (P[:, 1] < 0))) if P.size else np.nan,
        r_grid_pc=r_grid_pc, line=line,
        band_lo=np.nanpercentile(band, 16, axis=0),
        band_hi=np.nanpercentile(band, 84, axis=0),
        n_converged=int(ok.sum()),
        label=f"sigmoid: hi={hi0:.1f}, lo={lo0:.1f}",
    )
