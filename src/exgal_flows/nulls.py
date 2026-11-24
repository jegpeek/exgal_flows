"""Null models for the CO-fraction-vs-|velocity gradient| statistics.

The worry these address: CO is found where HI is bright, while beam smearing
and low signal-to-noise make the measured velocity field roughest -- and
therefore ``|shear|`` and ``|divergence|`` largest -- where HI is faint
(between clumps).  That alone would make the CO-bright fraction fall with
``|gradient|`` even if gas kinematics had no influence on CO at all.

Two complementary nulls are provided.

1. **HI-only CO null on the real maps** (:func:`hi_structure_features`,
   :func:`hi_only_co_null`, :func:`relative_drop_vs_null`).  A classifier
   predicts CO-brightness from HI structure alone (column density and
   Laplacian at several scales, plus radius), cross-fitted on held-out
   spatial tiles.  Its prediction contains *no* velocity information, so the
   CO-vs-``|gradient|`` curve it implies is exactly what HI structure (and
   any beam/noise coupling between HI and the gradients) produces by itself.
   The observed curve is compared to it with a spatial block bootstrap.
   This is conservative: any *real* pathway that runs through HI structure
   (e.g. converging flows compressing HI, which then forms CO) is absorbed
   into the null.

2. **Beam-smearing forward simulation** (:class:`SmearingSimParams`,
   :func:`simulate_observation`).  A lognormal, patchy HI disk and a velocity
   field drawn *independently* of it are "observed" through the real
   (deprojected, elliptical) beam, with mom0/mom1 noise matched to the data
   and the same faint-pixel blanking.  CO is placed by the true HI only.  Any
   CO-vs-``|gradient|`` trend in the simulation is therefore pure artefact,
   and its size can be compared to the observed one.  Moment-space
   observation is exact for the noise-free part: the intensity-weighted
   mean velocity of a beam-convolved cube is ``G*(Sigma v) / G*Sigma``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage as nd

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


# ----------------------------------------------------------------------------
# helpers shared by both nulls
# ----------------------------------------------------------------------------
def deprojected_beam_sigma_pix(fwhm_sky_pc: float, pixel_scale_pc: float, cos_i: float):
    """Gaussian ``(sigma_y, sigma_x)`` in pixels of a circular sky beam on the
    deprojected grid.

    :func:`exgal_flows.geometry.deproject` stretches the x (minor) axis by
    ``1/cos i`` so pixels are square in the galaxy plane; a circular sky beam
    of FWHM ``fwhm_sky_pc`` therefore becomes ``fwhm/cos i`` along x.
    """
    sy = fwhm_sky_pc / pixel_scale_pc * FWHM_TO_SIGMA
    return sy, sy / cos_i


def relative_drop(gradient, co_bright, weights=None, split: float = 0.1, upper: float = 0.19):
    """``1 - f(split <= g < upper) / f(0 <= g < split)`` -- the sensitivity
    used throughout (matches :func:`statistics.co_fraction_vs_gradient_boot`
    with the default ``shear_range=(0, 0.2)`` and 0.01 bins)."""
    g = np.asarray(gradient)
    c = np.asarray(co_bright, float)
    w = np.ones_like(c) if weights is None else np.asarray(weights, float)
    ok = np.isfinite(g) & np.isfinite(c)
    lo = ok & (g >= 0) & (g < split)
    hi = ok & (g >= split) & (g < upper)
    return 1.0 - (np.sum((w * c)[hi]) / np.sum(w[hi])) / (np.sum((w * c)[lo]) / np.sum(w[lo]))


def relative_drop_vs_null(
    gradient: np.ndarray,
    co_bright: np.ndarray,
    p_null: np.ndarray,
    tile_idx: np.ndarray,
    n_tiles: int,
    split: float = 0.1,
    upper: float = 0.19,
    n_boot: int = 1000,
    seed: int | None = 0,
) -> dict:
    """Observed vs null relative drop with a spatial-tile block bootstrap.

    The observed drop uses the measured CO-bright indicator; the null drop
    replaces it by the null probability ``p_null`` (pixel by pixel).  Both are
    recomputed on the *same* resampled tiles, so the bootstrap distribution of
    their difference (the "excess" not explained by the null) carries the
    spatial correlation of the data.

    Returns
    -------
    dict
        ``obs, null, excess`` (point estimates), their 16th/84th percentiles
        (``*_lo, *_hi``), ``p_excess_le0`` (one-sided bootstrap probability
        that the observed drop does not exceed the null) and
        ``frac_explained = null / obs``.
    """
    g = np.asarray(gradient, float)
    y = np.asarray(co_bright, float)
    p = np.asarray(p_null, float)
    ok = np.isfinite(g) & np.isfinite(y) & np.isfinite(p)
    lo = ok & (g >= 0) & (g < split)
    hi = ok & (g >= split) & (g < upper)

    def per_tile(w):
        return np.bincount(tile_idx, weights=w, minlength=n_tiles)

    A = np.vstack([per_tile(lo * 1.0), per_tile(hi * 1.0), per_tile(y * lo),
                   per_tile(y * hi), per_tile(p * lo), per_tile(p * hi)])

    def calc(S):
        nlo, nhi, ylo, yhi, plo, phi = S
        return 1 - (yhi / nhi) / (ylo / nlo), 1 - (phi / nhi) / (plo / nlo)

    o0, n0 = calc(A.sum(axis=1))
    rng = np.random.default_rng(seed)
    W = rng.multinomial(n_tiles, np.full(n_tiles, 1.0 / n_tiles), size=n_boot).astype(float)
    B = np.array(calc(A @ W.T)).T          # (n_boot, 2)
    ex = B[:, 0] - B[:, 1]
    q = lambda a, k: float(np.nanpercentile(a, k))
    return dict(obs=float(o0), obs_lo=q(B[:, 0], 16), obs_hi=q(B[:, 0], 84),
                null=float(n0), null_lo=q(B[:, 1], 16), null_hi=q(B[:, 1], 84),
                excess=float(o0 - n0), excess_lo=q(ex, 16), excess_hi=q(ex, 84),
                p_excess_le0=float(np.mean(ex <= 0)), frac_explained=float(n0 / o0))


# ----------------------------------------------------------------------------
# Null 1: HI-only CO model on the real maps
# ----------------------------------------------------------------------------
def hi_structure_features(mom0: np.ndarray, sigmas_pix=(0, 2, 4, 8, 16, 32)) -> dict:
    """Multi-scale HI structure maps: column and scale-normalised Laplacian.

    Smoothing uses normalised convolution so blanked pixels do not bias the
    smoothed column near mask edges.  ``Lap_s{s}`` is ``max(s,1)**2 *
    laplace(N_s{s})`` (negative on HI peaks, positive in troughs -- the sign
    convention of the original notebook's Laplacian split).

    Returns
    -------
    dict of str -> numpy.ndarray (float32), keys ``N_s{s}`` and ``Lap_s{s}``.
    """
    m = np.asarray(mom0, float)
    valid = np.isfinite(m).astype(float)
    m0 = np.where(np.isfinite(m), m, 0.0)
    out = {}
    for s in sigmas_pix:
        if s == 0:
            sm = m0
        else:
            sm = nd.gaussian_filter(m0, s) / np.maximum(nd.gaussian_filter(valid, s), 1e-3)
        out[f"N_s{s}"] = sm.astype(np.float32)
        out[f"Lap_s{s}"] = ((max(s, 1) ** 2) * nd.laplace(sm)).astype(np.float32)
    return out


def hi_only_co_null(
    features: np.ndarray,
    co_bright: np.ndarray,
    tile_idx: np.ndarray,
    n_tiles: int,
    n_folds: int = 5,
    seed: int = 0,
    **gbm_kwargs,
) -> np.ndarray:
    """Out-of-fold P(CO-bright | HI structure) by spatial cross-fitting.

    Whole spatial tiles are assigned to folds, so the classifier never sees the
    neighbourhood of a pixel it predicts (no memorising of cloud positions).

    Parameters
    ----------
    features : (N, F) array of HI-only predictors (no velocity information).
    co_bright : (N,) bool array.
    tile_idx, n_tiles : spatial tiles from :func:`statistics.assign_tiles`.

    Returns
    -------
    numpy.ndarray
        (N,) null probability for each pixel.
    """
    from sklearn.ensemble import HistGradientBoostingClassifier

    kw = dict(max_iter=300, learning_rate=0.08, max_leaf_nodes=63,
              min_samples_leaf=200, l2_regularization=1.0, random_state=seed)
    kw.update(gbm_kwargs)
    rng = np.random.default_rng(seed)
    fold = rng.integers(0, n_folds, size=n_tiles)[tile_idx]
    X = np.asarray(features, np.float32)
    y = np.asarray(co_bright, bool)
    p = np.full(y.shape, np.nan)
    for k in range(n_folds):
        tr, te = fold != k, fold == k
        clf = HistGradientBoostingClassifier(**kw).fit(X[tr], y[tr])
        p[te] = clf.predict_proba(X[te])[:, 1]
    return p


# ----------------------------------------------------------------------------
# Null 2: beam-smearing forward simulation
# ----------------------------------------------------------------------------
def gaussian_random_field(shape, slope: float, rng) -> np.ndarray:
    """Zero-mean, unit-variance periodic Gaussian random field with power
    spectrum ``P(k) ~ k**slope`` (``slope`` negative, e.g. -3.5)."""
    ny, nx = shape
    ky = np.fft.fftfreq(ny)[:, None]
    kx = np.fft.rfftfreq(nx)[None, :]
    k = np.hypot(kx, ky)
    k[0, 0] = np.inf
    f = np.fft.irfft2(np.fft.rfft2(rng.standard_normal(shape)) * k ** (slope / 2.0), s=shape)
    f -= f.mean()
    return f / f.std()


def correlated_noise(shape, sigma_yx, rng) -> np.ndarray:
    """Unit-variance white noise convolved with a Gaussian beam (periodic)."""
    n = nd.gaussian_filter(rng.standard_normal(shape), sigma_yx, mode="wrap")
    return n / n.std()


@dataclass
class SmearingSimParams:
    """Parameters of the beam-smearing forward model.

    Defaults are calibrated to the new-epoch M33 maps (see FINDINGS); every
    length is in galaxy-plane pixels of ``pixel_scale_pc``.
    """

    shape: tuple = (2048, 2048)
    pixel_scale_pc: float = 3.0543
    cos_i: float = 0.5724
    hi_beam_fwhm_sky_pc: float = 21.58
    co_beam_fwhm_sky_pc: float = 63.9
    density_slope: float = -3.0          # P(k) slope of ln(Sigma_true)
    density_lnsigma: float = 1.0         # std of ln(Sigma_true)
    hole_fraction: float = 0.0           # filling factor of HI holes (independent GRF, thresholded)
    hole_depth: float = 0.02             # Sigma multiplier at hole centres
    hole_slope: float = -3.5             # P(k) slope of the hole-defining field
    hole_edge: float = 0.25              # softness of hole edges (in units of that field's std)
    median_mom0: float = 216.4           # K km/s, median of observed mom0
    velocity_slope: float = -2.8         # P(k) slope of the true velocity field
    velocity_rms: float = 4.0            # km/s, rms of the true velocity field
    mom0_noise: float = 0.64             # K km/s, beam-correlated
    mom1_noise_mom0: np.ndarray = field(default_factory=lambda: np.array([1.7, 4.2, 7.7, 14.1, 28.3, 56.6, 110, 212, 424, 849]))
    mom1_noise_rms: np.ndarray = field(default_factory=lambda: np.array([2.6, 2.4, 2.1, 2.0, 1.45, 1.08, 0.78, 0.55, 0.43, 0.38]))
    blank_mom0: float = 5.0              # mom1 blanked where observed mom0 < this
    co_alpha: float = 2.0                # CO emissivity ~ Sigma_true**alpha * scatter
    co_lnscatter: float = 0.5            # lognormal scatter of CO emissivity
    co_scatter_fwhm_sky_pc: float = 63.9  # correlation scale of that scatter (cloud occupancy)
    co_fraction: float = 0.139           # fraction of valid pixels that are CO-bright
    inject_shear_scale: float | None = None  # positive control: CO emissivity *= exp(-|true shear|/scale)
    inject_shear_fwhm_pc: float = 0.0        # smooth the true velocity to this FWHM before taking that shear
    inject_relative: bool = False            # if True, scale is in units of the median |true shear|
    inject_component: str = "shear"          # "shear" (d/dy) or "div" (d/dx) for the positive control
    inject_signed: bool = False              # if True, CO *= exp(-grad/scale): converging (grad<0) gas gains CO


def simulate_observation(params: SmearingSimParams, rng) -> dict:
    """One realisation of the beam-smearing forward model.

    The true velocity field is independent of the true HI unless
    ``params.inject_shear_scale`` is set (positive control).

    Returns
    -------
    dict with ``mom0, mom1`` (observed, deprojected-grid maps; mom1 NaN where
    blanked), ``co_bright`` (bool), ``co_map``, ``sigma_true``, ``v_true``.
    """
    P = params
    beam_hi = deprojected_beam_sigma_pix(P.hi_beam_fwhm_sky_pc, P.pixel_scale_pc, P.cos_i)
    beam_co = deprojected_beam_sigma_pix(P.co_beam_fwhm_sky_pc, P.pixel_scale_pc, P.cos_i)
    G = lambda a, b: nd.gaussian_filter(a, b, mode="wrap")

    sig = np.exp(P.density_lnsigma * gaussian_random_field(P.shape, P.density_slope, rng))
    if P.hole_fraction > 0:
        gh = gaussian_random_field(P.shape, P.hole_slope, rng)
        t = np.quantile(gh, P.hole_fraction)
        sig *= P.hole_depth + (1.0 - P.hole_depth) / (1.0 + np.exp(-(gh - t) / P.hole_edge))
    v = P.velocity_rms * gaussian_random_field(P.shape, P.velocity_slope, rng)
    m0 = G(sig, beam_hi)
    scale = P.median_mom0 / np.median(m0)
    sig *= scale
    m0 *= scale
    m1 = G(sig * v, beam_hi) / m0

    mom0 = m0 + P.mom0_noise * correlated_noise(P.shape, beam_hi, rng)
    lm = np.log(np.clip(mom0, 1e-3, None))
    rms_v = np.interp(lm, np.log(P.mom1_noise_mom0), P.mom1_noise_rms)
    mom1 = m1 + rms_v * correlated_noise(P.shape, beam_hi, rng)
    mom1[mom0 < P.blank_mom0] = np.nan

    beam_sc = deprojected_beam_sigma_pix(P.co_scatter_fwhm_sky_pc, P.pixel_scale_pc, P.cos_i)
    emis = sig ** P.co_alpha * np.exp(P.co_lnscatter * correlated_noise(P.shape, beam_sc, rng))
    if P.inject_shear_scale:
        vs = v
        if P.inject_shear_fwhm_pc > 0:
            vs = G(v, P.inject_shear_fwhm_pc / P.pixel_scale_pc * FWHM_TO_SIGMA)
        gsig = np.gradient(vs, axis=0 if P.inject_component == "shear" else 1) / P.pixel_scale_pc / P.cos_i
        s0 = P.inject_shear_scale * (np.median(np.abs(gsig)) if P.inject_relative else 1.0)
        emis = emis * np.exp(-(gsig if P.inject_signed else np.abs(gsig)) / s0)
    co_map = G(emis, beam_co)
    ok = np.isfinite(mom1)
    thr = np.quantile(co_map[ok], 1.0 - P.co_fraction)
    return dict(mom0=mom0, mom1=mom1, co_map=co_map, co_bright=co_map > thr,
                sigma_true=sig, v_true=v)


def simulated_gradients(mom1: np.ndarray, pixel_scale_pc: float, cos_i: float):
    """``(divergence, shear)`` in km/s/pc from a simulated mom1, matching
    :func:`flows.velocity_gradients` (x = minor axis, both divided by cos i)."""
    gy, gx = np.gradient(mom1)
    return gx / pixel_scale_pc / cos_i, gy / pixel_scale_pc / cos_i


# ----------------------------------------------------------------------------
# Rank-based (tercile) drop statistic
# ----------------------------------------------------------------------------
def tercile_windows(gradient: np.ndarray):
    """Bottom- and top-third masks of ``gradient`` (finite values only).

    Thresholds are the 1/3 and 2/3 quantiles of the finite values, so each
    window holds a third of the pixels at every smoothing scale and nothing is
    excluded except the middle third.

    Returns
    -------
    lo, hi : numpy.ndarray of bool
    thresholds : (float, float)
    """
    g = np.asarray(gradient, float)
    ok = np.isfinite(g)
    t1, t2 = np.quantile(g[ok], [1 / 3, 2 / 3])
    return ok & (g <= t1), ok & (g > t2), (float(t1), float(t2))


def quantile_strata(x: np.ndarray, n_strata: int = 10, mask: np.ndarray | None = None) -> np.ndarray:
    """Equal-count stratum index (0..n_strata-1) of ``x`` (e.g. mom0 deciles);
    quantiles are computed over ``mask`` (default: finite ``x``)."""
    x = np.asarray(x, float)
    m = np.isfinite(x) if mask is None else (mask & np.isfinite(x))
    qs = np.quantile(x[m], np.linspace(0, 1, n_strata + 1))
    return np.clip(np.digitize(x, qs[1:-1]), 0, n_strata - 1)


def window_drop_boot(
    co_bright: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    tile_idx: np.ndarray,
    n_tiles: int,
    strata: np.ndarray | None = None,
    p_null: np.ndarray | None = None,
    min_count: int = 50,
    n_boot: int = 1000,
    seed: int | None = 0,
) -> dict:
    """``1 - f_hi / f_lo`` with optional stratification and a tile bootstrap.

    ``f`` is the CO-bright fraction in the ``lo`` / ``hi`` gradient windows.
    With ``strata`` (e.g. HI-column deciles from :func:`quantile_strata`) the
    drop is standardised to the high-window population::

        1 - sum_k n_hi,k f_hi,k / sum_k n_hi,k f_lo,k

    i.e. observed CO-bright count among high-gradient pixels vs the count
    expected if each had the low-gradient CO fraction of its own stratum.
    Strata with fewer than ``min_count`` pixels in either window are skipped.
    With ``p_null`` the same statistic is evaluated on the null probabilities
    and the observed-minus-null excess is bootstrapped on the same tiles.

    Returns
    -------
    dict with ``obs, obs_lo, obs_hi`` (16/84 percentiles) and, if ``p_null`` is
    given, ``null, null_lo, null_hi, excess, excess_lo, excess_hi,
    p_excess_le0``.
    """
    y = np.asarray(co_bright, float)
    lo = np.asarray(lo, bool) & np.isfinite(y)
    hi = np.asarray(hi, bool) & np.isfinite(y)
    y = np.nan_to_num(y)
    s = np.zeros(y.shape, int) if strata is None else np.asarray(strata, int)
    ns = int(s.max()) + 1
    key = np.asarray(tile_idx) * ns + s

    def cnt(w):
        return np.bincount(key, weights=w, minlength=n_tiles * ns).reshape(n_tiles, ns)

    A = [cnt(lo * 1.0), cnt(hi * 1.0), cnt(y * lo), cnt(y * hi)]
    if p_null is not None:
        p = np.nan_to_num(np.asarray(p_null, float))
        A += [cnt(p * lo), cnt(p * hi)]

    def stat(w):
        T = [w @ a for a in A]
        nlo, nhi = T[0], T[1]
        use = (nlo >= min_count) & (nhi >= min_count)
        def drop(ylo, yhi):
            flo = ylo[use] / nlo[use]
            return 1.0 - yhi[use].sum() / (nhi[use] * flo).sum()
        out = [drop(T[2], T[3])]
        if p_null is not None:
            out.append(drop(T[4], T[5]))
        return out

    point = stat(np.ones(n_tiles))
    rng = np.random.default_rng(seed)
    W = rng.multinomial(n_tiles, np.full(n_tiles, 1.0 / n_tiles), size=n_boot).astype(float)
    B = np.array([stat(w) for w in W])
    q = lambda a, k: float(np.nanpercentile(a, k))
    res = dict(obs=float(point[0]), obs_lo=q(B[:, 0], 16), obs_hi=q(B[:, 0], 84))
    if p_null is not None:
        ex = B[:, 0] - B[:, 1]
        res.update(null=float(point[1]), null_lo=q(B[:, 1], 16), null_hi=q(B[:, 1], 84),
                   excess=float(point[0] - point[1]), excess_lo=q(ex, 16), excess_hi=q(ex, 84),
                   p_excess_le0=float(np.mean(ex <= 0)))
    return res


# ----------------------------------------------------------------------------
# M33 calibration
# ----------------------------------------------------------------------------
_M33_MOM1_NOISE_MOM0 = np.array([1.732, 4.243, 7.746, 14.142, 28.284, 56.569, 109.545, 212.132, 424.264, 848.528])
_M33_MOM1_NOISE_WORST = np.array([4.27, 3.41, 2.93, 2.93, 1.85, 1.24, 0.93, 0.36, 0.03, 0.01])


def m33_calibrated_params(noise: str = "worst", **overrides) -> SmearingSimParams:
    """Simulation parameters calibrated to the new-epoch (5.3") M33 maps.

    Matched quantities: the mom0 power-spectrum slope and percentiles
    (patchy lognormal HI with ~25% hole filling factor), the velocity-field
    rms and spectrum, the median |shear| and |div|, the mom0 noise (0.64 K
    km/s), the HI-column -> CO predictability (AUC ~0.8 via
    ``co_lnscatter=2.4``) and the CO-bright fraction (13.9%).

    noise : {"worst", "radiometer"}
        ``"worst"`` sets the mom1 noise vs mom0 so the *observed* sub-beam
        velocity roughness is reproduced entirely by noise (an upper bound on
        the artefact).  ``"radiometer"`` uses ``sigma_v = 0.64 * 40/sqrt(12) /
        mom0`` (capped at 5 km/s), i.e. thermal noise only.
    """
    if noise == "worst":
        rms = _M33_MOM1_NOISE_WORST
    elif noise == "radiometer":
        rms = np.minimum(0.64 * 40 / np.sqrt(12) / _M33_MOM1_NOISE_MOM0, 5.0)
    else:
        raise ValueError("noise must be 'worst' or 'radiometer'")
    kw = dict(shape=(2048, 2048), pixel_scale_pc=3.0543261924, cos_i=0.572432125594591,
              density_slope=-3.2, density_lnsigma=0.8, hole_fraction=0.25, hole_edge=0.08,
              hole_depth=0.003, velocity_slope=-3.0, velocity_rms=4.5, co_lnscatter=2.4,
              mom1_noise_mom0=_M33_MOM1_NOISE_MOM0.copy(), mom1_noise_rms=np.asarray(rms, float))
    kw.update(overrides)
    return SmearingSimParams(**kw)


def window_drop_joint(
    co_bright: np.ndarray,
    windows: list,
    tile_idx: np.ndarray,
    n_tiles: int,
    strata: np.ndarray | None = None,
    min_count: int = 50,
    n_boot: int = 1000,
    seed: int | None = 0,
    resample: str = "bootstrap",
) -> dict:
    """:func:`window_drop_boot` for several (lo, hi) window pairs at once,
    with the *same* tile resampling for all of them, so the bootstrap
    returns the full covariance between the statistics (e.g. across smoothing
    scales).

    Parameters
    ----------
    windows : list of (lo, hi) boolean-mask pairs, one per statistic.

    resample : {"bootstrap", "jackknife"}
        ``"jackknife"`` drops one group (tile / wedge) at a time and adds
        ``jack_sd`` -- the delete-one jackknife standard error
        ``sqrt((G-1)/G * sum (t_g - mean)^2)`` -- which is better behaved than
        the bootstrap when there are only ~10-30 groups.

    Returns
    -------
    dict with ``point`` (n_stat,) and ``boot`` (n_resamples, n_stat); plus
    ``jack_sd`` for the jackknife.
    """
    y = np.nan_to_num(np.asarray(co_bright, float))
    s = np.zeros(y.shape, int) if strata is None else np.asarray(strata, int)
    ns = int(s.max()) + 1
    key = np.asarray(tile_idx) * ns + s
    cnt = lambda w: np.bincount(key, weights=w, minlength=n_tiles * ns).reshape(n_tiles, ns)
    As = [[cnt(lo * 1.0), cnt(hi * 1.0), cnt(y * lo), cnt(y * hi)] for lo, hi in windows]

    def stat(w):
        out = []
        for A in As:
            nlo, nhi, ylo, yhi = (w @ a for a in A)
            u = (nlo >= min_count) & (nhi >= min_count)
            out.append(1.0 - yhi[u].sum() / (nhi[u] * ylo[u] / nlo[u]).sum())
        return out

    point = np.array(stat(np.ones(n_tiles)))
    if resample == "jackknife":
        occ = np.nonzero(np.bincount(np.asarray(tile_idx), minlength=n_tiles))[0]
        J = np.array([stat(np.where(np.arange(n_tiles) == g, 0.0, 1.0)) for g in occ])
        G = len(occ)
        return dict(point=point, boot=J, jack_sd=np.sqrt((G - 1) / G * ((J - J.mean(0)) ** 2).sum(0)))
    rng = np.random.default_rng(seed)
    W = rng.multinomial(n_tiles, np.full(n_tiles, 1.0 / n_tiles), size=n_boot).astype(float)
    return dict(point=point, boot=np.array([stat(w) for w in W]))


def wedge_index(azimuth_rad: np.ndarray, n_wedges: int, offset_rad: float = 0.0) -> np.ndarray:
    """Fixed-angle azimuthal wedge index (0..n_wedges-1) around the galaxy centre."""
    a = np.mod(np.asarray(azimuth_rad, float) - offset_rad, 2 * np.pi)
    return np.minimum((a / (2 * np.pi) * n_wedges).astype(int), n_wedges - 1)


def gls_amplitude(c: np.ndarray, template: np.ndarray, cov: np.ndarray):
    """Generalised-least-squares amplitude ``A`` of ``c ~ A * template`` with
    covariance ``cov``; returns ``(A, sigma_A)``."""
    Ci = np.linalg.inv(cov)
    t = np.asarray(template, float)
    den = t @ Ci @ t
    return float(t @ Ci @ np.asarray(c, float) / den), float(1.0 / np.sqrt(den))


# ----------------------------------------------------------------------------
# In-plane *vector* velocity fields: calibrating the convergence amplitude
# ----------------------------------------------------------------------------
def vector_velocity_field(shape, slope: float, compressive_fraction: float, rng):
    """Gaussian random in-plane velocity field ``(u_x, u_y)`` (unit rms of u_x).

    Built from a scalar potential and a stream function,
    ``u = grad(phi) + curl(psi z)``, each with velocity power spectrum
    ``P(k) ~ k**slope``.  ``compressive_fraction`` is the fraction of the
    velocity variance in the curl-free (compressive) part; 0.5 is the natural
    2-D mix.  Axis convention matches the deprojected maps: ``x`` (axis 1) is
    the minor axis, ``y`` (axis 0) the major axis.
    """
    ny, nx = shape
    ky = np.fft.fftfreq(ny)[:, None] * 2 * np.pi
    kx = np.fft.rfftfreq(nx)[None, :] * 2 * np.pi
    k = np.hypot(kx, ky)
    k[0, 0] = np.inf
    amp = k ** ((slope - 2.0) / 2.0)            # potential spectrum -> velocity ~ k**slope

    def part(compressive):
        f = np.fft.rfft2(rng.standard_normal(shape)) * amp
        # sin(k): the Fourier symbol of a central difference, so that the
        # finite-difference divergence / curl of each part is exactly zero
        sx, sy = np.sin(kx), np.sin(ky)
        if compressive:                          # u = grad(phi)
            ux, uy = 1j * sx * f, 1j * sy * f
        else:                                    # u = (d psi/dy, -d psi/dx)
            ux, uy = 1j * sy * f, -1j * sx * f
        return np.fft.irfft2(ux, s=shape), np.fft.irfft2(uy, s=shape)

    cx, cy = part(True)
    sx, sy = part(False)
    vc = np.var(cx) + np.var(cy)
    vs = np.var(sx) + np.var(sy)
    a = np.sqrt(compressive_fraction / vc)
    b = np.sqrt((1.0 - compressive_fraction) / vs)
    ux, uy = a * cx + b * sx, a * cy + b * sy
    s = ux.std()
    return ux / s, uy / s


def gradient_tensor(ux, uy, pixel_scale_pc: float, fwhm_pc: float = 0.0, mode: str = "wrap") -> dict:
    """Velocity-gradient tensor ``d u_i / d x_j`` (km/s/pc) after Gaussian
    smoothing of ``u`` to ``fwhm_pc``.  Keys ``xx, xy, yx, yy``."""
    if fwhm_pc > 0:
        s = fwhm_pc / pixel_scale_pc * FWHM_TO_SIGMA
        ux = nd.gaussian_filter(ux, s, mode=mode)
        uy = nd.gaussian_filter(uy, s, mode=mode)
    uxy, uxx = np.gradient(ux)
    uyy, uyx = np.gradient(uy)
    d = pixel_scale_pc
    return dict(xx=uxx / d, xy=uxy / d, yx=uyx / d, yy=uyy / d)


def polar_tensor(T: dict, theta) -> dict:
    """Rotate a gradient tensor into the local (R, phi) frame at azimuth
    ``theta`` (from the +x / minor axis): ``RR``, ``pp``, ``Rp``, ``pR`` and
    ``div = RR + pp``."""
    c, s = np.cos(theta), np.sin(theta)
    eR, ep = (c, s), (-s, c)

    def comp(a, b):   # a . G . b  with G_ij = d u_i / d x_j
        return (a[0] * (T["xx"] * b[0] + T["xy"] * b[1]) + a[1] * (T["yx"] * b[0] + T["yy"] * b[1]))

    RR, pp = comp(eR, eR), comp(ep, ep)
    return dict(RR=RR, pp=pp, Rp=comp(eR, ep), pR=comp(ep, eR), div=RR + pp)


def simulate_vector_observation(params: SmearingSimParams, rng, compressive_fraction: float = 0.5,
                                inject: str | None = None, inject_fwhm_pc: float = 57.5,
                                inject_strength: float = 6.0, sin_i: float = 0.8199,
                                los_pattern: tuple = ()) -> dict:
    """Like :func:`simulate_observation`, but with an in-plane vector velocity field.

    The observed line-of-sight velocity is ``v_los = sin i * u_x`` (no vertical
    motions).  To match the data's units (``deproject`` multiplies velocities
    by cos i) the returned ``mom1`` is ``cos i * v_los``, scaled so its rms is
    ``params.velocity_rms``; the package gradient estimators then return
    d v_los / dx, and ``/ sin i`` gives the in-plane d u_x / dx.

    inject : {None, "div", "azimuthal", "radial"}
        CO emissivity is multiplied by ``exp(-T / (inject_strength * median|T|))``
        with ``T`` the true in-plane ``div u``, ``G_phiphi`` or ``G_RR`` of the
        velocity smoothed to ``inject_fwhm_pc`` (converging gas gains CO).
    los_pattern : sequence of (m, amplitude_kms, phase_rad)
        Large-scale residual patterns ``A cos(m theta + phase)`` added to the
        line-of-sight velocity, uncorrelated with CO (rotation-model errors,
        m = 1, or spiral streaming, m = 2), ramped up linearly to 1 kpc.

    Returns the :func:`simulate_observation` keys plus ``ux, uy`` (km/s,
    in-plane, true) and ``theta``, ``radius_pc`` (about the box centre).
    """
    P = params
    beam_hi = deprojected_beam_sigma_pix(P.hi_beam_fwhm_sky_pc, P.pixel_scale_pc, P.cos_i)
    beam_co = deprojected_beam_sigma_pix(P.co_beam_fwhm_sky_pc, P.pixel_scale_pc, P.cos_i)
    G = lambda a, b: nd.gaussian_filter(a, b, mode="wrap")
    ny, nx = P.shape
    yy, xx = np.indices(P.shape)
    theta = np.arctan2(yy - ny // 2, xx - nx // 2)
    radius_pc = np.hypot(yy - ny // 2, xx - nx // 2) * P.pixel_scale_pc

    sig = np.exp(P.density_lnsigma * gaussian_random_field(P.shape, P.density_slope, rng))
    if P.hole_fraction > 0:
        gh = gaussian_random_field(P.shape, P.hole_slope, rng)
        t = np.quantile(gh, P.hole_fraction)
        sig *= P.hole_depth + (1.0 - P.hole_depth) / (1.0 + np.exp(-(gh - t) / P.hole_edge))
    ux, uy = vector_velocity_field(P.shape, P.velocity_slope, compressive_fraction, rng)
    scale_u = P.velocity_rms / (P.cos_i * sin_i)            # so that rms(cos i sin i u_x) = velocity_rms
    ux, uy = ux * scale_u, uy * scale_u
    v = P.cos_i * sin_i * ux
    for m, amp, ph in los_pattern:
        v = v + P.cos_i * amp * np.cos(m * theta + ph) * np.clip(radius_pc / 1000.0, 0, 1)

    m0 = G(sig, beam_hi)
    scale = P.median_mom0 / np.median(m0)
    sig *= scale
    m0 *= scale
    m1 = G(sig * v, beam_hi) / m0
    mom0 = m0 + P.mom0_noise * correlated_noise(P.shape, beam_hi, rng)
    lm = np.log(np.clip(mom0, 1e-3, None))
    rms_v = np.interp(lm, np.log(P.mom1_noise_mom0), P.mom1_noise_rms)
    mom1 = m1 + rms_v * correlated_noise(P.shape, beam_hi, rng)
    mom1[mom0 < P.blank_mom0] = np.nan

    beam_sc = deprojected_beam_sigma_pix(P.co_scatter_fwhm_sky_pc, P.pixel_scale_pc, P.cos_i)
    emis = sig ** P.co_alpha * np.exp(P.co_lnscatter * correlated_noise(P.shape, beam_sc, rng))
    if inject is not None:
        T = polar_tensor(gradient_tensor(ux, uy, P.pixel_scale_pc, inject_fwhm_pc), theta)
        tr = {"div": T["div"], "azimuthal": T["pp"], "radial": T["RR"]}[inject] * 1e3
        emis = emis * np.exp(-tr / (inject_strength * np.median(np.abs(tr))))
    co_map = G(emis, beam_co)
    ok = np.isfinite(mom1)
    thr = np.quantile(co_map[ok], 1.0 - P.co_fraction)
    return dict(mom0=mom0, mom1=mom1, co_map=co_map, co_bright=co_map > thr, sigma_true=sig,
                v_true=v, ux=ux, uy=uy, theta=theta, radius_pc=radius_pc)


def azimuthal_excess(field_vals, co_bright, mom0, radius_pc, azimuth, ring_edges_pc=(1200, 2000, 3000, 4000),
                     n_wedges: int = 24, n_hi: int = 10):
    """Per-wedge excess of ``field_vals`` at CO-bright pixels over HI-matched pixels.

    Pixels are grouped into (wedge x radial ring x HI-column decile) strata;
    the excess of a CO-bright pixel is its value minus the mean over *all*
    pixels in its stratum.  Returns per-wedge summed excess ``E`` and CO-bright
    count ``N`` (so ``E / N`` is the mean excess in each wedge).
    """
    g = np.asarray(field_vals, float)
    y = np.asarray(co_bright, bool)
    edges = np.asarray(ring_edges_pc, float)
    nr = len(edges) - 1
    rb = np.clip(np.digitize(radius_pc, edges[1:-1]), 0, nr - 1)
    w = wedge_index(azimuth, n_wedges)
    ok = np.isfinite(g)
    hd = np.zeros(g.size, int)
    for k in range(nr):
        m = rb == k
        if m.sum() > n_hi * 10:
            hd[m] = quantile_strata(mom0[m], n_hi)
    key = (w * nr + rb) * n_hi + hd
    K = n_wedges * nr * n_hi
    n_all = np.bincount(key[ok], minlength=K)
    s_all = np.bincount(key[ok], weights=g[ok], minlength=K)
    mean_all = np.where(n_all > 0, s_all / np.maximum(n_all, 1), 0.0)
    yo = ok & y
    n_co = np.bincount(key[yo], minlength=K)
    s_co = np.bincount(key[yo], weights=g[yo], minlength=K)
    wk = np.arange(K) // (nr * n_hi)
    E = np.bincount(wk, weights=s_co - n_co * mean_all, minlength=n_wedges)
    N = np.bincount(wk, weights=n_co, minlength=n_wedges).astype(float)
    return E, N


def azimuthal_fit(E, N, jackknife: bool = True) -> dict:
    """CO-weighted fit of the per-wedge mean excess to ``A + B cos2t + C sin2t``.

    For ``field_vals = d u_x / d x`` (in-plane), ``2A`` is the in-plane
    divergence, ``A + B`` the radial term and ``A - B`` the azimuthal term.
    Returns ``A, B, C`` and (if ``jackknife``) delete-one-wedge ``*_sd``.
    """
    E, N = np.asarray(E, float), np.asarray(N, float)
    nw = E.size
    th = (np.arange(nw) + 0.5) * 2 * np.pi / nw - np.pi

    def fit(m):
        e = E[m] / np.maximum(N[m], 1)
        X = np.c_[np.ones(m.sum()), np.cos(2 * th[m]), np.sin(2 * th[m])]
        W = N[m]
        return np.linalg.solve((X * W[:, None]).T @ X, (X * W[:, None]).T @ e)

    allm = np.ones(nw, bool)
    b = fit(allm)
    out = dict(A=b[0], B=b[1], C=b[2])
    if jackknife:
        J = np.array([fit(np.arange(nw) != k) for k in range(nw)])
        sd = np.sqrt((nw - 1) / nw * ((J - J.mean(0)) ** 2).sum(0))
        out.update(A_sd=sd[0], B_sd=sd[1], C_sd=sd[2])
    return out
