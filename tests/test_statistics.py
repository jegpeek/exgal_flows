"""Tests for the statistics module -- the scientific claims in miniature.

These use synthetic data with a *known* injected signal, so they verify the
estimators recover the right sign and that the bootstrap behaves.
"""

import numpy as np
import pytest

from exgal_flows import statistics as st


def test_assign_wedges_covers_all_pixels():
    az = np.linspace(-np.pi, np.pi, 1000)
    w = st.assign_wedges(az, n_wedges=8)
    assert w.min() == 0
    assert w.max() == 7
    assert set(np.unique(w)) == set(range(8))


def test_wedge_bootstrap_recovers_median():
    rng = np.random.default_rng(0)
    vals = rng.normal(loc=-5.0, scale=1.0, size=8000)
    wedges = rng.integers(0, 8, size=8000)
    boot = st.wedge_bootstrap_median(vals, wedges, n_wedges=8, n_boot=500, rng=rng)
    # bootstrap median should sit near the true median (-5)
    assert abs(np.nanmedian(boot) - (-5.0)) < 0.3
    # and produce a non-degenerate spread
    assert np.nanstd(boot) > 0


def test_co_fraction_responds_to_threshold():
    # CO brightness anti-correlated with the gradient: high gradient -> low CO
    n = 20000
    rng = np.random.default_rng(1)
    grad = rng.uniform(0, 0.2, n)
    co = 2.0 * (grad < 0.05) + rng.normal(0, 0.01, n)  # bright only at low grad
    edges = np.arange(0, 0.2, 0.01)
    curve = st.co_fraction_vs_gradient(grad, co, co_threshold=0.6, bin_edges=edges)
    # first bin (low gradient) should have far more CO than the last
    assert curve.fraction[0] > 0.8
    assert curve.fraction[-1] < 0.2


def test_divergence_of_co_bright_sign():
    # inject: CO-bright pixels have negative divergence (colliding flows)
    n = 40000
    rng = np.random.default_rng(2)
    div = rng.normal(0, 0.05, n)              # km/s/pc
    az = rng.uniform(-np.pi, np.pi, n)
    # CO brightness higher where divergence is more negative
    co = 0.6 * (1.0 + np.clip(-div / 0.05, 0, None)) + rng.normal(0, 0.05, n)
    res = st.divergence_of_co_bright(
        div, co, az, co_threshold=0.6, factors=(2.0, 3.0), n_wedges=8, n_boot=200, seed=0
    )
    # brightest bin should have a negative median divergence
    assert res.medians[-1] < 0
    # error bars ordered
    assert np.all(res.err_low <= res.medians + 1e-6)
    assert np.all(res.err_high >= res.medians - 1e-6)


def test_bootstrap_is_reproducible():
    vals = np.random.default_rng(3).normal(size=2000)
    wedges = np.random.default_rng(4).integers(0, 8, size=2000)
    a = st.wedge_bootstrap_median(vals, wedges, 8, 100, rng=np.random.default_rng(42))
    b = st.wedge_bootstrap_median(vals, wedges, 8, 100, rng=np.random.default_rng(42))
    np.testing.assert_array_equal(a, b)


def test_weighted_median_matches_repeat():
    from exgal_flows.statistics import _weighted_median
    import numpy as np
    v = np.array([1.0, 2.0, 3.0, 4.0])
    w = np.array([1.0, 1.0, 1.0, 1.0])
    # equal weights -> ordinary median convention (lower-middle here)
    assert _weighted_median(v, w) in (2.0, 3.0)
    # heavy weight pins the median
    assert _weighted_median(v, np.array([0.0, 0.0, 10.0, 0.0])) == 3.0


def test_tiled_block_bootstrap_recovers_radial_sign_flip():
    """Synthetic map: inner pixels positive divergence, outer negative.

    The point-estimate profile must cross zero, and the block-bootstrap
    band must be at least as wide as the iid band.
    """
    import numpy as np
    from exgal_flows.statistics import tiled_block_bootstrap_radial_profile

    rng = np.random.default_rng(1)
    n = 4000
    i_pix = rng.integers(0, 400, size=n)
    j_pix = rng.integers(0, 400, size=n)
    # radius in pc from a corner, pixel scale 3 pc
    px = 3.0
    radius_pc = np.hypot(i_pix, j_pix) * px
    # inner positive, outer negative divergence
    div = np.where(radius_pc < 600, 5.0, -5.0) + rng.normal(0, 1, size=n)

    prof = tiled_block_bootstrap_radial_profile(
        i_pix, j_pix, div, radius_pc, px,
        r_edges=np.array([0, 300, 600, 900, 1200, 1800]),
        block_pc=150.0, n_boot=100, seed=0,
    )
    # crosses zero somewhere in the map
    assert np.isfinite(prof.turnover_pc)
    assert 300 < prof.turnover_pc < 1200
    # block band no narrower than iid band on average (correlation widens it)
    block_w = np.nanmean(prof.err_high - prof.err_low)
    iid_w = np.nanmean(prof.iid_high - prof.iid_low)
    assert block_w >= 0.5 * iid_w


def test_co_fraction_vs_gradient_boot_recovers_suppression():
    """Injected shear-suppression signal: CO-bright fraction should fall with
    gradient, the block-bootstrap band should be finite, and the scalar
    sensitivity (rel_drop) should be clearly positive."""
    n = 40000
    rng = np.random.default_rng(3)
    grad = rng.uniform(0.0, 0.2, n)
    az = rng.uniform(-np.pi, np.pi, n)
    # probability of being CO-bright falls linearly with gradient
    p_bright = np.clip(0.5 * (1.0 - grad / 0.2), 0, 1)
    co = np.where(rng.uniform(size=n) < p_bright, 1.0, 0.0)  # 1 K = bright, 0 = not
    edges = np.arange(0.0, 0.2, 0.01)
    res = st.co_fraction_vs_gradient_boot(
        grad, co, az, co_threshold=0.6, bin_edges=edges,
        n_wedges=8, n_boot=200, seed=0,
    )
    # fraction falls from first to last bin
    assert res.fraction[0] > res.fraction[-1]
    # sensitivity is positive and its band brackets the point estimate
    assert res.rel_drop > 0.2
    assert res.rel_drop_lo <= res.rel_drop <= res.rel_drop_hi
    # error band is finite and ordered
    good = np.isfinite(res.err_low) & np.isfinite(res.err_high)
    assert good.sum() > 5
    assert np.all(res.err_high[good] >= res.err_low[good])


def test_co_fraction_vs_gradient_boot_mask_restricts():
    """A radius-like mask should reduce the pixel counts per bin."""
    n = 10000
    rng = np.random.default_rng(4)
    grad = rng.uniform(0.0, 0.2, n)
    az = rng.uniform(-np.pi, np.pi, n)
    co = rng.uniform(0, 2, n)
    edges = np.arange(0.0, 0.2, 0.01)
    full = st.co_fraction_vs_gradient_boot(grad, co, az, 0.6, edges, n_boot=50, seed=0)
    half = st.co_fraction_vs_gradient_boot(
        grad, co, az, 0.6, edges, n_boot=50, seed=0,
        mask=(np.arange(n) < n // 2),
    )
    assert half.counts.sum() < full.counts.sum()


@pytest.mark.parametrize("method", ["ols", "l1"])
def test_block_bootstrap_line_fit_recovers_injected_line(method):
    """Both estimators recover a known intercept/slope and its sign."""
    rng = np.random.default_rng(0)
    # 200 tiles, ~50 pixels each, radius 0-4000 pc, div = 6 - 4*r_kpc + noise
    n_tiles = 200
    tile_idx = np.repeat(np.arange(n_tiles), 50)
    r_pc = rng.uniform(0, 4000, tile_idx.size)
    a_true, b_true = 6.0, -4.0
    div = a_true + b_true * (r_pc / 1000.0) + rng.normal(0, 3, r_pc.size)
    fit = st.block_bootstrap_line_fit(
        r_pc, div, tile_idx, n_tiles, n_boot=200, seed=1, method=method)
    assert fit.method == method
    # point estimates near truth
    assert abs(fit.intercept - a_true) < 1.5
    assert abs(fit.slope - b_true) < 1.5
    # bands ordered and bracket the point estimate
    assert fit.slope_lo <= fit.slope <= fit.slope_hi
    assert fit.intercept_lo <= fit.intercept <= fit.intercept_hi
    # slope is confidently negative; zero-crossing near 1500 pc
    assert fit.slope_neg_prob > 0.9
    assert 1000 < fit.zero_crossing_pc < 2000


def test_block_bootstrap_line_fit_rejects_bad_method():
    rng = np.random.default_rng(0)
    tile_idx = np.repeat(np.arange(10), 5)
    r = rng.uniform(0, 4000, tile_idx.size)
    with pytest.raises(ValueError):
        st.block_bootstrap_line_fit(r, r, tile_idx, 10, n_boot=5, method="bogus")


def test_block_bootstrap_sigmoid_fit_recovers_step():
    """Injected inner-positive / outer-negative step is recovered, straddles zero."""
    rng = np.random.default_rng(3)
    n_tiles = 200
    per = 50
    # tile centres spread across 0-4000 pc; pixels inherit tile radius + jitter
    tile_r = rng.uniform(0, 4000, size=n_tiles)
    r = np.repeat(tile_r, per) + rng.normal(0, 40, size=n_tiles * per)
    r = np.clip(r, 0, 4000)
    tile_idx = np.repeat(np.arange(n_tiles), per)
    hi_true, lo_true, r0_true, w_true = 6.0, -5.0, 1.2, 0.15  # kpc
    truth = lo_true + (hi_true - lo_true) / (1.0 + np.exp((r / 1000.0 - r0_true) / w_true))
    y = truth + rng.normal(0, 3.0, size=r.size)
    fit = st.block_bootstrap_sigmoid_fit(r, y, tile_idx, n_tiles, n_boot=150, seed=0)
    assert fit.hi > 0 and fit.lo < 0                    # correct signs
    assert 800 < fit.zero_crossing_pc < 1600            # crosses near truth
    assert fit.amplitude_lo < fit.amplitude < fit.amplitude_hi
    assert fit.straddle_zero_prob > 0.9                 # decisively straddles
    assert fit.method == "sigmoid"
