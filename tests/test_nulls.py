"""Tests for the CO-vs-|gradient| null models."""
import numpy as np
import pytest

from exgal_flows import nulls as nl
from exgal_flows import statistics as st


def test_grf_unit_variance_and_slope():
    rng = np.random.default_rng(0)
    f = nl.gaussian_random_field((256, 256), -3.0, rng)
    assert abs(f.mean()) < 1e-10 and abs(f.std() - 1) < 1e-10
    F = np.abs(np.fft.fft2(f)) ** 2
    k = np.hypot(*np.meshgrid(np.fft.fftfreq(256), np.fft.fftfreq(256)))
    sel = (k > 0.02) & (k < 0.3)
    slope = np.polyfit(np.log(k[sel]), np.log(F[sel]), 1)[0]
    assert slope == pytest.approx(-3.0, abs=0.15)


def test_deprojected_beam_is_stretched_along_x():
    sy, sx = nl.deprojected_beam_sigma_pix(21.58, 3.054, 0.5724)
    assert sx / sy == pytest.approx(1 / 0.5724)
    assert sy * 2.3548 == pytest.approx(21.58 / 3.054, rel=1e-3)


def test_relative_drop_vs_null_zero_excess_when_null_is_truth():
    rng = np.random.default_rng(1)
    n = 40000
    g = rng.uniform(0, 0.19, n)
    p = np.clip(0.3 - 1.0 * g, 0.01, 1)          # CO suppressed at high g ...
    y = rng.random(n) < p                         # ... and the null knows it exactly
    tiles = rng.integers(0, 200, n)
    d = nl.relative_drop_vs_null(g, y, p, tiles, 200, n_boot=200)
    assert d["obs"] > 0.3
    assert abs(d["excess"]) < 3 * (d["excess_hi"] - d["excess_lo"])
    assert 0.8 < d["frac_explained"] < 1.25


def test_relative_drop_vs_null_detects_excess():
    rng = np.random.default_rng(2)
    n = 40000
    g = rng.uniform(0, 0.19, n)
    y = rng.random(n) < np.clip(0.3 - 1.0 * g, 0.01, 1)
    p = np.full(n, y.mean())                      # null with no g dependence
    tiles = rng.integers(0, 200, n)
    d = nl.relative_drop_vs_null(g, y, p, tiles, 200, n_boot=200)
    assert d["null"] == pytest.approx(0.0, abs=1e-12)
    assert d["excess"] > 0.3 and d["p_excess_le0"] < 0.01


def test_simulation_shapes_and_co_fraction():
    P = nl.SmearingSimParams(shape=(128, 128), hole_fraction=0.2)
    o = nl.simulate_observation(P, np.random.default_rng(3))
    assert o["mom0"].shape == (128, 128)
    ok = np.isfinite(o["mom1"])
    assert o["co_bright"][ok].mean() == pytest.approx(P.co_fraction, abs=0.01)
    dv, sh = nl.simulated_gradients(o["mom1"], P.pixel_scale_pc, P.cos_i)
    assert dv.shape == sh.shape == (128, 128)


def test_hi_structure_features_sign_convention():
    yy, xx = np.indices((64, 64))
    bump = np.exp(-((xx - 32) ** 2 + (yy - 32) ** 2) / 50.0)
    F = nl.hi_structure_features(bump, sigmas_pix=(0, 2))
    assert F["Lap_s0"][32, 32] < 0 and F["Lap_s2"][32, 32] < 0   # peak -> negative Laplacian


def test_tercile_windows_hold_a_third_each():
    g = np.random.default_rng(4).exponential(size=9000)
    lo, hi, (t1, t2) = nl.tercile_windows(g)
    assert lo.mean() == pytest.approx(1 / 3, abs=1e-3) and hi.mean() == pytest.approx(1 / 3, abs=1e-3)
    assert t1 < t2


def test_stratified_window_drop_removes_confounder():
    """CO depends only on HI column; gradient anti-correlates with column.
    The raw tercile drop is large, the stratified one ~0."""
    rng = np.random.default_rng(5)
    n = 60000
    col = rng.lognormal(0, 1, n)
    g = 1.0 / col * rng.lognormal(0, 0.5, n)               # faint HI -> big gradients
    y = rng.random(n) < np.clip(0.05 * col, 0, 0.9)        # CO from HI only
    tiles = rng.integers(0, 300, n)
    lo, hi, _ = nl.tercile_windows(g)
    raw = nl.window_drop_boot(y, lo, hi, tiles, 300, n_boot=100)
    strat = nl.window_drop_boot(y, lo, hi, tiles, 300, strata=nl.quantile_strata(col, 20), n_boot=100)
    assert raw["obs"] > 0.5
    assert abs(strat["obs"]) < 0.1


def test_signed_injection_enriches_converging_gas():
    P = nl.m33_calibrated_params(shape=(256, 256), inject_shear_scale=1.0, inject_relative=True,
                                 inject_component="div", inject_signed=True, mom1_noise_rms=np.full(10, 0.01))
    o = nl.simulate_observation(P, np.random.default_rng(0))
    gx = np.gradient(o["v_true"], axis=1)
    conv, dvg = gx < np.quantile(gx, 1 / 3), gx > np.quantile(gx, 2 / 3)
    assert o["co_map"][conv].mean() > o["co_map"][dvg].mean()


def test_window_drop_joint_matches_single_and_gls():
    rng = np.random.default_rng(6)
    n = 30000
    g = rng.normal(size=n); y = rng.random(n) < 0.2 - 0.05 * np.tanh(g)
    tiles = rng.integers(0, 150, n)
    lo, hi, _ = nl.tercile_windows(g)
    j = nl.window_drop_joint(y, [(lo, hi), (lo, hi)], tiles, 150, n_boot=50)
    s = nl.window_drop_boot(y, lo, hi, tiles, 150, n_boot=50)
    assert j["point"][0] == pytest.approx(s["obs"]) and j["boot"].shape == (50, 2)
    A, sA = nl.gls_amplitude(np.array([1.0, 3.0]), np.ones(2), np.eye(2))
    assert A == pytest.approx(2.0) and sA == pytest.approx(1 / np.sqrt(2))


def test_jackknife_sd_matches_iid_for_independent_groups():
    rng = np.random.default_rng(8)
    n = 40000
    g = rng.normal(size=n); y = rng.random(n) < 0.2
    groups = rng.integers(0, 20, n)
    lo, hi, _ = nl.tercile_windows(g)
    jk = nl.window_drop_joint(y, [(lo, hi)], groups, 20, resample="jackknife")
    bs = nl.window_drop_joint(y, [(lo, hi)], groups, 20, n_boot=400, seed=1)
    assert jk["boot"].shape == (20, 1)
    assert 0.6 < jk["jack_sd"][0] / bs["boot"].std() < 1.6


def test_wedge_index_range():
    w = nl.wedge_index(np.linspace(-np.pi, np.pi, 1000), 12)
    assert w.min() == 0 and w.max() == 11 and len(np.unique(w)) == 12


def test_vector_field_compressive_fraction_and_tensor():
    rng = np.random.default_rng(9)
    ux, uy = nl.vector_velocity_field((256, 256), -3.0, 1.0, rng)          # purely compressive
    T = nl.gradient_tensor(ux, uy, 1.0, 0.0)
    sl = (slice(2, -2),) * 2                                                  # np.gradient is one-sided at the edges
    curl = (T["yx"] - T["xy"])[sl]
    assert np.std(curl) < 0.05 * np.std((T["xx"] + T["yy"])[sl])
    ux, uy = nl.vector_velocity_field((256, 256), -3.0, 0.0, rng)          # purely solenoidal
    T = nl.gradient_tensor(ux, uy, 1.0, 0.0)
    assert np.std((T["xx"] + T["yy"])[sl]) < 0.05 * np.std((T["yx"] - T["xy"])[sl])


def test_polar_tensor_radial_inflow():
    yy, xx = np.indices((101, 101)) - 50.0
    ux, uy = -0.1 * xx, -0.1 * yy                       # uniform radial inflow: G = -0.1 I
    T = nl.polar_tensor(nl.gradient_tensor(ux, uy, 1.0, 0.0, mode="nearest"), np.arctan2(yy, xx))
    sl = (slice(10, 90),) * 2
    assert np.allclose(T["RR"][sl], -0.1) and np.allclose(T["pp"][sl], -0.1) and np.allclose(T["div"][sl], -0.2)


def test_azimuthal_fit_recovers_isotropic_and_anisotropic_terms():
    rng = np.random.default_rng(10)
    n = 200000
    th = rng.uniform(-np.pi, np.pi, n); R = rng.uniform(1200, 4000, n); m0 = rng.lognormal(5, 1, n)
    y = rng.random(n) < 0.15
    A, B = -2.0, 1.5
    g = rng.normal(0, 5, n) + y * (A + B * np.cos(2 * th))     # excess only at CO-bright pixels
    E, N = nl.azimuthal_excess(g, y, m0, R, th)
    f = nl.azimuthal_fit(E, N)
    assert abs(f["A"] - A * 0.85) < 0.3 and abs(f["B"] - B * 0.85) < 0.3   # 0.85 = 1 - CO fraction (stratum mean includes CO px)
