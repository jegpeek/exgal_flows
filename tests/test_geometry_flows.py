"""Tests for deprojection, polar grid, and velocity-gradient physics."""

import numpy as np
import astropy.units as u

from exgal_flows import geometry, flows
from exgal_flows.io import pixel_scale_pc


def test_pixel_scale_pc():
    class H(dict):
        pass
    hdr = {"CDELT2": 1.0 / 3600.0}  # 1 arcsec/pix
    # at 840 kpc, 1 arcsec ~ 4.07 pc
    px = pixel_scale_pc(hdr, 840.0)
    assert np.isclose(px, np.radians(1 / 3600) * 840000, rtol=1e-6)
    assert 4.0 < px < 4.1


def test_deproject_footprint_finite(synthetic_flat, galaxy):
    out = geometry.deproject(synthetic_flat, galaxy)
    # some pixels finite (the galaxy footprint), some NaN (outside)
    assert np.isfinite(out).any()
    assert (~np.isfinite(out)).any()


def test_polar_grid_center():
    r, az = geometry.polar_grid((101, 101))
    # radius minimal at the centre
    cy, cx = np.unravel_index(np.argmin(r), r.shape)
    assert abs(cy - 50) <= 1 and abs(cx - 50) <= 1
    assert az.min() >= -np.pi - 1e-9
    assert az.max() <= np.pi + 1e-9


def test_velocity_gradients_recover_slope(synthetic_velocity, galaxy):
    proj, slope, cdelt = synthetic_velocity
    px = pixel_scale_pc(proj.header, galaxy.distance_kpc)
    # feed the raw field (no deprojection) so the gradient is exactly known
    ff = flows.velocity_gradients(proj.value, px, galaxy)
    # d v / d x = slope[km/s/pix] / px[pc/pix] / cos i, away from edges
    expected = slope / px / galaxy.cos_inc
    interior = ff.divergence.value[5:-5, 5:-5]
    assert np.allclose(interior, expected, rtol=1e-6)
    # pure x-gradient => shear (y-gradient) ~ 0
    assert np.allclose(ff.shear.value[5:-5, 5:-5], 0.0, atol=1e-9)


def test_hi_surface_density_units():
    sd = flows.hi_surface_density(np.ones((3, 3)) * 100.0, factor=0.0196)
    assert sd.unit == u.solMass / u.pc**2
    assert np.isclose(sd.value[0, 0], 1.96)


def test_smoothed_velocity_zero_sigma_identity():
    m0 = np.ones((10, 10))
    m1 = np.arange(100).reshape(10, 10).astype(float)
    out = flows.smoothed_velocity_field(m0, m1, 0)
    np.testing.assert_array_equal(out, m1)



def _holey_field(seed=0, n=160):
    """Offset velocity field (far from 0) with isolated blanks and one big hole."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:n, :n]
    truth = 250.0 + 0.2 * xx                       # km/s, gradient 0.2 km/s/pix
    m1 = truth + rng.normal(0, 1.0, (n, n))
    m0 = 50.0 + 10.0 * rng.random((n, n))
    bad = rng.random((n, n)) < 0.002               # isolated blank pixels
    bad[50:100, 50:100] = True                     # one large hole
    m1[bad] = np.nan
    m0[bad] = np.nan
    return m0, m1, truth, bad


def test_smoothed_velocity_blanks_get_zero_weight():
    """Blanks drop out of the average: no bias toward 0, no NaN spreading."""
    m0, m1, truth, bad = _holey_field()
    out = flows.smoothed_velocity_field(m0, m1, 6.0)
    good = ~bad
    kept = np.isfinite(out) & good
    # isolated blanks used to wipe out an ~8 sigma square each (the whole map
    # here); now only the large hole's rim is lost to the coverage cut
    assert kept.sum() / good.sum() > 0.95
    # values stay at the 250 km/s level instead of being pulled toward 0
    assert np.max(np.abs(out[kept] - truth[kept])) < 3.0
    # blank pixels are not in-painted
    assert not np.isfinite(out[bad]).any()
    # the gradient (what divergence / shear use) is unbiased
    gx = np.gradient(out, axis=1)
    assert abs(np.nanmedian(gx[good]) - 0.2) < 0.01


def test_smoothed_velocity_coverage_cut():
    """Pixels with a mostly-masked kernel are blanked; min_coverage=0 keeps them."""
    m0, m1, truth, bad = _holey_field()
    strict = flows.smoothed_velocity_field(m0, m1, 6.0, min_coverage=0.9)
    loose = flows.smoothed_velocity_field(m0, m1, 6.0, min_coverage=0.0)
    assert np.isfinite(loose[~bad]).all()
    lost = ~bad & np.isfinite(loose) & ~np.isfinite(strict)
    assert lost.any()
    # the pixels cut sit on the rim of the large hole
    ii, jj = np.nonzero(lost)
    rim = (ii >= 40) & (ii < 110) & (jj >= 40) & (jj < 110)
    assert rim.mean() > 0.95


def test_smoothed_velocity_matches_plain_smoothing_without_blanks():
    """With no blanks the result is the plain mom0-weighted Gaussian smooth."""
    from scipy import ndimage as nd
    rng = np.random.default_rng(1)
    m0 = 1.0 + rng.random((60, 60))
    m1 = rng.normal(0, 5, (60, 60))
    out = flows.smoothed_velocity_field(m0, m1, 3.0)
    ref = nd.gaussian_filter(m0 * m1, 3.0) / nd.gaussian_filter(m0, 3.0)
    np.testing.assert_allclose(out, ref, rtol=1e-12, atol=1e-12)

def test_fwhm_pc_conversion():
    assert np.isclose(flows.fwhm_pc(2.0, 3.0), 2.0 * 3.0 * 2.355)


def test_highpass_removes_large_scale_keeps_small_scale():
    ny, nx = 400, 400
    yy, xx = np.mgrid[:ny, :nx].astype(float)
    plane = 0.5 * xx - 0.2 * yy + 30.0           # large-scale "rotation"
    small = 2.0 * np.sin(2 * np.pi * xx / 12.0)  # 12-pixel ripple
    v = plane + small
    v[100:110, 200:215] = np.nan                 # a blanked hole
    hp = flows.highpass_velocity_field(v, sigma_pix=40.0)
    assert np.isnan(hp[105, 205])
    inner = (slice(150, 250), slice(150, 250))
    # plane (range ~200) gone, ripple preserved; the blanked hole skews the
    # normalised low-pass of a sloped field slightly, hence the 0.2 allowance
    resid = hp[inner] - small[inner]
    assert np.nanmax(np.abs(resid)) < 0.2
    gx = np.gradient(hp, axis=1)[inner]
    gx_true = np.gradient(small, axis=1)[inner]
    assert np.nanmax(np.abs(gx - gx_true)) < 0.02


def test_lowpass_block_matches_full_resolution():
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[:240, :240].astype(float)
    v = 10 * np.sin(xx / 60.0) + 5 * np.cos(yy / 45.0) + rng.normal(0, 1, xx.shape)
    v[50:60, 50:60] = np.nan
    full = flows.lowpass_velocity_field(v, 32.0, block=1)
    fast = flows.lowpass_velocity_field(v, 32.0, block=4)
    inner = (slice(40, 200), slice(40, 200))
    assert np.nanmax(np.abs(full[inner] - fast[inner])) < 0.05
