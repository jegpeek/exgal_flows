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


def test_fwhm_pc_conversion():
    assert np.isclose(flows.fwhm_pc(2.0, 3.0), 2.0 * 3.0 * 2.355)
