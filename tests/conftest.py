"""Shared fixtures: a small synthetic Projection with a real WCS.

Building a genuine :class:`~spectral_cube.Projection` (rather than a mock)
keeps the geometry tests honest -- the deprojection uses
``spatial_coordinate_map``, which needs a real WCS.
"""

import numpy as np
import pytest
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import Projection

from exgal_flows import GalaxyConfig


@pytest.fixture
def galaxy():
    return GalaxyConfig()


def _make_header(ny, nx, center_ra, center_dec, cdelt_deg):
    w = WCS(naxis=2)
    w.wcs.crpix = [nx / 2, ny / 2]
    w.wcs.cdelt = [-cdelt_deg, cdelt_deg]
    w.wcs.crval = [center_ra, center_dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    hdr = w.to_header()
    hdr["CDELT2"] = cdelt_deg
    return hdr


@pytest.fixture
def synthetic_velocity(galaxy):
    """A velocity field with a clean, known x-gradient (pure divergence).

    ``v = slope * x`` in pixels, so ``d v / d x`` is a constant we can check.
    """
    ny = nx = 101
    cdelt = 1.0 / 3600.0  # 1 arcsec/pix
    hdr = _make_header(ny, nx, galaxy.center_ra_deg, galaxy.center_dec_deg, cdelt)
    xx = np.arange(nx)[None, :] * np.ones((ny, 1))
    slope = 0.5  # km/s per pixel
    data = (slope * (xx - nx / 2)).astype(float)
    proj = Projection(data * u.km / u.s, wcs=WCS(hdr), header=hdr)
    return proj, slope, cdelt


@pytest.fixture
def synthetic_flat(galaxy):
    """A flat, all-ones intensity map for deprojection footprint tests."""
    ny = nx = 101
    cdelt = 1.0 / 3600.0
    hdr = _make_header(ny, nx, galaxy.center_ra_deg, galaxy.center_dec_deg, cdelt)
    proj = Projection(np.ones((ny, nx)) * u.K, wcs=WCS(hdr), header=hdr)
    return proj
