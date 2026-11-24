"""Velocity-gradient decomposition: divergence and shear fields.

Given a deprojected, model-subtracted velocity field in the galaxy plane,
we form the two plane-of-galaxy velocity gradients and identify:

* **divergence along the minor axis** -- the x-gradient
  (``d v / d x``).  Negative values are *converging* / colliding flows;
  positive values are *diverging* flows.
* **shear along the minor axis** -- the y-gradient (``d v / d y``).

Both are divided by ``cos i`` exactly as in the notebook.  As explained in
:mod:`exgal_flows.geometry`, for the velocity field this factor cancels the
``cos i`` applied during deprojection, so the net gradient is unaffected by
that convention; the ``1/cos i`` x-stretch during deprojection is the
operation that actually matters.

Units: gradients are km/s/pc.  Multiply by 1000 for the km/s/kpc used on
most of the notebook's axes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import astropy.units as u
from scipy import ndimage as nd

from .config import GalaxyConfig


@dataclass
class FlowFields:
    """Velocity-gradient fields derived from a deprojected velocity map.

    All arrays share the shape of the deprojected velocity field.

    Attributes
    ----------
    divergence : astropy.units.Quantity
        Minor-axis divergence ``d v / d x`` in km/s/pc.  Negative =
        converging / colliding.
    shear : astropy.units.Quantity
        Minor-axis shear ``d v / d y`` in km/s/pc.
    grad_x, grad_y : astropy.units.Quantity
        Raw x and y velocity gradients in km/s/pc (before the ``1/cos i``).
    """

    divergence: u.Quantity
    shear: u.Quantity
    grad_x: u.Quantity
    grad_y: u.Quantity


def velocity_gradients(
    velocity_deproj,
    pixel_scale_pc: float,
    galaxy: GalaxyConfig,
) -> FlowFields:
    """Compute divergence and shear from a deprojected velocity field.

    Parameters
    ----------
    velocity_deproj : array-like or Quantity
        Deprojected, model-subtracted velocity field (km/s).  A plain
        ndarray is interpreted as km/s.
    pixel_scale_pc : float
        Physical pixel size in parsecs.
    galaxy : GalaxyConfig
        Galaxy geometry (uses the inclination for the ``1/cos i`` factor).

    Returns
    -------
    FlowFields
    """
    v = u.Quantity(velocity_deproj, u.km / u.s)
    dx = pixel_scale_pc * u.pc

    grad_x = np.gradient(v, axis=1) / dx
    grad_y = np.gradient(v, axis=0) / dx

    cos_i = galaxy.cos_inc
    divergence = grad_x / cos_i
    shear = grad_y / cos_i

    return FlowFields(divergence=divergence, shear=shear, grad_x=grad_x, grad_y=grad_y)


def hi_surface_density(mom0_deproj, factor: float = 0.0196) -> u.Quantity:
    """HI surface density from a deprojected moment-0 map.

    Parameters
    ----------
    mom0_deproj : array-like or Quantity
        Deprojected HI moment-0 (K km/s).  ``.value`` is used if present.
    factor : float
        K km/s -> solar masses / pc**2 conversion (optically thin 21-cm).

    Returns
    -------
    astropy.units.Quantity
        HI surface density in solar masses per pc**2.
    """
    vals = getattr(mom0_deproj, "value", mom0_deproj)
    return factor * np.asarray(vals) * u.solMass / u.pc**2


def hi_laplacian(mom0_deproj) -> np.ndarray:
    """Laplacian of the deprojected HI moment-0 map (structure proxy)."""
    vals = getattr(mom0_deproj, "value", mom0_deproj)
    return nd.laplace(np.asarray(vals))


def smoothed_velocity_field(mom0_deproj, mom1_deproj, sigma_pix: float) -> np.ndarray:
    """Intensity-weighted smoothed velocity field.

    Reproduces the notebook's resolution-study smoothing:
    ``smooth(mom0 * mom1) / smooth(mom0)`` -- i.e. a moment-0-weighted
    Gaussian smooth of the velocity, which is more robust than smoothing
    the velocity directly where the HI is faint.

    Parameters
    ----------
    mom0_deproj, mom1_deproj : array-like or Quantity
        Deprojected moment-0 and moment-1 maps.
    sigma_pix : float
        Gaussian sigma in pixels.  ``0`` returns the unsmoothed velocity.

    Returns
    -------
    numpy.ndarray
        Smoothed velocity field (km/s, plain ndarray).
    """
    m0 = np.asarray(getattr(mom0_deproj, "value", mom0_deproj), dtype=float)
    m1 = np.asarray(getattr(mom1_deproj, "value", mom1_deproj), dtype=float)
    if sigma_pix == 0:
        return m1
    num = nd.gaussian_filter(m0 * m1, sigma_pix)
    den = nd.gaussian_filter(m0, sigma_pix)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def fwhm_pc(sigma_pix: float, pixel_scale_pc: float) -> float:
    """Convert a Gaussian sigma in pixels to a FWHM in parsecs."""
    return sigma_pix * pixel_scale_pc * 2.355


def normalized_smooth(values, weights, sigma_pix, min_coverage: float = 0.5):
    """Weighted Gaussian smoothing with NaNs given zero weight.

    ``G*(w v) / G*w``, with ``w = 0`` wherever ``values`` or ``weights`` is
    non-finite (so blanked pixels do not propagate NaN), negative weights
    clipped to zero, and the result blanked where the finite-pixel coverage
    ``G*valid`` is below ``min_coverage``.  ``sigma_pix`` may be a scalar or a
    ``(sigma_y, sigma_x)`` pair.
    """
    v = np.asarray(values, float)
    w = np.asarray(weights, float) * np.ones_like(v)
    valid = np.isfinite(v) & np.isfinite(w)
    w = np.where(valid, np.clip(w, 0, None), 0.0)
    vz = np.where(valid, v, 0.0)
    num = nd.gaussian_filter(w * vz, sigma_pix)
    den = nd.gaussian_filter(w, sigma_pix)
    cov = nd.gaussian_filter(valid.astype(float), sigma_pix)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out[(cov < min_coverage) | ~(den > 0)] = np.nan
    return out


def local_plane_gradient(values, weights, sigma_pix: float, min_coverage: float = 0.5):
    """Velocity gradient from a weighted local plane fit in a Gaussian aperture.

    At every pixel ``p`` this minimises
    ``sum_d K(d) w(p+d) [v(p+d) - a - b d_x - c d_y]^2`` with ``K`` a Gaussian of
    width ``sigma_pix`` and returns ``(b, c) = (dv/dx, dv/dy)`` in km/s per
    pixel.  Unlike differentiating a weighted-smoothed map, the plane fit is
    first-order unbiased when the weights vary across the aperture: a bright
    clump on one side of the aperture does not drag the estimate, so intensity
    structure is not converted into spurious velocity gradients.

    All weighted moments are computed with separable Gaussian-derivative
    filters (``d K = sigma^2 K'``, ``d^2 K = sigma^4 K'' + sigma^2 K``), so the
    cost is ~9 Gaussian filters.  Non-finite pixels get zero weight; pixels
    whose aperture coverage is below ``min_coverage`` are returned as NaN.

    Returns
    -------
    grad_x, grad_y : numpy.ndarray
        Gradient along axis 1 (x) and axis 0 (y), km/s per pixel.
    """
    v = np.asarray(values, float)
    w = np.asarray(weights, float) * np.ones_like(v)
    valid = np.isfinite(v) & np.isfinite(w)
    w = np.where(valid, np.clip(w, 0, None), 0.0)
    wv = w * np.where(valid, v, 0.0)
    s, s2 = float(sigma_pix), float(sigma_pix) ** 2
    g = lambda f, o: nd.gaussian_filter(f, s, order=o)
    S0, T0 = g(w, (0, 0)), g(wv, (0, 0))
    Sx, Sy = s2 * g(w, (0, 1)), s2 * g(w, (1, 0))
    Tx, Ty = s2 * g(wv, (0, 1)), s2 * g(wv, (1, 0))
    Sxx = s2 * s2 * g(w, (0, 2)) + s2 * S0
    Syy = s2 * s2 * g(w, (2, 0)) + s2 * S0
    Sxy = s2 * s2 * g(w, (1, 1))
    # solve [[S0,Sx,Sy],[Sx,Sxx,Sxy],[Sy,Sxy,Syy]] (a,b,c) = (T0,Tx,Ty) by Cramer's rule
    det = (S0 * (Sxx * Syy - Sxy ** 2) - Sx * (Sx * Syy - Sxy * Sy) + Sy * (Sx * Sxy - Sxx * Sy))
    db = (S0 * (Tx * Syy - Sxy * Ty) - T0 * (Sx * Syy - Sxy * Sy) + Sy * (Sx * Ty - Tx * Sy))
    dc = (S0 * (Sxx * Ty - Tx * Sxy) - Sx * (Sx * Ty - Tx * Sy) + T0 * (Sx * Sxy - Sxx * Sy))
    cov = nd.gaussian_filter(valid.astype(float), s)
    with np.errstate(invalid="ignore", divide="ignore"):
        bx, cy = db / det, dc / det
    bad = (cov < min_coverage) | ~(np.abs(det) > 1e-12 * np.nanmax(np.abs(det)))
    bx[bad] = np.nan
    cy[bad] = np.nan
    return bx, cy
