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


def lowpass_velocity_field(velocity, sigma_pix: float, block: int | None = None) -> np.ndarray:
    """NaN-aware Gaussian low-pass of a velocity field.

    A normalised convolution, ``smooth(v * w) / smooth(w)`` with ``w`` the
    finite-pixel mask, so blanked pixels neither contribute nor drag the
    average toward zero.  Because the output is smooth on the scale
    ``sigma_pix`` it is computed on a block-averaged grid and interpolated
    back, which makes kiloparsec-scale kernels affordable on the full map.

    Parameters
    ----------
    velocity : array-like or Quantity
        Deprojected velocity field (km/s); NaN where blanked.
    sigma_pix : float
        Gaussian sigma in pixels.
    block : int, optional
        Block-averaging factor.  Defaults to ``sigma_pix / 8`` (at least 1),
        which keeps the kernel well sampled on the coarse grid.

    Returns
    -------
    numpy.ndarray
        Low-passed field, same shape as the input, defined wherever the
        kernel reaches finite data (NaN elsewhere).
    """
    v = np.asarray(getattr(velocity, "value", velocity), dtype=float)
    if block is None:
        block = max(1, int(sigma_pix // 8))
    ny, nx = v.shape
    py, px = (-ny) % block, (-nx) % block
    finite = np.isfinite(v)
    vz = np.pad(np.where(finite, v, 0.0), ((0, py), (0, px)))
    wz = np.pad(finite.astype(float), ((0, py), (0, px)))
    sy, sx = vz.shape[0] // block, vz.shape[1] // block
    vsum = vz.reshape(sy, block, sx, block).sum(axis=(1, 3))
    wsum = wz.reshape(sy, block, sx, block).sum(axis=(1, 3))
    del vz, wz

    s = sigma_pix / block
    num = nd.gaussian_filter(vsum, s, mode="constant")
    den = nd.gaussian_filter(wsum, s, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        low = np.where(den > 1e-3 * block * block, num / den, np.nan)

    if block > 1:
        # fill NaN before the linear upsample so they do not spread, then
        # interpolate pixel-area aligned (grid_mode) back to the full grid
        bad = ~np.isfinite(low)
        low_f = np.where(bad, 0.0, low)
        up = nd.zoom(low_f, block, order=1, grid_mode=True, mode="nearest")
        upbad = nd.zoom(bad.astype(float), block, order=0, grid_mode=True, mode="nearest")
        up[upbad > 0.5] = np.nan
        low = up
    return low[:ny, :nx]


def highpass_velocity_field(velocity, sigma_pix: float, block: int | None = None) -> np.ndarray:
    """Velocity field with everything smoother than ``sigma_pix`` removed.

    ``v - lowpass(v)``: a model-free replacement for subtracting a fitted
    rotation model.  It removes rotation, warps / position-angle twists,
    axisymmetric radial flows and centre errors together, at the cost of
    also removing any real flow coherent on scales larger than the kernel.
    It is not valid within roughly one kernel width of the galaxy centre,
    where the rotation field itself varies on the kernel scale.

    Parameters
    ----------
    velocity : array-like or Quantity
        Deprojected velocity field (km/s); NaN where blanked.
    sigma_pix : float
        Gaussian sigma in pixels of the low-pass that is subtracted.
    block : int, optional
        Passed to :func:`lowpass_velocity_field`.

    Returns
    -------
    numpy.ndarray
        High-pass filtered velocity (km/s, plain ndarray).
    """
    v = np.asarray(getattr(velocity, "value", velocity), dtype=float)
    return v - lowpass_velocity_field(v, sigma_pix, block=block)
