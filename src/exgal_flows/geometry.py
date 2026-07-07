"""Deprojection to the galaxy plane and galaxy-frame coordinate grids.

The deprojection recipe is the one from the original notebook: shift the
image so the galaxy centre is at the array centre, rotate so the major
axis is vertical, then stretch the minor axis by ``1/cos(i)`` to undo the
inclination foreshortening.  It is a pixel-space (``scipy.ndimage``)
approximation -- fast and adequate for the statistical flow analysis here,
but not a full WCS-aware reprojection.

Behaviour is kept identical to the original notebook so that results
reproduce exactly.  In particular ``inc_correction`` (multiply the
deprojected values by ``cos i``) defaults to **True**, as it did in the
notebook, where *every* ``deproject`` call used the default.

There is a subtlety worth flagging (it does not change any number, only
your interpretation of them):

* For a **velocity** field, this ``cos i`` factor cancels exactly against
  the ``1 / cos i`` that :mod:`exgal_flows.flows` applies when forming the
  divergence/shear -- both are scalar multiplies, so the net effect on the
  velocity gradients is nil.  The physically meaningful operation on the
  velocity field is the ``1/cos i`` x-axis *stretch* (the ``nd.zoom``),
  which is always applied.
* For an **intensity** map (HI moment-0, CO peak intensity), the ``cos i``
  is a genuine surface-brightness (foreshortening) correction.  A
  consequence: the CO peak-intensity threshold ``cont_max`` is compared
  against CO values that have already been multiplied by ``cos i``, so the
  effective threshold carries that factor.  This is exactly what the
  notebook did; we keep it, and expose the flag so you can experiment.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as nd
from astropy.coordinates import SkyCoord

from .config import GalaxyConfig


def deproject(
    image,
    galaxy: GalaxyConfig,
    inc_correction: bool = True,
):
    """Deproject a 2-D map into the galaxy plane.

    Parameters
    ----------
    image : spectral_cube.Projection
        The map to deproject.  Must expose ``.shape``,
        ``.spatial_coordinate_map`` and behave like an ndarray.
    galaxy : GalaxyConfig
        Galaxy geometry (inclination, PA, centre).
    inc_correction : bool, optional
        If True (default, matching the notebook), multiply the deprojected
        values by ``cos(i)``.  See the module docstring for what this means
        for velocity vs. intensity maps.

    Returns
    -------
    numpy.ndarray
        The deprojected image (plain ndarray, NaN outside the valid
        footprint).
    """
    inc = galaxy.inclination
    pa = galaxy.position_angle
    gal_centre = galaxy.center

    image_copy = np.asarray(image).copy()
    mask = np.isfinite(np.asarray(image))
    image_copy[~mask] = 0.0

    pix_centre = (image.shape[0] // 2, image.shape[1] // 2)

    dec, ra = image.spatial_coordinate_map
    skycoord_map = SkyCoord(ra, dec)
    sep_dist = gal_centre.separation(skycoord_map)

    gal_centre_pix = np.unravel_index(np.argmin(sep_dist), sep_dist.shape)

    shift_pix = tuple(
        img_mid - gal_mid
        for img_mid, gal_mid in zip(pix_centre, gal_centre_pix)
    )

    shift_img = nd.shift(image_copy, shift_pix)
    shift_mask = nd.shift(mask.astype(float), shift_pix)

    # Rotate so the PA is along the y axis.
    rot = nd.rotate(shift_img, pa.to("deg").value - 180)
    rot_mask = nd.rotate(shift_mask, pa.to("deg").value - 180)

    # Stretch the x axis to undo inclination foreshortening.
    cos_i = galaxy.cos_inc
    deproj = nd.zoom(rot, (1.0, 1.0 / cos_i))
    deproj_mask = nd.zoom(rot_mask, (1.0, 1.0 / cos_i))

    if inc_correction:
        deproj = deproj * cos_i

    deproj[deproj_mask < 0.5] = np.nan
    return deproj


def polar_grid(shape: tuple[int, int]):
    """Radius and azimuth arrays for a deprojected image.

    The centre is taken at ``shape / 2`` (consistent with the deprojection,
    which places the galaxy centre at the array centre).

    Parameters
    ----------
    shape : tuple[int, int]
        (ny, nx) of the deprojected image.

    Returns
    -------
    radius : numpy.ndarray
        Radius from the centre in pixels.
    azimuth : numpy.ndarray
        Azimuthal angle in radians, in ``(-pi, pi]`` via ``arctan2``.
    """
    ny, nx = shape
    xx, yy = np.mgrid[:ny, :nx]
    azimuth = np.arctan2(yy - nx / 2, xx - ny / 2)
    radius = np.sqrt((yy - nx / 2) ** 2 + (xx - ny / 2) ** 2)
    return radius, azimuth
