"""Input / output: loading FITS maps and putting them on a common grid.

The analysis needs, on a single common pixel grid:

* the HI moment-1 (velocity) field,
* an axisymmetric rotation *model* velocity field to subtract,
* the HI moment-0 (integrated intensity) field,
* one or more CO peak-intensity maps,
* optionally GALEX FUV and a star-formation-history map.

This module loads them as :class:`spectral_cube.Projection` objects and
reprojects the CO/FUV maps onto the HI grid.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.io import fits
from spectral_cube import Projection


def load_projection(path: str | Path, to_unit: u.Unit | None = None) -> Projection:
    """Load a 2-D FITS image as a :class:`~spectral_cube.Projection`.

    Parameters
    ----------
    path : str or Path
        Path to a 2-D FITS file.
    to_unit : astropy Unit, optional
        If given, convert the projection to this unit.
    """
    proj = Projection.from_hdu(fits.open(str(path)))
    if to_unit is not None:
        proj = proj.to(to_unit)
    return proj


def pixel_scale_pc(header, distance_kpc: float) -> float:
    """Physical pixel size in parsecs from a FITS header + distance.

    Uses ``CDELT2`` (degrees/pixel) and the small-angle approximation,
    replacing the notebook's hardcoded ``12 pc``.

    Parameters
    ----------
    header : astropy.io.fits.Header
        Header carrying a ``CDELT2`` keyword.
    distance_kpc : float
        Distance to the galaxy in kpc.
    """
    pixel_in_rad = abs(header["CDELT2"]) * np.pi / 180.0
    distance_pc = distance_kpc * 1000.0
    return pixel_in_rad * distance_pc


@dataclass
class SkyMaps:
    """HI, CO and model maps all reprojected onto the HI moment-1 grid.

    Attributes
    ----------
    mom1 : Projection
        HI moment-1 (velocity) field, km/s.
    mom0 : Projection
        HI moment-0 (integrated intensity) field, K km/s.
    model1 : Projection
        Axisymmetric rotation model velocity field, km/s, on the HI grid.
    co : dict[str, Projection]
        CO peak-intensity maps (already reprojected), keyed by a short label.
    pixel_scale_pc : float
        Physical pixel size in parsecs (from the HI header).
    header : astropy.io.fits.Header
        The HI moment-1 header defining the common grid.
    extras : dict[str, Projection]
        Any optional maps (e.g. ``"fuv"``), already reprojected.
    """

    mom1: Projection
    mom0: Projection
    model1: Projection
    co: dict[str, Projection]
    pixel_scale_pc: float
    header: fits.Header
    extras: dict[str, Projection]


def load_maps(
    mom1_path: str | Path,
    mom0_path: str | Path,
    model_path: str | Path,
    co_paths: dict[str, str | Path],
    distance_kpc: float,
    extra_paths: dict[str, str | Path] | None = None,
) -> SkyMaps:
    """Load all maps and reproject CO / model / extras onto the HI grid.

    Parameters
    ----------
    mom1_path, mom0_path : str or Path
        HI moment-1 (velocity, km/s) and moment-0 (K km/s) FITS files.
    model_path : str or Path
        Rotation-model velocity field FITS file (km/s).
    co_paths : dict[str, str or Path]
        CO peak-intensity maps keyed by a short label (e.g. ``"hires"``).
    distance_kpc : float
        Galaxy distance, kpc, for the pixel-scale calculation.
    extra_paths : dict[str, str or Path], optional
        Any additional maps to load and reproject (e.g. ``{"fuv": ...}``).

    Returns
    -------
    SkyMaps
    """
    mom1_hdu = fits.open(str(mom1_path))[0]
    mom1 = Projection.from_hdu(mom1_hdu).to(u.km / u.s)
    header = mom1_hdu.header

    mom0 = Projection.from_hdu(fits.open(str(mom0_path))[0])

    model1 = Projection.from_hdu(fits.open(str(model_path))).to(u.km / u.s)
    model1 = model1.reproject(header)

    co = {}
    for label, path in co_paths.items():
        proj = Projection.from_hdu(fits.open(str(path)))
        co[label] = proj.reproject(header)

    extras = {}
    if extra_paths:
        for label, path in extra_paths.items():
            proj = Projection.from_hdu(fits.open(str(path)))
            extras[label] = proj.reproject(header)

    px = pixel_scale_pc(header, distance_kpc)

    return SkyMaps(
        mom1=mom1,
        mom0=mom0,
        model1=model1,
        co=co,
        pixel_scale_pc=px,
        header=header,
        extras=extras,
    )
