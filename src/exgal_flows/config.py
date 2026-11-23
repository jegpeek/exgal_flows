"""Configuration for the ExGal-Flows analysis.

All the physical constants and analysis parameters that were previously
scattered (and sometimes inconsistent) across notebook cells live here in
one place.  Two dataclasses:

* :class:`GalaxyConfig` -- geometry / distance of the target galaxy.
* :class:`AnalysisConfig` -- knobs of the flow-vs-CO analysis.

Both can be built from a YAML file (see ``config/m33.yaml``) via
:func:`load_config`, so a collaborator can reproduce or retune a run
without editing code.

Notes on values unified from the original notebook
--------------------------------------------------
* Inclination appeared as both ``55.08`` deg (deprojection) and ``55.1``
  deg (the ``1/cos(i)`` minor-axis correction).  These differ by 0.02 deg
  (a <0.03% effect on ``cos i``); we adopt a single ``inclination`` and use
  it in both places.
* Distance appeared as ``840`` kpc (pixel-scale calc) and ``800`` kpc
  (back-of-envelope cells).  We adopt ``840`` kpc, the value used for the
  headline results.
* The pixel scale was hardcoded as ``12 pc`` in some cells and computed
  from the FITS ``CDELT2`` in others.  We always compute it from the
  header (:func:`exgal_flows.io.pixel_scale_pc`); ``12 pc`` was simply an
  approximation of that computed value for this dataset.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord


@dataclass
class GalaxyConfig:
    """Geometry and distance of the target galaxy.

    Parameters
    ----------
    name : str
        Galaxy name (informational).
    inclination_deg : float
        Disk inclination in degrees (0 = face-on).
    position_angle_deg : float
        Position angle of the receding major axis, degrees East of North.
    center_ra_deg, center_dec_deg : float
        Kinematic/photometric centre, degrees (FK5).
    distance_kpc : float
        Distance to the galaxy in kpc.
    """

    name: str = "M33"
    inclination_deg: float = 55.08
    position_angle_deg: float = 201.12
    center_ra_deg: float = 23.4607
    center_dec_deg: float = 30.6583
    distance_kpc: float = 840.0

    # --- convenience accessors as astropy quantities ---
    @property
    def inclination(self) -> Angle:
        return Angle(self.inclination_deg * u.deg)

    @property
    def position_angle(self) -> Angle:
        return Angle(self.position_angle_deg * u.deg)

    @property
    def distance(self) -> u.Quantity:
        return self.distance_kpc * u.kpc

    @property
    def center(self) -> SkyCoord:
        return SkyCoord(
            self.center_ra_deg,
            self.center_dec_deg,
            unit=(u.deg, u.deg),
            frame="fk5",
        )

    @property
    def cos_inc(self) -> float:
        import numpy as np

        return float(np.cos(self.inclination.to(u.rad).value))


@dataclass
class AnalysisConfig:
    """Parameters of the flow-vs-CO correlation analysis.

    Parameters
    ----------
    co_threshold_K : float
        CO peak-intensity threshold (``cont_max`` in the notebook) used to
        define a "CO-bright" pixel, in Kelvin.
    hi_col_to_surfdens : float
        Multiplicative factor converting HI moment-0 (K km/s) to HI surface
        density in solar masses per pc**2.  0.0196 is the standard optically
        thin 21-cm conversion including helium is *not* applied here.
    radius_max_pc : float
        Outer radius (in parsecs) of the analysis region for the headline
        correlation plots.  Converted to deprojected pixels at runtime with
        the map's pixel scale, so the *physical* region is identical across
        datasets of different resolution.  (The notebook used ``rmax = 300``
        deprojected pixels, which at the new data's 3.05 pc/pix is ~916 pc.)
    smoothing_fwhm_pc : tuple[float, ...]
        Gaussian smoothing FWHMs *in parsecs* for the resolution study.
        Converted to a per-pixel sigma at runtime.  Expressing them in pc
        (rather than pixels) is what makes the "50-100 pc" claim comparable
        between the old and new maps.  The defaults are the physical
        equivalents of the notebook's ``sigma = [0, 2, 4, 8, 12]`` pixels on
        the new 3.05 pc/pix grid.
    grad_bin_width : float
        Bin width in km/s/pc for the shear / divergence axes (``ds``).
    shear_range : tuple[float, float]
        (min, max) absolute shear in km/s/pc for binning.
    div_range : tuple[float, float]
        (min, max) divergence in km/s/pc for binning.
    n_wedges : int
        Number of azimuthal wedges for the block bootstrap.
    n_boot : int
        Number of bootstrap resamples.
    co_factor_bins : tuple[float, ...]
        Multiples of ``co_threshold_K`` defining CO-brightness bins for the
        stacked divergence histogram.
    random_seed : int
        Seed for reproducible bootstrap resampling.
    radial_edges_pc : tuple[float, ...]
        Galactocentric annulus edges (pc) for the tiled block-bootstrap
        radial divergence profile.
    block_pc : float
        Block side length (pc) for the spatial (tiled) block bootstrap.
    smoothing_min_coverage : float
        For the smoothed velocity fields (resolution study, smoothed shear
        sensitivity): blank pixels where less than this fraction of the
        Gaussian kernel lands on valid pixels.  Blanks get zero weight in
        the smoothing, so this only removes strongly one-sided edge values.
    """

    co_threshold_K: float = 0.6
    hi_col_to_surfdens: float = 0.0196
    radius_max_pc: float = 4000.0
    smoothing_fwhm_pc: tuple[float, ...] = (0.0, 14.4, 28.8, 57.5, 86.3)
    grad_bin_width: float = 0.01
    shear_range: tuple[float, float] = (0.0, 0.2)
    div_range: tuple[float, float] = (-0.2, 0.2)
    n_wedges: int = 8
    n_boot: int = 1000
    co_factor_bins: tuple[float, ...] = (2.0, 3.0, 4.0, 7.0)
    random_seed: int = 0
    # Spatial (tiled) block bootstrap of the radial divergence profile.
    radial_edges_pc: tuple[float, ...] = (
        0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 2000.0, 2500.0, 3200.0, 4000.0
    )
    block_pc: float = 300.0
    smoothing_min_coverage: float = 0.5


@dataclass
class Config:
    """Top-level container bundling galaxy + analysis config."""

    galaxy: GalaxyConfig = field(default_factory=GalaxyConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)

    def to_dict(self) -> dict[str, Any]:
        return {"galaxy": asdict(self.galaxy), "analysis": asdict(self.analysis)}


def load_config(path: str | Path | None = None) -> Config:
    """Load a :class:`Config` from a YAML file.

    If ``path`` is None, returns the built-in defaults (M33).
    """
    if path is None:
        return Config()

    import yaml

    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}

    gal_raw = raw.get("galaxy", {})
    ana_raw = raw.get("analysis", {})

    # tuples come back from YAML as lists; coerce the fields that want tuples
    for key in ("smoothing_fwhm_pc", "co_factor_bins", "radial_edges_pc"):
        if key in ana_raw and ana_raw[key] is not None:
            ana_raw[key] = tuple(ana_raw[key])
    for key in ("shear_range", "div_range"):
        if key in ana_raw and ana_raw[key] is not None:
            ana_raw[key] = tuple(ana_raw[key])

    return Config(
        galaxy=GalaxyConfig(**gal_raw),
        analysis=AnalysisConfig(**ana_raw),
    )
