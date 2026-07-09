"""End-to-end pipeline tying the modules together.

:class:`FlowAnalysis` runs the whole chain a collaborator cares about:

1. load HI + CO + model maps onto a common grid (:mod:`exgal_flows.io`),
2. subtract the rotation model, deproject to the galaxy plane
   (:mod:`exgal_flows.geometry`),
3. form divergence / shear fields (:mod:`exgal_flows.flows`),
4. hold the raveled arrays ready for the statistics
   (:mod:`exgal_flows.statistics`).

The class exposes the intermediate maps as attributes so nothing is
hidden, and thin ``result_*`` methods return the headline numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import astropy.units as u

from .config import Config, load_config
from . import io as _io
from .geometry import deproject, polar_grid
from . import flows as _flows
from . import statistics as _stats


@dataclass
class FlowAnalysis:
    """Run and hold the full flow-vs-CO analysis for one galaxy.

    Construct with :meth:`from_files` (or :meth:`from_config`), then access
    the deprojected maps as attributes and call the ``result_*`` methods.
    """

    config: Config
    maps: _io.SkyMaps
    co_label: str = "hires"

    # populated by run()
    velocity_deproj: np.ndarray = field(default=None, repr=False)
    mom0_deproj: np.ndarray = field(default=None, repr=False)
    co_deproj: np.ndarray = field(default=None, repr=False)
    flow: _flows.FlowFields = field(default=None, repr=False)
    radius: np.ndarray = field(default=None, repr=False)
    azimuth: np.ndarray = field(default=None, repr=False)

    # ---- constructors ---------------------------------------------------
    @classmethod
    def from_config(
        cls,
        config: Config,
        mom1_path,
        mom0_path,
        model_path,
        co_paths: dict,
        co_label: str = "hires",
        extra_paths: dict | None = None,
    ) -> "FlowAnalysis":
        maps = _io.load_maps(
            mom1_path=mom1_path,
            mom0_path=mom0_path,
            model_path=model_path,
            co_paths=co_paths,
            distance_kpc=config.galaxy.distance_kpc,
            extra_paths=extra_paths,
        )
        obj = cls(config=config, maps=maps, co_label=co_label)
        obj.run()
        return obj

    @classmethod
    def from_files(
        cls,
        mom1_path,
        mom0_path,
        model_path,
        co_paths: dict,
        config_path=None,
        co_label: str = "hires",
        extra_paths: dict | None = None,
    ) -> "FlowAnalysis":
        """Convenience constructor loading config from YAML (or defaults)."""
        config = load_config(config_path)
        return cls.from_config(
            config, mom1_path, mom0_path, model_path, co_paths, co_label, extra_paths
        )

    # ---- the pipeline ---------------------------------------------------
    def run(self) -> "FlowAnalysis":
        """Execute steps 2-4 (maps are already loaded)."""
        gal = self.config.galaxy
        maps = self.maps

        # model subtraction + deprojection of the residual velocity field
        residual = maps.mom1 - maps.model1
        self.velocity_deproj = deproject(residual, gal) * u.km / u.s

        # deproject moment-0 and the chosen CO map
        self.mom0_deproj = deproject(maps.mom0, gal)
        self.co_deproj = deproject(maps.co[self.co_label], gal)

        # velocity-gradient decomposition
        self.flow = _flows.velocity_gradients(
            self.velocity_deproj, maps.pixel_scale_pc, gal
        )

        # galaxy-frame polar grid on the divergence footprint
        self.radius, self.azimuth = polar_grid(self.flow.divergence.shape)
        return self

    # ---- physical <-> pixel conversions ---------------------------------
    @property
    def pixel_scale_pc(self) -> float:
        """Physical pixel size (pc); convenience alias for ``maps.pixel_scale_pc``."""
        return self.maps.pixel_scale_pc

    @property
    def _radius_max_pix(self) -> float:
        """Analysis radius in deprojected pixels for this map's scale."""
        return self.config.analysis.radius_max_pc / self.maps.pixel_scale_pc

    def _sigma_pix(self, fwhm_pc: float) -> float:
        """Gaussian sigma in pixels for a smoothing FWHM in parsecs."""
        if fwhm_pc == 0:
            return 0.0
        return fwhm_pc / self.maps.pixel_scale_pc / 2.355

    # ---- raveled views for statistics -----------------------------------
    def _ravel(self, decim: int = 1):
        sl = (slice(None, None, decim), slice(None, None, decim))
        div = self.flow.divergence.value[sl].ravel()
        shear = np.abs(self.flow.shear.value[sl].ravel())
        co = np.asarray(self.co_deproj)[sl].ravel()
        rad = self.radius[sl].ravel()
        az = self.azimuth[sl].ravel()
        return div, shear, co, rad, az

    # ---- shear sensitivity ----------------------------------------------
    def shear_sensitivity(
        self,
        decim: int = 1,
        co_threshold: float | None = None,
        fwhm_pc: float = 0.0,
        radius_window_pc: tuple[float, float] | None = None,
        n_boot: int = 300,
        seed: int = 0,
    ):
        """CO-bright fraction as a function of |shear|, with a block bootstrap.

        This is the shear analog of the divergence results: it measures how
        strongly the fraction of CO-bright pixels is suppressed as the local
        line-of-sight velocity shear increases, with an honest azimuthal
        wedge block-bootstrap error band (the per-pixel iid error is invalid
        because the maps are spatially correlated).

        The same call supports every breakdown of the shear study by varying
        one argument:

        * **smoothing** -- set ``fwhm_pc`` to re-derive the shear from a
          Gaussian-smoothed velocity field (as in the resolution study);
        * **CO-intensity cut** -- set ``co_threshold`` to change the peak-Tb
          level that defines "CO-bright" (defaults to the config value);
        * **radial trend** -- set ``radius_window_pc=(r_lo, r_hi)`` to
          restrict to a galactocentric annulus (image-plane radius);
        * **epoch** -- build two :class:`FlowAnalysis` objects (new / old
          config) and call this on each.

        Parameters
        ----------
        decim : int
            Pixel decimation for speed.
        co_threshold : float, optional
            Peak-Tb threshold (K) defining a CO-bright pixel.  Defaults to
            ``config.analysis.co_threshold_K``.
        fwhm_pc : float
            Gaussian smoothing FWHM (pc) applied to the velocity field
            before the gradient decomposition.  ``0`` uses the native field.
        radius_window_pc : (float, float), optional
            Inclusive galactocentric radius window (pc).  ``None`` uses the
            full analysis aperture ``(0, radius_max_pc)``.
        n_boot : int
            Number of wedge bootstrap resamples.
        seed : int
            RNG seed.

        Returns
        -------
        statistics.GradientCurveCI
        """
        an = self.config.analysis
        gal = self.config.galaxy
        px = self.maps.pixel_scale_pc
        thr = an.co_threshold_K if co_threshold is None else co_threshold

        # shear field (optionally from a smoothed velocity field)
        if fwhm_pc and fwhm_pc > 0:
            vsm = _flows.smoothed_velocity_field(
                self.mom0_deproj, self.velocity_deproj, self._sigma_pix(fwhm_pc)
            )
            flow = _flows.velocity_gradients(vsm, px, gal)
        else:
            flow = self.flow

        sl = (slice(None, None, decim), slice(None, None, decim))
        shear = np.abs(flow.shear.value[sl].ravel())
        co = np.asarray(self.co_deproj)[sl].ravel()
        az = self.azimuth[sl].ravel()
        rad_pc = self.radius[sl].ravel() * px

        # aperture + optional annulus mask
        r_lo, r_hi = (0.0, an.radius_max_pc) if radius_window_pc is None else radius_window_pc
        mask = (
            np.isfinite(shear) & np.isfinite(co) & np.isfinite(az)
            & (rad_pc >= r_lo) & (rad_pc < r_hi)
        )

        edges = np.arange(an.shear_range[0], an.shear_range[1], an.grad_bin_width)
        return _stats.co_fraction_vs_gradient_boot(
            shear, co, az, thr, edges,
            n_wedges=an.n_wedges, n_boot=n_boot, seed=seed, mask=mask,
        )

    # ---- headline results -----------------------------------------------
    def result_divergence_of_co(self, decim: int = 1) -> _stats.DivergenceResult:
        """Median divergence of CO-bright pixels (colliding-flow result)."""
        div, shear, co, rad, az = self._ravel(decim)
        ana = self.config.analysis
        rmask = rad < self._radius_max_pix
        return _stats.divergence_of_co_bright(
            divergence=div[rmask],
            co=co[rmask],
            azimuth=az[rmask],
            co_threshold=ana.co_threshold_K,
            factors=ana.co_factor_bins,
            n_wedges=ana.n_wedges,
            n_boot=ana.n_boot,
            seed=ana.random_seed,
        )

    def result_co_vs_shear(self, decim: int = 1) -> _stats.GradientCurve:
        """CO-bright fraction vs. |shear| (shear-suppresses-CO result)."""
        div, shear, co, rad, az = self._ravel(decim)
        ana = self.config.analysis
        edges = np.arange(ana.shear_range[0], ana.shear_range[1], ana.grad_bin_width)
        return _stats.co_fraction_vs_gradient(
            gradient=shear,
            co=co,
            co_threshold=ana.co_threshold_K,
            bin_edges=edges,
            mask=rad < self._radius_max_pix,
        )

    def result_co_vs_divergence(self, decim: int = 1) -> _stats.GradientCurve:
        """CO-bright fraction vs. divergence."""
        div, shear, co, rad, az = self._ravel(decim)
        ana = self.config.analysis
        edges = np.arange(ana.div_range[0], ana.div_range[1], ana.grad_bin_width)
        return _stats.co_fraction_vs_gradient(
            gradient=div,
            co=co,
            co_threshold=ana.co_threshold_K,
            bin_edges=edges,
            mask=rad < self._radius_max_pix,
        )

    def radial_divergence_profile(
        self, decim: int = 1, co_factor: float = 2.0
    ) -> _stats.RadialProfile:
        """Radial profile of CO-bright divergence with a spatial block bootstrap.

        Selects CO-bright pixels (CO > ``co_factor`` x threshold) inside the
        analysis aperture and computes the median line-of-sight divergence in
        galactocentric annuli, with the tiled block bootstrap that respects
        the maps' spatial correlation.  This is the honest-error version of
        the "does the CO-divergence sign flip with radius?" test; the
        point-estimate zero-crossing is returned as ``turnover_pc``.

        Parameters
        ----------
        decim : int
            Decimation stride applied to the maps before selection.
        co_factor : float
            CO-brightness cut as a multiple of ``co_threshold_K``.

        Returns
        -------
        exgal_flows.statistics.RadialProfile
        """
        ana = self.config.analysis
        px = self.maps.pixel_scale_pc
        sl = (slice(None, None, decim), slice(None, None, decim))

        div = self.flow.divergence.value[sl]
        co = np.asarray(self.co_deproj)[sl]
        rad = self.radius[sl]

        sel = (
            (co > ana.co_threshold_K * co_factor)
            & np.isfinite(div)
            & (rad < self._radius_max_pix)
        )
        ii, jj = np.nonzero(sel)
        return _stats.tiled_block_bootstrap_radial_profile(
            i_pix=ii,
            j_pix=jj,
            divergence=div[sel] * 1000.0,  # km/s/kpc
            radius_pc=rad[sel] * px,
            pixel_scale_pc=px * decim,
            r_edges=np.asarray(ana.radial_edges_pc, float),
            block_pc=ana.block_pc,
            n_boot=ana.n_boot,
            seed=ana.random_seed,
        )

    def radial_divergence_fit(
        self, decim: int = 1, co_factor: float = 2.0, n_boot: int | None = None,
        method: str = "ols",
    ) -> _stats.LinearRadialFit:
        """Fit the radial divergence trend with a straight line.

        Fits ``div(r) = a + b r`` directly to the CO-bright pixels (the same
        selection as :meth:`radial_divergence_profile`): ``a`` is the overall
        sense of divergence extrapolated to the centre, ``b`` the rate of
        change with radius.  Errors are a spatial block bootstrap over tiles,
        so the fit -- not just the binned error bars -- carries an honest
        uncertainty band.

        Parameters
        ----------
        decim : int
            Decimation stride applied to the maps before selection.
        co_factor : float
            CO-brightness cut as a multiple of ``co_threshold_K``.
        n_boot : int, optional
            Number of bootstrap resamples (defaults to the config value).
        method : str
            ``"ols"`` (least squares) or ``"l1"`` (robust least absolute
            deviations, the median-like fit).

        Returns
        -------
        exgal_flows.statistics.LinearRadialFit
        """
        t = self.block_tiling(decim=decim, co_factor=co_factor)
        ana = self.config.analysis
        return _stats.block_bootstrap_line_fit(
            radius_pc=t["radius_pc"],
            divergence=t["divergence"],
            tile_idx=t["tile_idx"],
            n_tiles=t["n_tiles"],
            n_boot=ana.n_boot if n_boot is None else n_boot,
            seed=ana.random_seed,
            method=method,
        )

    def radial_divergence_sigmoid(
        self, decim: int = 1, co_factor: float = 2.0, n_boot: int | None = None
    ) -> _stats.SigmoidRadialFit:
        """Fit the radial divergence trend with a saturating sigmoid.

        A bounded alternative to :meth:`radial_divergence_fit`: the divergence
        saturates to an inner plateau and an outer plateau across a transition
        radius, instead of a line that grows without limit.  Same CO-bright
        selection and spatial block bootstrap.  The headline number is
        ``straddle_zero_prob`` -- the probability the inner plateau is positive
        and the outer negative, i.e. that the reversal rules out zero.

        Parameters
        ----------
        decim : int
            Decimation stride applied to the maps before selection.
        co_factor : float
            CO-brightness cut as a multiple of ``co_threshold_K``.
        n_boot : int, optional
            Number of bootstrap resamples (defaults to the config value).

        Returns
        -------
        exgal_flows.statistics.SigmoidRadialFit
        """
        t = self.block_tiling(decim=decim, co_factor=co_factor)
        ana = self.config.analysis
        return _stats.block_bootstrap_sigmoid_fit(
            radius_pc=t["radius_pc"],
            divergence=t["divergence"],
            tile_idx=t["tile_idx"],
            n_tiles=t["n_tiles"],
            n_boot=ana.n_boot if n_boot is None else n_boot,
            seed=ana.random_seed,
        )

    def block_tiling(self, decim: int = 1, co_factor: float = 2.0) -> dict:
        """Return the CO-bright pixels and their block-bootstrap tile IDs.

        Exposes the exact spatial tiling that
        :meth:`radial_divergence_profile` resamples, so a figure can show
        *what the blocks look like*.  The selection (CO-bright pixels inside
        the aperture) and the tile keying are identical to that method.

        Parameters
        ----------
        decim : int
            Decimation stride (use the same value you pass to
            :meth:`radial_divergence_profile` for a faithful picture).
        co_factor : float
            CO-brightness cut as a multiple of ``co_threshold_K``.

        Returns
        -------
        dict
            Keys: ``i`` / ``j`` (pixel row/col of each CO-bright pixel, in
            decimated-pixel units), ``divergence`` (km/s/kpc), ``radius_pc``,
            ``tile_idx`` (tile index per pixel), ``n_tiles``, ``block_pix``
            (tile side, decimated pixels), ``block_pc``, and
            ``pixel_scale_pc`` (physical pc per decimated pixel).
        """
        ana = self.config.analysis
        px = self.maps.pixel_scale_pc
        px_dec = px * decim
        sl = (slice(None, None, decim), slice(None, None, decim))

        div = self.flow.divergence.value[sl]
        co = np.asarray(self.co_deproj)[sl]
        rad = self.radius[sl]
        sel = (
            (co > ana.co_threshold_K * co_factor)
            & np.isfinite(div)
            & (rad < self._radius_max_pix)
        )
        ii, jj = np.nonzero(sel)
        block_pix = max(1, int(round(ana.block_pc / px_dec)))
        tile_idx, n_tiles, block_pix = _stats.assign_tiles(ii, jj, block_pix)
        return {
            "i": ii,
            "j": jj,
            "divergence": div[sel] * 1000.0,
            "radius_pc": rad[sel] * px,
            "tile_idx": tile_idx,
            "n_tiles": n_tiles,
            "block_pix": block_pix,
            "block_pc": ana.block_pc,
            "pixel_scale_pc": px_dec,
        }

    def resolution_study(self, decim: int = 1):
        """Repeat the shear/divergence curves at each smoothing scale.

        Returns
        -------
        dict[float, dict]
            Keyed by FWHM in parsecs; each value has ``"shear"`` and
            ``"divergence"`` :class:`~exgal_flows.statistics.GradientCurve`
            plus the raw sigma in pixels.
        """
        gal = self.config.galaxy
        ana = self.config.analysis
        px = self.maps.pixel_scale_pc

        out = {}
        for fwhm in ana.smoothing_fwhm_pc:
            sigma = self._sigma_pix(fwhm)
            vsm = _flows.smoothed_velocity_field(
                self.mom0_deproj, self.velocity_deproj, sigma
            )
            flow = _flows.velocity_gradients(vsm, px, gal)
            sl = (slice(None, None, decim), slice(None, None, decim))
            shear = np.abs(flow.shear.value[sl].ravel())
            div = flow.divergence.value[sl].ravel()
            co = np.asarray(self.co_deproj)[sl].ravel()
            rad = self.radius[sl].ravel()
            rmask = rad < self._radius_max_pix

            sh_edges = np.arange(*ana.shear_range, ana.grad_bin_width)
            dv_edges = np.arange(*ana.div_range, ana.grad_bin_width)
            out[fwhm] = {
                "sigma_pix": sigma,
                "shear": _stats.co_fraction_vs_gradient(
                    shear, co, ana.co_threshold_K, sh_edges, rmask
                ),
                "divergence": _stats.co_fraction_vs_gradient(
                    div, co, ana.co_threshold_K, dv_edges, rmask
                ),
            }
        return out
