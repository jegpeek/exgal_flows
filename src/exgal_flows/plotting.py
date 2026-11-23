"""Publication-oriented plots for the three headline results.

Each function takes results from :mod:`exgal_flows.statistics` (or a
:class:`~exgal_flows.pipeline.FlowAnalysis`) and returns a Matplotlib
``Figure`` so the caller controls saving/showing.  No global ``plt`` state
is mutated beyond creating the figure.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .statistics import GradientCurve, GradientCurveCI, DivergenceResult


def plot_co_fraction_vs_gradient(
    curve: GradientCurve,
    xlabel: str = "Shear [km/s/kpc]",
    ax=None,
    kms_per_kpc: bool = True,
    **plot_kw,
):
    """Plot a CO-bright-fraction curve vs. a velocity-gradient axis.

    Parameters
    ----------
    curve : GradientCurve
        From :func:`exgal_flows.statistics.co_fraction_vs_gradient`.
    xlabel : str
        X-axis label.
    ax : matplotlib Axes, optional
        Axis to draw on; created if omitted.
    kms_per_kpc : bool
        If True, multiply the (km/s/pc) bin centres by 1000 for the axis.
    **plot_kw
        Passed to ``ax.plot``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = ax.figure

    x = curve.centers * (1000.0 if kms_per_kpc else 1.0)
    ax.plot(x, curve.fraction, **plot_kw)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fraction of pixels with CO")
    return fig


def plot_divergence_of_co_bright(result: DivergenceResult, ax=None):
    """Plot median divergence vs. CO brightness with bootstrap errors.

    A negative trend that steepens with CO brightness is the
    "CO sits at colliding flows" result.

    Parameters
    ----------
    result : DivergenceResult
        From :func:`exgal_flows.statistics.divergence_of_co_bright`.
    ax : matplotlib Axes, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = ax.figure

    x = np.arange(len(result.labels))
    yerr = np.vstack(
        [result.medians - result.err_low, result.err_high - result.medians]
    )
    ax.errorbar(x, result.medians, yerr=yerr, fmt="o-", capsize=4)
    ax.axhline(0.0, color="k", linestyle="--", alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(result.labels)
    ax.set_xlabel("CO brightness bin (multiples of threshold)")
    ax.set_ylabel("Median LOS divergence [km/s/kpc]")
    return fig


def plot_resolution_study(study: dict, which: str = "shear", axes=None,
                          label: str | None = None, **plot_kw):
    """Small-multiples of a gradient curve across smoothing scales.

    Parameters
    ----------
    study : dict
        From :meth:`exgal_flows.pipeline.FlowAnalysis.resolution_study`.
    which : {"shear", "divergence"}
        Which curve to show in each panel.
    axes : sequence of matplotlib Axes, optional
        One axis per smoothing scale.  Pass the axes from a previous call
        to overlay several studies (e.g. different pixel selections).
    label : str, optional
        Legend label for this study's line (legend drawn on the first axis).
    **plot_kw
        Passed to ``ax.plot`` (default ``color="C1"``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    fwhms = sorted(study)
    n = len(fwhms)
    if axes is None:
        fig, axes = plt.subplots(1, n, figsize=(4 * n, 3.5), sharey=True)
        if n == 1:
            axes = [axes]
    else:
        axes = list(np.ravel(axes))
        fig = axes[0].figure
    plot_kw.setdefault("color", "C1")

    xlabel = "Shear [km/s/kpc]" if which == "shear" else "Divergence [km/s/kpc]"
    for ax, fwhm in zip(axes, fwhms):
        curve = study[fwhm][which]
        ax.plot(curve.centers * 1000.0, curve.fraction, label=label, **plot_kw)
        title = "no smoothing" if study[fwhm]["sigma_pix"] == 0 else f"FWHM = {fwhm:.0f} pc"
        ax.set_title(title)
        ax.set_xlabel(xlabel)
    axes[0].set_ylabel("Fraction of pixels with CO")
    if label is not None:
        axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def plot_radial_divergence_profile(profile, ax=None, fit=None):
    """Radial profile of CO-bright divergence with block-bootstrap errors.

    Draws the point-estimate median divergence per annulus with the spatial
    block-bootstrap band, the zero line, and the turnover radius.  If a
    linear fit is supplied, the fitted line ``div(r) = a + b r`` and its
    block-bootstrap band are overlaid, with the fitted parameters annotated.

    Parameters
    ----------
    profile : exgal_flows.statistics.RadialProfile
        From :meth:`exgal_flows.pipeline.FlowAnalysis.radial_divergence_profile`
        or :func:`exgal_flows.statistics.tiled_block_bootstrap_radial_profile`.
    ax : matplotlib Axes, optional
    fit : LinearRadialFit or list of LinearRadialFit, optional
        From :meth:`exgal_flows.pipeline.FlowAnalysis.radial_divergence_fit`.
        When given, each fitted line and its block-bootstrap band are drawn
        (e.g. pass ``[ols_fit, l1_fit]`` to overlay both estimators).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 4.2))
    else:
        fig = ax.figure

    r = profile.r_centers
    ax.fill_between(r, profile.err_low, profile.err_high,
                    color="#1f6feb", alpha=0.25,
                    label="block bootstrap (68%)")
    ax.plot(r, profile.medians, "o", color="#1f6feb", label="median divergence",
            zorder=5)
    ax.axhline(0.0, color="k", ls="--", alpha=0.6)

    if fit is not None:
        fits = [fit] if not isinstance(fit, (list, tuple)) else list(fit)
        _fit_color = {"ols": "#d62728", "l1": "#2ca02c", "sigmoid": "#9467bd"}
        for f in fits:
            color = _fit_color.get(f.method, "#d62728")
            ax.fill_between(f.r_grid_pc, f.band_lo, f.band_hi,
                            color=color, alpha=0.18)
            lbl = getattr(f, "label", "") or f.method
            ax.plot(f.r_grid_pc, f.line, "-", color=color, lw=2.0,
                    label=lbl, zorder=6)
        for f in fits:
            color = _fit_color.get(f.method, "#d62728")
            if np.isfinite(f.zero_crossing_pc):
                ax.axvline(f.zero_crossing_pc, color=color, ls=":", lw=1.4)
    elif np.isfinite(profile.turnover_pc):
        ax.axvline(profile.turnover_pc, color="#d62728", ls=":", lw=1.8)
        ymin, ymax = ax.get_ylim()
        ax.text(profile.turnover_pc, ymax,
                f"  turnover ≈ {profile.turnover_pc:.0f} pc",
                color="#d62728", va="top", ha="left", fontsize=8)

    ax.set_xlabel("Galactocentric radius [pc]")
    ax.set_ylabel("LOS divergence of CO-bright gas [km/s/kpc]")
    ax.legend(fontsize=8, loc="lower right")
    return fig


def plot_shear_sensitivity(curves: dict, ax=None, normalize: bool = True,
                           colors=None, band: bool = True):
    """Overlay CO-bright-fraction-vs-|shear| curves for a set of conditions.

    Parameters
    ----------
    curves : dict
        Mapping ``label -> GradientCurveCI`` (as returned by
        :meth:`pipeline.FlowAnalysis.shear_sensitivity`).  One line is drawn
        per entry, in dict order.
    ax : matplotlib.axes.Axes, optional
        Target axes; a new figure is made if omitted.
    normalize : bool
        If True (default) each curve is divided by its own low-shear
        (first-bin) value, so the plotted quantity is the *relative*
        suppression -- the sensitivity -- and curves with very different
        absolute CO fractions are comparable.  If False, raw fractions.
    colors : list, optional
        One colour per curve (defaults to the current prop-cycle).
    band : bool
        Draw the block-bootstrap 68% band as a shaded region.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(4.5, 3.4))
    labels = list(curves)
    if colors is None:
        cyc = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#1f4e79"])
        colors = [cyc[i % len(cyc)] for i in range(len(labels))]
    for (label, cur), c in zip(curves.items(), colors):
        x = cur.centers * 1000.0            # km/s/pc -> km/s/kpc
        f = cur.fraction.astype(float)
        lo = cur.err_low.astype(float)
        hi = cur.err_high.astype(float)
        if normalize:
            ref = f[0]
            f, lo, hi = f / ref, lo / ref, hi / ref
        if band:
            ax.fill_between(x, lo, hi, color=c, alpha=0.15, lw=0)
        ax.plot(x, f, "-o", color=c, ms=3, label=str(label))
    if normalize:
        ax.axhline(1.0, color="0.6", lw=0.8, ls=":", zorder=0)
        ax.set_ylabel("CO-bright fraction (normalized to low-shear)")
    else:
        ax.set_ylabel("CO-bright fraction")
    ax.set_xlabel("|shear|  [km/s/kpc]")
    ax.legend(frameon=False)
    return ax


def plot_block_bootstrap_tiles(tiling: dict, seed: int = 0, axes=None, vlim: float = 100.0):
    """Illustrate the spatial tiled block bootstrap.

    Two panels: (left) the CO-bright pixels coloured by line-of-sight
    divergence with the block grid drawn over them; (right) one bootstrap
    resample, with each occupied tile shaded by how many times it was drawn
    (0 = dropped, 1, 2, 3+).  The point is visual: whole tiles -- and the
    correlated pixels inside them -- are kept or dropped together, which is
    why the block bootstrap gives wider, honest error bars than the iid one.

    Parameters
    ----------
    tiling : dict
        Output of :meth:`pipeline.FlowAnalysis.block_tiling`.
    seed : int
        Seed for the single illustrative resample drawn on the right panel.
    axes : (Axes, Axes), optional
        Target axes pair; a new 1x2 figure is made if omitted.
    vlim : float
        Divergence colour limit (km/s/kpc).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import numpy as np
    from matplotlib.patches import Rectangle, Patch
    from matplotlib.collections import PatchCollection

    ii = np.asarray(tiling["i"]); jj = np.asarray(tiling["j"])
    dv = np.asarray(tiling["divergence"])
    tile_idx = np.asarray(tiling["tile_idx"]); n_tiles = int(tiling["n_tiles"])
    block_pix = int(tiling["block_pix"]); px_dec = float(tiling["pixel_scale_pc"])
    block_pc = float(tiling["block_pc"])

    if axes is None:
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 5.0), sharex=True, sharey=True)
    else:
        axL, axR = axes
        fig = axL.figure

    # local origin at the lower-left occupied tile; work in kpc
    i0 = (ii.min() // block_pix) * block_pix
    j0 = (jj.min() // block_pix) * block_pix
    x = (jj - j0) * px_dec / 1000.0
    y = (ii - i0) * px_dec / 1000.0
    tile_kpc = block_pix * px_dec / 1000.0

    # tile lower-left corners (kpc), recovered from the pixel members
    tx = np.empty(n_tiles); ty = np.empty(n_tiles)
    for t in range(n_tiles):
        m = tile_idx == t
        tx[t] = (jj[m].min() // block_pix) * block_pix
        ty[t] = (ii[m].min() // block_pix) * block_pix
    tx = (tx - j0) * px_dec / 1000.0
    ty = (ty - i0) * px_dec / 1000.0

    # LEFT: pixels + grid
    scL = axL.scatter(x, y, c=np.clip(dv, -vlim, vlim), cmap="RdBu_r",
                      s=2, vmin=-vlim, vmax=vlim, rasterized=True)
    for xv in np.arange(0, tx.max() + 2 * tile_kpc, tile_kpc):
        axL.axvline(xv, color="0.55", lw=0.4, zorder=1)
    for yv in np.arange(0, ty.max() + 2 * tile_kpc, tile_kpc):
        axL.axhline(yv, color="0.55", lw=0.4, zorder=1)
    axL.set_title(f"CO-bright pixels + {block_pc:.0f} pc tiles")
    cb = fig.colorbar(scL, ax=axL, fraction=0.046, pad=0.02)
    cb.set_label("LOS divergence [km/s/kpc]")

    # RIGHT: one bootstrap draw
    rng = np.random.default_rng(seed)
    counts = np.bincount(rng.integers(0, n_tiles, size=n_tiles), minlength=n_tiles)
    palette = {0: "#f0f0f0", 1: "#9ecae1", 2: "#3182bd"}
    faces = [palette.get(int(c), "#08306b") for c in counts]
    rects = [Rectangle((tx[t], ty[t]), tile_kpc, tile_kpc) for t in range(n_tiles)]
    axR.add_collection(PatchCollection(rects, facecolor=faces, edgecolor="white",
                                       lw=0.3, zorder=1))
    axR.scatter(x, y, c="0.25", s=0.5, alpha=0.35, rasterized=True, zorder=2)
    axR.set_title("One bootstrap draw: whole tiles resampled")
    axR.legend(handles=[Patch(fc="#f0f0f0", ec="0.6", label="0 (dropped)"),
                        Patch(fc="#9ecae1", ec="0.6", label="1"),
                        Patch(fc="#3182bd", ec="0.6", label="2"),
                        Patch(fc="#08306b", ec="0.6", label="3+")],
               title="times drawn", frameon=False, loc="upper right")

    for ax in (axL, axR):
        ax.set_aspect("equal"); ax.set_xlabel("x [kpc]"); ax.margins(0.02)
    axL.set_ylabel("y [kpc]")
    n_drop = int(np.sum(counts == 0)); n_dup = int(np.sum(counts >= 2))
    fig.suptitle(f"Spatial tiled block bootstrap  |  {n_tiles} occupied tiles "
                 f"({block_pix}x{block_pix} px = {block_pc:.0f} pc);  "
                 f"this draw: {n_drop} dropped, {n_dup} duplicated")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig
