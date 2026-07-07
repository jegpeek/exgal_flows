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


def plot_resolution_study(study: dict, which: str = "shear"):
    """Small-multiples of a gradient curve across smoothing scales.

    Parameters
    ----------
    study : dict
        From :meth:`exgal_flows.pipeline.FlowAnalysis.resolution_study`.
    which : {"shear", "divergence"}
        Which curve to show in each panel.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fwhms = sorted(study)
    n = len(fwhms)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 3.5), sharey=True)
    if n == 1:
        axes = [axes]

    xlabel = "Shear [km/s/kpc]" if which == "shear" else "Divergence [km/s/kpc]"
    for ax, fwhm in zip(axes, fwhms):
        curve = study[fwhm][which]
        ax.plot(curve.centers * 1000.0, curve.fraction, color="C1")
        title = "no smoothing" if study[fwhm]["sigma_pix"] == 0 else f"FWHM = {fwhm:.0f} pc"
        ax.set_title(title)
        ax.set_xlabel(xlabel)
    axes[0].set_ylabel("Fraction of pixels with CO")
    fig.tight_layout()
    return fig


def plot_radial_divergence_profile(profile, ax=None):
    """Radial profile of CO-bright divergence: iid vs block-bootstrap errors.

    Draws the point-estimate median divergence per annulus with the wide
    spatial block-bootstrap band and, for contrast, the narrow (over-
    optimistic) iid band.  The zero line and the turnover radius are marked.

    Parameters
    ----------
    profile : exgal_flows.statistics.RadialProfile
        From :meth:`exgal_flows.pipeline.FlowAnalysis.radial_divergence_profile`
        or :func:`exgal_flows.statistics.tiled_block_bootstrap_radial_profile`.
    ax : matplotlib Axes, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 4.2))
    else:
        fig = ax.figure

    r = profile.r_centers
    ax.fill_between(r, profile.iid_low, profile.iid_high,
                    color="0.6", alpha=0.5, label="iid bootstrap (68%)")
    ax.fill_between(r, profile.err_low, profile.err_high,
                    color="#1f6feb", alpha=0.30, label="block bootstrap (68%)")
    ax.plot(r, profile.medians, "o-", color="#1f6feb", label="median divergence")
    ax.axhline(0.0, color="k", ls="--", alpha=0.6)

    if np.isfinite(profile.turnover_pc):
        ax.axvline(profile.turnover_pc, color="#d62728", ls=":", lw=1.8)
        ymin, ymax = ax.get_ylim()
        ax.text(profile.turnover_pc, ymax, f"  turnover ≈ {profile.turnover_pc:.0f} pc",
                color="#d62728", va="top", ha="left", fontsize=8)

    ax.set_xlabel("Galactocentric radius [pc]")
    ax.set_ylabel("Median LOS divergence of CO-bright gas [km/s/kpc]")
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
