"""Run the colliding-flow analysis on both HI epochs and compare.

Produces a side-by-side of the "median divergence of CO-bright pixels"
result for the new (A+B+C+D+GBT, ~3 pc/pix) and old (14B+17B 8", ~4 pc/pix)
HI data, so you can see how much the higher-resolution A-array data changes
the headline number.

Usage
-----
    python scripts/compare_epochs.py --data data --outdir outputs --decim 2
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path


def _isolate_astropy_dirs():
    # keep astropy's config/cache out of a possibly read-only HOME
    tmp = tempfile.mkdtemp(prefix="astropy_")
    os.environ.setdefault("XDG_CONFIG_HOME", tmp)
    os.environ.setdefault("XDG_CACHE_HOME", tmp)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="data", help="Directory with the FITS files")
    p.add_argument("--outdir", default="outputs", help="Where to write the figure")
    p.add_argument("--decim", type=int, default=2, help="Pixel decimation for speed")
    args = p.parse_args(argv)

    _isolate_astropy_dirs()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from exgal_flows import FlowAnalysis
    from exgal_flows import plotting

    data = Path(args.data)
    co_hires = data / "m33_7m_co21_2p5kms_pbcorr_round.image.max.fits"
    model = data / "rad.mod.fits"

    epochs = {
        "new (A+B+C+D+GBT, 5.3\")": {
            "mom1": data / "m33_A+B+C+D+GBT_hi21cm_strict_mom1.fits",
            "mom0": data / "m33_A+B+C+D+GBT_hi21cm_strict_mom0.fits",
            "config": "config/m33.yaml",
        },
        "old (14B+17B, 8.2\")": {
            "mom1": data / "M33_14B_17B_HI_contsub_width_1kms.image.pbcor.GBT_feathered_K.mom1.fits",
            "mom0": data / "M33_14B_17B_HI_contsub_width_1kms.image.pbcor.GBT_feathered_K.mom0.fits",
            "config": "config/m33_old.yaml",
        },
    }

    results = {}
    for label, paths in epochs.items():
        print(f"[compare_epochs] running {label} ...")
        ana = FlowAnalysis.from_files(
            mom1_path=str(paths["mom1"]),
            mom0_path=str(paths["mom0"]),
            model_path=str(model),
            co_paths={"hires": str(co_hires)},
            config_path=paths["config"],
            co_label="hires",
        )
        results[label] = ana.result_divergence_of_co(decim=args.decim)

    fig, axes = plt.subplots(1, len(results), figsize=(6 * len(results), 4), sharey=True)
    if len(results) == 1:
        axes = [axes]
    for ax, (label, res) in zip(axes, results.items()):
        plotting.plot_divergence_of_co_bright(res, ax=ax)
        ax.set_title(label)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "compare_epochs_divergence.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"[compare_epochs] wrote {out}")

    # also print the brightest-bin numbers for a quick text comparison
    for label, res in results.items():
        print(f"  {label}: brightest-CO median divergence = "
              f"{res.medians[-1]:+.2f} [{res.err_low[-1]:+.2f}, {res.err_high[-1]:+.2f}] km/s/kpc")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
