"""Command-line entry point: run the analysis and write figures + a summary.

Example
-------
::

    exgal-flows \
        --mom1 data/..._mom1.fits --mom0 data/..._mom0.fits \
        --model data/rad.mod.fits --co-hires data/..._max.fits \
        --config config/m33.yaml --outdir outputs/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="exgal-flows", description=__doc__)
    p.add_argument("--mom1", required=True, help="HI moment-1 (velocity) FITS")
    p.add_argument("--mom0", required=True, help="HI moment-0 (intensity) FITS")
    p.add_argument("--model", required=True, help="Rotation-model velocity FITS")
    p.add_argument("--co-hires", required=True, help="CO peak-intensity FITS (hi-res)")
    p.add_argument("--co-lores", default=None, help="CO peak-intensity FITS (lo-res)")
    p.add_argument("--config", default=None, help="YAML config (defaults to M33 built-ins)")
    p.add_argument("--outdir", default="outputs", help="Directory for figures + summary")
    p.add_argument("--decim", type=int, default=1, help="Pixel decimation for speed")
    args = p.parse_args(argv)

    # imports here so `--help` is fast and import errors are actionable
    import matplotlib
    matplotlib.use("Agg")
    from .pipeline import FlowAnalysis
    from . import plotting

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    co_paths = {"hires": args.co_hires}
    if args.co_lores:
        co_paths["lores"] = args.co_lores

    ana = FlowAnalysis.from_files(
        mom1_path=args.mom1,
        mom0_path=args.mom0,
        model_path=args.model,
        co_paths=co_paths,
        config_path=args.config,
        co_label="hires",
    )

    div = ana.result_divergence_of_co(decim=args.decim)
    shear = ana.result_co_vs_shear(decim=args.decim)
    study = ana.resolution_study(decim=args.decim)

    fig1 = plotting.plot_divergence_of_co_bright(div)
    fig1.savefig(outdir / "divergence_of_co.png", dpi=200, bbox_inches="tight")

    fig2 = plotting.plot_co_fraction_vs_gradient(shear, xlabel="Shear [km/s/kpc]")
    fig2.savefig(outdir / "co_vs_shear.png", dpi=200, bbox_inches="tight")

    fig3 = plotting.plot_resolution_study(study, which="divergence")
    fig3.savefig(outdir / "resolution_study_divergence.png", dpi=200, bbox_inches="tight")

    summary = {
        "divergence_of_co": {
            "labels": div.labels,
            "median_kms_per_kpc": [float(x) for x in div.medians],
            "err_low": [float(x) for x in div.err_low],
            "err_high": [float(x) for x in div.err_high],
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Wrote figures and summary.json to {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
