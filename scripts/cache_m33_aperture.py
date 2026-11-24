"""Cache the deprojected M33 (new epoch) maps inside the analysis aperture.

Writes a compact float32 .npz (aperture bounding box + margin) so null
simulations do not re-read the ~5 GB HI products.
"""
import sys, numpy as np
from exgal_flows.pipeline import FlowAnalysis

D = "data"
out = sys.argv[1] if len(sys.argv) > 1 else "/tmp/m33_aperture_cache.npz"
margin_pix = 256
ana = FlowAnalysis.from_files(
    mom1_path=f"{D}/m33_A+B+C+D+GBT_hi21cm_strict_mom1.fits",
    mom0_path=f"{D}/m33_A+B+C+D+GBT_hi21cm_strict_mom0.fits",
    model_path=f"{D}/rad.mod.fits",
    co_paths={"hires": f"{D}/m33_7m_co21_2p5kms_pbcorr_round.image.max.fits"},
    config_path="config/m33.yaml",
)
px = ana.maps.pixel_scale_pc
rmax_pix = ana.config.analysis.radius_max_pc / px
rr = ana.radius
ii, jj = np.nonzero(rr < rmax_pix)
i0, i1 = max(ii.min() - margin_pix, 0), min(ii.max() + margin_pix + 1, rr.shape[0])
j0, j1 = max(jj.min() - margin_pix, 0), min(jj.max() + margin_pix + 1, rr.shape[1])
sl = (slice(i0, i1), slice(j0, j1))
g = lambda a: np.asarray(getattr(a, "value", a))[sl].astype(np.float32)
np.savez(out,
         mom0=g(ana.mom0_deproj), vel=g(ana.velocity_deproj), co=g(ana.co_deproj),
         shear=g(ana.flow.shear), div=g(ana.flow.divergence),
         radius_pix=g(ana.radius), azimuth=g(ana.azimuth),
         pixel_scale_pc=px, rmax_pix=rmax_pix, bbox=np.array([i0, i1, j0, j1]),
         co_threshold_K=ana.config.analysis.co_threshold_K,
         inclination_deg=ana.config.galaxy.inclination_deg)
print("cached", out, "shape", (i1 - i0, j1 - j0), "px", px)
