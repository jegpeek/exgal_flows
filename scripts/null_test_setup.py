"""Shared state for the beam-smearing null tests (exec'd or imported).

Builds, for the new-epoch M33 data inside the 4 kpc aperture (every 2nd pixel):
``use`` mask, CO-bright flags ``y``, HI column ``m0``, spatial ``tile_idx`` /
``n_tiles`` (300 pc), the HI-only CO null probability ``p_null`` (spatially
cross-fitted), and ``gradients(mom0, mom1, estimator, fwhm_pc)``.

Requires the aperture cache written by ``scripts/cache_m33_aperture.py``.
"""
import numpy as np
from exgal_flows import nulls as nl, statistics as st, flows as fl

CACHE = "/tmp/m33_aperture_cache.npz"
Z = np.load(CACHE)
mom0, vel, co, shear, div = (Z[k] for k in ["mom0", "vel", "co", "shear", "div"])
px, thr = float(Z["pixel_scale_pc"]), float(Z["co_threshold_K"])
cos_i = float(np.cos(np.deg2rad(float(Z["inclination_deg"]))))
DECIM = 2
sl = (slice(None, None, DECIM),) * 2
rav = lambda a: a[sl].ravel()
R = rav(Z["radius_pix"]) * px
F = nl.hi_structure_features(mom0)
X = np.column_stack([rav(F[k]) for k in F] + [R]).astype(np.float32)
use = (np.isfinite(rav(shear)) & np.isfinite(rav(co)) & np.isfinite(rav(Z["azimuth"])) & (R < 4000)
       & np.all(np.isfinite(X), axis=1))
y = rav(co)[use] > thr
m0 = rav(mom0)[use]
ii, jj = np.indices(mom0.shape)
tile_idx, n_tiles, _ = st.assign_tiles(rav(ii)[use] // DECIM, rav(jj)[use] // DECIM, 300.0 / (px * DECIM))
p_null = nl.hi_only_co_null(X[use], y, tile_idx, n_tiles, seed=0)
strata = nl.quantile_strata(m0, 10)
SCALES = [14.4, 28.8, 57.5, 86.3]


def gradients(m0_map, m1_map, estimator="plane", fwhm_pc=0.0):
    """Signed (div, shear) in km/s/pc. estimator: 'native' | 'current' | 'plane'."""
    sg = fwhm_pc / px / 2.355
    if estimator == "native" or fwhm_pc == 0:
        gy, gx = np.gradient(m1_map)
    elif estimator == "current":
        gy, gx = np.gradient(fl.smoothed_velocity_field(m0_map, m1_map, sg))
    else:
        gx, gy = fl.local_plane_gradient(m1_map, np.clip(np.nan_to_num(m0_map), 0, None) ** 2, sg)
    return gx / px / cos_i, gy / px / cos_i
