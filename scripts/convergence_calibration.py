"""Measure, and calibrate, the convergence of HI flows onto CO-bright gas.

Data (new-epoch M33, R = 1.2-4 kpc): the in-plane gradient d u_x / d x =
divergence / sin i is compared between CO-bright pixels and HI-matched pixels
(same wedge x ring x HI decile).  The excess is fit around the disk as
A + B cos 2theta + C sin 2theta (theta from the minor axis), giving
div u = 2A, radial term A + B, azimuthal term A - B.

Simulations: the same estimator run on mock observations with an in-plane
*vector* velocity field (``nulls.simulate_vector_observation``), where the
true convergence at CO-bright gas is known.  The recovered / true ratio at
each scale calibrates the measurement.

Usage (from the repo root; needs the aperture cache from
``scripts/cache_m33_aperture.py``)::

    python scripts/convergence_calibration.py

Writes ``outputs/convergence_calibration_sims.csv`` and
``outputs/convergence_amplitude_data.csv``.
"""
import time

import numpy as np
import pandas as pd

if "use" not in globals():                      # allow exec() from a session that already ran the setup
    exec(open("scripts/null_test_setup.py").read())
from exgal_flows import nulls as nl

SIN_I = np.sin(np.deg2rad(55.08))
RING = (1200, 2000, 3000, 4000)
CASES = [("native", 0.0)] + [("plane", f) for f in SCALES]
OUT = "outputs"


def measure(m0_map, m1_map, co_b, R_map, az_map, mask):
    """Azimuthal fit of the in-plane excess d u_x / d x at each scale (km/s/kpc)."""
    R_, az_, y_, m0_ = (rav(a)[mask] for a in (R_map, az_map, co_b, m0_map))
    rows = []
    for est, fw in CASES:
        gD, _ = gradients(m0_map, m1_map, est, fw)
        g = rav(gD)[mask] * 1e3 / SIN_I
        E, N = nl.azimuthal_excess(g, y_, m0_, R_, az_, RING)
        rows.append(dict(fwhm=fw, **nl.azimuthal_fit(E, N)))
    return rows


def truth(o, mask, fwhm_list):
    """True excess of div u, G_RR, G_phiphi (km/s/kpc) at CO-bright pixels, HI-matched."""
    R_, az_, y_, m0_ = (rav(a)[mask] for a in (o["radius_pc"], o["theta"], o["co_bright"], o["mom0"]))
    out = {}
    for fw in fwhm_list:
        T = nl.polar_tensor(nl.gradient_tensor(o["ux"], o["uy"], px, fw), o["theta"])
        for k in ("div", "RR", "pp"):
            E, N = nl.azimuthal_excess(rav(T[k])[mask] * 1e3, y_, m0_, R_, az_, RING)
            out[(fw, k)] = E.sum() / N.sum()
    return out


def run_data():
    sel = (R >= 1200) & (R < 4000) & use
    az = Z["azimuth"]
    rows = measure(mom0, vel, (co > thr), Z["radius_pix"] * px, az, sel)
    return pd.DataFrame(rows)


SIM_CASES = [
    # name, seeds, kwargs for simulate_vector_observation, noise model
    ("null", range(1001, 1009), dict(), "worst"),
    ("null + m=1,2 residual patterns", range(1101, 1105),
     dict(los_pattern=((1, 5.0, 0.3), (2, 5.0, 1.1))), "worst"),
    ("div -> CO, 57 pc, strong", range(1201, 1207), dict(inject="div", inject_fwhm_pc=57.5, inject_strength=6.0), "worst"),
    ("div -> CO, 57 pc, weak", range(1301, 1307), dict(inject="div", inject_fwhm_pc=57.5, inject_strength=12.0), "worst"),
    ("div -> CO, 20 pc (cloud scale)", range(1401, 1407), dict(inject="div", inject_fwhm_pc=20.0, inject_strength=6.0), "worst"),
    ("azimuthal -> CO, 57 pc", range(1501, 1505), dict(inject="azimuthal", inject_fwhm_pc=57.5, inject_strength=6.0), "worst"),
    ("radial -> CO, 57 pc", range(1601, 1605), dict(inject="radial", inject_fwhm_pc=57.5, inject_strength=6.0), "worst"),
    ("div -> CO, 57 pc, radiometer noise", range(1701, 1705), dict(inject="div", inject_fwhm_pc=57.5, inject_strength=6.0), "radiometer"),
    ("div -> CO, 57 pc, compressive 0.25", range(1801, 1804),
     dict(inject="div", inject_fwhm_pc=57.5, inject_strength=6.0, compressive_fraction=0.25), "worst"),
    ("div -> CO, 57 pc, compressive 0.75", range(1901, 1904),
     dict(inject="div", inject_fwhm_pc=57.5, inject_strength=6.0, compressive_fraction=0.75), "worst"),
]


def run_sims(cases=SIM_CASES, log=print):
    rows, t0 = [], time.time()
    for name, seeds, kw, noise in cases:
        P = nl.m33_calibrated_params(noise)
        for sd in seeds:
            o = nl.simulate_vector_observation(P, np.random.default_rng(sd), sin_i=SIN_I, **kw)
            mask = np.isfinite(rav(o["mom1"])) & (rav(o["radius_pc"]) >= 1200) & (rav(o["radius_pc"]) < 4000)
            meas = measure(o["mom0"], o["mom1"], o["co_bright"], o["radius_pc"], o["theta"], mask)
            inj_fw = kw.get("inject_fwhm_pc", 57.5)
            tr = truth(o, mask, sorted({0.0, inj_fw, *SCALES}))
            for r in meas:
                fw = r["fwhm"]
                rows.append(dict(case=name, seed=sd, **r,
                                 true_div=tr[(fw, "div")], true_RR=tr[(fw, "RR")], true_pp=tr[(fw, "pp")],
                                 true_div_inj=tr[(inj_fw, "div")], true_div_grid=tr[(0.0, "div")]))
        log(f"{name}: {len(list(seeds))} sims [{time.time() - t0:.0f}s]")
    return pd.DataFrame(rows)


if __name__ == "__main__" or globals().get("RUN_CALIBRATION", False):
    D = run_data()
    D.to_csv(f"{OUT}/convergence_amplitude_data.csv", index=False, float_format="%.5g")
    S = run_sims()
    S.to_csv(f"{OUT}/convergence_calibration_sims.csv", index=False, float_format="%.5g")
