"""Selection cuts (``sel``) and displaced CO (``co_mode``) for the shear results.

Stubbed fields as in ``test_pipeline_null.py``: ``FlowAnalysis`` is a plain
dataclass, so the attributes the shear methods read are filled in directly.
The velocity field is built so that |shear| is low in the left half of the
grid and high in the right half, with bright CO only in a left-half blob.
"""

from types import SimpleNamespace

import numpy as np
import pytest
import astropy.units as u

from exgal_flows import FlowAnalysis, Config
from exgal_flows.config import AnalysisConfig
from exgal_flows.flows import velocity_gradients
from exgal_flows.geometry import deproject, polar_grid

N = 201
PX = 10.0                 # pc / pixel
LO_SHEAR, HI_SHEAR = 0.025, 0.155   # km/s/pc, mid-bin for 0.01 bins
BLOB = (slice(60, 140), slice(40, 80))


@pytest.fixture
def analysis():
    cfg = Config(analysis=AnalysisConfig(n_boot=20))
    ana = FlowAnalysis(config=cfg, maps=SimpleNamespace(pixel_scale_pc=PX, extras={}))

    yy, xx = np.mgrid[:N, :N].astype(float)
    g = np.where(xx < N // 2, LO_SHEAR, HI_SHEAR)        # target |shear| per column
    # v = y * g * cos(i) * px  ->  shear = dv/dy / px / cos(i) = g
    v = (yy - N / 2) * g * cfg.galaxy.cos_inc * PX
    ana.velocity_deproj = v * u.km / u.s
    ana.flow = velocity_gradients(ana.velocity_deproj, PX, cfg.galaxy)

    r2 = (yy - N / 2) ** 2 + (xx - N / 2) ** 2
    ana.mom0_deproj = 1.0 + 10.0 * np.exp(-r2 / (2 * 15.0**2))   # HI bump at centre

    co = np.full((N, N), 0.1)
    co[BLOB] = 5.0
    ana.co_deproj = co
    ana.radius, ana.azimuth = polar_grid((N, N))
    return ana


def _bin_of(curve, value):
    return int(np.argmin(np.abs(curve.centers - value)))


def test_shear_field_is_as_built(analysis):
    sh = np.abs(analysis.flow.shear.value)
    assert np.allclose(sh[:, : N // 2], LO_SHEAR)
    assert np.allclose(sh[:, N // 2 + 1:], HI_SHEAR)


def test_all_true_sel_reproduces_default(analysis):
    sel = np.ones((N, N), bool)
    a, b = analysis.result_co_vs_shear(), analysis.result_co_vs_shear(sel=sel)
    np.testing.assert_array_equal(a.fraction, b.fraction)

    a, b = analysis.shear_sensitivity(n_boot=20), analysis.shear_sensitivity(n_boot=20, sel=sel)
    np.testing.assert_array_equal(a.fraction, b.fraction)
    np.testing.assert_array_equal(a.err_low, b.err_low)
    assert a.rel_drop == b.rel_drop

    a, b = analysis.resolution_study(decim=2), analysis.resolution_study(decim=2, sel=sel)
    for fw in a:
        np.testing.assert_array_equal(a[fw]["shear"].fraction, b[fw]["shear"].fraction)


def test_sel_restricts_pixels(analysis):
    base = analysis.result_co_vs_shear()
    sel = np.zeros((N, N), bool)
    sel[BLOB] = True
    cut = analysis.result_co_vs_shear(sel=sel)
    i_lo, i_hi = _bin_of(base, LO_SHEAR), _bin_of(base, HI_SHEAR)
    assert 0 < base.fraction[i_lo] < 1
    assert cut.fraction[i_lo] == 1.0                 # only the CO blob is left
    assert np.isnan(cut.fraction[i_hi])              # no high-shear pixels selected

    sens = analysis.shear_sensitivity(n_boot=20, sel=sel, decim=2)
    assert sens.lo_frac == 1.0 and np.isnan(sens.hi_frac)


def test_co_mode_reverses_trend_and_restores_co(analysis):
    co_before = analysis.co_deproj
    real = analysis.shear_sensitivity(n_boot=20)
    rot = analysis.shear_sensitivity(n_boot=20, co_mode="rot180")
    assert real.lo_frac > 0 and real.hi_frac == 0    # CO only at low shear
    assert rot.lo_frac == 0 and rot.hi_frac > 0      # rotated blob lands at high shear

    c = analysis.result_co_vs_shear(co_mode="mirror_minor")   # x -> x: blob stays left
    np.testing.assert_array_equal(c.fraction, analysis.result_co_vs_shear().fraction)
    c = analysis.result_co_vs_shear(co_mode="mirror_major")   # x -> -x: blob moves right
    assert c.fraction[_bin_of(c, LO_SHEAR)] == 0

    study = analysis.resolution_study(decim=2, co_mode="rot180")
    assert study[0.0]["shear"].fraction[_bin_of(c, LO_SHEAR)] == 0
    assert analysis.co_deproj is co_before


def test_bad_inputs(analysis):
    with pytest.raises(ValueError):
        analysis.shear_sensitivity(co_mode="shuffle")
    with pytest.raises(ValueError):
        analysis.result_co_vs_shear(sel=np.ones((N - 1, N), bool))
    with pytest.raises(ValueError):
        analysis.selection_mask(laplacian_sign=2)
    with pytest.raises(ValueError):
        analysis.selection_mask(dispersion_range_kms=(0, 10))
    with pytest.raises(KeyError):
        analysis.selection_mask(dispersion_label="hi_mom2", dispersion_range_kms=(0, 10))


def test_laplacian_sign_mask(analysis):
    neg = analysis.selection_mask(laplacian_sign=-1)
    pos = analysis.selection_mask(laplacian_sign=+1)
    c = N // 2
    assert neg[c, c] and not pos[c, c]               # HI peak: del^2 < 0
    assert pos[c, c + 35] and not neg[c, c + 35]     # flank (r > sqrt(2) sigma): del^2 > 0
    assert not np.any(neg & pos)
    # smoothing keeps the sign structure of a smooth bump
    neg_s = analysis.selection_mask(laplacian_sign=-1, laplacian_fwhm_pc=50.0)
    assert neg_s[c, c] and not neg_s[c, c + 35]


def test_dispersion_mask_converts_units(analysis):
    yy, xx = np.mgrid[:N, :N]
    disp_ms = np.where(xx < N // 2, 5000.0, 15000.0)          # 5 / 15 km/s, in m/s
    analysis.maps.extras["hi_mom2"] = SimpleNamespace(unit=u.m / u.s)
    analysis._extras_deproj[("hi_mom2", False)] = disp_ms      # skip the real deprojection
    sel = analysis.selection_mask(dispersion_label="hi_mom2", dispersion_range_kms=(0, 10))
    assert sel[:, : N // 2].all() and not sel[:, N // 2:].any()
    both = analysis.selection_mask(
        laplacian_sign=-1, dispersion_label="hi_mom2", dispersion_range_kms=(10, np.inf)
    )
    assert both.sum() > 0 and not both[:, : N // 2].any()


def test_extra_deproj_caches_real_projection(synthetic_flat):
    ana = FlowAnalysis(config=Config(), maps=SimpleNamespace(extras={"m": synthetic_flat}))
    a = ana.extra_deproj("m")
    assert a is ana.extra_deproj("m")
    # same grid as every other deprojected map, and no cos i factor by default
    np.testing.assert_array_equal(a, deproject(synthetic_flat, ana.config.galaxy,
                                               inc_correction=False))
    assert not np.allclose(a[np.isfinite(a)],
                           ana.extra_deproj("m", inc_correction=True)[np.isfinite(a)])
