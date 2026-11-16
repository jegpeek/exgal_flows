"""Chance-alignment null test on FlowAnalysis, with stubbed fields.

No FITS or deprojection: ``FlowAnalysis`` is a plain dataclass, so the
attributes the statistics methods read are filled in directly.
"""

from types import SimpleNamespace

import numpy as np
import pytest
import astropy.units as u

from exgal_flows import FlowAnalysis, Config
from exgal_flows.config import AnalysisConfig
from exgal_flows.geometry import polar_grid


@pytest.fixture
def analysis():
    n = 201
    rng = np.random.default_rng(0)
    div = rng.normal(0, 0.002, (n, n))           # km/s/pc
    co = np.full((n, n), 0.1)
    div[40:60, 40:60] = -0.02                    # converging blob...
    co[40:60, 40:60] = 5.0                       # ...with bright CO on it
    cfg = Config(analysis=AnalysisConfig(n_boot=50, block_pc=100.0))
    ana = FlowAnalysis(config=cfg, maps=SimpleNamespace(pixel_scale_pc=10.0))
    ana.flow = SimpleNamespace(
        divergence=div * u.km / u.s / u.pc,
        shear=np.zeros_like(div) * u.km / u.s / u.pc,
    )
    ana.co_deproj = co
    ana.radius, ana.azimuth = polar_grid(div.shape)
    return ana


def test_null_falls_to_control_while_real_does_not(analysis):
    co_before = analysis.co_deproj
    res = analysis.null_test_chance_alignment("rot180")
    assert res.real_bins.medians[-1] < -15          # km/s/kpc
    assert abs(res.null_bins.medians[-1]) < 5
    # the blob sits at ~560 pc: annuli 300-600 pc and 600-900 pc
    assert np.nanmin(res.real_excess) < -15
    assert np.nanmax(np.abs(res.null_excess[np.isfinite(res.null_excess)])) < 5
    assert abs(res.control_overall) < 5
    # the real CO map is restored and nothing was copied into it
    assert analysis.co_deproj is co_before


def test_co_restored_after_error(analysis, monkeypatch):
    co_before = analysis.co_deproj

    def boom(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(analysis, "radial_divergence_profile", boom)
    with pytest.raises(RuntimeError):
        analysis.null_test_chance_alignment()
    assert analysis.co_deproj is co_before


def test_null_views_do_not_copy(analysis):
    for fn in FlowAnalysis._NULL_VIEWS.values():
        assert np.shares_memory(fn(analysis.co_deproj), analysis.co_deproj)


def test_unknown_mode(analysis):
    with pytest.raises(ValueError):
        analysis.null_test_chance_alignment("shuffle")
