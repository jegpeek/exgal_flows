"""Tests for the config layer: defaults, YAML round-trip, unit accessors."""

import numpy as np
import astropy.units as u

from exgal_flows import load_config, Config, GalaxyConfig


def test_defaults_are_m33():
    cfg = Config()
    assert cfg.galaxy.name == "M33"
    assert cfg.galaxy.distance_kpc == 840.0
    # single inclination used everywhere (no more 55.08 vs 55.1 split)
    assert cfg.galaxy.inclination_deg == 55.08


def test_quantity_accessors():
    g = GalaxyConfig()
    assert g.distance.unit == u.kpc
    assert np.isclose(g.cos_inc, np.cos(np.radians(g.inclination_deg)))
    assert g.center.ra.deg == g.center_ra_deg


def test_yaml_roundtrip(tmp_path):
    import yaml

    cfg = Config()
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(cfg.to_dict()))
    loaded = load_config(p)
    assert loaded.galaxy.distance_kpc == cfg.galaxy.distance_kpc
    # tuple fields survive the list<->tuple coercion
    assert isinstance(loaded.analysis.smoothing_fwhm_pc, tuple)
    assert isinstance(loaded.analysis.co_factor_bins, tuple)


def test_load_none_returns_defaults():
    assert load_config(None).galaxy.name == "M33"
