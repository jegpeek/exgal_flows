"""ExGal-Flows: CO abundance vs. colliding flows and shear in galaxy disks.

A small, tested package refactored from an exploratory M33 notebook.  It
tests whether cold molecular gas (CO) preferentially forms where the
neutral-gas velocity field is *converging* (colliding flows) and is
suppressed where the flow is *shearing*, on 50-100 pc scales.

Typical use::

    from exgal_flows import FlowAnalysis, load_config

    ana = FlowAnalysis.from_files(
        mom1_path="data/m33_A+B+C+D+GBT_hi21cm_strict_mom1.fits",
        mom0_path="data/m33_A+B+C+D+GBT_hi21cm_strict_mom0.fits",
        model_path="data/rad.mod.fits",
        co_paths={"hires": "data/m33_7m_co21_..._max.fits"},
        config_path="config/m33.yaml",
    )
    div_result = ana.result_divergence_of_co()   # colliding-flow result
    shear_curve = ana.result_co_vs_shear()        # shear-suppresses-CO
    study = ana.resolution_study()                # the 50-100 pc claim
"""

from .config import (
    Config,
    GalaxyConfig,
    AnalysisConfig,
    load_config,
)
from .pipeline import FlowAnalysis
from . import io, geometry, flows, statistics, plotting

__all__ = [
    "Config",
    "GalaxyConfig",
    "AnalysisConfig",
    "load_config",
    "FlowAnalysis",
    "io",
    "geometry",
    "flows",
    "statistics",
    "plotting",
]

__version__ = "0.1.0"
