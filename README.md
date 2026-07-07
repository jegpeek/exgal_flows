# ExGal-Flows

**Does cold molecular gas form where interstellar flows collide, and is it
suppressed where they shear?**

This package tests, in the resolved disk of M33, three linked claims that
came out of an exploratory notebook:

1. **CO sits at colliding flows.** Bright CO peaks preferentially at
   *negative* line-of-sight velocity divergence (converging / colliding
   gas), not at diverging flows.
2. **Shear suppresses CO.** The fraction of pixels that are CO-bright falls
   as the local velocity shear increases.
3. **This happens on 50-100 pc scales.** Both trends persist across a range
   of smoothing scales, and are strongest at ~50-100 pc.

It is a cleaned-up, tested, importable refactor of the original
`M33 HI Velocity fields and CO.ipynb`.

## Install

```bash
pip install -e .            # from the repo root
pip install -e '.[test]'    # with pytest
```

Requires Python >= 3.10 and (numpy, scipy, matplotlib, astropy,
spectral-cube, reproject, pyyaml).

## Quick start

```python
from exgal_flows import FlowAnalysis

ana = FlowAnalysis.from_files(
    mom1_path="data/m33_A+B+C+D+GBT_hi21cm_strict_mom1.fits",
    mom0_path="data/m33_A+B+C+D+GBT_hi21cm_strict_mom0.fits",
    model_path="data/rad.mod.fits",
    co_paths={"hires": "data/m33_7m_co21_2p5kms_pbcorr_round.image.max.fits"},
    config_path="config/m33.yaml",
)

div   = ana.result_divergence_of_co()   # result 1: colliding flows
shear = ana.result_co_vs_shear()        # result 2: shear suppresses CO
study = ana.resolution_study()          # result 3: the 50-100 pc scan
prof  = ana.radial_divergence_profile() # radial trend w/ spatial block bootstrap
```

Or from the command line:

```bash
exgal-flows \
  --mom1 data/..._mom1.fits --mom0 data/..._mom0.fits \
  --model data/rad.mod.fits --co-hires data/..._max.fits \
  --config config/m33.yaml --outdir outputs/
```

## What the pipeline does

| Step | Module | What happens |
|------|--------|--------------|
| Load & regrid | `io` | Load HI mom0/mom1, rotation model, CO peak maps; reproject everything onto the HI grid. |
| Deproject | `geometry` | Subtract the rotation model; shift/rotate/stretch into the galaxy plane. |
| Decompose flow | `flows` | Velocity gradients -> **divergence** (colliding vs diverging) and **shear**. |
| Statistics | `statistics` | CO-bright fraction vs gradient; median divergence of CO-bright pixels with a **wedge block-bootstrap** error; **tiled spatial block bootstrap** of the radial profile (honest errors under map correlation). |
| Figures | `plotting` | The three headline plots. |
| Orchestration | `pipeline` | `FlowAnalysis` ties it together and exposes intermediate maps. |
| Config | `config` | All physical constants and knobs, from YAML. |

See [`docs/architecture.md`](docs/architecture.md) for a one-page module map
and data-flow diagram.

## Configuration

Everything the notebook hardcoded lives in `config/m33.yaml`. Analysis
scales are expressed **physically (in parsecs)**, not pixels, so a run on
the new 3 pc/pixel data and the old 4 pc/pixel data probes the *same*
physical scales -- essential for the new-vs-old comparison.

Provide a different YAML (or edit this one) to retune the CO threshold,
analysis radius, smoothing scales, or bootstrap settings.

## Datasets

Place FITS files in `data/` (git-ignored). Two HI epochs are supported:

* **New** -- VLA A+B+C+D+GBT, 5.3" beam (~3.05 pc/pix): the default.
* **Old** -- VLA 14B+17B 8" feathered, 8.2" beam (~4.07 pc/pix):
  `config/m33_old.yaml` points the HI paths here for comparison.

## Reproducibility notes

The refactor reproduces the notebook's recipe. In doing so it fixes several
inconsistencies that were scattered across cells (none change the intended
science; they make it consistent and reviewable):

* **Inclination** was `55.08` deg in one place and `55.1` in another --
  unified to a single value.
* **Distance** was `840` kpc (pixel scale) and `800` kpc (rough cells) --
  adopts `840` kpc, the value behind the headline results.
* **Pixel scale** was hardcoded `12 pc` in some cells and computed from the
  header in others -- always computed from `CDELT2` + distance now.
* **`rmax`** took four different values across cells -- one config value.
* **Analysis scales** moved from pixels to parsecs (see above).
* **`np.NaN`** (removed in NumPy 2) replaced by `np.nan`.
* The `cos i` inclination correction inside `deproject` is kept **on** (as
  in the notebook) so numbers reproduce; its meaning for velocity vs
  intensity maps is documented in `geometry.py`.

## Tests

```bash
pytest -q
```

The statistics tests inject a *known* signal and check the estimators
recover the right sign and that the bootstrap is reproducible.

## Layout

```
src/exgal_flows/   config, io, geometry, flows, statistics, plotting, pipeline, cli
config/            m33.yaml (new data), m33_old.yaml (old data)
tests/             pytest suite (synthetic data, no FITS needed)
notebooks/         reproduce_results.ipynb -- thin driver over the package
scripts/           compare_epochs.py -- new-vs-old runner
data/              FITS files (git-ignored)
```
