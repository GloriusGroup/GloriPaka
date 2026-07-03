<p align="center">
  <img src="docs/assets/GloriPaka.png" alt="GloriPaka" width="340">
</p>

<p align="center">
  Sequence-based hydrophobic-interaction-chromatography (HIC) retention prediction for VHH nanobodies.
</p>

<p align="center">
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.11%2B-blue.svg" alt="Python 3.11+"></a>
  <a href="https://docs.astral.sh/uv/"><img src="https://img.shields.io/badge/env-uv-6340ac.svg" alt="uv"></a>
  <a href="https://github.com/astral-sh/ruff"><img src="https://img.shields.io/badge/style-ruff%20(google)-blueviolet.svg" alt="Code style: Ruff / Google"></a>
  <a href="docs/datasets.md"><img src="https://img.shields.io/badge/docs-datasets-blue.svg" alt="Docs"></a>
</p>

---

## Contents
- [Overview](#overview)
- [Method](#method)
- [How it works](#how-it-works)
- [Installation](#installation)
- [Usage](#usage)
- [Package layout](#package-layout)
- [Data and reproducibility](#data-and-reproducibility)
- [Scientific background](#scientific-background)
- [Citation](#citation)
- [Acknowledgements](#acknowledgements)
- [License](#license)

## Overview

GloriPaka predicts hydrophobic-interaction-chromatography (HIC) retention, a common
surface-hydrophobicity and developability proxy, for VHH single-domain antibodies
(nanobodies) directly from their amino-acid sequence. The user provides a sequence;
the package produces a predicted 3-D structure, biophysical features, and a
leakage-free HIC prediction.

The public HIC labels for part of the challenge test set are used only as a held-out
validation set. They never enter training or feature selection, so reported metrics
estimate genuine generalization rather than a fit to the evaluation data.

## Method

The predictor is a rank-mean ensemble of two complementary estimators:

1. An unfit physics composite: the equal-weight average of four solvent-exposed
   hydrophobic-surface descriptor ranks, defined a priori with no fitted coefficients.
   It requires no labels, cannot overfit, and is the reference estimator in the
   low-data regime.
2. An interpretable Ridge regression trained on native VHH data, dominated by the same
   exposed-hydrophobic-surface, aromatic-exposure, and charge/pI biophysics.

On a held-out set of 34 VHH with public HIC labels, the ensemble reaches a Spearman
correlation of about 0.68 and a Pearson correlation of about 0.79, compared with about
0.55 Spearman for a strong mAb-trained baseline. With n = 34 the confidence intervals
are wide; the lab journal reports them in full.

## How it works

```
VHH heavy-chain sequence
        |
        v
3-D structure          NanoBodyBuilder2 (4-model conformational ensemble)
        |
        v
biophysical features   FreeSASA surface and patches, exposed aromatic content,
        |               charge/pI, IMGT CDR statistics, DeepSP SAP/SCM, ESM2
        v
two estimators         (a) unfit physics composite  +  (b) native-trained Ridge
        |
        v
rank-mean ensemble ->  calibrated HIC value and per-sequence confidence tier
```

## Installation

GloriPaka uses [uv](https://docs.astral.sh/uv/) for environment and package management.

```bash
git clone https://github.com/GloriusGroup/GloriPaka.git
cd GloriPaka

# core library
uv sync

# optional: featurization and structure-prediction stack
uv sync --extra featurize

# optional: developer tooling and documentation
uv sync --extra dev --extra docs
```

Run commands inside the environment with `uv run`, for example:

```bash
uv run python -c "import gloripaka; print(gloripaka.__version__)"
```

## Usage

The core library is available now. It provides the honest metric suite and the unfit
physics composite:

```python
import pandas as pd
from gloripaka.physics import physics_composite, PHYSICS_COMPOSITE
from gloripaka.metrics import compute_metrics

# `features` holds the four physics descriptors as columns; see the `featurize`
# extra for how they are computed from a structure.
score = physics_composite(features)       # equal-weight rank composite, no training
print(compute_metrics(true_hic, score))   # Pearson, Spearman, Kendall, RMSE, ...
```

The high-level predictor runs the full ensemble, from precomputed features or, with the
`featurize` extra installed, directly from raw sequences:

```python
from gloripaka import GloriPakaPredictor

predictor = GloriPakaPredictor()
predictions = predictor.predict_from_features(features)   # from precomputed features
# predictions = predictor.predict(sequences)              # from sequence (uv sync --extra featurize)
print(predictions[["predicted_hic", "confidence"]])
```

Reproduce the challenge predictions and score the held-out proxy from the command line:

```bash
gloripaka reproduce
```

## Package layout

| Module | Purpose |
|:--|:--|
| `gloripaka.config` | Central paths, feature-group prefixes, and shared constants. |
| `gloripaka.physics` | The unfit physics composite (transferable hydrophobic-surface baseline). |
| `gloripaka.models` | Rank-target regression with in-fold feature selection, open backends, and scale calibration. |
| `gloripaka.features` | Feature-table loading, group slicing, and label-joined dataset assembly. |
| `gloripaka.labels` | Label loading and the known/unknown split of the challenge VHH. |
| `gloripaka.validation` | Clone-family clustering and leave-clone-family-out cross-validation. |
| `gloripaka.metrics` | The full metric suite with bootstrap confidence intervals and residual-outlier flags. |
| `gloripaka.report` | Lightweight experiment reporting to CSV and markdown. |
| `gloripaka.secrets` | Out-of-repo access to the optional TabPFN API token. |
| `gloripaka.viz` | A single central publication figure style. |
| `gloripaka.pipeline` | The high-level `GloriPakaPredictor` (sequence or features to HIC). |

## Data and reproducibility

The featurization stack is installed, not vendored. `uv sync --extra featurize` pulls
the published tools we build on (NanoBodyBuilder2 via ImmuneBuilder, ANARCI/anarcii,
FreeSASA, and DeepSP; ESM2 is a separate `esm` extra), each of which downloads its own
model weights from its original source at first use. This keeps the repository free of
third-party model files while staying fully reproducible.

Predicted structures and feature tables are not shipped because they are deterministic:
NanoBodyBuilder2 is a feed-forward predictor and GloriPaka uses its unrefined output, so
the same sequence regenerates the same structure and the same features on any machine.
Nothing beyond the sequences needs to be hosted or downloaded.

Third-party datasets are cited, not redistributed. The Ginkgo / PROPHET-Ab (GDPa)
measurements in particular are not included here; `data/README.md` and `docs/datasets.md`
record their provenance and how to obtain them from the original source.

## Scientific background

The full experimental record, including negative results, is kept in the research lab
journal (`notes/lab_journal.md` in the source workspace) and summarized in the
documentation here. In brief: training on native VHH data improves on mAb-only training;
transferable surface hydrophobicity is the signal that generalizes; the physics-plus-ML
ensemble improves on either component alone; and the main limitation is the quantity of
data available on a shared assay rather than the modelling method.

## Citation

If you use GloriPaka, please cite it (see [CITATION.cff](CITATION.cff)):

```bibtex
@software{gloripaka2026,
  title  = {GloriPaka: sequence-based HIC-retention prediction for VHH nanobodies},
  author = {Boser, F. M. and Dobberstein, N. D. and Niedzella, M. G. and Spies, J. C.},
  year   = {2026},
  url    = {https://github.com/GloriusGroup/GloriPaka},
  note   = {Version 0.1.0}
}
```

## Acknowledgements

Developed in the Glorius Group (University of Münster). Built on open tools including
NanoBodyBuilder2, FreeSASA, ANARCI/IMGT, DeepSP, ESM2, and scikit-learn, and on public
antibody and HIC datasets acknowledged in the documentation. See [NOTICE](NOTICE) for
third-party asset attributions, including the logo.

## License

Released under the [MIT License](LICENSE).
