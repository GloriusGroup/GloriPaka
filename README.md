<!-- Badge/DOI/docs URLs marked (placeholder) switch to live targets in Phase 5/6 (see ROADMAP.md). -->
<p align="center">
  <img src="docs/assets/GloriPaka.png" alt="GloriPaka" width="360">
</p>

<p align="center">
  <em>Sequence in, HIC out — hydrophobic-interaction-chromatography retention prediction for VHH nanobodies.</em>
</p>

<p align="center">
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.11%2B-blue.svg" alt="Python 3.11+"></a>
  <a href="https://docs.astral.sh/uv/"><img src="https://img.shields.io/badge/env-uv-purple.svg" alt="uv"></a>
  <a href="https://github.com/astral-sh/ruff"><img src="https://img.shields.io/badge/code%20style-ruff%20(google)-blueviolet.svg" alt="Code style: Ruff / Google"></a>
  <a href="https://GloriusGroup.github.io/GloriPaka/"><img src="https://img.shields.io/badge/docs-building-lightgrey.svg" alt="Docs"></a>
  <a href="tests/"><img src="https://img.shields.io/badge/tests-CI%20pending-lightgrey.svg" alt="Tests"></a>
  <a href="#citation"><img src="https://img.shields.io/badge/Zenodo-DOI%20pending-lightgrey.svg" alt="Zenodo DOI"></a>
</p>

---

## Contents
- [Overview](#overview)
- [Key features](#key-features)
- [How it works](#how-it-works)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Package layout](#package-layout)
- [Data & reproducibility](#data--reproducibility)
- [Scientific background](#scientific-background)
- [Roadmap](#roadmap)
- [Citation](#citation)
- [Authors & acknowledgements](#authors--acknowledgements)
- [License](#license)

## Overview

**GloriPaka** predicts **hydrophobic-interaction-chromatography (HIC) retention** — a
widely-used surface-hydrophobicity / developability proxy — for **VHH single-domain
antibodies (nanobodies)** directly from their amino-acid sequence. You provide a
sequence; the package handles the rest: predicted 3-D structure, biophysical
featurization, and a leakage-free predictor.

The model is a **rank-mean ensemble of two complementary estimators**:

1. an **unfit "physics" composite** — the equal-weight average of four solvent-exposed
   hydrophobic-surface descriptor ranks, defined a priori with **no fitted coefficients**.
   It needs no labels, cannot overfit, and is the honest workhorse in the low-data regime; and
2. an **interpretable Ridge** trained on native VHH, dominated by the same exposed
   hydrophobic-surface / aromatic-exposure / charge biophysics.

On a **held-out validation proxy of 34 VHH** with public HIC labels (never used in
training or feature selection), the ensemble reaches **Spearman ≈ 0.68 / Pearson ≈ 0.79**,
versus **≈ 0.55 Spearman** for a strong mAb-trained baseline. *(n = 34: confidence
intervals are wide — see the lab journal. Honesty about uncertainty is a project value.)*

> **Why "honest"?** The public labels for part of the challenge test set are treated
> strictly as a *held-out* proxy. GloriPaka never injects them into training or feature
> selection — so the reported numbers estimate genuine generalization, not leaderboard fit.

## Key features

- 🧬 **Sequence → prediction.** One input (a VHH heavy-chain sequence); the package does
  structure prediction, featurization, and inference. *(End-to-end one-call API lands in
  v0.2 — see [Roadmap](#roadmap); the core library and honest evaluation are usable now.)*
- 🧠 **Transferable surface physics.** The most predictive, most transferable signal is
  solvent-exposed hydrophobic surface — captured by descriptors that work across the
  mAb → VHH domain gap, not by memorized sequence identity.
- ⚖️ **Leakage-free by construction.** Clone-family grouped cross-validation, in-fold
  feature selection, in-fold calibration; the validation proxy is quarantined.
- 🔬 **Interpretable + honest.** Every feature has a biophysical meaning; the full metric
  suite (Pearson, Spearman, Kendall, RMSE, MAE, R², outliers) is always recorded.
- 🔁 **Reproducible.** Fully-open default backends (ExtraTrees / HistGB / CatBoost /
  LightGBM / Ridge). `uv`-locked environment. TabPFN is an optional, clearly-marked backend.

## How it works

```
   VHH heavy-chain sequence
            │
            ▼
   3-D structure          NanoBodyBuilder2 (4-model conformational ensemble)
            │
            ▼
   biophysical features   FreeSASA surface + patches · exposed aromatic ·
            │             charge/pI · IMGT CDR stats · DeepSP SAP/SCM · ESM2
            ▼
   two predictors         (a) unfit physics composite  +  (b) native-trained Ridge
            │
            ▼
   rank-mean ensemble  →  calibrated HIC value + per-sequence confidence tier
```

## Installation

GloriPaka uses [**uv**](https://docs.astral.sh/uv/) for environment and package management.

```bash
# 1. Clone
git clone https://github.com/GloriusGroup/GloriPaka.git
cd GloriPaka

# 2. Create the environment and install (core library)
uv sync

# 3. (optional) add the featurization / structure-prediction stack
uv sync --extra featurize

# 4. (optional) developer tooling and docs
uv sync --extra dev --extra docs
```

Run anything inside the environment with `uv run`, e.g. `uv run python -c "import gloripaka; print(gloripaka.__version__)"`.

## Quick start

The core library is available today — honest metrics and the unfit physics composite:

```python
import pandas as pd
from gloripaka.physics import physics_composite, PHYSICS_COMPOSITE
from gloripaka.metrics import compute_metrics

# `features` is a table with the four physics descriptors as columns
# (see gloripaka.config.PHYSICS_COMPOSITE / the featurization extra to compute them).
score = physics_composite(features)          # equal-weight rank composite, no training
print(compute_metrics(true_hic, score))      # Pearson / Spearman / Kendall / RMSE / ...
```

The one-call, sequence-to-prediction API — **the headline goal, arriving in v0.2**:

```python
# COMING IN v0.2 (Phase 2, see ROADMAP.md)
from gloripaka import GloriPakaPredictor

predictor = GloriPakaPredictor()             # structure + featurize + ensemble
result = predictor.predict([
    "QVQLVESGGGLVQPGGSLRLSCAAS...",           # VHH heavy-chain sequence(s)
])
print(result[["sequence_id", "predicted_hic", "confidence"]])
```

## Package layout

| Module | Purpose |
|:--|:--|
| `gloripaka.config` | Central paths, feature-group prefixes, shared constants (self-contained). |
| `gloripaka.physics` | The unfit, interpretable physics composite (transferable hydrophobic-surface baseline). |
| `gloripaka.models` | Rank-target regression + in-fold feature selection; open backends + optional TabPFN; scale calibration. |
| `gloripaka.features` | Feature-table loading, group slicing, label-joined dataset assembly. |
| `gloripaka.labels` | Label loading + the honest known/unknown split of the challenge VHH. |
| `gloripaka.validation` | Clone-family clustering + leave-clone-family-out cross-validation. |
| `gloripaka.metrics` | The full metric suite + bootstrap confidence intervals + residual-outlier flags. |
| `gloripaka.report` | Lightweight experiment reporting (CSV + markdown tables). |
| `gloripaka.secrets` | Secure, out-of-repo access to the optional TabPFN API token. |
| `gloripaka.viz` | Single central publication figure style. |
| `gloripaka.featurize` *(v0.2)* | Sequence → structure → biophysical features. |
| `gloripaka.pipeline` *(v0.2)* | High-level `GloriPakaPredictor` (sequence → HIC). |

## Data & reproducibility

Small canonical inputs (challenge sequences, labels, and the feature tables the models
consume) ship in `data/`. The **heavy, regenerable artifacts** — predicted PDB
structures and full feature matrices — are hosted on **Zenodo** (DOI pending) and fetched
with a helper script, keeping the repository lean. A small in-repo sample (the 34 known
VHH) supports examples and tests without any download. See `data/README.md` (added in
Phase 3) for provenance and fetch instructions.

## Scientific background

The full experimental record — every experiment, what worked, and the honest negatives —
lives in the research lab journal (`modeling_strategyB/notes/lab_journal.md` in the source
workspace) and is being curated into `docs/` here. Core findings: native-VHH training beats
mAb-only training; transferable **surface hydrophobicity** is the signal that generalizes; a
**physics + ML ensemble** beats either alone; and the ceiling is **data quantity on a shared
assay**, not method. *(The accompanying manuscript is kept out of this repository pending
figure/publisher licensing; it will be linked once cleared.)*

## Roadmap

GloriPaka is built minimal-first, then expanded. See **[`ROADMAP.md`](ROADMAP.md)** for the
full phased plan (pipeline → data/Zenodo → literature → tests/docs/CI → release).

## Citation

If you use GloriPaka, please cite it (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@software{gloripaka2026,
  title  = {GloriPaka: sequence-based HIC-retention prediction for VHH nanobodies},
  author = {Boser, Florian and {Glorius Group}},
  year   = {2026},
  url    = {https://github.com/GloriusGroup/GloriPaka},
  note   = {Version 0.1.0}
}
```

## Authors & acknowledgements

Developed in the **Glorius Group** (University of Münster). Built on open tools —
NanoBodyBuilder2, FreeSASA, ANARCI/IMGT, DeepSP, ESM2, scikit-learn — and on public
antibody/HIC datasets acknowledged in the documentation.

## License

Released under the [MIT License](LICENSE).
