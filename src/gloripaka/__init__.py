"""GloriPaka: sequence-based HIC-retention prediction for VHH nanobodies.

GloriPaka predicts hydrophobic-interaction-chromatography (HIC) retention - a
surface-hydrophobicity / developability proxy - for VHH single-domain antibodies
(nanobodies) directly from their amino-acid sequence. Give it a sequence; the
package handles the rest: 3-D structure prediction (NanoBodyBuilder2), biophysical
featurization (solvent-exposed hydrophobic surface, patches, aromatic exposure,
charge/pI), and a leakage-free predictor.

The predictor is a rank-mean ensemble of two complementary estimators:

* an **unfit physics composite** - an equal-weight average of four solvent-exposed
  hydrophobic-surface descriptor ranks, defined a priori with no fitted coefficients,
  which is the honest low-data workhorse; and
* an **interpretable Ridge** trained on native VHH, dominated by the same exposed
  hydrophobic-surface / aromatic / charge biophysics.

Design commitments (see ``README.md``):

* **No leakage.** The 34 public ARRAS VHH labels are a *held-out validation proxy*;
  they never enter training or feature selection. Any absolute calibration is fit
  strictly in-fold.
* **Reproducible.** The default backends are fully open (ExtraTrees / HistGB /
  CatBoost / LightGBM / Ridge). TabPFN is an optional, clearly-marked backend.
* **Honest metrics.** Pearson + Spearman are the objective, but the full suite
  (Kendall, RMSE, MAE, R2, outliers) is always recorded to surface quirks.
"""

from __future__ import annotations

from gloripaka.pipeline import GloriPakaPredictor

__version__ = "0.1.0"

__all__ = ["GloriPakaPredictor", "__version__"]
