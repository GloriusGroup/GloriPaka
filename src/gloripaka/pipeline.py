"""High-level prediction pipeline: biophysical features (or sequences) to HIC.

This module ties the individual components together into one predictor. The shipped
model is a rank-mean ensemble of two estimators (see :mod:`gloripaka.physics` and
:mod:`gloripaka.models`):

1. an unfit physics composite (four exposed-hydrophobic-surface descriptor ranks,
   no training), and
2. an interpretable Ridge regression trained on the native-VHH pool.

Two entry points are provided:

* :meth:`GloriPakaPredictor.predict_from_features` runs on a table of precomputed
  features. This is fully self-contained and needs only the core dependencies.
* :meth:`GloriPakaPredictor.predict` starts from raw sequences. It requires the
  optional ``featurize`` extra (structure prediction and biophysical featurization)
  and raises an informative error if that stack is not installed.

No test labels are used to fit anything; the held-out labels are only ever read to
score predictions after the fact.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from gloripaka.config import PATHS
from gloripaka.features import feature_columns, load_feature_table
from gloripaka.labels import known_test_labels, load_test_sequences
from gloripaka.metrics import RegressionMetrics, bootstrap_ci, compute_metrics
from gloripaka.models import ModelRegistry, within_group_quantiles
from gloripaka.physics import physics_composite
from gloripaka.validation import sequence_identity

# Public ARRAS HIC range (minutes). Used only to place rank predictions on an
# interpretable absolute axis; the map is monotone, so Spearman is unaffected.
ARRAS_MIN, ARRAS_MAX = 4.5, 7.3

# Interpretable feature families the Ridge is allowed to select from (ESM2
# embeddings are deliberately excluded to keep the model interpretable).
INTERP_GROUPS = ["desc", "cdr", "deepsp", "struct", "vhh"]


def _clean_seq(sequence: str) -> str:
    """Return an upper-cased, whitespace-stripped sequence."""
    return str(sequence).strip().upper()


def _rank01(values: np.ndarray) -> np.ndarray:
    """Map values to ranks scaled to the closed unit interval ``[0, 1]``."""
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.zeros(len(values), dtype=float)
    return (rankdata(values) - 1.0) / (len(values) - 1.0)


def _confidence_tier(nn_identity: float, is_known: bool) -> str:
    """Classify a prediction's trustworthiness from its nearest training neighbour.

    Args:
        nn_identity: Maximum sequence identity to any training sequence, in ``[0, 1]``.
        is_known: Whether the sequence has a public HIC label (part of the proxy set).

    Returns:
        One of ``"known"``, ``"high"``, ``"medium"`` or ``"low"``.
    """
    if is_known:
        return "known"
    if nn_identity >= 0.90:
        return "high"
    if nn_identity >= 0.80:
        return "medium"
    return "low"


@dataclass
class _NativePool:
    """The native-VHH training pool: features, rank-normalized target, sequences."""

    X: pd.DataFrame
    y_quantile: np.ndarray
    sequences: list[str]


class GloriPakaPredictor:
    """Predict VHH HIC retention with the leakage-free physics-plus-Ridge ensemble.

    The predictor trains its Ridge component on the bundled native-VHH pool the first
    time it is needed and caches it. The physics component needs no training.

    Args:
        w_physics: Weight of the physics component in the rank-mean ensemble; the
            Ridge component gets ``1 - w_physics``. Defaults to ``0.5`` (validated
            optimum; tuning was found to be saturated).
        k: Number of features the Ridge selects in-fold via univariate F-test.
    """

    def __init__(self, w_physics: float = 0.5, k: int = 40) -> None:
        if not 0.0 <= w_physics <= 1.0:
            raise ValueError("w_physics must be in [0, 1].")
        self.w_physics = w_physics
        self.k = k
        self._pool: _NativePool | None = None
        self._ridge = None
        self._interp: list[str] | None = None

    # -- training ------------------------------------------------------------ #
    def _load_pool(self) -> _NativePool:
        """Load and cache the native-VHH training pool from the bundled data."""
        if self._pool is not None:
            return self._pool
        features = load_feature_table("native")
        labels = pd.read_csv(PATHS.data / "native_vhh_labels.csv")
        labels["Name"] = labels["Name"].astype(str)
        labels["hic"] = pd.to_numeric(labels["hic"], errors="coerce")
        labels = labels[labels["hic"].notna()].set_index("Name")
        common = features.index.intersection(labels.index)
        features, labels = features.loc[common], labels.loc[common]
        y_quantile = within_group_quantiles(
            labels["hic"].to_numpy(dtype=float), labels["scale_group"].to_numpy()
        )
        sequences = [_clean_seq(s) for s in labels["sequence"].astype(str)]
        self._pool = _NativePool(X=features, y_quantile=y_quantile, sequences=sequences)
        return self._pool

    def _fit_ridge(self) -> None:
        """Fit and cache the interpretable Ridge on the native pool."""
        if self._ridge is not None:
            return
        pool = self._load_pool()
        self._interp = feature_columns(pool.X, INTERP_GROUPS)
        k = min(self.k, len(self._interp))
        model = ModelRegistry.factory("ridge", rank_select=True, k=k)()
        model.fit(pool.X[self._interp].to_numpy(dtype=float), pool.y_quantile)
        self._ridge = model

    # -- prediction ---------------------------------------------------------- #
    def predict_from_features(
        self, features: pd.DataFrame, sequences: pd.Series | None = None
    ) -> pd.DataFrame:
        """Predict HIC for a table of precomputed features.

        Args:
            features: Feature table indexed by name, with at least the physics
                composite columns and the interpretable feature families.
            sequences: Optional per-row sequences (aligned to ``features.index``),
                used to attach a nearest-neighbour identity and confidence tier.

        Returns:
            A DataFrame indexed like ``features`` with columns ``predicted_rank`` and
            ``predicted_hic`` (and ``nn_identity_to_train`` / ``confidence`` when
            ``sequences`` is given), sorted by descending predicted HIC.
        """
        self._fit_ridge()
        assert self._interp is not None and self._ridge is not None

        physics = _rank01(physics_composite(features))
        ridge = _rank01(self._ridge.predict(features[self._interp].to_numpy(dtype=float)))
        score = self.w_physics * physics + (1.0 - self.w_physics) * ridge
        hic = ARRAS_MIN + _rank01(score) * (ARRAS_MAX - ARRAS_MIN)

        out = pd.DataFrame(
            {"predicted_rank": np.round(score, 4), "predicted_hic": np.round(hic, 4)},
            index=features.index,
        )
        if sequences is not None:
            out = out.join(self._confidence(features.index, sequences))
        return out.sort_values("predicted_hic", ascending=False)

    def predict(self, sequences: list[str]) -> pd.DataFrame:
        """Predict HIC directly from raw VHH heavy-chain sequences.

        This is the headline interface: structure prediction and biophysical
        featurization happen internally before the ensemble runs. It requires the
        optional ``featurize`` extra.

        Args:
            sequences: VHH heavy-chain amino-acid sequences.

        Returns:
            The prediction table from :meth:`predict_from_features`.

        Raises:
            RuntimeError: If the ``featurize`` extra (structure prediction and
                biophysical featurization) is not installed.
        """
        try:
            from gloripaka.featurize import featurize_sequences
        except ImportError as exc:  # pragma: no cover - depends on optional extra
            raise RuntimeError(
                "Predicting from raw sequences needs the featurization stack "
                "(NanoBodyBuilder2, FreeSASA, ANARCI, DeepSP). Install it with "
                "`uv sync --extra featurize`. Alternatively, call "
                "`predict_from_features` on a precomputed feature table."
            ) from exc
        seqs = [_clean_seq(s) for s in sequences]
        features = featurize_sequences(seqs)
        return self.predict_from_features(features, sequences=pd.Series(seqs, index=features.index))

    # -- evaluation ---------------------------------------------------------- #
    def evaluate(self, features: pd.DataFrame, y_true: np.ndarray) -> RegressionMetrics:
        """Score predictions on labelled features (labels are used only here)."""
        preds = self.predict_from_features(features).reindex(features.index)
        return compute_metrics(np.asarray(y_true, dtype=float), preds["predicted_rank"].to_numpy())

    # -- helpers ------------------------------------------------------------- #
    def _confidence(self, index: pd.Index, sequences: pd.Series) -> pd.DataFrame:
        """Attach nearest-neighbour identity and a confidence tier per row."""
        pool = self._load_pool()
        try:
            known_names = set(known_test_labels()["Name"])
        except Exception:  # pragma: no cover - harvest data optional
            known_names = set()
        train = pool.sequences
        rows = []
        for name in index:
            seq = _clean_seq(sequences.get(name, ""))
            nn = max((sequence_identity(seq, t) for t in train if t), default=0.0)
            rows.append(
                {
                    "nn_identity_to_train": round(float(nn), 3),
                    "confidence": _confidence_tier(nn, name in known_names),
                }
            )
        return pd.DataFrame(rows, index=index)


def predict_challenge_test_set() -> pd.DataFrame:
    """Predict all 81 challenge VHH from the bundled feature table.

    Returns:
        The 81-row prediction table with confidence tiers.
    """
    features = load_feature_table("vhh")
    seqs = load_test_sequences().set_index("Name")["sequence"]
    predictor = GloriPakaPredictor()
    return predictor.predict_from_features(features, sequences=seqs.reindex(features.index))


def evaluate_known_proxy() -> tuple[RegressionMetrics, tuple[float, float, float], tuple[float, float, float]]:
    """Score the ensemble on the 34 held-out VHH with public HIC labels.

    Returns:
        The metric suite plus bootstrap confidence intervals for Spearman and
        Pearson, computed on the held-out proxy (never used in training).
    """
    known = known_test_labels().set_index("Name")
    features = load_feature_table("vhh").reindex(known.index)
    predictor = GloriPakaPredictor()
    preds = predictor.predict_from_features(features).reindex(known.index)
    y = known["hic"].to_numpy(dtype=float)
    score = preds["predicted_rank"].to_numpy()
    metrics = compute_metrics(y, score)
    return metrics, bootstrap_ci(y, score, "spearman"), bootstrap_ci(y, score, "pearson")
