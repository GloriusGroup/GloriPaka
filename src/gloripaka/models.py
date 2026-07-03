"""Models: rank-target regression with in-fold feature selection, plus calibration.

All estimators here follow scikit-learn's ``fit``/``predict`` protocol so they drop
into :class:`gloripaka.validation.GroupedCrossValidator` unchanged.

Two ideas do most of the work on this small, assay-noisy problem:

  1. **Rank-target learning.** HIC retention has assay outliers and is only
     rank-comparable across scales; training against ``rankdata(y)`` (and reporting
     Spearman) is both more robust and better matched to the objective.
  2. **In-fold univariate selection.** With 613 features and ~240 rows, blindly
     feeding everything invites overfitting. ``SelectKBest(f_regression)`` fit *only*
     on the training fold trims to the most target-correlated columns without leakage.

The default backends are fully open (ExtraTrees / HistGB / CatBoost / LightGBM).
TabPFN is reproducible but license-gated and download-dependent, so it is an opt-in
backend that fails loudly with guidance if its weights are unavailable.
"""

from __future__ import annotations

from typing import Callable, Protocol

import numpy as np
from scipy.stats import rankdata
from sklearn.feature_selection import SelectKBest, f_regression


class Regressor(Protocol):
    """Structural type for anything usable as a base estimator."""

    def fit(self, X: np.ndarray, y: np.ndarray) -> "Regressor": ...
    def predict(self, X: np.ndarray) -> np.ndarray: ...


class RankTransformSelectRegressor:
    """Base regressor + rank-target transform + in-fold top-k feature selection.

    Model-agnostic reimplementation of the modeling team's ``RankSelectRegressor``
    (which was TabPFN-specific). All three transforms are fit strictly on the data
    passed to :meth:`fit`, so they add no leakage in cross-validation.

    Predictions are mapped back onto the training target's empirical distribution
    (monotone), keeping outputs interpretable as HIC-like values while leaving
    Spearman unchanged.
    """

    def __init__(self, base_factory: Callable[[], Regressor], k: int = 60) -> None:
        self.base_factory = base_factory
        self.k = k

    def fit(self, X: np.ndarray, y: np.ndarray) -> "RankTransformSelectRegressor":
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        y_rank = rankdata(y)

        k = min(self.k, X.shape[1])
        self.selector_ = SelectKBest(f_regression, k=k)
        # f_regression is undefined on zero-variance columns; guard with nan_to_num.
        Xs = self.selector_.fit_transform(np.nan_to_num(X), y_rank)

        self.model_ = self.base_factory()
        self.model_.fit(Xs, y_rank)

        # empirical quantile function: predicted rank -> training-target value
        self._ranks = np.arange(1, len(y) + 1)
        self._sorted_y = np.sort(y)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=float)
        Xs = self.selector_.transform(np.nan_to_num(X))
        r = self.model_.predict(Xs)
        return np.interp(r, self._ranks, self._sorted_y)


# --------------------------------------------------------------------------- #
# Base-model factories (open, reproducible, no license or network required).   #
# --------------------------------------------------------------------------- #
def _extratrees() -> Regressor:
    from sklearn.ensemble import ExtraTreesRegressor

    return ExtraTreesRegressor(
        n_estimators=500, min_samples_leaf=3, max_features=0.8, random_state=0, n_jobs=-1
    )


def _histgb() -> Regressor:
    from sklearn.ensemble import HistGradientBoostingRegressor

    return HistGradientBoostingRegressor(
        max_depth=3, learning_rate=0.05, max_iter=400, l2_regularization=1.0,
        min_samples_leaf=15, random_state=0,
    )


def _catboost() -> Regressor:
    from catboost import CatBoostRegressor

    return CatBoostRegressor(
        iterations=400, learning_rate=0.025, depth=3, l2_leaf_reg=10.0,
        loss_function="RMSE", random_seed=0, verbose=False,
        allow_writing_files=False, thread_count=-1,
    )


def _lightgbm() -> Regressor:
    from lightgbm import LGBMRegressor

    return LGBMRegressor(
        n_estimators=300, learning_rate=0.02, num_leaves=7, min_child_samples=12,
        subsample=0.85, colsample_bytree=0.75, reg_lambda=3.0, random_state=0,
        verbosity=-1, force_col_wise=True,
    )


def _ridge() -> Regressor:
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    return make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25)))


# TabPFN's default model version (v2.5/v3) is license-gated and prompts for a
# browser token - impossible in a headless run. The original v2 regressor weights,
# however, are a PUBLIC HuggingFace repo: fetching that checkpoint and passing it as
# ``model_path`` bypasses the license prompt while staying fully reproducible
# (public weights + open inference code). This is how we run TabPFN in a closed env.
_TABPFN_V2_REPO = "Prior-Labs/TabPFN-v2-reg"
_TABPFN_V2_CKPT = "tabpfn-v2-regressor.ckpt"


def _resolve_tabpfn_v2_weights() -> str:
    """Return a local path to the public TabPFN-v2 regressor checkpoint.

    Downloads once into the HuggingFace cache (no auth token required) and reuses it
    thereafter, so repeated runs are offline and deterministic.
    """
    from huggingface_hub import hf_hub_download

    return hf_hub_download(repo_id=_TABPFN_V2_REPO, filename=_TABPFN_V2_CKPT)


def _tabpfn() -> Regressor:
    """Reproducible TabPFN-v2 backend using public weights (no license token)."""
    try:
        from tabpfn import TabPFNRegressor
    except Exception as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "TabPFN is not installed. Add it with `uv add tabpfn` (extra 'tabpfn')."
        ) from exc
    return TabPFNRegressor(device="cpu", random_state=0, model_path=_resolve_tabpfn_v2_weights())


def _tabpfn_api() -> Regressor:
    """Hosted TabPFN-v3 API backend (thinking mode).

    NOT reproducible in the strict sense: it egresses features to Prior Labs servers
    and the server-side model version can change. Kept for a strong *cross-check* of
    the local models - never as the reproducible primary. Requires a token
    (``gloripaka.secrets``); raises with guidance if absent.
    """
    from gloripaka.secrets import get_tabpfn_token

    token = get_tabpfn_token()
    if not token:
        raise RuntimeError(
            "No TabPFN API token found (env TABPFN_TOKEN or ~/.config/priorlabs/tabpfn_token)."
        )
    import tabpfn_client

    tabpfn_client.set_access_token(token)
    from tabpfn_client import TabPFNRegressor as APIRegressor

    try:
        return APIRegressor(thinking_mode=True)  # v3-Plus reasoning, if supported
    except TypeError:
        return APIRegressor()  # older client without thinking_mode kwarg


_BASE_FACTORIES: dict[str, Callable[[], Regressor]] = {
    "extratrees": _extratrees,
    "histgb": _histgb,
    "catboost": _catboost,
    "lightgbm": _lightgbm,
    "ridge": _ridge,
    "tabpfn": _tabpfn,  # local v2 public weights - reproducible, offline
    "tabpfn_api": _tabpfn_api,  # hosted v3 thinking-mode - cross-check only (data egress)
}


class ModelRegistry:
    """Build model *factories* by name, optionally wrapped in rank-transform+select.

    A *factory* (zero-arg callable returning a fresh estimator) is what the CV runner
    needs, so it can rebuild an untrained model per fold.
    """

    @staticmethod
    def available() -> list[str]:
        return sorted(_BASE_FACTORIES)

    @staticmethod
    def factory(
        name: str, rank_select: bool = True, k: int = 60
    ) -> Callable[[], Regressor]:
        if name not in _BASE_FACTORIES:
            raise KeyError(f"Unknown model '{name}'. Available: {ModelRegistry.available()}")
        base = _BASE_FACTORIES[name]
        if not rank_select:
            return base
        return lambda: RankTransformSelectRegressor(base, k=k)  # type: ignore[return-value]


def within_group_quantiles(y: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Map targets to per-group quantiles in [0, 1] so different assay scales pool.

    Absolute HIC minutes are not comparable across assays (World-1 mAb 2.4-4.5 vs
    World-2 ARRAS VHH 4.7-6.7). Ranking each scale group internally and expressing it
    as a quantile puts every measurement on a shared 0..1 "relative retention" axis,
    which is the only honest way to train one model across scales. Singleton groups
    map to 0.5.
    """
    y = np.asarray(y, dtype=float)
    groups = np.asarray(groups)
    q = np.full(len(y), 0.5, dtype=float)
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        if len(idx) == 1:
            continue
        r = rankdata(y[idx])
        q[idx] = (r - 1) / (len(idx) - 1)
    return q


class MonotoneCalibrator:
    """Isotonic map from model scores to a target (absolute) scale.

    Spearman is invariant to any monotone transform, so calibration only affects
    Pearson/RMSE. It MUST be fit in-fold on training-scale pairs - never on the
    validation labels - to stay leak-free. Used to place rank-unit predictions back
    onto the ARRAS minute scale for the absolute-Pearson objective.
    """

    def __init__(self, increasing: bool = True) -> None:
        self.increasing = increasing

    def fit(self, scores: np.ndarray, target: np.ndarray) -> "MonotoneCalibrator":
        from sklearn.isotonic import IsotonicRegression

        self.iso_ = IsotonicRegression(increasing=self.increasing, out_of_bounds="clip")
        self.iso_.fit(np.asarray(scores, dtype=float), np.asarray(target, dtype=float))
        return self

    def transform(self, scores: np.ndarray) -> np.ndarray:
        return self.iso_.predict(np.asarray(scores, dtype=float))
