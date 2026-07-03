"""Evaluation metrics.

The challenge is scored by **Pearson + Spearman**, but we always record the full
suite so that quirks (rank/linear disagreement, heavy-tailed residuals, single-point
leverage) are visible in the lab journal rather than hidden behind one number.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.stats import kendalltau, pearsonr, spearmanr


@dataclass(frozen=True)
class RegressionMetrics:
    """A complete, comparable snapshot of predictive performance.

    ``pearson`` and ``spearman`` are the challenge objective; the rest are diagnostic.
    ``challenge_score`` is their mean (the leaderboard aggregate) for quick ranking.
    """

    n: int
    pearson: float
    spearman: float
    kendall: float
    rmse: float
    mae: float
    r2: float
    challenge_score: float

    def as_dict(self) -> dict[str, float]:
        return asdict(self)

    def __str__(self) -> str:  # compact one-liner for logs / journal
        return (
            f"n={self.n} | P={self.pearson:.3f} S={self.spearman:.3f} "
            f"K={self.kendall:.3f} | RMSE={self.rmse:.3f} MAE={self.mae:.3f} "
            f"R2={self.r2:.3f} | score={self.challenge_score:.3f}"
        )


def _safe_corr(fn, y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Correlation that returns NaN instead of raising on degenerate input."""
    if len(y_true) < 3 or np.std(y_true) == 0 or np.std(y_pred) == 0:
        return float("nan")
    return float(fn(y_true, y_pred)[0])


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> RegressionMetrics:
    """Compute the full metric suite for one set of predictions.

    ``r2`` is the coefficient of determination on the raw values; it can be very
    negative when predictions live on a different absolute scale than the target
    (e.g. rank-unit outputs vs. HIC minutes) - that is expected and informative,
    which is precisely why we keep Pearson/Spearman as the primary objective.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true, y_pred = y_true[mask], y_pred[mask]

    pearson = _safe_corr(pearsonr, y_true, y_pred)
    spearman = _safe_corr(spearmanr, y_true, y_pred)
    kendall = _safe_corr(kendalltau, y_true, y_pred)
    resid = y_pred - y_true
    rmse = float(np.sqrt(np.mean(resid**2))) if len(resid) else float("nan")
    mae = float(np.mean(np.abs(resid))) if len(resid) else float("nan")
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    challenge = np.nanmean([pearson, spearman])

    return RegressionMetrics(
        n=int(len(y_true)),
        pearson=pearson,
        spearman=spearman,
        kendall=kendall,
        rmse=rmse,
        mae=mae,
        r2=r2,
        challenge_score=float(challenge),
    )


def bootstrap_ci(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    metric: str = "spearman",
    n_boot: int = 2000,
    seed: int = 0,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """Percentile bootstrap confidence interval for a correlation metric.

    Returns ``(point_estimate, lo, hi)``. With only ~34 validation points the CI is
    wide; reporting it keeps us honest about how much any single-number improvement
    actually means.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    fn = {"spearman": spearmanr, "pearson": pearsonr, "kendall": kendalltau}[metric]
    point = _safe_corr(fn, y_true, y_pred)

    rng = np.random.default_rng(seed)
    n = len(y_true)
    stats = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        stats[b] = _safe_corr(fn, y_true[idx], y_pred[idx])
    lo, hi = np.nanpercentile(stats, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return point, float(lo), float(hi)


def residual_outliers(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    labels: list[str] | None = None,
    z_thresh: float = 2.5,
) -> list[dict[str, object]]:
    """Flag points whose standardized residual exceeds ``z_thresh``.

    Uses a rank-based standardization of residuals (robust MAD scaling) so a couple
    of extreme misses don't inflate the scale and hide the rest. Returned records
    are sorted by absolute z, most extreme first - the shortlist for the journal.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    resid = y_pred - y_true
    med = np.median(resid)
    mad = np.median(np.abs(resid - med)) or 1e-9
    z = 0.6745 * (resid - med) / mad  # 0.6745 => MAD-to-sigma for normal data

    out: list[dict[str, object]] = []
    for i, zi in enumerate(z):
        if abs(zi) >= z_thresh:
            out.append(
                {
                    "label": labels[i] if labels is not None else str(i),
                    "y_true": float(y_true[i]),
                    "y_pred": float(y_pred[i]),
                    "residual": float(resid[i]),
                    "z": float(zi),
                }
            )
    out.sort(key=lambda r: -abs(r["z"]))  # type: ignore[arg-type]
    return out
