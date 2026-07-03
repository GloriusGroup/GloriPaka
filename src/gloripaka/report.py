"""Lightweight experiment reporting: accumulate metric rows, save CSV + markdown.

Every experiment produces one :class:`ExperimentReport`. It writes a machine-readable
CSV into ``results/`` and returns a markdown table ready to paste into the lab journal,
so the journal and the results stay in lock-step without manual transcription.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from gloripaka.config import PATHS
from gloripaka.metrics import RegressionMetrics


@dataclass
class ExperimentReport:
    """Collects labelled metric rows for one experiment."""

    name: str
    rows: list[dict[str, object]] = field(default_factory=list)

    def add(self, label: str, metrics: RegressionMetrics, **extra: object) -> None:
        """Append one result row (a model/config + its metrics + any extra columns)."""
        self.rows.append({"label": label, **metrics.as_dict(), **extra})

    def to_frame(self) -> pd.DataFrame:
        df = pd.DataFrame(self.rows)
        if "challenge_score" in df.columns:
            df = df.sort_values("challenge_score", ascending=False, ignore_index=True)
        return df

    def save_csv(self) -> Path:
        PATHS.results_dir.mkdir(parents=True, exist_ok=True)
        path = PATHS.results_dir / f"{self.name}.csv"
        self.to_frame().to_csv(path, index=False)
        return path

    def to_markdown(self, cols: list[str] | None = None, floatfmt: str = ".3f") -> str:
        df = self.to_frame()
        if cols is not None:
            df = df[[c for c in cols if c in df.columns]]
        return df.to_markdown(index=False, floatfmt=floatfmt)
