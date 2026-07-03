"""Feature-matrix loading, group slicing, and label-joined dataset assembly.

We reuse the modeling team's precomputed 613-feature tables (ESM2 + structural +
DeepSP + descriptors + CDR stats + VHH-specific) rather than recomputing them. This
module joins those features to the honest labels from :mod:`gloripaka.labels`, attaches
clone-family groups for leak-free CV, and exposes clean numpy matrices to the models.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from gloripaka.config import FEATURE_GROUPS, PATHS
from gloripaka.labels import known_test_labels, load_gdpa1_labels, load_test_sequences
from gloripaka.validation import CloneFamilyClusterer

# Tables from the shared featurization era (read-only) + tables WE produce (workspace).
_FEATURE_FILES: dict[str, str] = {
    "gdpa1": "gdpa1_all.csv",
    "vhh": "vhh_all.csv",
    "arras37": "arras37_all.csv",
    "native": "native_vhh_all.csv",  # GDPa2 + TNP + Siegmund (we featurized)
    "world1": "world1_all.csv",  # GDPa3 + Shehata + Jain VH (we featurized)
}


def load_feature_table(key: str) -> pd.DataFrame:
    """Load a 613-feature table indexed by its ``Name`` key.

    Searches the shared read-only features dir first, then our writable workspace
    dir, so tables we generate coexist with the inherited ones transparently.
    """
    if key not in _FEATURE_FILES:
        raise KeyError(f"Unknown feature table '{key}'. Available: {list(_FEATURE_FILES)}")
    fname = _FEATURE_FILES[key]
    for base in (PATHS.features_dir, PATHS.workspace_features_dir):
        path = base / fname
        if path.exists():
            df = pd.read_csv(path)
            df["Name"] = df["Name"].astype(str)
            return df.set_index("Name")
    raise FileNotFoundError(f"Feature table '{key}' ({fname}) not found in shared or workspace dir.")


def feature_columns(df: pd.DataFrame, groups: list[str] | None = None) -> list[str]:
    """Return feature column names, optionally restricted to feature ``groups``.

    ``groups`` is a subset of :data:`gloripaka.config.FEATURE_GROUPS` keys (e.g.
    ``["struct", "vhh"]``). ``None`` returns every feature column.
    """
    if groups is None:
        prefixes = tuple(FEATURE_GROUPS.values())
    else:
        unknown = set(groups) - set(FEATURE_GROUPS)
        if unknown:
            raise KeyError(f"Unknown feature groups {unknown}. Known: {list(FEATURE_GROUPS)}")
        prefixes = tuple(FEATURE_GROUPS[g] for g in groups)
    return [c for c in df.columns if c.startswith(prefixes)]


@dataclass
class Dataset:
    """A label-joined, CV-ready dataset.

    Attributes:
        key: short identifier (for logs / figures).
        X: feature DataFrame (index = Name, columns = features only).
        y: aligned target vector (HIC), or ``None`` for unlabeled prediction sets.
        groups: aligned clone-family / fold ids for grouped CV, or ``None``.
        meta: aligned side-table (sequence, scale, source, etc.).
    """

    key: str
    X: pd.DataFrame
    y: np.ndarray | None
    groups: np.ndarray | None
    meta: pd.DataFrame

    def matrix(self, groups: list[str] | None = None) -> tuple[np.ndarray, list[str]]:
        """Return ``(values, column_names)`` for the requested feature groups."""
        cols = feature_columns(self.X, groups)
        return self.X[cols].to_numpy(dtype=float), cols

    @property
    def names(self) -> list[str]:
        return list(self.X.index)


class DatasetFactory:
    """Builds the standard Strategy-B datasets with consistent joins and grouping.

    Clone families are computed **once** on all 81 test VHH so the family id of a
    known VHH is stable whether it appears in the validation proxy or the full set.
    """

    def __init__(self, clone_threshold: float = 0.95) -> None:
        self._test_seqs = load_test_sequences()
        families = CloneFamilyClusterer(threshold=clone_threshold).fit_predict(
            list(self._test_seqs["sequence"])
        )
        self._test_seqs = self._test_seqs.assign(family=families)

    def gdpa1(self) -> Dataset:
        """GDPa1 training set: 242 mAb with HIC and predefined hierarchical folds."""
        feats = load_feature_table("gdpa1")
        labels = load_gdpa1_labels().set_index("Name")
        common = feats.index.intersection(labels.index)
        feats, labels = feats.loc[common], labels.loc[common]
        return Dataset(
            key="gdpa1",
            X=feats,
            y=labels["hic"].to_numpy(dtype=float),
            groups=labels["fold"].to_numpy(),
            meta=labels[["hic", "fold", "vh_sequence"]].rename(columns={"vh_sequence": "sequence"}),
        )

    def test_all(self) -> Dataset:
        """All 81 challenge VHH (unlabeled), with clone-family groups and known flag."""
        feats = load_feature_table("vhh")
        meta = self._test_seqs.set_index("Name")
        known = set(known_test_labels()["Name"])
        meta = meta.assign(is_known=[n in known for n in meta.index])
        common = feats.index.intersection(meta.index)
        feats, meta = feats.loc[common], meta.loc[common]
        return Dataset(
            key="test81",
            X=feats,
            y=None,
            groups=meta["family"].to_numpy(),
            meta=meta,
        )

    def native_pool(self) -> Dataset:
        """Native VHH we featurized ourselves (GDPa2 + TNP + Siegmund), labeled.

        These are disjoint by source from the 34 ARRAS validation VHH, so they can
        train a model that is then evaluated on the held-out 34 with no leakage. Each
        source sits on its own HIC assay scale (``groups`` = scale_group) -> rank-
        normalize per group before pooling (see :func:`gloripaka.models.within_group_quantiles`).
        """
        feats = load_feature_table("native")
        labels = pd.read_csv(PATHS.workspace / "data" / "native_vhh_labels.csv")
        labels["Name"] = labels["Name"].astype(str)
        labels["hic"] = pd.to_numeric(labels["hic"], errors="coerce")
        labels = labels[labels["hic"].notna()].set_index("Name")
        common = feats.index.intersection(labels.index)
        feats, labels = feats.loc[common], labels.loc[common]
        return Dataset(
            key="native_pool",
            X=feats,
            y=labels["hic"].to_numpy(dtype=float),
            groups=labels["scale_group"].to_numpy(),
            meta=labels[["hic", "scale_group", "source_key", "format"]],
        )

    def known_vhh(self) -> Dataset:
        """The 34 known ARRAS test VHH = the held-out validation proxy (labeled)."""
        feats = load_feature_table("vhh")
        known = known_test_labels().set_index("Name")
        fam = self._test_seqs.set_index("Name")["family"]
        known = known.assign(family=fam.reindex(known.index).to_numpy())
        common = feats.index.intersection(known.index)
        feats, known = feats.loc[common], known.loc[common]
        return Dataset(
            key="known_vhh34",
            X=feats,
            y=known["hic"].to_numpy(dtype=float),
            groups=known["family"].to_numpy(),
            meta=known[["sequence", "hic", "family", "scale_group"]],
        )
