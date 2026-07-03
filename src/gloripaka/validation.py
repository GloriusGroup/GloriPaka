"""Honest cross-validation: clone-family grouping + leak-free fold running.

The 81 challenge VHH contain engineered near-duplicates (point-mutant series of the
same parent). If two variants of one clone straddle a train/test split, the model
memorizes the clone rather than learning transferable physics — silent leakage that
inflates CV scores. We defend against it by clustering sequences into **clone
families** and only ever splitting *between* families (leave-clone-family-out).
"""

from __future__ import annotations

from dataclasses import dataclass
from difflib import SequenceMatcher
from typing import Callable, Iterator

import numpy as np


def sequence_identity(a: str, b: str) -> float:
    """Approximate global sequence identity in [0, 1].

    Uses :class:`difflib.SequenceMatcher` (Ratcliff/Obershelp), which is symmetric,
    dependency-free, and reliable for near-duplicate detection among VHH domains of
    similar length. It is a *clustering* tool, not a substitute for a true alignment
    score, and is used only to decide which sequences belong to the same clone family.
    """
    return SequenceMatcher(None, a, b).ratio()


@dataclass
class CloneFamilyClusterer:
    """Single-linkage clustering of sequences into clone families.

    Two sequences join the same family if their identity is >= ``threshold``.
    Single linkage is deliberate: an engineered mutation ladder A->B->C where A~C
    are more diverged than adjacent steps should still form one family, because any
    two adjacent members exceed the threshold.

    The default 0.95 was calibrated on the 81 test VHH: their nearest-neighbour
    identity distribution has ~31 genuinely distinct clones at NN-id < 0.93 and a
    tail of engineered point-mutant twins at NN-id >= 0.97. A naive 0.90 threshold
    merges *different* clones through their shared ~100-residue framework (one family
    of 54); 0.95 groups true variants (differing by <~6 residues) while keeping
    distinct clones apart. See the lab journal E1 entry for the full analysis.
    """

    threshold: float = 0.95

    def fit_predict(self, sequences: list[str]) -> np.ndarray:
        """Return an integer family label per input sequence (contiguous from 0)."""
        n = len(sequences)
        parent = list(range(n))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x: int, y: int) -> None:
            parent[find(x)] = find(y)

        for i in range(n):
            for j in range(i + 1, n):
                if sequence_identity(sequences[i], sequences[j]) >= self.threshold:
                    union(i, j)

        roots = [find(i) for i in range(n)]
        remap = {r: k for k, r in enumerate(sorted(set(roots)))}
        return np.array([remap[r] for r in roots], dtype=int)


@dataclass
class FoldResult:
    """Out-of-fold predictions plus the group membership that produced them."""

    y_true: np.ndarray
    y_pred: np.ndarray
    groups: np.ndarray  # group (family/fold) id each row was held out in


class GroupedCrossValidator:
    """Run a model factory over leave-one-group-out folds, collecting OOF predictions.

    The model is rebuilt from ``model_factory`` for every fold, fit on the training
    rows only, and used to predict the held-out group. Nothing about the held-out
    group (labels, target statistics, feature selection) touches the training fit —
    that discipline lives inside the model, and this runner simply never shows it the
    test fold.
    """

    def __init__(self, model_factory: Callable[[], object]) -> None:
        self.model_factory = model_factory

    @staticmethod
    def _splits(groups: np.ndarray) -> Iterator[tuple[np.ndarray, np.ndarray]]:
        for g in np.unique(groups):
            test = np.where(groups == g)[0]
            train = np.where(groups != g)[0]
            if len(train) and len(test):
                yield train, test

    def run(self, X: np.ndarray, y: np.ndarray, groups: np.ndarray) -> FoldResult:
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        groups = np.asarray(groups)
        oof = np.full(len(y), np.nan)
        for train, test in self._splits(groups):
            model = self.model_factory()
            model.fit(X[train], y[train])  # type: ignore[attr-defined]
            oof[test] = model.predict(X[test])  # type: ignore[attr-defined]
        return FoldResult(y_true=y, y_pred=oof, groups=groups)
