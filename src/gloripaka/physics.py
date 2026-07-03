"""Transferable-physics baselines: unfit, interpretable, leakage-free.

These are the domain-knowledge features that should transfer across the mAb->VHH gap
(hydrophobic surface / patch), assembled a priori — NOT fit to any ARRAS label. They
are the honest yardstick every learned model must beat to justify its complexity.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata

# Hydrophobic-surface descriptors, all "higher => more hydrophobic/retained".
# Chosen from domain knowledge (memory: reference-hic-features), not from ARRAS fit.
PHYSICS_COMPOSITE: list[str] = [
    "vhh_A_positive_sasa",  # WW positive-only hydrophobic SASA (most transferable single feature)
    "vhh_A_patch_max",  # largest contiguous hydrophobic surface patch
    "struct_psh_cdr",  # CDR patches of surface hydrophobicity (TAP-style)
    "vhh_A_arom_exposed",  # solvent-exposed aromatic content
]


def physics_composite(df: pd.DataFrame, cols: list[str] | None = None) -> np.ndarray:
    """Equal-weight mean of per-feature ranks over the available composite columns."""
    cols = cols or PHYSICS_COMPOSITE
    present = [c for c in cols if c in df.columns]
    if not present:
        raise ValueError("None of the physics-composite columns are present.")
    ranks = np.column_stack([rankdata(df[c].to_numpy(dtype=float)) for c in present])
    return ranks.mean(axis=1)
