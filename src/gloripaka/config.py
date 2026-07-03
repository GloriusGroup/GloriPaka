"""Central configuration: filesystem paths, dataset keys, and shared constants.

Everything a script might otherwise hard-code as a string literal lives here, so
paths and conventions are defined once and imported everywhere. Unlike a research
scratch workspace, this package is **self-contained**: every path resolves relative
to the installed repository root, so a fresh ``git clone`` (plus the data bundle from
Zenodo, see ``README.md``) is fully runnable with no reference to any external tree.

Typical use:
    >>> from gloripaka.config import PATHS, FEATURE_GROUPS
    >>> PATHS.data.exists()  # doctest: +SKIP
    True
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Paths:
    """Resolved absolute paths for the whole package.

    This file lives at ``gloripaka/src/gloripaka/config.py``; counting parents:
    ``[0]`` = the ``gloripaka`` package dir, ``[1]`` = ``src``, ``[2]`` = the
    repository root. All data, results, and figures live under that root, so the
    package never reaches outside its own checkout.
    """

    @property
    def root(self) -> Path:
        """The repository root (the ``git clone`` directory)."""
        return Path(__file__).resolve().parents[2]

    @property
    def workspace(self) -> Path:
        """Alias for :attr:`root`.

        Kept so figure/report code that historically referred to a separate
        "workspace" keeps working after the package was made self-contained.
        """
        return self.root

    # --- data ------------------------------------------------------------- #
    @property
    def data(self) -> Path:
        """The bundled/downloaded data directory (``<root>/data``)."""
        return self.root / "data"

    @property
    def features_dir(self) -> Path:
        """Precomputed feature tables (``<root>/data/features``)."""
        return self.data / "features"

    @property
    def workspace_features_dir(self) -> Path:
        """Second search location for feature tables.

        Identical to :attr:`features_dir` in the standalone package; retained so
        :func:`gloripaka.features.load_feature_table` can search both without a
        code change if a split layout is ever reintroduced.
        """
        return self.data / "features"

    @property
    def test_csv(self) -> Path:
        """The 81 challenge VHH (space-delimited: ``ID VHHsequence``)."""
        return self.data / "official_test_set.csv"

    @property
    def gdpa1_csv(self) -> Path:
        """GDPa1 (clinical mAb with HIC + predefined CV folds)."""
        return self.data / "training" / "GDPa1_v1.2_20250814.csv"

    @property
    def vhh_harvest_csv(self) -> Path:
        """Native-VHH HIC harvest (Arras / GDPa2 / TNP / Siegmund)."""
        return self.data / "hic_harvest" / "hic_records_vhh.csv"

    # --- generated outputs ------------------------------------------------ #
    @property
    def results_dir(self) -> Path:
        """Where experiment result tables are written (``<root>/results``)."""
        return self.root / "results"

    @property
    def figures_dir(self) -> Path:
        """Where figures are written (``<root>/figures``)."""
        return self.root / "figures"

    @property
    def cache_dir(self) -> Path:
        """Scratch cache for intermediate artifacts (``<root>/cache``)."""
        return self.root / "cache"

    def ensure_outputs(self) -> None:
        """Create the output directories if missing (safe to call repeatedly)."""
        for directory in (self.results_dir, self.figures_dir, self.cache_dir):
            directory.mkdir(parents=True, exist_ok=True)


PATHS = Paths()


# --- feature-group prefixes in the ~613-feature tables ----------------------- #
# Used to slice the feature matrix into interpretable blocks for ablations.
FEATURE_GROUPS: dict[str, str] = {
    "desc": "desc_",  # 29  whole-domain composition / physico-chemical descriptors
    "cdr": "cdr_",  # 27  per-CDR length / hydrophobicity / aromatic stats (IMGT)
    "deepsp": "deepsp_",  # 12  DeepSP SAP/SCM CNN surface-patch heads
    "struct": "struct_",  # 58  NanoBodyBuilder2 + FreeSASA structural descriptors
    "esm": "esm_",  # 480 ESM2 (esm2_t12_35M_UR50D) mean-pooled embedding
    "vhh": "vhh_",  # 7   VHH-specific (positive-SASA, aromatic-exposed, patch, pI, CDR3, FR2)
}

# The single most transferable engineered feature (Wimley-White positive-only
# hydrophobic SASA). See the manuscript / lab journal for why the lysine-heavy
# "physics4" descriptor does NOT transfer and this one does.
TRANSFERABLE_PHYSICS_FEATURE: str = "vhh_A_positive_sasa"

# Scale groups (assay identities). Absolute HIC minutes are NOT comparable across
# assays; only within-group rank is. ``SCALE_TEST`` is the challenge target scale.
SCALE_TEST: str = "ARRAS_TSKButylNPR_MeOH"

# Amino-acid alphabet used throughout.
AA: str = "ACDEFGHIKLMNPQRSTVWY"
