"""Label loading and the honest known/unknown split of the 81 challenge VHH.

The ARRAS test labels for 34 of the 81 VHH are public (they are Merck's own Arras
2023 data). Strategy B uses those 34 **only as a held-out validation proxy**. This
module is the single source of truth for:

  * the 81 test sequences,
  * which 34 have known HIC and their values (by exact sequence match), and
  * the GDPa1 mAb training labels + predefined CV folds.

Keeping the join logic here (not in experiment scripts) guarantees every experiment
uses the exact same definition of "known" — no accidental drift into leakage.
"""

from __future__ import annotations

import pandas as pd

from gloripaka.config import PATHS, SCALE_TEST


def _clean_seq(s: str) -> str:
    return str(s).strip().upper()


def load_test_sequences() -> pd.DataFrame:
    """The 81 challenge VHH as ``[ID, sequence]`` (file is space-delimited)."""
    rows: list[tuple[int, str]] = []
    with open(PATHS.test_csv) as fh:
        next(fh)  # header: "ID VHHsequence"
        for line in fh:
            line = line.strip()
            if not line:
                continue
            sid, seq = line.split(None, 1)
            rows.append((int(sid), _clean_seq(seq)))
    df = pd.DataFrame(rows, columns=["ID", "sequence"])
    df["Name"] = "VHH_" + df["ID"].astype(str)  # matches the feature-table Name key
    return df


def load_vhh_harvest() -> pd.DataFrame:
    """Native-VHH HIC records (Arras / GDPa2 / TNP / Siegmund), sequences cleaned."""
    df = pd.read_csv(PATHS.vhh_harvest_csv)
    df["sequence"] = df["sequence"].map(_clean_seq)
    return df


def known_test_labels() -> pd.DataFrame:
    """The subset of the 81 test VHH with public HIC, by **exact** sequence match.

    Returns ``[ID, Name, sequence, hic, scale_group, source_key]`` for the 34 ARRAS
    VHH. Exact match (not fuzzy) is deliberate: engineered point-mutant variants of a
    known VHH are treated as *unknown*, so they cannot leak a near-label.
    """
    test = load_test_sequences()
    harvest = load_vhh_harvest()
    arras = harvest[harvest["scale_group"] == SCALE_TEST]

    merged = test.merge(
        arras[["sequence", "hic_value", "scale_group", "source_key"]].drop_duplicates("sequence"),
        on="sequence",
        how="inner",
    )
    merged = merged.rename(columns={"hic_value": "hic"})
    return merged[["ID", "Name", "sequence", "hic", "scale_group", "source_key"]].reset_index(
        drop=True
    )


def load_gdpa1_labels() -> pd.DataFrame:
    """GDPa1 mAb training labels: ``[Name, hic, fold, vh_sequence]`` (HIC present only)."""
    df = pd.read_csv(PATHS.gdpa1_csv)
    df = df[df["HIC"].notna() & df["vh_protein_sequence"].notna()].copy()
    out = pd.DataFrame(
        {
            "Name": df["antibody_name"].astype(str),
            "hic": df["HIC"].astype(float),
            "fold": df["hierarchical_cluster_fold"].fillna(-1).astype(int),
            "vh_sequence": df["vh_protein_sequence"].map(_clean_seq),
        }
    )
    return out.reset_index(drop=True)
