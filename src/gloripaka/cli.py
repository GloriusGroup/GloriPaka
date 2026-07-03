"""Command-line interface for GloriPaka.

Examples:
    Reproduce the challenge predictions and score the held-out proxy::

        gloripaka reproduce

    Predict from a precomputed feature table::

        gloripaka predict --features my_features.csv --out predictions.csv

    Predict from raw sequences (requires the ``featurize`` extra)::

        gloripaka predict --sequences my_sequences.csv --out predictions.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from gloripaka import __version__


def _read_indexed(path: str) -> pd.DataFrame:
    """Read a CSV and index it by its ``Name`` column when present."""
    frame = pd.read_csv(path)
    if "Name" in frame.columns:
        frame["Name"] = frame["Name"].astype(str)
        frame = frame.set_index("Name")
    return frame


def _write(frame: pd.DataFrame, out: str | None) -> None:
    """Write predictions to ``out`` or print a preview to stdout."""
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out)
        print(f"Wrote {len(frame)} predictions to {out}")
    else:
        print(frame.to_string())


def _cmd_predict(args: argparse.Namespace) -> int:
    """Run the ``predict`` subcommand."""
    from gloripaka import GloriPakaPredictor

    predictor = GloriPakaPredictor()
    if args.sequences:
        table = _read_indexed(args.sequences)
        column = "sequence" if "sequence" in table.columns else table.columns[0]
        try:
            result = predictor.predict(list(table[column].astype(str)))
        except RuntimeError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
        _write(result, args.out)
        return 0
    if args.features:
        features = _read_indexed(args.features)
        _write(predictor.predict_from_features(features), args.out)
        return 0
    print("error: provide --features or --sequences", file=sys.stderr)
    return 2


def _cmd_reproduce(args: argparse.Namespace) -> int:
    """Run the ``reproduce`` subcommand: predict the 81 VHH and score the proxy."""
    from gloripaka.pipeline import evaluate_known_proxy, predict_challenge_test_set

    preds = predict_challenge_test_set()
    _write(preds, args.out or "results/predictions_81.csv")

    metrics, (sp, sp_lo, sp_hi), (pe, pe_lo, pe_hi) = evaluate_known_proxy()
    print(
        "\nHeld-out proxy (34 VHH with public HIC labels; never trained on):\n"
        f"  Spearman {sp:.3f} (95% CI {sp_lo:.2f}..{sp_hi:.2f})\n"
        f"  Pearson  {pe:.3f} (95% CI {pe_lo:.2f}..{pe_hi:.2f})\n"
        f"  n = {metrics.n}"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point for the ``gloripaka`` console script.

    Args:
        argv: Optional argument list (defaults to ``sys.argv``).

    Returns:
        A process exit code.
    """
    parser = argparse.ArgumentParser(
        prog="gloripaka",
        description="Sequence-based HIC-retention prediction for VHH nanobodies.",
    )
    parser.add_argument("--version", action="version", version=f"gloripaka {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    predict = sub.add_parser("predict", help="Predict HIC from features or sequences.")
    predict.add_argument("--features", help="CSV of precomputed features (indexed by Name).")
    predict.add_argument("--sequences", help="CSV with a sequence column (needs the featurize extra).")
    predict.add_argument("--out", help="Output CSV path (default: print to stdout).")
    predict.set_defaults(func=_cmd_predict)

    reproduce = sub.add_parser(
        "reproduce", help="Predict the 81 challenge VHH and score the held-out proxy."
    )
    reproduce.add_argument("--out", help="Output CSV path for the 81 predictions.")
    reproduce.set_defaults(func=_cmd_reproduce)

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
