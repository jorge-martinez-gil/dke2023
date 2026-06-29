"""Backward-compatible thin wrapper around the unified benchmark runner.

The full, maintained entry point is ``python -m neurofuzzy.benchmark`` (see
``docs/benchmark.md``). This script preserves the historical invocation
``python experiments/run_experiment.py --dataset mc --runs 10`` by translating it
to the benchmark CLI, then running and aggregating.

Examples
--------
    python experiments/run_experiment.py --dataset mc --runs 10
    python experiments/run_experiment.py --dataset mc --dataset geresid --runs 10
"""

from __future__ import annotations

import argparse
from pathlib import Path

from neurofuzzy.benchmark import main as benchmark_main


def parse_args() -> argparse.Namespace:
    """Parse the legacy command-line options."""
    parser = argparse.ArgumentParser(description="Run neurofuzzy experiments (compatibility shim).")
    parser.add_argument(
        "--dataset",
        action="append",
        choices=["mc", "geresid"],
        required=True,
        help="Dataset name. Can be repeated.",
    )
    parser.add_argument("--runs", type=int, default=10, help="Number of seeds (0..N-1).")
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--maxiter", type=int, default=40)
    parser.add_argument("--popsize", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    """Translate legacy options to the benchmark CLI and run it."""
    args = parse_args()
    argv = [
        "--all",
        "--datasets",
        *args.dataset,
        "--seeds",
        str(args.runs),
        "--maxiter",
        str(args.maxiter),
        "--popsize",
        str(args.popsize),
        "--out",
        str(args.output),
    ]
    benchmark_main(argv)


if __name__ == "__main__":
    main()
