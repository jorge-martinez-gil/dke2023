"""Tests for the benchmark runner and aggregation."""

import numpy as np
import pandas as pd

from neurofuzzy.benchmark import (
    BenchmarkConfig,
    _bootstrap_ci,
    _significance,
    aggregate,
    run,
)

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]


def test_bootstrap_ci_orders_and_brackets_mean():
    vals = np.array([0.4, 0.5, 0.6, 0.55, 0.45])
    lo, hi = _bootstrap_ci(vals, n_boot=2000, seed=0)
    assert lo <= vals.mean() <= hi
    assert _bootstrap_ci(np.array([0.5])) == (0.5, 0.5)


def test_run_and_aggregate_end_to_end(tmp_path):
    config = BenchmarkConfig(
        methods=["mean-of-features", "best-single-feature", "linear-regression"],
        seeds=[0, 1, 2],
    )
    completed, remaining = run(config, ["mc"], tmp_path, ROOT / "datasets")
    assert remaining == 0
    assert completed == 9  # 3 methods x 3 seeds

    # Resumable: a second run does nothing new.
    completed2, remaining2 = run(config, ["mc"], tmp_path, ROOT / "datasets")
    assert completed2 == 0 and remaining2 == 0

    summary = aggregate(tmp_path)
    assert {"mean", "std", "ci95_low", "ci95_high", "n"}.issubset(summary.columns)
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "results_table.tex").exists()
    assert (tmp_path / "manifest.json").exists()


def test_significance_columns():
    # Build a minimal raw frame with neurofuzzy + one baseline across 4 seeds.
    rows = []
    for seed in range(4):
        rows.append(
            {
                "dataset": "mc",
                "method": "neurofuzzy",
                "seed": seed,
                "split": "test",
                "pearson": 0.7 + 0.01 * seed,
                "spearman": 0.6,
            }
        )
        rows.append(
            {
                "dataset": "mc",
                "method": "mean-of-features",
                "seed": seed,
                "split": "test",
                "pearson": 0.8 + 0.01 * seed,
                "spearman": 0.65,
            }
        )
    sig = _significance(pd.DataFrame(rows))
    assert {"dataset", "metric", "comparison", "p_value"}.issubset(sig.columns)
    assert (sig["comparison"] == "neurofuzzy vs mean-of-features").any()
