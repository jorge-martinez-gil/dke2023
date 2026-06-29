"""Unified, resumable benchmark runner for neurofuzzy semantic similarity.

The runner evaluates the :class:`~neurofuzzy.model.NeurofuzzyModel` and a set of
reference baselines across multiple datasets and random seeds under an identical
train/test protocol, then aggregates per-seed metrics into means, standard
deviations, bootstrap confidence intervals, and paired significance tests.

Design goals
------------
* **One command** runs everything (``--run`` then ``--aggregate``, or ``--all``).
* **Resumable**: results are appended row-by-row to ``raw_results.csv`` and
  already-computed ``(dataset, method, seed)`` combinations are skipped, so a
  long benchmark can be completed across several invocations using
  ``--time-budget``.
* **Traceable**: every run writes a ``manifest.json`` capturing package
  versions, the git commit, the exact configuration, and timestamps.
* **No fabricated numbers**: all metrics are computed from data at run time.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .baselines import BASELINES, get_baseline
from .data_loader import load_dataset, train_test_split
from .metrics import evaluate_all
from .model import NeurofuzzyModel

#: Built-in datasets shipped with the repository (paths relative to data dir).
DATASET_FILES: dict[str, str] = {
    "mc": "mc.txt",
    "geresid": "geresid.txt",
}

#: Name used for the proposed method in reports and figures.
NEUROFUZZY = "neurofuzzy"

METRIC_KEYS = [
    "pearson",
    "spearman",
    "cosine_similarity",
    "mean_absolute_error",
    "root_mean_squared_error",
]


@dataclass
class BenchmarkConfig:
    """Configuration for a benchmark run."""

    maxiter: int = 40
    popsize: int = 8
    train_ratio: float = 0.6
    n_features: int = 4
    methods: list[str] = field(default_factory=lambda: [NEUROFUZZY, *BASELINES.keys()])
    seeds: list[int] = field(default_factory=lambda: list(range(10)))


def _git_hash() -> str:
    """Return the current git commit hash, or ``"unknown"`` if unavailable."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip()
    except Exception:  # pragma: no cover - environment dependent
        return "unknown"


def _environment() -> dict[str, str]:
    """Capture a minimal, reproducibility-relevant snapshot of the environment."""
    import scipy

    return {
        "neurofuzzy": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "git_commit": _git_hash(),
    }


def _evaluate_method(
    method: str,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    seed: int,
    config: BenchmarkConfig,
) -> tuple[dict[str, float], dict[str, float]]:
    """Fit one method on the training split and score it on train and test.

    Returns
    -------
    tuple of dict
        ``(train_metrics, test_metrics)`` dictionaries.
    """
    if method == NEUROFUZZY:
        model = NeurofuzzyModel(
            n_features=config.n_features,
            train_ratio=1.0,
            maxiter=config.maxiter,
            popsize=config.popsize,
            random_state=seed,
            verbose=False,
        ).fit(X_train, y_train)
        train_pred = model.predict(X_train)
        test_pred = model.predict(X_test)
    else:
        baseline = get_baseline(method).fit(X_train, y_train)
        train_pred = baseline.predict(X_train)
        test_pred = baseline.predict(X_test)

    return (
        evaluate_all(y_train, train_pred),
        evaluate_all(y_test, test_pred),
    )


def _load_raw(raw_csv: Path) -> pd.DataFrame:
    if raw_csv.exists():
        return pd.read_csv(raw_csv)
    return pd.DataFrame()


def run(
    config: BenchmarkConfig,
    datasets: list[str],
    out_dir: Path,
    data_dir: Path,
    time_budget: float | None = None,
) -> tuple[int, int]:
    """Run (or resume) the benchmark, appending per-seed rows to disk.

    Parameters
    ----------
    config : BenchmarkConfig
        Methods, seeds, and optimizer settings.
    datasets : list of str
        Dataset names (keys of :data:`DATASET_FILES`).
    out_dir : pathlib.Path
        Directory where ``raw_results.csv`` and ``manifest.json`` are written.
    data_dir : pathlib.Path
        Directory containing the dataset text files.
    time_budget : float or None
        Soft wall-clock budget in seconds. New fits are not started once the
        budget is exceeded, enabling chunked/resumable execution.

    Returns
    -------
    tuple of int
        ``(n_completed_this_call, n_remaining)``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_csv = out_dir / "raw_results.csv"
    existing = _load_raw(raw_csv)
    done = set()
    if not existing.empty:
        done = {(r.dataset, r.method, int(r.seed)) for r in existing.itertuples(index=False)}

    todo = [
        (ds, method, seed)
        for ds in datasets
        for seed in config.seeds
        for method in config.methods
        if (ds, method, seed) not in done
    ]

    start = time.time()
    completed = 0
    cache: dict[tuple[str, int], tuple] = {}
    for ds, method, seed in todo:
        if time_budget is not None and (time.time() - start) > time_budget:
            break
        key = (ds, seed)
        if key not in cache:
            X, y = load_dataset(data_dir / DATASET_FILES[ds], n_features=config.n_features)
            X_tr, X_te, y_tr, y_te = train_test_split(
                X, y, train_ratio=config.train_ratio, random_state=seed
            )
            cache[key] = (X_tr, X_te, y_tr, y_te)
        X_tr, X_te, y_tr, y_te = cache[key]

        t0 = time.time()
        train_metrics, test_metrics = _evaluate_method(method, X_tr, y_tr, X_te, y_te, seed, config)
        fit_seconds = time.time() - t0

        rows = []
        for split, metrics in (("train", train_metrics), ("test", test_metrics)):
            row = {
                "dataset": ds,
                "method": method,
                "seed": seed,
                "split": split,
                "n_train": int(X_tr.shape[0]),
                "n_test": int(X_te.shape[0]),
                "fit_seconds": round(fit_seconds, 4),
            }
            row.update({k: metrics[k] for k in METRIC_KEYS})
            rows.append(row)
        pd.DataFrame(rows).to_csv(raw_csv, mode="a", header=not raw_csv.exists(), index=False)
        completed += 1

    total = len({(ds, m, s) for ds in datasets for s in config.seeds for m in config.methods})
    refreshed = _load_raw(raw_csv)
    if refreshed.empty:
        done_now = 0
    else:
        done_now = refreshed[["dataset", "method", "seed"]].drop_duplicates().shape[0]
    remaining = max(total - done_now, 0)

    manifest = {
        "environment": _environment(),
        "config": asdict(config),
        "datasets": datasets,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "completed_combinations": done_now,
        "total_combinations": total,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return completed, remaining


def _bootstrap_ci(values: np.ndarray, n_boot: int = 10000, seed: int = 0) -> tuple[float, float]:
    """Percentile bootstrap 95% confidence interval for the mean."""
    values = np.asarray(values, dtype=float)
    if values.size <= 1:
        v = float(values[0]) if values.size else float("nan")
        return v, v
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, values.size, size=(n_boot, values.size))
    means = values[idx].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def aggregate(out_dir: Path) -> pd.DataFrame:
    """Aggregate raw per-seed results into summary statistics and reports."""
    raw_csv = out_dir / "raw_results.csv"
    if not raw_csv.exists():
        raise FileNotFoundError(f"No raw results at {raw_csv}; run the benchmark first.")
    df = pd.read_csv(raw_csv)

    records = []
    for (ds, method, split), group in df.groupby(["dataset", "method", "split"]):
        for metric in METRIC_KEYS:
            vals = group[metric].to_numpy(dtype=float)
            lo, hi = _bootstrap_ci(vals)
            records.append(
                {
                    "dataset": ds,
                    "method": method,
                    "split": split,
                    "metric": metric,
                    "n": int(vals.size),
                    "mean": float(np.mean(vals)),
                    "std": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
                    "ci95_low": lo,
                    "ci95_high": hi,
                }
            )
    summary = pd.DataFrame.from_records(records)
    summary.to_csv(out_dir / "summary.csv", index=False)
    (out_dir / "summary.json").write_text(summary.to_json(orient="records", indent=2))

    significance = _significance(df)
    if not significance.empty:
        significance.to_csv(out_dir / "significance.csv", index=False)

    _write_markdown_report(out_dir, summary, significance, df)
    _write_latex_table(out_dir, summary)
    return summary


def _significance(df: pd.DataFrame) -> pd.DataFrame:
    """Paired Wilcoxon signed-rank tests: neurofuzzy vs each baseline (test split)."""
    from scipy.stats import wilcoxon

    rows = []
    test = df[df["split"] == "test"]
    for ds in sorted(test["dataset"].unique()):
        sub = test[test["dataset"] == ds]
        if NEUROFUZZY not in sub["method"].values:
            continue
        nf = sub[sub["method"] == NEUROFUZZY].set_index("seed")
        for method in sorted(sub["method"].unique()):
            if method == NEUROFUZZY:
                continue
            other = sub[sub["method"] == method].set_index("seed")
            common = sorted(set(nf.index) & set(other.index))
            if len(common) < 3:
                continue
            for metric in ("pearson", "spearman"):
                a = nf.loc[common, metric].to_numpy(dtype=float)
                b = other.loc[common, metric].to_numpy(dtype=float)
                try:
                    stat, p = wilcoxon(a, b)
                    p = float(p)
                except ValueError:
                    stat, p = float("nan"), float("nan")
                rows.append(
                    {
                        "dataset": ds,
                        "metric": metric,
                        "comparison": f"{NEUROFUZZY} vs {method}",
                        "n_pairs": len(common),
                        "mean_diff": float(np.mean(a - b)),
                        "wilcoxon_stat": float(stat),
                        "p_value": p,
                    }
                )
    return pd.DataFrame.from_records(rows)


def _fmt(mean: float, lo: float, hi: float) -> str:
    return f"{mean:.4f} [{lo:.4f}, {hi:.4f}]"


def _write_markdown_report(out_dir, summary, significance, raw) -> None:
    """Write a human-readable Markdown benchmark report."""
    lines = ["# Neurofuzzy Benchmark Report", ""]
    manifest_path = out_dir / "manifest.json"
    if manifest_path.exists():
        man = json.loads(manifest_path.read_text())
        env = man.get("environment", {})
        cfg = man.get("config", {})
        lines += [
            f"- Generated: {man.get('timestamp_utc', 'n/a')}",
            f"- neurofuzzy {env.get('neurofuzzy', '?')} | numpy {env.get('numpy', '?')} | "
            f"scipy {env.get('scipy', '?')} | python {env.get('python', '?')}",
            f"- git commit: `{env.get('git_commit', 'unknown')}`",
            f"- optimizer: maxiter={cfg.get('maxiter')}, popsize={cfg.get('popsize')}, "
            f"train_ratio={cfg.get('train_ratio')}",
            f"- seeds: {len(cfg.get('seeds', []))} | "
            f"datasets: {', '.join(man.get('datasets', []))}",
            "",
            "Values are **mean [95% bootstrap CI]** over seeds on the held-out test split.",
            "",
        ]
    for ds in sorted(summary["dataset"].unique()):
        lines += [f"## Dataset: {ds.upper()}", ""]
        test = summary[(summary["dataset"] == ds) & (summary["split"] == "test")]
        methods = sorted(test["method"].unique())
        lines.append("| Method | " + " | ".join(METRIC_KEYS) + " |")
        lines.append("|" + "---|" * (len(METRIC_KEYS) + 1))
        for method in methods:
            cells = [method]
            for metric in METRIC_KEYS:
                row = test[(test["method"] == method) & (test["metric"] == metric)]
                if row.empty:
                    cells.append("-")
                else:
                    r = row.iloc[0]
                    cells.append(_fmt(r["mean"], r["ci95_low"], r["ci95_high"]))
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
        if not significance.empty:
            sig = significance[significance["dataset"] == ds]
            if not sig.empty:
                lines += ["**Paired Wilcoxon signed-rank vs baselines (test split):**", ""]
                lines.append("| Comparison | Metric | n | mean diff | p-value |")
                lines.append("|---|---|---|---|---|")
                for r in sig.itertuples(index=False):
                    lines.append(
                        f"| {r.comparison} | {r.metric} | {r.n_pairs} | "
                        f"{r.mean_diff:+.4f} | {r.p_value:.4f} |"
                    )
                lines.append("")
    lines += [
        "> Results are produced by `python -m neurofuzzy.benchmark`. "
        "They depend on the optimizer configuration above and the small benchmark "
        "sample sizes; report the configuration alongside any cited number.",
        "",
    ]
    (out_dir / "report.md").write_text("\n".join(lines))


def _write_latex_table(out_dir, summary) -> None:
    """Write a publication-ready LaTeX (booktabs) table of test-split results."""
    test = summary[summary["split"] == "test"]
    lines = [
        "% Auto-generated by neurofuzzy.benchmark -- do not edit by hand.",
        "\\begin{table}[t]",
        "\\centering",
        "\\caption{Semantic similarity benchmark (mean over seeds on the test "
        "split; 95\\% bootstrap CI in brackets).}",
        "\\label{tab:neurofuzzy-benchmark}",
        "\\begin{tabular}{llrr}",
        "\\toprule",
        "Dataset & Method & Pearson $r$ & Spearman $\\rho$ \\\\",
        "\\midrule",
    ]
    for ds in sorted(test["dataset"].unique()):
        for method in sorted(test[test["dataset"] == ds]["method"].unique()):

            def cell(metric, ds=ds, method=method):
                row = test[
                    (test["dataset"] == ds)
                    & (test["method"] == method)
                    & (test["metric"] == metric)
                ]
                if row.empty:
                    return "--"
                r = row.iloc[0]
                return f"{r['mean']:.3f} [{r['ci95_low']:.3f}, {r['ci95_high']:.3f}]"

            lines.append(f"{ds.upper()} & {method} & {cell('pearson')} & {cell('spearman')} \\\\")
        lines.append("\\midrule")
    if lines[-1] == "\\midrule":
        lines.pop()
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    (out_dir / "results_table.tex").write_text("\n".join(lines))


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="neurofuzzy-benchmark",
        description="Resumable benchmark runner for neurofuzzy semantic similarity.",
    )
    parser.add_argument("--run", action="store_true", help="Run/resume fitting and scoring.")
    parser.add_argument("--aggregate", action="store_true", help="Aggregate results into reports.")
    parser.add_argument("--all", action="store_true", help="Run then aggregate.")
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=sorted(DATASET_FILES),
        choices=sorted(DATASET_FILES),
        help="Datasets to evaluate.",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=[NEUROFUZZY, *BASELINES.keys()],
        help="Methods to evaluate (neurofuzzy and/or baselines).",
    )
    parser.add_argument("--seeds", type=int, default=10, help="Number of seeds (0..N-1).")
    parser.add_argument("--maxiter", type=int, default=40, help="DE iterations.")
    parser.add_argument("--popsize", type=int, default=8, help="DE population multiplier.")
    parser.add_argument("--train-ratio", type=float, default=0.6, help="Train split fraction.")
    parser.add_argument(
        "--time-budget",
        type=float,
        default=None,
        help="Soft wall-clock budget (s) for resumable chunking.",
    )
    parser.add_argument("--out", type=Path, default=Path("results"), help="Output directory.")
    parser.add_argument(
        "--data-dir", type=Path, default=Path("datasets"), help="Dataset directory."
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Command-line entry point for the benchmark runner."""
    args = _parse_args(argv)
    do_run = args.run or args.all or not args.aggregate
    do_aggregate = args.aggregate or args.all

    config = BenchmarkConfig(
        maxiter=args.maxiter,
        popsize=args.popsize,
        train_ratio=args.train_ratio,
        methods=list(args.methods),
        seeds=list(range(args.seeds)),
    )

    if do_run:
        completed, remaining = run(
            config, list(args.datasets), args.out, args.data_dir, args.time_budget
        )
        print(f"[run] completed {completed} new combination(s); {remaining} remaining.")
        if remaining and not args.aggregate:
            print("[run] re-invoke with the same arguments to continue (resumable).")

    if do_aggregate:
        summary = aggregate(args.out)
        print(
            f"[aggregate] wrote summary for {summary['dataset'].nunique()} dataset(s) "
            f"to {args.out}/ (summary.csv, report.md, results_table.tex)."
        )


if __name__ == "__main__":
    main()
