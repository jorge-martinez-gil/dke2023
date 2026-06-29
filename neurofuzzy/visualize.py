"""Publication-quality figure generation for neurofuzzy benchmarks.

All figures are produced deterministically from the benchmark artifacts written
by :mod:`neurofuzzy.benchmark` (plus, for diagnostic scatter/membership plots, a
single re-fit under the documented configuration). Run as a module::

    python -m neurofuzzy.visualize --results benchmarks/results --out benchmarks/figures

A non-interactive backend is used so the script runs headless in CI.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Colorblind-friendly qualitative palette (Wong, 2011).
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]

plt.rcParams.update(
    {
        "figure.dpi": 200,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.axisbelow": True,
    }
)


def _load_summary(results_dir: Path) -> pd.DataFrame:
    path = Path(results_dir) / "summary.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}; run `neurofuzzy.benchmark --aggregate` first.")
    return pd.read_csv(path)


def plot_metric_comparison(
    summary: pd.DataFrame, out_path: Path, metric: str = "pearson", split: str = "test"
) -> Path:
    """Grouped bar chart of one metric per method per dataset, with 95% CIs."""
    data = summary[(summary["metric"] == metric) & (summary["split"] == split)]
    datasets = sorted(data["dataset"].unique())
    methods = sorted(data["method"].unique())

    fig, ax = plt.subplots(figsize=(1.6 * len(methods) + 2, 4.2))
    width = 0.8 / max(len(methods), 1)
    x = np.arange(len(datasets))
    for i, method in enumerate(methods):
        means, lo, hi = [], [], []
        for ds in datasets:
            row = data[(data["dataset"] == ds) & (data["method"] == method)]
            m = float(row["mean"].iloc[0]) if not row.empty else np.nan
            means.append(m)
            lo.append(m - float(row["ci95_low"].iloc[0]) if not row.empty else 0.0)
            hi.append(float(row["ci95_high"].iloc[0]) - m if not row.empty else 0.0)
        ax.bar(
            x + i * width,
            means,
            width,
            yerr=[lo, hi],
            capsize=4,
            label=method,
            color=PALETTE[i % len(PALETTE)],
            edgecolor="white",
        )
    ax.set_xticks(x + width * (len(methods) - 1) / 2)
    ax.set_xticklabels([d.upper() for d in datasets])
    ax.set_ylabel(f"{metric} ({split})")
    ax.set_title(f"Held-out {metric} by method (mean ± 95% CI)")
    ax.legend(frameon=False, fontsize=9, ncol=2)
    fig.tight_layout()
    out_path = Path(out_path)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_train_test_gap(summary: pd.DataFrame, out_path: Path, method: str = "neurofuzzy") -> Path:
    """Bar chart contrasting train vs test Pearson to expose overfitting."""
    data = summary[(summary["method"] == method) & (summary["metric"] == "pearson")]
    datasets = sorted(data["dataset"].unique())
    fig, ax = plt.subplots(figsize=(1.6 * len(datasets) + 2, 4.2))
    x = np.arange(len(datasets))
    for i, split in enumerate(("train", "test")):
        means = [
            float(data[(data["dataset"] == ds) & (data["split"] == split)]["mean"].iloc[0])
            for ds in datasets
        ]
        ax.bar(x + i * 0.4, means, 0.4, label=split, color=PALETTE[i], edgecolor="white")
    ax.set_xticks(x + 0.2)
    ax.set_xticklabels([d.upper() for d in datasets])
    ax.set_ylabel("Pearson r")
    ax.set_title(f"Train vs test Pearson ({method}) — generalization gap")
    ax.legend(frameon=False)
    fig.tight_layout()
    out_path = Path(out_path)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_prediction_scatter(
    dataset: str, out_path: Path, data_dir: Path, maxiter: int, popsize: int, seed: int = 0
) -> Path:
    """Scatter of predicted vs. true similarity for a single re-fit model."""
    from .baselines import LinearRegressionBaseline
    from .benchmark import DATASET_FILES
    from .data_loader import load_dataset, train_test_split
    from .model import NeurofuzzyModel

    X, y = load_dataset(Path(data_dir) / DATASET_FILES[dataset], n_features=4)
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, train_ratio=0.6, random_state=seed)
    model = NeurofuzzyModel(
        train_ratio=1.0, maxiter=maxiter, popsize=popsize, random_state=seed, verbose=False
    ).fit(X_tr, y_tr)
    base = LinearRegressionBaseline().fit(X_tr, y_tr)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot([0, 1], [0, 1], color="gray", ls="--", lw=1, label="ideal")
    ax.scatter(y_te, model.predict(X_te), color=PALETTE[0], s=50, label="neurofuzzy", zorder=3)
    ax.scatter(
        y_te,
        base.predict(X_te),
        color=PALETTE[1],
        s=50,
        marker="^",
        label="linear-regression",
        zorder=3,
        alpha=0.8,
    )
    ax.set_xlabel("True similarity")
    ax.set_ylabel("Predicted similarity")
    ax.set_title(f"{dataset.upper()} test predictions (seed {seed})")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False)
    fig.tight_layout()
    out_path = Path(out_path)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_membership_functions(controller, out_path: Path) -> Path:
    """Visualize the learned input membership functions and output sets."""
    grid = np.linspace(0.0, 1.0, 201)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    centers = controller._centers
    widths = controller._widths
    from .fuzzy_controller import trimf

    for j in range(centers.shape[1]):
        # Show membership terms of the first input as a representative example.
        c, w = centers[0, j], widths[0, j]
        axes[0].plot(
            grid,
            trimf(grid, c - w, c, c + w),
            color=PALETTE[j % len(PALETTE)],
            label=f"term {j + 1}",
        )
    axes[0].set_title("Input #1 membership functions")
    axes[0].set_xlabel("feature value")
    axes[0].set_ylabel("membership")
    axes[0].legend(frameon=False, fontsize=9)

    low, mid, high = controller._output_sets()
    for arr, name, col in zip((low, mid, high), ("low", "mid", "high"), PALETTE):
        axes[1].plot(grid, arr, color=col, label=name)
    axes[1].set_title("Output fuzzy sets")
    axes[1].set_xlabel("similarity score")
    axes[1].set_ylabel("membership")
    axes[1].legend(frameon=False, fontsize=9)
    fig.suptitle("Learned fuzzy sets")
    fig.tight_layout()
    out_path = Path(out_path)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def generate_all(
    results_dir: Path, out_dir: Path, data_dir: Path, with_fit: bool = True
) -> list[Path]:
    """Generate every benchmark figure and return the written paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = _load_summary(results_dir)
    manifest = Path(results_dir) / "manifest.json"
    maxiter, popsize = 40, 8
    if manifest.exists():
        import json

        cfg = json.loads(manifest.read_text()).get("config", {})
        maxiter = int(cfg.get("maxiter", maxiter))
        popsize = int(cfg.get("popsize", popsize))

    written = [
        plot_metric_comparison(summary, out_dir / "pearson_comparison.png", "pearson"),
        plot_metric_comparison(summary, out_dir / "spearman_comparison.png", "spearman"),
        plot_train_test_gap(summary, out_dir / "train_test_gap.png"),
    ]
    if with_fit:
        from .benchmark import DATASET_FILES
        from .data_loader import load_dataset, train_test_split
        from .model import NeurofuzzyModel

        for ds in sorted(summary["dataset"].unique()):
            written.append(
                plot_prediction_scatter(
                    ds, out_dir / f"{ds}_scatter.png", data_dir, maxiter, popsize
                )
            )
        # Membership functions from a representative fit on the first dataset.
        ds0 = sorted(summary["dataset"].unique())[0]
        X, y = load_dataset(Path(data_dir) / DATASET_FILES[ds0], n_features=4)
        X_tr, _, y_tr, _ = train_test_split(X, y, train_ratio=0.6, random_state=0)
        model = NeurofuzzyModel(
            train_ratio=1.0, maxiter=maxiter, popsize=popsize, random_state=0, verbose=False
        ).fit(X_tr, y_tr)
        written.append(plot_membership_functions(model.controller_, out_dir / "fuzzy_sets.png"))
    return written


def main(argv: list[str] | None = None) -> None:
    """Command-line entry point for figure generation."""
    parser = argparse.ArgumentParser(
        prog="neurofuzzy-figures", description="Generate benchmark figures."
    )
    parser.add_argument("--results", type=Path, default=Path("benchmarks/results"))
    parser.add_argument("--out", type=Path, default=Path("benchmarks/figures"))
    parser.add_argument("--data-dir", type=Path, default=Path("datasets"))
    parser.add_argument(
        "--no-fit", action="store_true", help="Skip figures that require re-fitting."
    )
    args = parser.parse_args(argv)
    paths = generate_all(args.results, args.out, args.data_dir, with_fit=not args.no_fit)
    print(f"Wrote {len(paths)} figure(s) to {args.out}/:")
    for p in paths:
        print(f"  - {p.name}")


if __name__ == "__main__":
    main()
