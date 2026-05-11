"""CLI runner for neurofuzzy similarity experiments."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from neurofuzzy.data_loader import load_dataset, train_test_split
from neurofuzzy.model import NeurofuzzyModel


DATASET_FILES = {
    'mc': Path('datasets/mc.txt'),
    'geresid': Path('datasets/geresid.txt'),
}


def _render_table(dataset_name: str, runs: int, summary: dict[str, tuple[float, float]]) -> str:
    title = 'Neurofuzzy Similarity Results'
    header = f'Dataset: {dataset_name.upper()} ({runs} runs)'
    lines = [
        '╔══════════════════════════════════════════════════╗',
        f'║ {title:<48} ║',
        f'║ {header:<48} ║',
        '╠════════════════╦═════════════╦═════════════╣',
        '║ Metric         ║ Mean        ║ Std         ║',
        '╠════════════════╬═════════════╬═════════════╣',
    ]
    labels = [
        ('pearson', 'Pearson r'),
        ('spearman', 'Spearman ρ'),
        ('cosine_similarity', 'Cosine sim'),
        ('mean_absolute_error', 'MAE'),
        ('root_mean_squared_error', 'RMSE'),
    ]
    for key, label in labels:
        mean, std = summary[key]
        lines.append(f'║ {label:<14} ║ {mean:>10.4f}  ║ {std:>10.4f}  ║')
    lines.append('╚════════════════╩═════════════╩═════════════╝')
    return '\n'.join(lines)


def _run_dataset(dataset: str, runs: int, output_dir: Path, maxiter: int, popsize: int) -> None:
    X, y = load_dataset(DATASET_FILES[dataset], n_features=4)
    rows: list[dict[str, float | int]] = []

    for run_idx in range(runs):
        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            train_ratio=0.6,
            random_state=run_idx,
        )
        model = NeurofuzzyModel(
            n_features=4,
            train_ratio=1.0,
            maxiter=maxiter,
            popsize=popsize,
            random_state=run_idx,
            verbose=False,
        )
        model.fit(X_train, y_train)
        metrics = model.evaluate(X_test, y_test)
        rows.append({'run': run_idx + 1, **metrics})

    df = pd.DataFrame(rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'{dataset}_results.csv'
    df.to_csv(output_file, index=False)

    summary = {
        column: (float(df[column].mean()), float(df[column].std(ddof=1 if len(df) > 1 else 0)))
        for column in [
            'pearson',
            'spearman',
            'cosine_similarity',
            'mean_absolute_error',
            'root_mean_squared_error',
        ]
    }
    print(_render_table(dataset, runs, summary))
    print(f'Saved run-level metrics to {output_file}')


def parse_args() -> argparse.Namespace:
    """Parse command-line options.

    Returns
    -------
    argparse.Namespace
        Parsed CLI arguments.
    """
    parser = argparse.ArgumentParser(description='Run neurofuzzy experiments')
    parser.add_argument(
        '--dataset',
        action='append',
        choices=sorted(DATASET_FILES),
        required=True,
        help='Dataset name. Can be repeated.',
    )
    parser.add_argument(
        '--runs',
        action='append',
        type=int,
        required=False,
        help='Number of repeated runs for corresponding --dataset entry.',
    )
    parser.add_argument('--output', type=Path, default=Path('results'))
    parser.add_argument('--maxiter', type=int, default=200)
    parser.add_argument('--popsize', type=int, default=15)
    return parser.parse_args()


def main() -> None:
    """Execute experiment CLI entrypoint."""
    args = parse_args()

    datasets = args.dataset
    runs_list = args.runs if args.runs is not None else [10]
    if len(runs_list) == 1 and len(datasets) > 1:
        runs_list = runs_list * len(datasets)
    if len(datasets) != len(runs_list):
        raise ValueError('The number of --dataset and --runs entries must match')

    for dataset, runs in zip(datasets, runs_list):
        _run_dataset(dataset, runs, args.output, args.maxiter, args.popsize)


if __name__ == '__main__':
    main()
