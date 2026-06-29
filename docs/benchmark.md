# Benchmark protocol

The benchmark harness (`neurofuzzy.benchmark`) evaluates the proposed model and
all baselines under one protocol so results are directly comparable.

## Protocol

1. **Split.** For each random seed, a 60/40 train/test split is drawn
   (`train_test_split`, seeded).
2. **Fit.** Each method is fit on the training split only. The fuzzy controller
   is trained by Differential Evolution maximizing the training-set Pearson
   correlation.
3. **Score.** Each method is scored on both the train and held-out test splits
   with Pearson, Spearman, cosine similarity, MAE, and RMSE.
4. **Repeat & aggregate.** Steps 1–3 run over many seeds. Per-seed metrics are
   aggregated into mean, standard deviation, and a **95% percentile bootstrap
   confidence interval**.
5. **Significance.** A **paired Wilcoxon signed-rank test** compares the fuzzy
   method against each baseline across seeds on the test split.

## Outputs

Running `--aggregate` (or `--all`) writes to the output directory:

| File | Contents |
|---|---|
| `raw_results.csv` | One row per (dataset, method, seed, split). |
| `summary.csv` / `summary.json` | Mean, std, and 95% CI per metric. |
| `significance.csv` | Paired Wilcoxon tests vs. each baseline. |
| `report.md` | Human-readable report. |
| `results_table.tex` | LaTeX (booktabs) table for papers. |
| `manifest.json` | Versions, git commit, configuration, timestamp. |

## Resumability

Results are appended row-by-row and completed `(dataset, method, seed)` triples
are skipped on re-run. Combine with `--time-budget` to run in fixed-length
chunks:

```bash
# Repeat until it reports "0 remaining":
python -m neurofuzzy.benchmark --run --seeds 10 --time-budget 30 --out benchmarks/results
python -m neurofuzzy.benchmark --aggregate --out benchmarks/results
```

## Reproducibility caveats

All metrics depend on the optimizer configuration (`--maxiter`, `--popsize`,
`--train-ratio`) and on the small benchmark sample sizes. Always report the
configuration alongside any cited number; the `manifest.json` captures it for you.
