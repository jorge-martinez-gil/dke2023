# Results

The committed results under `benchmarks/` are fully regenerable:

```bash
make repro     # or: python -m neurofuzzy.benchmark --all && python -m neurofuzzy.visualize
```

- Full report: [`benchmarks/results/report.md`](https://github.com/jorge-martinez-gil/dke2023/blob/main/benchmarks/results/report.md)
- LaTeX table: [`benchmarks/results/results_table.tex`](https://github.com/jorge-martinez-gil/dke2023/blob/main/benchmarks/results/results_table.tex)
- Figures: [`benchmarks/figures/`](https://github.com/jorge-martinez-gil/dke2023/tree/main/benchmarks/figures)

## Headline (test split, mean over 10 seeds)

| Dataset | Method | Pearson r | Spearman ρ |
|---|---|---:|---:|
| MC | best-single-feature | 0.826 | 0.750 |
| MC | mean-of-features | 0.825 | 0.756 |
| MC | linear-regression | 0.797 | 0.700 |
| MC | neurofuzzy | 0.741 | 0.679 |
| Geresid | linear-regression | 0.590 | 0.440 |
| Geresid | best-single-feature | 0.576 | 0.532 |
| Geresid | mean-of-features | 0.551 | 0.518 |
| Geresid | neurofuzzy | 0.368 | 0.375 |

!!! note "Honest interpretation"
    On these small datasets the fuzzy controller overfits the training split and
    does not beat the simple baselines on held-out data. This is reported
    transparently; see the README and [roadmap](https://github.com/jorge-martinez-gil/dke2023/blob/main/ROADMAP.md).
