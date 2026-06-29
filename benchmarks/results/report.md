# Neurofuzzy Benchmark Report

- Generated: 2026-06-29T04:59:09Z
- neurofuzzy 1.0.0 | numpy 2.2.6 | scipy 1.15.3 | python 3.10.12
- git commit: `526041114c15226e5bb9b9bf9d2161a50d2a5d2d`
- optimizer: maxiter=40, popsize=8, train_ratio=0.6
- seeds: 10 | datasets: mc, geresid

Values are **mean [95% bootstrap CI]** over seeds on the held-out test split.

## Dataset: GERESID

| Method | pearson | spearman | cosine_similarity | mean_absolute_error | root_mean_squared_error |
|---|---|---|---|---|---|
| best-single-feature | 0.5762 [0.4972, 0.6426] | 0.5318 [0.4714, 0.5779] | 0.7184 [0.6782, 0.7518] | 0.2635 [0.2468, 0.2803] | 0.3325 [0.3126, 0.3523] |
| linear-regression | 0.5898 [0.4906, 0.6695] | 0.4395 [0.3492, 0.5080] | 0.8766 [0.8496, 0.8994] | 0.1854 [0.1728, 0.2000] | 0.2253 [0.2084, 0.2451] |
| mean-of-features | 0.5510 [0.4548, 0.6242] | 0.5179 [0.4209, 0.5873] | 0.7854 [0.7465, 0.8180] | 0.2469 [0.2297, 0.2635] | 0.3243 [0.3031, 0.3442] |
| neurofuzzy | 0.3684 [0.2014, 0.5316] | 0.3749 [0.2231, 0.5182] | 0.8331 [0.7940, 0.8682] | 0.2116 [0.1914, 0.2333] | 0.2579 [0.2385, 0.2815] |

**Paired Wilcoxon signed-rank vs baselines (test split):**

| Comparison | Metric | n | mean diff | p-value |
|---|---|---|---|---|
| neurofuzzy vs best-single-feature | pearson | 10 | -0.2078 | 0.0488 |
| neurofuzzy vs best-single-feature | spearman | 10 | -0.1569 | 0.1055 |
| neurofuzzy vs linear-regression | pearson | 10 | -0.2214 | 0.0645 |
| neurofuzzy vs linear-regression | spearman | 10 | -0.0646 | 0.4316 |
| neurofuzzy vs mean-of-features | pearson | 10 | -0.1827 | 0.1934 |
| neurofuzzy vs mean-of-features | spearman | 10 | -0.1430 | 0.2324 |

## Dataset: MC

| Method | pearson | spearman | cosine_similarity | mean_absolute_error | root_mean_squared_error |
|---|---|---|---|---|---|
| best-single-feature | 0.8256 [0.7757, 0.8793] | 0.7501 [0.7106, 0.7810] | 0.9461 [0.9332, 0.9595] | 0.1518 [0.1338, 0.1687] | 0.2121 [0.1796, 0.2414] |
| linear-regression | 0.7966 [0.7400, 0.8549] | 0.6995 [0.6239, 0.7676] | 0.9411 [0.9251, 0.9568] | 0.1654 [0.1438, 0.1875] | 0.2253 [0.1935, 0.2562] |
| mean-of-features | 0.8246 [0.7834, 0.8679] | 0.7560 [0.7124, 0.7889] | 0.9479 [0.9380, 0.9578] | 0.1606 [0.1474, 0.1746] | 0.2143 [0.1900, 0.2383] |
| neurofuzzy | 0.7408 [0.6701, 0.8056] | 0.6791 [0.6367, 0.7149] | 0.9290 [0.9134, 0.9436] | 0.2509 [0.2284, 0.2761] | 0.3015 [0.2790, 0.3253] |

**Paired Wilcoxon signed-rank vs baselines (test split):**

| Comparison | Metric | n | mean diff | p-value |
|---|---|---|---|---|
| neurofuzzy vs best-single-feature | pearson | 10 | -0.0848 | 0.0020 |
| neurofuzzy vs best-single-feature | spearman | 10 | -0.0710 | 0.0039 |
| neurofuzzy vs linear-regression | pearson | 10 | -0.0558 | 0.1309 |
| neurofuzzy vs linear-regression | spearman | 10 | -0.0204 | 0.6250 |
| neurofuzzy vs mean-of-features | pearson | 10 | -0.0838 | 0.0098 |
| neurofuzzy vs mean-of-features | spearman | 10 | -0.0768 | 0.0020 |

> Results are produced by `python -m neurofuzzy.benchmark`. They depend on the optimizer configuration above and the small benchmark sample sizes; report the configuration alongside any cited number.
