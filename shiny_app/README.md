# BayeSPeer Shiny App

Companion app to **"Bias Meets Bayes: Bayesian Aggregation for Peer and Self
Ranking"** (Ruan, Wang, Lim & Wang). Upload a peer + self ranking dataset and
get the final ranking estimates from BayeSPeer and the **13 baseline methods**
compared in the paper, plus the real-world data analyses of Section 4.

## Quick start

```r
install.packages(c("shiny", "ggplot2", "GGally",
                   "truncnorm", "TopKLists", "RobustRankAggreg"))   # R >= 4.0
shiny::runApp("shiny_app")         # from the repository root
```

## Input format

Upload an `n x n` CSV, type the matrix directly into an editable grid in the
app (set the group size and click "Type it in"; supported for n ≤ 20), or load
one of the built-in examples. Format: **columns = rankers, rows = individuals,
entry (i, j) = rank given by ranker j to individual i (1 = best), diagonal =
self-ranks.** Every
column must be a complete permutation of 1..n (complete rankings, no ties —
as assumed in the paper). A header row and/or an ID column are auto-detected.
`example_classI.csv` (8 teams) and `example_classII.csv` (16 students) are the
two classroom datasets from Section 4; both are also loadable inside the app.

Optionally, on the **Bias vs performance** tab (after running the analysis),
type each individual's exam grade or other independent score into the score
form to reproduce the Fig. 9 correlation analysis against an external measure.
To try it without real scores, press **Use made-up example scores**: it fills
in simulated scores (each person's average peer rank plus random noise),
labels them "made up" on the plot, and notes that they are not real grades.

## What you get

| Tab | Paper analogue | Content |
|---|---|---|
| Aggregated rankings | Fig. 7 | Final ranks from BayeSPeer, BIRRA, RRA, Stuart, min, MC1–MC3, CEMC.s, CEMC.k, Mean, Geomean, L2norm, Median (+ optional `.selfrm` variants with self-ranks removed) |
| Rank probabilities | Table 4 | BayeSPeer posterior probability of each individual occupying each rank |
| Posteriors & bias | Fig. 8 | Posterior densities of true performance μᵢ and self-bias βⱼ; posterior summary with 95% credible intervals |
| Bias vs performance | Fig. 9 / Sec. 4.2 | Pairwise correlations among an external measure, μ̂, β̂, the observed rank difference (rank from the average peer rank, self excluded, minus the self-rank), and ranker variance σ̂²ⱼ |

Geweke convergence z-scores are included in the posterior summary CSV
download.

## Implementation notes

- **BayeSPeer** is the paper's two-step empirical-Bayes Gibbs sampler
  (Sec. 2.3–2.4). Step 1 uses vague priors: μ_β ~ U(−10, 10) and, for σ²_β
  and σ²ⱼ, the variance prior chosen in the app. The default is IG(0.1, 0.1),
  as in the paper; the alternatives are the priors of the paper's
  sensitivity analysis (Section S2): IG(0.001, 0.001), IG(0.01, 0.01),
  IG(1, 1) and half-Cauchy(0, A) on the standard deviations with
  A = 0.5, 1, 2 or 4. Step 2 sets the μ_β range to β̄ ± 3 SD of the step-1 β
  draws and keeps the same prior family, with its parameters set from the
  step-1 draws: IG shapes and rates from the median and SD of the variance
  draws, or half-Cauchy scales equal to the square root of their medians.
  The latent scores of the first and last ranked individuals are kept
  within ±100, as in the analyses reported in the paper. Both steps run the
  chosen number of iterations (default 4,000 with seed 123) with the first
  half as burn-in, and all results use the step-2 draws after burn-in. The
  sampler in `methods.R` follows `BayeSPeer.pre()` and `bayespeer_2step()`
  in `analysis/bayespeer_core.r` (and `BayeSPeer.pre.HC()` in
  `analysis/sampler_half_cauchy.R` for half-Cauchy priors) step by step, so
  with the same seed it returns the same draws.
- **Baselines** are called as in the paper's analysis scripts: Mean, Median,
  Geomean, L2norm (`TopKLists::Borda`), MC1–MC3 (`TopKLists::MC`), CEMC.s and
  CEMC.k (`TopKLists::CEMC`), RRA, Stuart and min
  (`RobustRankAggreg::aggregateRanks`), and BIRRA (the Badgeley et al. code,
  `num.iter = 50`). The `.selfrm` variants set the self-ranks to missing and
  pass the shortened lists to the same functions. CEMC is a stochastic search,
  so its ordering of nearly tied individuals can change with the seed.
- Runtime with 2,000 iterations per step: about 1 second for n = 5 and
  17 seconds for n = 25. BayeSPeer takes about 6 seconds of that at n = 25
  (the rest is mainly CEMC), and its time grows in proportion to the
  number of iterations.

## Figures and tables

The figures reuse the paper's own ggplot2/GGally code, so they look like the
paper: the aggregated rankings use the paper's method colors, method order,
grey panel headers and rotated method names (Fig. 7); the posterior densities
use filled curves with the paper's line types, in panels A and B (Fig. 8); and
the bias vs performance analysis is the paper's GGally pairs plot (Fig. 9).
The `.selfrm` variants, which the paper does not plot, use a lighter shade of
their method's color. Tables are ruled like the paper's (Tables 3 and 4). The
rank probability heatmap has no paper counterpart and follows the Fig. 8 style.
Figures keep the paper's proportions and are sized to fit a laptop screen
(`FIG` in `app.R`): at most 640 px wide for each posterior density panel,
600 px for the pairs plot and 820 px for the rank probability heatmap, and the
aggregated rankings figure puts as many teams per row as fit.

## Text sizes

All text in the app (page text, tables and figures) takes its size from one
scale defined near the top of `app.R`. `BASE_PX` (default 14) sets body and
table text, and every other role is a step of a 9:8 scale around it: notes,
axis ticks and bar or cell values 12, axis titles and legends 14, panel titles
16, section headings 18, correlation values 20. The figures keep the paper's
design and take only their text sizes from this scale. Figures render at 72 dpi, so
one point equals one CSS pixel and figure text matches the page text beside
it. To make the whole app larger or smaller, change `BASE_PX` only. The
aggregated rankings figure chooses how many panels to place per row from the
window width, so bar value labels never overlap at these sizes.

## Files

- `app.R` — Shiny UI + server
- `methods.R` — BayeSPeer sampler, the 13 baseline calls and the Fig. 9
  quantities (usable standalone:
  `source("methods.R"); bayespeer_fit(example_class1())`)
- `example_classI.csv`, `example_classII.csv` — the paper's two datasets
