# BayeSPeer

**BayeSPeer** is an R implementation of a bias-aware Bayesian rank-aggregation
framework for settings with peer- and self-evaluations. The method explicitly
models self-evaluation bias and uncertainty through a latent-utility
formulation, producing aggregated rankings with posterior uncertainty
quantification.

It accompanies the manuscript *"Bias Meets Bayes: Bayesian Aggregation for Peer
and Self Ranking."* This repository provides the Gibbs samplers, the
data-generating mechanisms, the evaluation metrics, and the full function
library (including the competing aggregators) used in the study.

---

## Repository structure

```
BayeSPeer/
├── analysis/                   # Full function library (modular), used in the study
│   ├── dependencies.R          #   package loads (source this first)
│   ├── utils.R                 #   rank conversions, calculation & Inverse-Gamma helpers
│   ├── data_generation.R       #   gen_data_dir / mix / t / gamma / tri
│   ├── metrics.R               #   Top-1 / Top-3 + rank-distance metrics
│   ├── bayespeer_core.r        #   Gibbs samplers + two-step EB wrapper (IG priors)
│   ├── sampler_half_cauchy.R   #   half-Cauchy prior variant
│   └── competing_methods.R     #   Borda family, BIRRA, BiGER
├── examples/
│   └── example_synthetic.R     # runnable end-to-end demo on simulated data
├── BayeSPeer.Rproj
├── .gitignore
└── README.md
```

> **Note on file names.** The core file is `analysis/bayespeer_core.r`
> (lower-case name, lower-case `.r` extension). On case-sensitive systems use
> that exact name in `source()`.

---

## Core API (`analysis/bayespeer_core.r`)

### Bayesian rank aggregation
- **`BayeSPeer()`** — main Gibbs sampler (Inverse-Gamma priors); returns
  posterior summaries of the latent performance scores, self-evaluation bias,
  and variance components.
- **`BayeSPeer.pre()`** — same sampler, but additionally returns the full
  MCMC draws (useful for diagnostics and hyper-prior tuning).
- **`bayespeer_2step()`** — two-stage procedure that first runs a pilot chain
  of `BayeSPeer.pre()` with diffuse hyper-priors to estimate scale parameters,
  then fits the final model with empirically tuned hyper-priors.
- **`BayeSPeer.HC()` / `BayeSPeer.pre.HC()`** (`sampler_half_cauchy.R`) —
  half-Cauchy prior variant of the variance components.

### Data generation (simulations, `data_generation.R`)
- **`gen_data_dir()`** — Normal latent effects
- **`gen_data_mix()`** — Normal latent effects with ranker-specific noise levels
- **`gen_data_t()`** — heavy-tailed (t-distributed) latent effects, rescaled to
  unit variance (requires the `df` argument)
- **`gen_data_gamma()`** — skewed positive latent effects
- **`gen_data_tri()`** — tri-modal self-evaluation bias
  (over-/un-/under-confident ranker types)

Each generator simulates latent utilities, converts them to observed ranks, and
regenerates the data in the degenerate case where all rows of the resulting
rank matrix coincide.

### Rank utilities & evaluation
- Conversion between rank matrices, ordinal lists, and numeric rankings
- Self-rank vs. peer-rank diagnostics; latent-utility initialization
- Top-1 / Top-3 accuracy and rank-distance metrics

---

## Quick start

```r
source("analysis/dependencies.R")     # package loads — source this first
source("analysis/utils.R")
source("analysis/data_generation.R")
source("analysis/metrics.R")
source("analysis/bayespeer_core.r")

set.seed(42)

# Simulate a small peer/self ranking data set
dat <- gen_data_dir(n = 8, mu_beta = 1.0, sigma_beta = sqrt(0.5), sigma_epsilon = 1.0)

# Fit BayeSPeer with the two-step tuned procedure
res <- bayespeer_2step(rank_matrix = dat$rank, M_itrns = 2000, m_burn = 1000, seed = 42)

res$post_mean_mus                     # posterior mean latent scores
agg_rank <- rank(-res$post_mean_mus)  # aggregated ranking (rank 1 = best)

calculate_top1(agg_rank, dat$true_rank)
calculate_top3(agg_rank, dat$true_rank)
```

A fuller, commented version is in
[`examples/example_synthetic.R`](examples/example_synthetic.R).

---

## Function library (`analysis/`)

The library is split into modules by concern:

- `dependencies.R` — loads the packages used by the library (source first);
- `utils.R`, `data_generation.R`, `metrics.R` — shared helpers, the five data
  generators, and the accuracy/rank-distance metrics;
- `bayespeer_core.r` — the Gibbs samplers under Inverse-Gamma priors
  (`BayeSPeer.pre`, `BayeSPeer`) together with the two-step empirical-Bayes
  wrapper `bayespeer_2step()`;
- `sampler_half_cauchy.R` — the half-Cauchy prior variant of the samplers;
- `competing_methods.R` — the competing aggregators: the Borda-score family,
  BIRRA, and BiGER.

The simulation, processing, and figure scripts used to produce the manuscript's
results are not distributed here.

---

## Dependencies

`MCMCpack`, `extraDistr`, `matrixStats`, `truncnorm`, `TopKLists`, `expm`,
`tidyverse`, `reshape2`.

```r
install.packages(c(
  "MCMCpack", "extraDistr", "matrixStats", "truncnorm",
  "TopKLists", "expm", "tidyverse", "reshape2"
))
```

`truncnorm` is called via `truncnorm::` inside the samplers, so it must be
installed even though `dependencies.R` does not attach it.

---

## Reproducibility note

BayeSPeer is a Monte Carlo method: posterior summaries carry Monte Carlo
error and depend on the random-number seed, R version, and platform. The
manuscript does not fix a seed, so small numerical differences from the
reported results are expected when re-running the code; differences of this
size do not affect the aggregated rankings or the substantive conclusions.

The samplers keep the latent scores of the first and last ranked
individuals within ±100, as in the analyses reported in the manuscript.
With the vague step 1 priors and very small groups, this bound keeps the
latent scale from drifting.

---

## Data & privacy

No real data is distributed with this repository; all code shipped here runs
on fully synthetic, simulated data. The classroom data analysed in the
manuscript were anonymized, and the resulting de-identified ranking matrices
are reported in full in the manuscript. The underlying identifiable records
(student identities and individual exam grades) are **not** available for
privacy reasons.
