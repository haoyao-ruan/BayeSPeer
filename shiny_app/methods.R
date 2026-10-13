# ============================================================================
# BayeSPeer and the 13 baseline rank aggregators, run as in the paper
# Companion code to:
#   Ruan, Wang, Lim & Wang, "Bias Meets Bayes: Bayesian Aggregation for
#   Peer and Self Ranking" (BayeSPeer)
#
# BayeSPeer: the two-step empirical-Bayes Gibbs sampler of Sections 2.3-2.4.
#   .bayespeer_gibbs() follows BayeSPeer.pre() in analysis/bayespeer_core.r
#   step by step (same full conditionals, update order and random-number
#   calls; only the bookkeeping is faster), and bayespeer_fit() follows
#   bayespeer_2step(), so the same seed gives the same draws.
# Baselines: the packages and calls used in the paper's analysis scripts:
#   TopKLists (Borda Mean/Median/Geomean/L2norm, MC1-MC3, CEMC.s, CEMC.k),
#   RobustRankAggreg (RRA, Stuart, min) and the BIRRA code of Badgeley et al.
#
# Input convention throughout:
#   R  : n x n integer matrix. Column j = ranker j's complete ranking of all
#        n individuals (rows). r_ij = rank given by ranker j to individual i,
#        1 = best. Diagonal r_jj = self-rank.
# Requires: truncnorm, TopKLists, RobustRankAggreg
# ============================================================================

# ---------------------------------------------------------------- validation
validate_rank_matrix <- function(R) {
  if (!is.matrix(R)) R <- as.matrix(R)
  if (nrow(R) != ncol(R))
    stop(sprintf("Matrix must be square (individuals x rankers). Got %d x %d.",
                 nrow(R), ncol(R)))
  n <- nrow(R)
  if (n < 3) stop("Need at least 3 individuals.")
  if (any(is.na(R))) stop("Missing values found. Complete rankings required.")
  if (any(abs(R - round(R)) > 1e-8)) stop("Ranks must be integers.")
  bad <- which(vapply(seq_len(n),
                      function(j) !all(sort(as.integer(R[, j])) == seq_len(n)),
                      logical(1)))
  if (length(bad))
    stop(sprintf("Column(s) %s are not complete permutations of 1..%d (ties/partial lists are not supported).",
                 paste(bad, collapse = ", "), n))
  storage.mode(R) <- "integer"
  if (is.null(rownames(R))) rownames(R) <- paste0("Item", seq_len(n))
  if (is.null(colnames(R))) colnames(R) <- paste0("Ranker", seq_len(n))
  R
}

# ============================================================================
# BayeSPeer
# ============================================================================

# Inverse-Gamma shape and rate from a mean and SD (find_shape_rate() in
# analysis/utils.R; the paper passes the median of the draws as "mean")
find_shape_rate <- function(mean, sd) {
  alpha <- (mean^2) / sd^2 + 2
  list(shape = alpha, rate = mean * (alpha - 1))
}

# One Gibbs chain: BayeSPeer.pre() of analysis/bayespeer_core.r.
# Full conditionals (paper Sec. 2.3), updated in this order each iteration:
#   mu_i, beta_j, mu_beta ~ TN(L, U), sigma2_beta, sigma2_j, then the latent
#   scores omega_ij column by column, each column in a fresh random order
#   (Li et al., 2018). Returns the full chains (burn-in included), as
#   BayeSPeer.pre() does; summaries use the draws after burn-in.
# Variance priors: inverse gamma IG(a, b) for sigma2_beta and IG(c, d) for
#   sigma2_j, or, with hc = c(A_beta, A_epsilon), half-Cauchy priors with these
#   scales on the standard deviations through the Makalic-Schmidt auxiliary
#   variables, as BayeSPeer.pre.HC() in analysis/sampler_half_cauchy.R (same
#   updates and order of random draws).
.bayespeer_gibbs <- function(n, r, M, burnin, W, mu0, mubeta0, Var_beta0,
                             Var_epsilon0, minbeta, maxbeta, a, b, c, d,
                             progress = NULL, hc = NULL) {
  N <- J <- n
  rtn <- truncnorm::rtruncnorm
  # the C routine behind truncnorm::rtruncnorm(), called directly in the
  # latent-score loop to skip the wrapper's argument checks (same draws)
  rtn_c <- tryCatch(get("C_do_rtruncnorm", envir = asNamespace("truncnorm")),
                    error = function(e) NULL)
  Mu <- Beta <- matrix(1, N, M)
  Var_epsilon <- matrix(1, N, M)
  Mu_beta <- Var_beta <- matrix(1, 1, M)
  Mu[, 1] <- mu0
  Mu_beta[1, 1] <- mubeta0
  Var_beta[1, 1] <- Var_beta0
  Beta[, 1] <- rnorm(J, mubeta0, sd = sqrt(Var_beta0))
  Var_epsilon[, 1] <- Var_epsilon0
  if (!is.null(hc)) { A_b <- hc[1]; A_e <- hc[2]; xi_b <- 1; xi_e <- rep(1, N) }
  # latent scores kept in a vector with two sentinels (-100, 100) after the
  # N*N cells; lo_at/hi_at hold, for every cell, the position of the item
  # ranked one place below (lower bound) and one place above (upper bound) in
  # the same ranking, or a sentinel for the end bounds of ranks n and 1
  # (+/-100, as in BayeSPeer.pre() and the analyses reported in the paper)
  NN <- N * J
  ord <- apply(r, 2, order)                       # ord[k, j]: item ranked k by j
  lo_at <- hi_at <- integer(NN)
  for (j in seq_len(J)) for (i in seq_len(N)) {
    k <- r[i, j]; cell <- (j - 1L) * N + i
    lo_at[cell] <- if (k < N) ord[k + 1, j] + (j - 1L) * N else NN + 1L
    hi_at[cell] <- if (k > 1) ord[k - 1, j] + (j - 1L) * N else NN + 2L
  }
  Wx <- c(as.vector(W), -100, 100)
  dg <- cbind(seq_len(N), seq_len(N))
  tick <- max(1L, M %/% 50L)
  for (m in 2:M) {
    W <- matrix(Wx[seq_len(NN)], N, J)
    ## 1. mu_i
    ve <- Var_epsilon[, m - 1]
    Wb <- W; Wb[dg] <- W[dg] - Beta[, m - 1]
    g <- vapply(seq_len(N), function(i) sum(Wb[i, ] / ve), numeric(1))
    h <- sum(1 / ve) + 1
    Mu[, m] <- rnorm(N, mean = g / h, sd = sqrt(1 / h))
    ## 2. beta_j
    g <- ((W[dg] - Mu[, m]) / ve) + (Mu_beta[m - 1] / Var_beta[m - 1])
    h <- (1 / ve) + (1 / Var_beta[m - 1])
    Beta[, m] <- rnorm(J, g / h, sqrt(1 / h))
    ## 3. mu_beta ~ TN(mean(beta), sigma2_beta / J; L, U)
    Mu_beta[m] <- rtn(1, a = minbeta, b = maxbeta, mean = mean(Beta[, m]),
                      sd = sqrt(Var_beta[m - 1] / J))
    E <- W - Mu[, m]; E[dg] <- E[dg] - Beta[, m]
    if (is.null(hc)) {
      ## 4. sigma2_beta
      Var_beta[m] <- 1 / rgamma(1, shape = a + (J / 2),
                                rate = b + sum((Beta[, m] - Mu_beta[, m])^2) / 2)
      ## 5. sigma2_j
      rate_e <- vapply(seq_len(J), function(j) d + sum(E[, j]^2) / 2, numeric(1))
      Var_epsilon[, m] <- 1 / rgamma(J, shape = c + (J / 2), rate = rate_e)
    } else {
      ## 4. sigma2_beta and its auxiliary variable (half-Cauchy)
      Var_beta[m] <- 1 / rgamma(1, shape = (J + 1) / 2,
                                rate = 1 / xi_b + sum((Beta[, m] - Mu_beta[m])^2) / 2)
      xi_b <- 1 / rgamma(1, shape = 1, rate = 1 / A_b^2 + 1 / Var_beta[m])
      ## 5. sigma2_j and their auxiliary variables, ranker by ranker
      for (j in seq_len(J)) {
        Var_epsilon[j, m] <- 1 / rgamma(1, shape = (N + 1) / 2,
                                        rate = 1 / xi_e[j] + sum(E[, j]^2) / 2)
        xi_e[j] <- 1 / rgamma(1, shape = 1, rate = 1 / A_e^2 + 1 / Var_epsilon[j, m])
      }
    }
    ## 6. latent scores omega_ij: truncated normal between the rank
    ##    neighbours, column by column, each column in a random order
    mu_m <- Mu[, m]; sd_m <- sqrt(Var_epsilon[, m])
    for (j in seq_len(J)) {
      mj <- mu_m; mj[j] <- mu_m[j] + Beta[j, m]    # self-bias on the diagonal
      sdj <- sd_m[j]; off <- (j - 1L) * N
      if (is.null(rtn_c)) {
        for (i in sample(N)) {
          cell <- off + i
          Wx[cell] <- rtn(1, a = Wx[lo_at[cell]], b = Wx[hi_at[cell]], mean = mj[i], sd = sdj)
        }
      } else {
        for (i in sample(N)) {
          cell <- off + i
          Wx[cell] <- .Call(rtn_c, 1L, Wx[lo_at[cell]], Wx[hi_at[cell]], mj[i], sdj)
        }
      }
    }
    if (!is.null(progress) && m %% tick == 0L) progress(m / M)
  }
  list(mu = Mu, Beta_draws = Beta, muBeta_draws = Mu_beta,
       varBeta_draws = Var_beta, varEpsilon_draws = Var_epsilon)
}

# Two-step empirical-Bayes fit: bayespeer_2step() of analysis/bayespeer_core.r
# (both steps run M iterations with burn-in M / 2; step 1 vague priors
# mu_beta ~ U(-10, 10) and the chosen variance prior; step 2 mu_beta range =
# mean +/- 3 SD of the step-1 beta draws, and the same variance prior family
# with its parameters set from the step-1 draws).
# prior = "IG": IG(prior_value, prior_value) for all variances in step 1
#   (0.1 in Section 2.4 of the paper); step 2 IG shapes and rates from the
#   median and SD of the step-1 variance draws (full chains, as
#   bayespeer_2step()).
# prior = "HC": half-Cauchy(0, prior_value) on the standard deviations in
#   step 1; step 2 half-Cauchy with scales sqrt(median) of the step-1 variance
#   draws after burn-in, as in the paper's sensitivity analysis (Section S2).
bayespeer_fit <- function(R, M = 40000, burnin = M / 2, seed = 42,
                          prior = c("IG", "HC"), prior_value = NULL,
                          unit_sd = 3, progress = NULL) {
  prior <- match.arg(prior)
  if (is.null(prior_value)) prior_value <- if (prior == "IG") 0.1 else 1
  hc1 <- if (prior == "HC") c(prior_value, prior_value) else NULL
  dlt <- if (prior == "IG") prior_value else NA_real_
  R <- validate_rank_matrix(R)
  n <- nrow(R)
  if (!is.null(seed)) set.seed(seed)
  ## initial latent scores: N(0, 1) draws ordered to agree with each ranking
  init_W <- matrix(rnorm(n^2, mean = 0, sd = 1), nrow = n)
  sorted_W <- apply(init_W, 2, sort, decreasing = TRUE)
  initW <- matrix(unlist(mapply(function(s, rk) s[rk], as.data.frame(sorted_W),
                                as.data.frame(R))), nrow = n)
  ## step 1: vague hyper-priors
  p1 <- if (!is.null(progress)) function(f) progress(0.5 * f) else NULL
  mu0 <- rnorm(n, 0, 1); mubeta0 <- rnorm(1, 0, 1)
  pre <- .bayespeer_gibbs(n, R, M, burnin, W = initW, mu0 = mu0, mubeta0 = mubeta0,
                          Var_beta0 = 1, Var_epsilon0 = 1, minbeta = -10, maxbeta = 10,
                          a = dlt, b = dlt, c = dlt, d = dlt, progress = p1, hc = hc1)
  ## empirical-Bayes hyper-priors from the step-1 draws (BayeSPeer.pre.HC()
  ## returns the draws after burn-in, BayeSPeer.pre() the full chains)
  if (prior == "HC") {
    k1 <- -(seq_len(burnin))
    pre$Beta_draws <- pre$Beta_draws[, k1, drop = FALSE]
    pre$varBeta_draws <- pre$varBeta_draws[, k1, drop = FALSE]
    pre$varEpsilon_draws <- pre$varEpsilon_draws[, k1, drop = FALSE]
  }
  Beta.bar <- mean(pre$Beta_draws); Beta.sd <- sd(pre$Beta_draws)
  L <- Beta.bar - unit_sd * Beta.sd; U <- Beta.bar + unit_sd * Beta.sd
  varBeta_est <- median(pre$varBeta_draws)
  varEpsilon_est <- median(pre$varEpsilon_draws)
  if (prior == "IG") {
    IG_b <- find_shape_rate(varBeta_est, sd(pre$varBeta_draws))
    IG_e <- find_shape_rate(varEpsilon_est, sd(pre$varEpsilon_draws))
    ab <- c(IG_b$shape, IG_b$rate, IG_e$shape, IG_e$rate); hc2 <- NULL
    hyper2 <- list(ig_beta = unlist(IG_b), ig_eps = unlist(IG_e))
  } else {
    ab <- rep(NA_real_, 4); hc2 <- c(sqrt(varBeta_est), sqrt(varEpsilon_est))
    hyper2 <- list(hc_scale_beta = hc2[1], hc_scale_eps = hc2[2])
  }
  rm(pre)
  ## step 2: main chain under the data-informed hyper-priors
  p2 <- if (!is.null(progress)) function(f) progress(0.5 + 0.5 * f) else NULL
  mu0 <- rnorm(n, 0, 1); mubeta0 <- rnorm(1, 0, 1)
  res <- .bayespeer_gibbs(n, R, M, burnin, W = initW, mu0 = mu0, mubeta0 = mubeta0,
                          Var_beta0 = varBeta_est, Var_epsilon0 = varEpsilon_est,
                          minbeta = L, maxbeta = U, a = ab[1], b = ab[2],
                          c = ab[3], d = ab[4], progress = p2, hc = hc2)
  ## posterior summaries from the draws after burn-in (paper Sec. 2.3)
  keep <- -(seq_len(burnin))
  mu_d <- t(res$mu[, keep, drop = FALSE])           # draws x individuals
  be_d <- t(res$Beta_draws[, keep, drop = FALSE])   # draws x rankers
  s2_d <- t(res$varEpsilon_draws[, keep, drop = FALSE])
  ids <- rownames(R)
  colnames(mu_d) <- colnames(be_d) <- colnames(s2_d) <- ids
  mu_hat <- apply(mu_d, 2, mean)                    # post_mean_mus in BayeSPeer.pre()
  final_rank <- rank(-mu_hat)                       # post_rank in BayeSPeer.pre()
  ## rank probability profile: share of draws in which individual i's mu has
  ## each rank among all mu's (paper Sec. 2.3, Table 4)
  rk <- t(apply(mu_d, 1, function(x) rank(-x, ties.method = "first")))
  rank_prob <- vapply(seq_len(n), function(i) tabulate(rk[, i], nbins = n) / nrow(rk),
                      numeric(n))
  dimnames(rank_prob) <- list(paste0("Rank", seq_len(n)), ids)
  qs <- function(X, p) apply(X, 2, quantile, p, names = FALSE)
  summ <- data.frame(
    id = ids, final_rank = final_rank,
    mu_mean = mu_hat, mu_lo = qs(mu_d, .025), mu_hi = qs(mu_d, .975),
    beta_mean = apply(be_d, 2, mean), beta_lo = qs(be_d, .025), beta_hi = qs(be_d, .975),
    sig2_med = apply(s2_d, 2, median),
    geweke_mu = apply(mu_d, 2, geweke_z), geweke_beta = apply(be_d, 2, geweke_z),
    row.names = NULL, check.names = FALSE)
  list(draws = list(mu = mu_d, beta = be_d, sig2 = s2_d,
                    mu_beta = as.numeric(res$muBeta_draws[, keep]),
                    sig2_beta = as.numeric(res$varBeta_draws[, keep])),
       summary = summ, rank_prob = rank_prob,
       final_rank = stats::setNames(final_rank, ids),
       hyper = c(list(mu_beta_range = c(L, U)), hyper2),
       ids = ids, n = n,
       settings = list(M = M, burnin = burnin, seed = seed, prior = prior,
                       prior_value = prior_value,
                       prior_label = if (prior == "IG") sprintf("IG(%s, %s)", prior_value, prior_value)
                                     else sprintf("half-Cauchy(0, %s)", prior_value)))
}

# Geweke z-score (first 10% vs last 50%) with AR spectral variance, as in coda
.spectrum0 <- function(x) {
  x <- as.numeric(x)
  if (var(x) == 0) return(0)
  a <- try(ar(x, aic = TRUE), silent = TRUE)
  if (inherits(a, "try-error")) return(var(x))
  if (length(a$ar) == 0) return(a$var.pred)
  a$var.pred / (1 - sum(a$ar))^2
}

geweke_z <- function(x, frac1 = 0.1, frac2 = 0.5) {
  nx <- length(x)
  x1 <- x[seq_len(max(2, floor(frac1 * nx)))]
  x2 <- x[(nx - floor(frac2 * nx) + 1):nx]
  (mean(x1) - mean(x2)) / sqrt(.spectrum0(x1) / length(x1) + .spectrum0(x2) / length(x2))
}

# ============================================================================
# Baselines, called as in the paper's scripts (05a_RWD_2020.R for the real
# data; 02A2_compare_selfrm_parallel.R for the .selfrm variants, which set
# the self-ranks to NA and pass the shortened lists to the same functions)
# ============================================================================

# convert_to_rank_lists() of analysis/utils.R: each ranker's list of names,
# best first (NA self-ranks drop out)
convert_to_rank_lists <- function(rank_matrix) {
  sorted_mat <- apply(rank_matrix, 2, function(col) names(sort(col)))
  if (is.matrix(sorted_mat))
    lapply(seq_len(ncol(sorted_mat)), function(j) sorted_mat[, j])
  else unname(sorted_mat)
}

# BIRRA (Badgeley et al., 2015): the code the paper sourced (BIRRA.R, also in
# analysis/competing_methods.R), unchanged
BIRRA <- function(data, prior=0.05, num.bin=50, num.iter=10, return.all=F, plot.legend=F, grp=NULL, cor.stop=1, ...){
  nr=nrow(data)
  nrp=floor(nrow(data)*prior)
  data=apply(data,2,rank)/nr
  nc=ncol(data)
  TPR=FPR=Bayes.factors=matrix(ncol=nc, nrow=num.bin)
  binned.data=ceiling(data*num.bin)
  bayes.data=matrix(nrow=nrow(data), ncol=ncol(data))
  guess=apply(data,1,mean)
  cprev=0
  for ( iter in 1:num.iter){
    if((cor.stop-cprev)>1e-15| is.na(cprev)){
      guesslast=guess
      oo=order(guess)
      guess[oo[1:nrp]]=1
      guess[oo[(nrp+1):nr]]=0
      for (i in 1:nc){
        for (bin in 1:num.bin){
          frac=bin/num.bin
          TPR=sum(guess[binned.data[,i]<=bin])
          FPR=sum((!guess)[binned.data[,i]<=bin])
          Bayes.factors[bin,i]=log((TPR+1)/(FPR+1)/(prior/(1-prior)))
        }
      }
      Bayes.factors=apply(Bayes.factors,2,smooth)
      Bayes.factors=apply(Bayes.factors,2,function(x){rev(cummax(rev(x)))})
      if(is.null(grp)){
        matplot(1:num.bin, Bayes.factors, type="l", lwd=2, ...)
      }
      else{
        matplot(1:num.bin, Bayes.factors, type="l", lwd=2, lty=grp, col=grp)
      }
      title(paste("Iteration", iter))
      if (iter==1&plot.legend){
        legend("topright", col=1:5, lty=1:4, legend=colnames(data), lwd=2, ncol=2)
      }
      for (bin in 1:num.bin){
        oo=order(Bayes.factors[bin,], decreasing=T)
        Bayes.factors[bin, oo[1]]=Bayes.factors[bin, oo[2]]
      }
      for (i in 1:nc){
        bayes.data[,i]=Bayes.factors[binned.data[,i],i]
      }
      bb=exp(apply(bayes.data,1, sum))
      f=prior/(1-prior)
      prob=bb*f/(1+bb*f)
      exp=sort(prob, decreasing=F)[nrp]
      guess=rank(-apply(bayes.data,1, sum))
      cprev=cor(guess, guesslast)
      message("correlation with pervious iteration=",cprev)
    }
    else{
      message("Converged");
      break
    }
  }
  if(return.all){
    return(list(result=guess, data=bayes.data, BF=Bayes.factors))
  }
  else{
    guess
  }
}

METHOD_INFO <- data.frame(
  method = c("BayeSPeer", "BIRRA", "RRA", "Stuart", "min",
             "MC1", "MC2", "MC3", "CEMC.s", "CEMC.k",
             "Mean", "Geomean", "L2norm", "Median"),
  category = c("BayeSPeer", "Bayesian", "Distribution-based", "Distribution-based",
               "Distribution-based", "Markov chain", "Markov chain", "Markov chain",
               "Optimization", "Optimization", "Summary statistic",
               "Summary statistic", "Summary statistic", "Summary statistic"),
  stringsAsFactors = FALSE)

# The 13 baselines on one rank matrix (self-ranks set to NA when selfrm = TRUE).
# Each aggregated ranking is the position of every individual in the method's
# ordered output list (match(itemnames, list) in the paper's scripts).
run_baselines <- function(R, selfrm = FALSE, seed = 42) {
  Rm <- R
  if (selfrm) diag(Rm) <- NA
  itemnames <- rownames(R)
  rank_list <- convert_to_rank_lists(Rm)
  pos <- function(lst) match(itemnames, lst)
  quiet <- function(expr) {              # these functions print and plot
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
    out <- NULL
    utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
    out
  }
  set.seed(seed)                          # CEMC is stochastic
  borda  <- quiet(TopKLists::Borda(rank_list)$TopK)
  mc     <- quiet(TopKLists::MC(rank_list))
  cemc_s <- quiet(TopKLists::CEMC(rank_list, dm = "s")$TopK)
  cemc_k <- quiet(TopKLists::CEMC(rank_list, dm = "k")$TopK)
  stuart <- RobustRankAggreg::aggregateRanks(glist = rank_list, method = "stuart")$Name
  rra    <- RobustRankAggreg::aggregateRanks(glist = rank_list, method = "RRA")$Name
  minr   <- RobustRankAggreg::aggregateRanks(glist = rank_list, method = "min")$Name
  birra  <- quiet(BIRRA(Rm, num.iter = 50))
  out <- list(BIRRA = as.numeric(birra),
              RRA = pos(rra), Stuart = pos(stuart), min = pos(minr),
              MC1 = pos(mc$MC1.TopK), MC2 = pos(mc$MC2.TopK), MC3 = pos(mc$MC3.TopK),
              CEMC.s = pos(cemc_s), CEMC.k = pos(cemc_k),
              Mean = pos(borda$mean), Geomean = pos(borda$geo.mean),
              L2norm = pos(borda$l2norm), Median = pos(borda$median))
  if (selfrm) names(out) <- paste0(names(out), ".selfrm")
  out
}

aggregate_all <- function(R, fit = NULL, include_selfrm = FALSE, seed = 42, ...) {
  R <- validate_rank_matrix(R)
  if (is.null(fit)) fit <- bayespeer_fit(R, seed = seed, ...)
  res <- c(list(BayeSPeer = as.numeric(fit$final_rank)),
           run_baselines(R, FALSE, seed))
  if (include_selfrm) res <- c(res, run_baselines(R, TRUE, seed))
  tab <- do.call(rbind, res)
  colnames(tab) <- rownames(R)
  list(table = tab, fit = fit)
}

# ---------------------------------------------------------------------------
# Quantities of the paper's Fig. 9 (Sec. 4.2): posterior mean of mu and beta,
# posterior median of sigma2_j, and the observed rank difference
# (rank_diff_self_vs_peers() in analysis/utils.R: the rank of each
# individual's average peer rank, self excluded, minus the self-rank),
# plus an optional external measure (e.g., an exam grade)
# ---------------------------------------------------------------------------
rank_diff_self_vs_peers <- function(rank_mat) {
  peer_mat <- rank_mat
  diag(peer_mat) <- NA
  peer_means <- rowMeans(peer_mat, na.rm = TRUE)
  peer_ranks <- rank(peer_means)
  list(meandiff = peer_means - diag(rank_mat), rankdiff = peer_ranks - diag(rank_mat),
       peer_means = peer_means, peer_ranks = peer_ranks)
}

analysis_frame <- function(R, fit, external = NULL) {
  df <- data.frame(
    id = fit$ids,
    Performance_mu  = fit$summary$mu_mean,
    SelfBias_beta   = fit$summary$beta_mean,
    RankDifference  = unname(rank_diff_self_vs_peers(R)$rankdiff),
    Variability_s2  = fit$summary$sig2_med,
    check.names = FALSE)
  if (!is.null(external)) df <- cbind(df[, 1, drop = FALSE],
                                      External = external,
                                      df[, -1, drop = FALSE])
  df
}

# ---------------------------------------------------------------- example data
example_class1 <- function() {
  R <- matrix(c(
    1,1,4,2,5,1,1,1,
    2,2,7,3,3,2,7,2,
    3,4,5,4,7,3,2,3,
    4,6,2,1,4,6,3,4,
    8,3,1,6,1,4,8,7,
    7,7,3,5,2,5,6,6,
    5,8,6,7,6,7,5,8,
    6,5,8,8,8,8,4,5), nrow = 8, byrow = TRUE)
  dimnames(R) <- list(paste0("Team", 1:8), paste0("Ranker", 1:8))
  R
}

example_class2 <- function() {
  R <- matrix(c(
     1,14,16,15, 2,16,13, 7,16,13,14,15,13,15,14,16,
     4, 4,14,10, 4,12, 8, 4,13, 2,11,13,10, 6,11,14,
    11, 1, 8, 3,11,13,16, 1, 5, 6, 5,10, 2, 8, 9, 1,
    10, 6, 3, 2, 6,10,11,11,12,12, 8, 8, 6, 2, 2, 6,
     2,13,15,14, 1,15,14, 8,15,14,13,16,14,16,13,15,
    14,10, 1, 9,14, 2, 4,16, 1,10, 4, 3, 3,10, 6, 4,
     5,15, 5,13, 8, 3, 5,13, 3,16, 1, 1,15,13,16,11,
     3, 3,13, 6, 3,11, 7, 3,14, 1,12,14, 9, 5,12,13,
    16,11,12, 7, 9, 5,10,10, 8, 8,10, 5, 8,11, 7,10,
     6,16, 6,16, 7, 4, 6,14, 4,15, 2, 2,16,14,15,12,
    13, 9, 2, 5,13, 1, 3,15, 2, 9, 3, 4, 4, 9, 5, 3,
    15,12,11, 8,10, 6, 9, 9, 7, 7, 9, 6, 7,12, 8, 9,
     9, 5, 4, 1, 5, 9,12,12,11,11, 7, 7, 5, 1, 1, 5,
    12, 2, 7, 4,12,14,15, 2, 6, 5, 6, 9, 1, 7,10, 2,
     8, 8, 9,11,16, 7, 2, 6,10, 4,16,12,12, 4, 4, 7,
     7, 7,10,12,15, 8, 1, 5, 9, 3,15,11,11, 3, 3, 8), nrow = 16, byrow = TRUE)
  dimnames(R) <- list(paste0("Student", 1:16), paste0("Ranker", 1:16))
  R
}
