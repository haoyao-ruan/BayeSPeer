# ============================================================================
# sampler_half_cauchy.R — Gibbs samplers with half-Cauchy priors
# Half-Cauchy prior variant of the variance components.
# Part of the BayeSPeer function library.
# ============================================================================

#------------------Gibbs Sampler Preliminary-------------------
BayeSPeer.pre.HC <- function(n, r, M, burnin, 
                             W,
                             mu0,
                             mubeta0,
                             beta0,
                             Var_beta0,
                             Var_epsilon0,
                             minbeta, maxbeta,
                             A_beta = 2.5,
                             A_epsilon = 2.5) {
  
  N <- n
  J <- n
  
  Mu <- Beta <- matrix(1, N, M)
  Var_epsilon <- matrix(1, N, M)
  Mu_beta <- Var_beta <- matrix(1, 1, M)
  
  ## Makalic-Schmidt auxiliary variables
  Xi_beta <- matrix(1, 1, M)
  Xi_epsilon <- matrix(1, N, M)
  
  Mu[, 1] <- mu0
  Mu_beta[1, 1] <- mubeta0
  Var_beta[1, 1] <- Var_beta0
  
  beta0 <- rnorm(J, mubeta0, sd = sqrt(Var_beta0))
  Beta[, 1] <- beta0
  
  Var_epsilon[, 1] <- Var_epsilon0
  
  Xi_beta[1, 1] <- 1
  Xi_epsilon[, 1] <- 1
  
  for (m in 2:M) {
    
    ## 1. Update Mu
    for (i in 1:N) {
      indBeta <- Beta[, m - 1]
      indBeta[-i] <- 0 
      
      g <- sum((W[i, ] - indBeta) / Var_epsilon[, m - 1])
      h <- sum(1 / Var_epsilon[, m - 1]) + 1
      
      Mu[i, m] <- rnorm(1, mean = g / h, sd = sqrt(1 / h)) 
    }
    
    
    ## 2. Update Beta
    for (j in 1:J) {
      g <- ((W[j, j] - Mu[j, m]) / Var_epsilon[j, m - 1]) +
        (Mu_beta[m - 1] / Var_beta[m - 1])
      
      h <- (1 / Var_epsilon[j, m - 1]) + (1 / Var_beta[m - 1])
      
      Beta[j, m] <- rnorm(1, mean = g / h, sd = sqrt(1 / h))
    } 
    
    ## 3. Update mu_beta
    Mu_beta[m] <- truncnorm::rtruncnorm(
      1,
      a = minbeta,
      b = maxbeta,
      mean = mean(Beta[, m]),
      sd = sqrt(Var_beta[m - 1] / J)
    )
    
    ## 4. Update Var_beta using Makalic-Schmidt half-Cauchy
    ## sigma_beta ~ C+(0, A_beta)
    Var_beta[m] <- 1 / rgamma(
      1,
      shape = (J + 1) / 2,
      rate = 1 / Xi_beta[m - 1] + 
        sum((Beta[, m] - Mu_beta[m])^2) / 2
    )
    
    Xi_beta[m] <- 1 / rgamma(
      1,
      shape = 1,
      rate = 1 / A_beta^2 + 1 / Var_beta[m]
    )
    
    ## 5. Update Var_epsilon using Makalic-Schmidt half-Cauchy
    ## sigma_epsilon_j ~ C+(0, A_epsilon)
    for (j in 1:J) {
      indBeta <- Beta[, m]
      indBeta[-j] <- 0 
      
      resid_j <- W[, j] - Mu[, m] - indBeta
      
      Var_epsilon[j, m] <- 1 / rgamma(
        1,
        shape = (N + 1) / 2,
        rate = 1 / Xi_epsilon[j, m - 1] + sum(resid_j^2) / 2
      )
      
      Xi_epsilon[j, m] <- 1 / rgamma(
        1,
        shape = 1,
        rate = 1 / A_epsilon^2 + 1 / Var_epsilon[j, m]
      )
    }
    
    ## 6. Update W
    for (j in 1:J) {
      for (i in sample(N)) {   # random update order per sweep (Li et al., 2018)
        
        if (r[i, j] == 1) {
          upper <- 100    # end bounds +/-100, as in the reported analyses
          lower <- W[which(r[, j] == 2), j]
          
        } else if (r[i, j] == n) {
          upper <- W[which(r[, j] == n - 1), j]
          lower <- -100
          
        } else {
          lower <- W[which(r[, j] == r[i, j] + 1), j]
          upper <- W[which(r[, j] == r[i, j] - 1), j]
        }
        
        W[i, j] <- truncnorm::rtruncnorm(
          1,
          a = lower,
          b = upper,
          mean = Mu[i, m] + Beta[j, m] * (i == j),
          sd = sqrt(Var_epsilon[j, m])
        )
      }
    }
  }
  
  drop <- 1:burnin
  
  post.mean.mus <- apply(Mu[, -drop], 1, mean)
  post.mean.betas <- apply(Beta[, -drop], 1, mean)
  post.median.Var_epsilon <- apply(Var_epsilon[, -drop], 1, median)
  
  post.mean.mu_beta <- mean(Mu_beta[, -drop])
  post.median.Var_beta <- median(Var_beta[, -drop])
  
  out <- list(
    mu = Mu[, -drop], 
    Beta_draws = Beta[, -drop],
    muBeta_draws = Mu_beta[, -drop],
    varBeta_draws = Var_beta[, -drop],
    varEpsilon_draws = Var_epsilon[, -drop],
    xiBeta_draws = Xi_beta[, -drop],
    xiEpsilon_draws = Xi_epsilon[, -drop],
    W = W,
    post_mean_mus = post.mean.mus,
    post_mean_betas = post.mean.betas, 
    post_mean_mubeta = post.mean.mu_beta, 
    post_median_Varbeta = post.median.Var_beta,
    post_median_VarEps = post.median.Var_epsilon,
    post_rank = rank(-post.mean.mus)
  )
  
  return(out)
}

#------------------Gibbs Sampler Simulation only report estimates-------------
BayeSPeer.HC <- function(n, r, M, burnin, 
                             W,
                             mu0,
                             mubeta0,
                             beta0,
                             Var_beta0,
                             Var_epsilon0,
                             minbeta, maxbeta,
                             A_beta = 2.5,
                             A_epsilon = 2.5) {
  
  N <- n
  J <- n
  
  Mu <- Beta <- matrix(1, N, M)
  Var_epsilon <- matrix(1, N, M)
  Mu_beta <- Var_beta <- matrix(1, 1, M)
  
  ## Makalic-Schmidt auxiliary variables
  Xi_beta <- matrix(1, 1, M)
  Xi_epsilon <- matrix(1, N, M)
  
  Mu[, 1] <- mu0
  Mu_beta[1, 1] <- mubeta0
  Var_beta[1, 1] <- Var_beta0
  
  beta0 <- rnorm(J, mubeta0, sd = sqrt(Var_beta0))
  Beta[, 1] <- beta0
  
  Var_epsilon[, 1] <- Var_epsilon0
  
  Xi_beta[1, 1] <- 1
  Xi_epsilon[, 1] <- 1
  
  for (m in 2:M) {
    
    ## 1. Update Mu
    for (i in 1:N) {
      indBeta <- Beta[, m - 1]
      indBeta[-i] <- 0 
      
      g <- sum((W[i, ] - indBeta) / Var_epsilon[, m - 1])
      h <- sum(1 / Var_epsilon[, m - 1]) + 1
      
      Mu[i, m] <- rnorm(1, mean = g / h, sd = sqrt(1 / h)) 
    }
    
    
    ## 2. Update Beta
    for (j in 1:J) {
      g <- ((W[j, j] - Mu[j, m]) / Var_epsilon[j, m - 1]) +
        (Mu_beta[m - 1] / Var_beta[m - 1])
      
      h <- (1 / Var_epsilon[j, m - 1]) + (1 / Var_beta[m - 1])
      
      Beta[j, m] <- rnorm(1, mean = g / h, sd = sqrt(1 / h))
    } 
    
    ## 3. Update mu_beta
    Mu_beta[m] <- truncnorm::rtruncnorm(
      1,
      a = minbeta,
      b = maxbeta,
      mean = mean(Beta[, m]),
      sd = sqrt(Var_beta[m - 1] / J)
    )
    
    ## 4. Update Var_beta using Makalic-Schmidt half-Cauchy
    ## sigma_beta ~ C+(0, A_beta)
    Var_beta[m] <- 1 / rgamma(
      1,
      shape = (J + 1) / 2,
      rate = 1 / Xi_beta[m - 1] + 
        sum((Beta[, m] - Mu_beta[m])^2) / 2
    )
    
    Xi_beta[m] <- 1 / rgamma(
      1,
      shape = 1,
      rate = 1 / A_beta^2 + 1 / Var_beta[m]
    )
    
    ## 5. Update Var_epsilon using Makalic-Schmidt half-Cauchy
    ## sigma_epsilon_j ~ C+(0, A_epsilon)
    for (j in 1:J) {
      indBeta <- Beta[, m]
      indBeta[-j] <- 0 
      
      resid_j <- W[, j] - Mu[, m] - indBeta
      
      Var_epsilon[j, m] <- 1 / rgamma(
        1,
        shape = (N + 1) / 2,
        rate = 1 / Xi_epsilon[j, m - 1] + sum(resid_j^2) / 2
      )
      
      Xi_epsilon[j, m] <- 1 / rgamma(
        1,
        shape = 1,
        rate = 1 / A_epsilon^2 + 1 / Var_epsilon[j, m]
      )
    }
    
    ## 6. Update W
    for (j in 1:J) {
      for (i in sample(N)) {   # random update order per sweep (Li et al., 2018)
        
        if (r[i, j] == 1) {
          upper <- 100    # end bounds +/-100, as in the reported analyses
          lower <- W[which(r[, j] == 2), j]
          
        } else if (r[i, j] == n) {
          upper <- W[which(r[, j] == n - 1), j]
          lower <- -100
          
        } else {
          lower <- W[which(r[, j] == r[i, j] + 1), j]
          upper <- W[which(r[, j] == r[i, j] - 1), j]
        }
        
        W[i, j] <- truncnorm::rtruncnorm(
          1,
          a = lower,
          b = upper,
          mean = Mu[i, m] + Beta[j, m] * (i == j),
          sd = sqrt(Var_epsilon[j, m])
        )
      }
    }
  }
  
  drop <- 1:burnin
  
  post.mean.mus <- apply(Mu[, -drop], 1, mean)
  post.mean.betas <- apply(Beta[, -drop], 1, mean)
  post.median.Var_epsilon <- apply(Var_epsilon[, -drop], 1, median)
  
  post.mean.mu_beta <- mean(Mu_beta[, -drop])
  post.median.Var_beta <- median(Var_beta[, -drop])
  
  out <- list(
    post_mean_mus = post.mean.mus,
    post_mean_betas = post.mean.betas, 
    post_mean_mubeta = post.mean.mu_beta, 
    post_median_Varbeta = post.median.Var_beta,
    post_median_VarEps = post.median.Var_epsilon,
    post_rank = rank(-post.mean.mus)
  )
  
  return(out)
}