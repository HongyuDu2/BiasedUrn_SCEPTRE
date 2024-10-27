# Packages
library(dplyr)
library(MASS)
library(ggplot2)
library(BiasedUrn)
library(sn)

# p-value in SCEOTRE
p2 <- 0
# p-value in BiasedUrn
p3 <- 0

# record all p-values in SCEPTRE and BiasedUrn
Big_P <- c()
Big_P2 <- c()

for (m in 1:20){ ## 20 is the number of simulations, change to any number you want
  # create the simulation dataset
  num_cells <- 23 ## length of RNA sequence, change to any number you want
  grna_avg_prob = 0.1
  size = 1
  reps = 1
  mean_expression = 5
  grna_effect = 0.6 ## coefficient of guide-rna
  zero_inflation = 0
  B = 500 ## number of permutation
  All_Pvalues = numeric(reps)
  Z <- rnorm(num_cells) ## generate confounding 
  pi_true <-  1 / (1 + exp(-(grna_avg_prob + 0.9*Z)))
  X <- rbinom(n = num_cells, size = 1, prob = pi_true) ## generate gudie-rna dataset
  Y <- rnbinom(n = num_cells,
               size = size,
               mu = exp(log(mean_expression) + 1*Z - grna_effect * X))
  Y = Y * rbinom(n = num_cells,
                 size = 1,
                 prob = 1 - zero_inflation) ## generate gene expression
  
  # generate the z_score or coefficient of guide-rna for SCEPTRE
  nb_loco_fit = glm(Y ~ Z,
                    family = negative.binomial(theta=1),
                    data = tibble(Y, Z))
  offsets = log(nb_loco_fit$fitted.values)

  model <- glm(Y ~ X + Z,
               family = negative.binomial(theta=1),)
  d_odds <- exp (model$linear.predictors)
  z_score_2 <- model$coefficients[2]
  
  # number of permutation
  Permutation <- 5500
  
  # X is a vector consist of 0 or 1
  # calculate the number of 1s in X
  n_cases <- sum(X)
  
  # generate N x n.perm matrix of permuted data sets
  m1 <- c(rep(1, length(X)))
  
  # generate multivariate Fisher's nonCentral hypergeometric distribution
  perm_hg <- rMFNCHypergeo(Permutation, m1, n_cases, d_odds)
  
  
  # generate the z_score or coefficient of guide-rna for BiasedUrn
  Y_new <- Y[X==1]
  nb_distilled_fit = glm(
    Y_new ~ 1,
    offset = offsets[X == 1],
    family = negative.binomial(theta=1),
    data = tibble(Y = Y[X == 1])
  )
  z_score <- coefficients(summary(nb_distilled_fit))[1]
  
  # regenerate the pi_hat parameter based on the algorithm of SCEPTRE
  log_reg_fit = glm(X ~ Z, family = "binomial", data = tibble(X, Z))
  Pi_hat <- 1 / (1 + exp(-(coefficients(summary(log_reg_fit))[1] 
                           + coefficients(summary(log_reg_fit))[2]*Z)))
  
  # generate new z_score
  resampled_zvalues = numeric(B)
  for (b in 1:B) { ## B is the number of new z_score
    X <- rbinom(n = num_cells, size = 1, prob = Pi_hat) ## regenerate the guide_rna dataset based on new pi_hat
    Y_old <- Y[X==1]
    nb_distilled_fit = glm(
      Y_old ~ 1,
      offset = offsets[X == 1],
      family = negative.binomial(theta=1),
      data = tibble(Y = Y[X == 1])
    )
    resampled_zvalues[b] = coefficients(summary(nb_distilled_fit))[1] ## store new z_score datasets in a vector
    #["(Intercept)", "z value"]
  }
  
  # p-value calculation
  j <- 0
  for (i in 1:B) {
    if (abs(resampled_zvalues[i]) >= abs(z_score)){ ## record the number of new z_score which is larger than the old one
      j = j+1
    }
  }
  P_values = (1+j)/(B+1) ## p-value is calculated through the percentage of new z_score which is larger than the old on in the new z_score
  
  Big_P2[m] <- P_values

  # record the p-value if the test is significant
  if (P_values <= 0.05){
    p2 <- p2 + 1
  }
  
  
  # generate new z_score
  resampled_zvalues_2 = numeric(Permutation)
  for (b in 1:Permutation) { ## permutation is the number of new z_score
    X <- perm_hg[,b] ## regenerate the guide_rna dataset based on BiasedUrn method
    Y_old <- Y[X==1]
    nb_distilled_fit_2 = glm(
      Y_old ~ 1,
      offset = offsets[X == 1],
      family = negative.binomial(theta=1),
      data = tibble(Y = Y[X == 1])
    )
    resampled_zvalues_2[b] = coefficients(summary(nb_distilled_fit_2))[1] ## store new z_score datasets in a vector
    #["(Intercept)", "z value"]
  }
  
  # p-value calculation
  j <- 0
  for (i in 1:Permutation) { ## record the number of new z_score which is larger than the old one
    if (abs(resampled_zvalues_2[i]) >= abs(z_score)){
      j = j+1
    }
  }
  P_values_2 = (1+j)/(Permutation+1) ## p-value is calculated through the percentage of new z_score which is larger than the old on in the new z_score
  
  Big_P[m] <- P_values_2

  # record the p-value if the test is significant
  if (P_values_2 <= 0.05){
    p3 <- p3 + 1
  }
}

# power calculation
# the percentage of p-values which is smaller than 0.05 or the percentage of tests significant
# in general, the power of BiasedUrn is larger than the power of SCEPTRE
p2/5000
p3/5000

# plot the power curves
expected_p <- c()
for (i in 1:20){
  expected_p[i] <- -log((i/21),base = 20)
}
plot(sort(expected_p), sort(Big_P), type='l',col='red')
lines(sort(expected_p), sort(Big_P2),type='l',col='blue')
lines(x=y)


t_nulls <- resampled_zvalues[5001:5500]
t_star <- z_score

# kkew-T distribution in BiasedUrn
fit_skew_t <- function(t_nulls, t_star, side) {
  p_value_skew_t <- NA
  skew_t_fit <- tryCatch(sn::selm(t_nulls ~ 1, family = "ST"), error = function(e) return(NA))
  if (class(skew_t_fit) == "selm") { # If the fit worked,
    dp <- skew_t_fit@param$dp # then extract the parameters.
    if (!any(is.na(dp))) { # If all the fitted parameters are numbers,
      p_value_skew_t <- switch(side,
                               'left' = pmax(.Machine$double.eps, sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)), # then compute the skew t-based p-value. pst(x = t_star, dp = dp)
                               'right' = pmax(.Machine$double.eps, 1 - sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
                               'both' = pmax(.Machine$double.eps, sn::pst(x = -abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps) +
                                               (1 - sn::pst(x = abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps)))
      )
    }
  }
  # check if the skew-t fit worked
  skew_t_fit_success <- !is.na(p_value_skew_t)
  if (skew_t_fit_success) {
    out_p <- p_value_skew_t
    skew_t_mle <- dp
  } else {
    out_p <- switch(side,
                    'left' = mean(c(-Inf, t_nulls) <= t_star),
                    'right' = mean(c(Inf, t_nulls) >= t_star),
                    'both' = mean(c(Inf, abs(t_nulls)) >= abs(t_star)))
    skew_t_mle <- c(xi = NA, omega = NA, alpha = NA, nu = NA)
  }
  return(list(skew_t_fit_success = skew_t_fit_success, out_p = out_p, skew_t_mle = skew_t_mle))
}

fit_skew_t(t_nulls, t_star, 'left')
fit_skew_t(t_nulls, t_star, 'right')
SCEPTRE_t <- fit_skew_t(t_nulls, t_star, 'both')


t_nulls <- resampled_zvalues_2[5001:5500]
t_star <- z_score

# skew-T distribution in SCEPTRE
fit_skew_t <- function(t_nulls, t_star, side) {
  p_value_skew_t <- NA
  skew_t_fit <- tryCatch(sn::selm(t_nulls ~ 1, family = "ST"), error = function(e) return(NA))
  if (class(skew_t_fit) == "selm") { # If the fit worked,
    dp <- skew_t_fit@param$dp # then extract the parameters.
    if (!any(is.na(dp))) { # If all the fitted parameters are numbers,
      p_value_skew_t <- switch(side,
                               'left' = pmax(.Machine$double.eps, sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)), # then compute the skew t-based p-value. pst(x = t_star, dp = dp)
                               'right' = pmax(.Machine$double.eps, 1 - sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
                               'both' = pmax(.Machine$double.eps, sn::pst(x = -abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps) +
                                               (1 - sn::pst(x = abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps)))
      )
    }
  }
  # check if the skew-t fit worked
  skew_t_fit_success <- !is.na(p_value_skew_t)
  if (skew_t_fit_success) {
    out_p <- p_value_skew_t
    skew_t_mle <- dp
  } else {
    out_p <- switch(side,
                    'left' = mean(c(-Inf, t_nulls) <= t_star),
                    'right' = mean(c(Inf, t_nulls) >= t_star),
                    'both' = mean(c(Inf, abs(t_nulls)) >= abs(t_star)))
    skew_t_mle <- c(xi = NA, omega = NA, alpha = NA, nu = NA)
  }
  return(list(skew_t_fit_success = skew_t_fit_success, out_p = out_p, skew_t_mle = skew_t_mle))
}

fit_skew_t(t_nulls, t_star, 'left')
fit_skew_t(t_nulls, t_star, 'right')
fit_skew_t(t_nulls, t_star, 'both')

t_nulls <- rnorm(500,0,0.2)
t_star <- resampled_zvalues_2

                         
# check uniform distribution of p-values in BiasedUrn
km <- c()

for (i in 1:500){
  km[i] <- (i)/(500+1)
}

plot(km, sort(Big_P))
hist(Big_P)

simpleQQPlot = function (observedPValues) {
  plot(-log10(1:length(observedPValues)/length(observedPValues)), 
       -log10(sort(observedPValues)))
  abline(0, 1, col = "red")
}
