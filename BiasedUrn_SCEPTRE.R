#Packages
library(dplyr)
library(MASS)
library(ggplot2)
library(BiasedUrn)
library(sn)

p2 <- 0
p3 <- 0

Big_P <- c()
Big_P2 <- c()

for (m in 1:20){
  #Creation of parameters and data
  num_cells <- 23
  grna_avg_prob = 0.1
  size = 1
  reps = 1
  mean_expression = 5
  grna_effect = 0.6
  zero_inflation = 0
  B = 500
  All_Pvalues = numeric(reps)
  Z <- rnorm(num_cells)
  pi_true <-  1 / (1 + exp(-(grna_avg_prob + 0.9*Z)))
  X <- rbinom(n = num_cells, size = 1, prob = pi_true)
  Y <- rnbinom(n = num_cells,
               size = size,
               mu = exp(log(mean_expression) + 1*Z - grna_effect * X))
  Y = Y * rbinom(n = num_cells,
                 size = 1,
                 prob = 1 - zero_inflation)
  
  #Step One/
  nb_loco_fit = glm(Y ~ Z,
                    family = negative.binomial(theta=1),
                    data = tibble(Y, Z))
  offsets = log(nb_loco_fit$fitted.values)
  
  #Permutation Number
  Permutation <- 5500
  
  #Odds
  n_cases <- sum(X)
  
  # step 3: generate N x n.perm matrix of permuted data sets
  m1 <- c(rep(1, length(X)))
  
  #d_odds
  model <- glm(Y ~ X + Z,
               family = negative.binomial(theta=1),)
  d_odds <- exp (model$linear.predictors)
  z_score_2 <- model$coefficients[2]
  
  #Multivariate Fisher's NonCentral Hypergeometric distributi
  perm_hg <- rMFNCHypergeo(Permutation, m1, n_cases, d_odds)
  
  
  #Step Two
  Y_new <- Y[X==1]
  nb_distilled_fit = glm(
    Y_new ~ 1,
    offset = offsets[X == 1],
    family = negative.binomial(theta=1),
    data = tibble(Y = Y[X == 1])
  )
  z_score <- coefficients(summary(nb_distilled_fit))[1]
  
  #Step Three
  log_reg_fit = glm(X ~ Z, family = "binomial", data = tibble(X, Z))
  Pi_hat <- 1 / (1 + exp(-(coefficients(summary(log_reg_fit))[1] 
                           + coefficients(summary(log_reg_fit))[2]*Z)))
  
  #Step Four
  resampled_zvalues = numeric(B)
  for (b in 1:B) {

    X <- rbinom(n = num_cells, size = 1, prob = Pi_hat)
    Y_old <- Y[X==1]
    nb_distilled_fit = glm(
      Y_old ~ 1,
      offset = offsets[X == 1],
      family = negative.binomial(theta=1),
      data = tibble(Y = Y[X == 1])
    )
    resampled_zvalues[b] = coefficients(summary(nb_distilled_fit))[1]
    #["(Intercept)", "z value"]
  }
  
  #P-value Calculation
  j <- 0
  for (i in 1:B) {
    if (abs(resampled_zvalues[i]) >= abs(z_score)){
      j = j+1
    }
  }
  P_values = (1+j)/(B+1)
  
  Big_P2[m] <- P_values
  
  if (P_values <= 0.05){
    p2 <- p2 + 1
  }
  
  
  #Step Four
  resampled_zvalues_2 = numeric(Permutation)
  for (b in 1:Permutation) {
    
    X <- perm_hg[,b]
    Y_old <- Y[X==1]
    nb_distilled_fit_2 = glm(
      Y_old ~ 1,
      offset = offsets[X == 1],
      family = negative.binomial(theta=1),
      data = tibble(Y = Y[X == 1])
    )
    resampled_zvalues_2[b] = coefficients(summary(nb_distilled_fit_2))[1]
    #["(Intercept)", "z value"]
  }
  
  #P-value Calculation
  j <- 0
  for (i in 1:Permutation) {
    if (abs(resampled_zvalues_2[i]) >= abs(z_score)){
      j = j+1
    }
  }
  P_values_2 = (1+j)/(Permutation+1)
  
  Big_P[m] <- P_values_2
  
  if (P_values_2 <= 0.05){
    p3 <- p3 + 1
  }
}
p2/5000
p3/5000

expected_p <- c()
for (i in 1:20){
  expected_p[i] <- -log((i/21),base = 20)
}
plot(sort(expected_p), sort(Big_P), type='l',col='red')
lines(sort(expected_p), sort(Big_P2),type='l',col='blue')
lines(x=y)



t_nulls <- resampled_zvalues[5001:5500]
t_star <- z_score

#Skew-T Distribution in BiasedUrn
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

#Skew-T Distribution in SCEPTRE
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



#Check uniform distribution of p-values in BiasedUrn
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
