library(tidyverse)
# library(nimble) # MCMC
library(nlme) # building linear mixed effects model for initial
library(nnet) # multinomial logistic regression for initial
library(survival)
library(eha)

# Function to simulation data
# n sample size
# X1  independent variable matrix for latent class
# X2  independent variable matrix for survival outcome
# K number of clusters
# gamma coefficient matrix for Multinominal regression
# theta and alpha parameter in informative cluster size
# alpha2 shape parameter for AFT survival model
# phi coefficient matrix for AFT regression
# perc percentage of censoring
# sigma_b0 random effect of subject
# signa_M random error
# beta coefficient matrix for longitindinal regression
sim_data = function(n, X1, X2, K, gamma, theta, alpha, alpha2, phi,
                    perc, sigma_b0, sigma_M, beta) {
  ID = 1:n
  # Class membership
  logits = exp(X1  %*% t(gamma))
  prob = logits/rowSums(logits)
  ## Multinominal regression
  lc = as.matrix(t(apply(prob, 1, function(p) rmultinom(1, size = 1, prob = p))))
  colnames(lc) = 1:K
  lc_cate = lc %>% 
    as_tibble() %>% 
    pivot_longer(cols = "1":as.character(K), 
                 names_to = "lc", values_to = "lc_10") %>% 
    filter(lc_10 != 0) %>% 
    dplyr::select(lc)
  # Informative cluster size
  obs_n = apply(alpha/(alpha + exp(lc %*% theta)), 1, function(prob) rnbinom(n = 1, size = alpha, prob = prob)) + 1 
  max_k = max(obs_n)
  # Longitudinal Mediator time from 1 to max_k (ensure # of time points is greater than cluster size)
  Time_cov = matrix(rep(1:max_k, n),
                    ncol = max_k, 
                    nrow = n, 
                    byrow = T)
  
  M = matrix(rep(lc %*% beta[,1], max_k),ncol = max_k,nrow = n) + 
    Time_cov*matrix(rep(lc %*% beta[,2], max_k),ncol = max_k,nrow = n) +
    # Random effect
    matrix(rep(rnorm(n)*sigma_b0, max_k), ncol = max_k,nrow = n) + 
    # Random error
    matrix(rnorm(n*max_k)*sigma_M, ncol = max_k,nrow = n)
  
  # Outcome 
  mu = rowSums(X2 * (lc %*% phi))
  ## scale for survival time
  lambda = exp(mu)/((log(2))^(1/alpha2))
  ## scale for censor time
  lambda2 = (1/perc-1)^(1/alpha2)*lambda
  
  t = rweibull(n, shape=alpha2, scale = lambda)
  cent = rweibull(n, shape=alpha2, scale = lambda2)
  
  censured = t > cent
  delta = as.logical(1-censured)
  ob_t = t
  ob_t[censured==1] = cent[censured==1]
  t[censured==1] = NA
  cent[censured==0] = Inf
  
  
  
  # Create missing pattern based on the cluster size
  mask = matrix(NA,nrow = n, ncol = max_k)
  for (i in 1:n) {
    mask[i,] = sample(c(rep(F, obs_n[i]), 
                        rep(T, max_k-obs_n[i])), replace = F)
  }
  M[mask] = NA
  colnames(M) = 1:max_k
  
  # Observed data
  d_long =  cbind(ID, lc_cate, A, t, cent, ob_t, censured, delta, obs_n, M) %>%
    as_tibble() %>% 
    pivot_longer(cols = "1":as.character(max_k), 
                 names_to = "time",
                 values_to = "M") %>% 
    filter(!is.na(M)) %>% 
    mutate(M = as.numeric(M),
           time = as.numeric(time)) 
  
  return(d_long)
}

# log likelihood of negative binomial
loglk = function(parameters, obs_n_vec) {
  alpha = parameters[1]
  theta = parameters[2]
  prob = alpha/(alpha+exp(theta))
  sum(dnbinom(obs_n_vec - 1, alpha, prob, log = T))
}


initial_values = function(d_long, K) {
  d_subject = d_long %>% dplyr::select(-time, -M) %>% distinct()
  
  # cluster membership based on K means on coefficients (random + fixed) in lme
  fit.nmle = lme(fixed = M ~ time, random = ~ (1) | ID, 
                 control = list(msMaxIter = 1000, msMaxEval = 1000),
                 data = d_long)
  B = coef(fit.nmle)
  lc_init = kmeans(B,K)$cluster
  lc_init = match(lc_init,unique(lc_init))
  # reorder class based on the percentage
  freq = table(lc_init)
  sorted_f = sort(freq, decreasing = F)
  value_order = as.numeric(names(sorted_f))
  lc_init = match(lc_init, value_order)
  
  # add latent class to the dataset
  ID = unique(d_long$ID)
  id_lc = cbind(ID, lc_init) %>% as_tibble()
  d_long_int = d_long %>% left_join(id_lc)
  d_subject_int = d_long_int %>% dplyr::select(-time, -M) %>% distinct()
  
  # ggplot(d_long_int,
  #        aes(x = time, y = M, group = ID, color = as.factor(lc_init))) +
  #   geom_line(alpha = 0.2) +
  #   labs(x = "Time", y = "Measurement") +
  #   theme_minimal()
  # 
  # ggplot(d_long, 
  #        aes(x = time, y = M, group = ID)) +
  #   geom_line(alpha = 0.2) +
  #   labs(x = "Time", y = "Measurement") +
  #   theme_minimal()
  
  # initial value for gamma
  fit.mult = multinom(as.factor(lc_init) ~ A, data = d_subject_int)
  gamma_init = as.matrix(coef(fit.mult))
  
  # initial value for longitidinal
  beta_init = matrix(0, ncol = 2, nrow = K)
  sigma_b0_init = rep(0,K)
  sigma_M_init = rep(0,K)
  for(i in 1:K) {
    temp_data = d_long_int %>% filter(lc_init == i)
    temp_model = lme4::lmer(M ~ time + (1|ID), data = temp_data)
    sum = summary(temp_model)
    beta_init[i, 1:2] = sum$coefficients[1:2,1]
    temp_cov = as.data.frame(VarCorr(temp_model))
    sigma_b0_init[i] = temp_cov[1,'sdcor']
    sigma_M_init[i] = temp_cov[2,'sdcor']
  }
  sigma_M_init = mean(sigma_M_init)
  sigma_b0_init = mean(sigma_b0_init)
  
  phi_init = matrix(0, ncol = 2, nrow = K)
  sigma_Y_init = rep(0,K)
  
  # random effect
  b0_init = rep(0,length(ID))
  # Outcome
  for(i in 1:K) {
    temp_data = d_subject_int %>% filter(lc_init == i)
    temp_model = aftreg(Surv(ob_t, delta) ~ A,
                        data = temp_data,
                        dist = "weibull")
    phi_init[i, 2] = - coef(temp_model) [1]
    phi_init[i, 1] = coef(temp_model) [2]
  }
  alpha2_init = 1
  
  # Initial value for cluster size model
  alpha_init = rep(0,K)
  theta_init = rep(0,K)
  for (i in 1:K) {
    temp_data = d_subject_int %>% filter(lc_init == i)
    optim_re = optim(
      par = c(1, 1),  
      fn = loglk,  
      obs_n_vec = temp_data$obs_n,  
      control = list(fnscale = -1)
    )
    alpha_init[i] = optim_re$par[1]
    theta_init[i] = optim_re$par[2]
  }
  alpha_init = 30
  
  return(list(lc_init, gamma_init, 
              b0_init, sigma_b0_init, sigma_M_init, beta_init, phi_init, alpha2_init,
              alpha_init, theta_init))
}


# for no confounder (for logT) 
mediation_effect = function(gamma, phi, X1, X2, K) {
  # P(U = c | A = 0)
  X1_A0 = X1
  X1_A1 = X1
  X1_A0[, 2] = 0
  X1_A1[, 2] = 1
  
  X2_A0 = X2
  X2_A1 = X2
  X2_A0[, 2] = 0
  X2_A1[, 2] = 1
  
  logits0 = exp(X1_A0  %*% t(gamma))
  prob0 = logits0/rowSums(logits0)
  
  # P(U = c | A = 1)
  logits1 = exp(X1_A1  %*% t(gamma))
  prob1 = logits1/rowSums(logits1)
  
  # NDE (no interaction)
  NDE = mean(rowSums(phi[,2]*prob0))
  
  # NIE
  NIE = mean(rowSums((X2_A1 %*% t(phi))*(prob1-prob0)))
  
  # TE
  TE = NDE + NIE
  prop = NIE/TE
  return(c(NDE, NIE, TE, prop))
}

mediation_summary = function(gamma_mat, phi_mat, X1, X2, K) {
  zeros = rep(0, ncol(X1))
  for (i in 1:nrow(gamma_mat)) {
    gamma = matrix(gamma_mat[i, ], nrow = K - 1, byrow = F)
    gamma = rbind(zeros, gamma)
    phi = matrix(phi_mat[i, ], nrow = K, byrow = F)
    if (!exists("re")) {
      re = mediation_effect(gamma, phi, X1, X2, K) 
    } else {
      re = rbind(re, mediation_effect(gamma, phi, X1, X2, K))
    }
  }
  return(re)
}











initial_values2 = function(d_long, K) {
  d_subject = d_long %>% dplyr::select(-time, -M) %>% distinct()
  
  # cluster membership based on K means on coefficients (random + fixed) in lme
  fit.nmle = lme(fixed = M ~ time, random = ~ (1 + time) | ID, data = d_long)
  B = coef(fit.nmle)
  lc_init = kmeans(B,K)$cluster
  lc_init = match(lc_init,unique(lc_init))
  # reorder class based on the percentage
  freq = table(lc_init)
  sorted_f = sort(freq, decreasing = F)
  value_order = as.numeric(names(sorted_f))
  lc_init = match(lc_init, value_order)
  
  # add latent class to the dataset
  ID = unique(d_long$ID)
  id_lc = cbind(ID, lc_init) %>% as_tibble()
  d_long_int = d_long %>% left_join(id_lc)
  d_subject_int = d_long_int %>% dplyr::select(-time, -M) %>% distinct()
  
  ggplot(d_long_int, 
         aes(x = time, y = M, group = ID, color = as.factor(lc_init))) +
    geom_line(alpha = 0.2) +
    geom_point(alpha = 0.2) +
    labs(x = "Time", y = "Measurement") +
    theme_minimal() +
    facet_grid(.~delta)
  
  
  # initial value for gamma
  fit.mult = multinom(as.factor(lc_init) ~ A + age_baseline + EDUC + SEX, data = d_subject_int)
  gamma_init = as.matrix(coef(fit.mult))
  
  # initial value for longitidinal
  beta_init = matrix(0, ncol = 2, nrow = K)
  sigma_b0_init = rep(0,K)
  sigma_M_init = rep(0,K)
  for(i in 1:K) {
    temp_data = d_long_int %>% filter(lc_init == i)
    temp_model = lme4::lmer(M ~ time + (1|ID), data = temp_data)
    sum = summary(temp_model)
    beta_init[i, 1:2] = sum$coefficients[1:2,1]
    temp_cov = as.data.frame(VarCorr(temp_model))
    sigma_b0_init[i] = temp_cov[1,'sdcor']
    sigma_M_init[i] = temp_cov[2,'sdcor']
  }
  sigma_M_init = mean(sigma_M_init)
  sigma_b0_init = mean(sigma_b0_init)
  # random effect
  b0_init = rep(0,length(ID))
  
  
  phi_init = matrix(0, ncol = 5, nrow = K)
  
  # Outcome
  for(i in 1:K) {
    temp_data = d_subject_int %>% filter(lc_init == i)
    temp_model = aftreg(Surv(ob_t, delta) ~ A + age_baseline + EDUC + SEX,
                        data = temp_data,
                        dist = "weibull")
    phi_init[i, 2:5] = - coef(temp_model) [1:4]
    phi_init[i, 1] = coef(temp_model) [5]
  }
  alpha2_init = 1
  
  # Initial value for cluster size model
  alpha_init = rep(0,K)
  theta_init = rep(0,K)
  for (i in 1:K) {
    temp_data = d_subject_int %>% filter(lc_init == i)
    optim_re = optim(
      par = c(1, 1),  
      fn = loglk,  
      obs_n_vec = temp_data$obs_n,  
      control = list(fnscale = -1)
    )
    alpha_init[i] = optim_re$par[1]
    theta_init[i] = optim_re$par[2]
  }
  alpha_init = mean(alpha_init)
  
  return(list(lc_init, gamma_init, 
              b0_init, sigma_b0_init, sigma_M_init, beta_init, phi_init, alpha2_init,
              alpha_init, theta_init))
}


