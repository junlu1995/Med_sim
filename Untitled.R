# Rjags code
library(rjags)

model_code <- "
model {
  # fixed constant 
  # N1 number of subjects
  # N2 number of observations
  # K number of latent category
  # P number of parameters in gamma
  # Q number of parameters in beta
  # logits for Latent  
  
  # Prior membership
  for (p in 1:P) {
    for (k in 1:(K-1)) {
      gamma[k,p] ~ dnorm(0, 0.0001)
    }
  }
  
  # Prior for latent class longitudinal 
  for (k in 1:K) {
    for (q in 1:Q) {
      beta[k, q] ~ dnorm(0, 0.0001)
    }
  }
  sigma_b ~ dunif(0, 10000)
  sigma_M ~ dunif(0, 10000)
  
  # Prior for Survival Outcome
  alpha2 ~ dunif(0.01,100)
  
  for (k in 1:K) {
    for (p in 1:P) {
      phi[k, p] ~ dnorm(0, 0.0001)
    }
  }
  
  # Random effect
  for (i in 1:N1){
    b0[i] ~ dnorm(0, 1)
  }
  
  # Prior for N
  for (i in 1:K) {
    theta[i] ~ dnorm(0, 0.0001)
  }
  alpha ~ dunif(0, 50000)
  
  for (i in 1:N1) {
    for (k in 2:K) {
      logits[i,k] <- inprod(X1[i,], gamma[k-1,])
    }
    logits[i,1] <- 0
  }
  
  # Softmax transformation
  for (i in 1:N1) {
    for (k in 1:K) {
      soft[i, k] <- exp(logits[i,k])
    }
    for (k in 1:K) {
      Prob[i, k] <- soft[i,k]/sum(soft[i,])
    }
  }
  
  # Sampling procedure
  for (i in 1:N1) {
    # Sampling latent cluster
    L[i] ~ dcat(Prob[i, ])
    
    # Sampling survival outcome
    mu[i] <- inprod(X2[i,], phi[L[i],])
    censured[i] ~ dinterval(t[i], cent[i])
    lambda[i] <- log(2)*exp(-mu[i]*alpha2) 
    t[i] ~ dweib(alpha2, lambda[i])
    
    # Sampling observation number
    O[i] ~ dnegbin(alpha / (alpha + exp(theta[L[i]])), alpha)
    
    # Sampling longitudinal outcome
    for (j in s[i]:e[i]) {
      EM[j] <- beta[L[i], 1] + Time[j]*beta[L[i], 2] + b0[i] * sigma_b
      M[j] ~ dnorm(EM[j], sigma_M^-2)
    }
  }
}
"

writeLines(model_code, con = "model.jags")
d_long = sim_data(n, X1, X2, K, gamma, theta, alpha, alpha2, phi,
                  perc, sigma_b0, sigma_M, beta) 

intital = initial_values(d_long, 3)
d_subject = d_long %>% select(-time, -M) %>% distinct()
n = nrow(d_subject)
X1 = cbind(rep(1, n), d_subject$A)
obs_n = d_subject$obs_n
time = d_long$time
s = rep(0, n)
e = rep(0, n)
# Some issue with writing loop through longitidinal M in NIMBLE
s[1] = 1
for (i in 1:n) {
  e[i] = s[i] + obs_n[i]-1
  if(i != n) {
    s[i+1] = e[i]+1
  }
}

Data <- list(N1 = n,
             X1 = X1,
             X2 = X1,
             Q = 2,
             P = 2,
             K = K,
             s = s,
             e = e,
             Time = d_long$time,
             M = d_long$M, O = d_subject$obs_n - 1, 
             censured = d_subject$censured,
             t = d_subject$t, cent = d_subject$cent)

Inits <- list(L = intital[[1]],
              gamma = intital[[2]],
              b0 = intital[[3]],
              sigma_b = intital[[4]],
              sigma_M = intital[[5]], 
              beta = intital[[6]],
              phi = intital[[7]],
              alpha2 = intital[[8]],
              alpha = intital[[9]],
              theta = intital[[10]])


jags_model <- jags.model(file = "model.jags",
                         data = Data, 
                         inits = Inits,
                         n.chains = 1, 
                         n.adapt = 1000)

update(jags_model, 10000)

samples <- coda.samples(model = jags_model, variable.names = c("gamma", "beta", "sigma_b", "sigma_M", "alpha2", "phi", "theta", "alpha"), n.iter = 1000)
summary(samples)