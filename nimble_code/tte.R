sumLogPostDens <- nimbleFunction(
  name = 'sumLogPostDens',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    stochNodes <- setdiff(model$getNodeNames(stochOnly = TRUE), target)
  },
  run = function() {
    model[[target]] <<- model$getLogProb(stochNodes)
  },
  methods = list( reset = function() {} )
)


code <- nimbleCode({
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
      gamma[k,p] ~ dnorm(0, sd = 100)
    }
  }
  
  # Prior for latent class longitudinal 
  for (k in 1:K) {
    for (q in 1:Q) {
      beta[k, q] ~ dnorm(0, sd = 100)
    }
  }
  sigma_b ~ dunif(0, 10000)
  sigma_M ~ dunif(0, 10000)
  
  #one ~ dconstraint(beta[1,1] > beta[2,1])
  #two ~ dconstraint(beta[2,1] > beta[3,1])
  
  # Prior for Survival Outcome
  alpha2 ~ dunif(0.01,100)
  
  for (k in 1:K) {
    for (p in 1:P) {
      phi[k, p] ~ dnorm(0, sd = 100)
    }
  }
  
  # Random effect
  for (i in 1:N1){
    b0[i] ~ dnorm(0, sd = 1)
  }
  
  # Prior for N
  for (i in 1:K) {
    theta[i] ~ dnorm(0, sd = 100)
  }
  alpha ~ dunif(0, 50000)
  
  logits[1:N1,2:K] <- X1[1:N1, 1:P] %*% t(gamma[1:(K-1),1:P])
  # Reference set at 0
  logits[1:N1,1] <- rep(0, N1)
  
  # Softmax transformation
  for (i in 1:N1) {
    for (k in 1:K) {
      soft[i, k] <- exp(logits[i,k])
    }
  }
  for (i in 1:N1) {
    for (k in 1:K) {
      Prob[i, k] <- soft[i,k]/sum(soft[i, 1:K]) 
    }
  }
  
  # Sampling procedure
  for (i in 1:N1) {
    # Sampling latent cluster
    L[i] ~ dcat(Prob[i, 1:K])
    
    # Sampling survival outcome
    mu[i] <- (t(X2[i, 1:P]) %*% phi[L[i], 1:P])[1,1]
    censured[i] ~ dinterval(t[i], cent[i])
    lambda[i] <- log(2)*exp(-mu[i]*alpha2) 
    t[i] ~ dweib(alpha2, lambda[i])
    
    # Sampling observation number
    O[i] ~ dnegbin(alpha / (alpha + exp(theta[L[i]])), alpha)
    
    # Sampling longitindinal outcome
    for (j in s[i]:e[i]) {
      EM[j] <- (beta[L[i], 1])  + Time[j]*beta[L[i], 2] + b0[i] * sigma_b
      M[j] ~ dnorm(EM[j], sd = sigma_M)
    }
  }
  # logDens ~ dnorm(0, 1)
})
