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
  for (i in 1:P) {
    for (k in 1:(K-1)) {
      gamma[i,k] ~ dnorm(0, sd = 100)
    }
  }
  
  # Prior for latent class longitudinal 
  for (i in 1:K) {
    for (j in 1:Q) {
      beta[i,j] ~ dnorm(0, sd = 100)
    }
  }
  sigma_b ~ dunif(0, 10000)
  sigma_M ~ dunif(0, 10000)
  
  #one ~ dconstraint(beta[1,1] > beta[2,1])
  #two ~ dconstraint(beta[2,1] > beta[3,1])
  
  # Prior for Y
  for (i in 1:K) {
    for (j in 1:P) {
      phi[i,j] ~ dnorm(0, sd = 5)
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
  
  logits[1:N1,2:K] <- X[1:N1, 1:P] %*% t(gamma[1:(K-1),1:P])
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
    L[i] ~ dcat(Prob[i, 1:K])
    logit(p[i]) <- (t(X[i, 1:P]) %*% phi[L[i], 1:P])[1,1]
    Y[i] ~ dbin(size = 1, p[i])
    O[i] ~ dnegbin(alpha / (alpha + exp(theta[L[i]])), alpha)
    for (j in s[i]:e[i]) {
      EM[j] <- (beta[L[i], 1])  + Time[j]*beta[L[i], 2] + b0[i] * sigma_b
      M[j] ~ dnorm(EM[j], sd = sigma_M)
    }
  }
  logDens ~ dnorm(0, 1)
})
