library(tidyverse)
library(MCMCvis)
library(nimble) 

code <- nimbleCode({
  mu ~ dnorm(0, sd = 100)
  alpha2 ~ dunif(0.01,100)
  lambda <- log(2)*exp(-mu*alpha2) 
  for (i in 1:N) {
    censured[i] ~ dinterval(t[i], cent[i])
    t[i] ~ dweib(alpha2, lambda)
  }
})




n = 1000
alpha2 = 0.8
mu = 2
perc = 0.1
## scale for survival time
lambda = exp(mu)/((log(2))^(1/alpha2))
## scale for censor time
lambda2 = (1/perc-1)^(1/alpha2)*lambda

t = rweibull(n, shape = alpha2, scale = lambda)
cent = rweibull(n, shape = alpha2, scale = lambda2)

censured = t > cent
delta = as.logical(1-censured)
ob_t = t
ob_t[censured==1] = cent[censured==1]
t[censured==1] = NA
cent[censured==0] = Inf
d = cbind(ob_t, delta)

Consts <- list(N = n)

Data <- list(censured = censured,
             t = t, cent = cent)

Inits <- list(mu = 0,
              alpha2 = 1)

n.burnin = 50000


mcmc.out <- nimbleMCMC(code = code, 
                       constants = Consts,
                       data = Data, 
                       inits = Inits,
                       monitors = c("mu", "alpha2"), nchains = 1, setSeed = TRUE)

MCMCsummary(object = mcmc.out, round = 2)









