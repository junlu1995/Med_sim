code <- nimbleCode({
  mu[i] <- (t(X[i, 1:P]) %*% phi[L[i], 1:P])[1,1]
  censured[i] ~ dinterval(t[i], cent[i])
  lambda[i] <- log(2)*exp(-mu[i]*alpha2) 
  t[i] ~ dweib(alpha2, lambda[i])
  
})