# This model modifies the original BUGS example (Example Vol I): 
# Dogs: loglinear model for binary data.
# Original model is modified to allow for dog-specific parameters (alpha and beta),
# and these parameters are modeled hierarchically, producing the HB model.
# We'll assume that this model is the focal model and that we wish to make 
# inferences about dogs such that we computed Ndog*(Ndog-1)/2 pairwise comparisons
# for each of the two parameters (alpha and beta)

data{
  for(i in 1:Ndog){
    for(t in 1:Ntrial){
      # For the probability model of shock, redefine Y=0 as y=1, dog gets shocked
      y[i,t] <- 1-Y[i,t]
    }
  }
}
model{
  for(i in 1:Ndog){
    p[i,1] <- 0
    lpd[i,1] <- 0
    for(t in 2:Ntrial){
      # Specify model for p, probability of a shock at trial t on log scale
      # (xa is the number of avoidances before trial t, and xs is the number of previous shocks.)
      log(p[i,t]) <- alpha[i]*xa[i,t] + beta[i]*xs[i,t]
      # Bernoulli likelihood for failure (got shocked, y = 1)
      y[i,t] ~ dbern(p[i,t])
      # log pointwise density
      lpd[i,t] <- log(p[i,t])*y[i,t] + log(1-p[i,t])*(1-y[i,t])
    }
    
    # Heirarchical priors for log-scale parameters, and back-transform
    # to original parameter scale, thus obeying domain for alpha and beta
    log.alpha[i] ~ dnorm(mu.log.alpha, tau.log.alpha)
    log.beta[i] ~ dnorm(mu.log.beta, tau.log.beta)
    alpha[i] <- -exp(log.alpha[i])
    beta[i] <- -exp(log.beta[i])
  }
  
  # Vague priors, log scale
  mu.log.alpha ~ dnorm(0, 0.0001)
  mu.log.beta ~ dnorm(0, 0.0001)
  tau.log.alpha ~ dgamma(0.1, 0.1)
  tau.log.beta ~ dgamma(0.1, 0.1)
  mu.alpha <- -exp(mu.log.alpha)
  mu.beta <- -exp(mu.log.beta)
  
  # Pairwise comparisons among dogs of alpha and beta:
  for(j in 1:Ncomp){
    alpha.diff[j] <- alpha[pairID[j,1]] - alpha[pairID[j,2]]
    beta.diff[j] <- beta[pairID[j,1]] - beta[pairID[j,2]]
    # Posterior means of the quantities below give Bayesian 
    # p-values for each pairwise comparison
    p.alpha[j] <- step(alpha.diff[j])
    p.beta[j] <- step(beta.diff[j])
  }
}