# Non-hierarchical version of Dog learning model with dog-specific parameters
# See the HB model version for additional comments/explanations.

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
      log(p[i,t]) <- alpha[i]*xa[i,t] + beta[i]*xs[i,t]
      # Bernoulli likelihood for failure (got shocked, y = 1)
      y[i,t] ~ dbern(p[i,t])   
      # log pointwise density
      lpd[i,t] <- log(p[i,t])*y[i,t] + log(1-p[i,t])*(1-y[i,t])      
    }
    
    # Vague, non-hierarchical priors for log-scale, dog-level parameters; back-
    # transform to original parameter scale, thus obeying domain for alpha and beta
    log.alpha[i] ~ dnorm(0, 0.0001)
    log.beta[i] ~ dnorm(0, 0.0001)
    alpha[i] <- -exp(log.alpha[i])
    beta[i] <- -exp(log.beta[i])
  }
}