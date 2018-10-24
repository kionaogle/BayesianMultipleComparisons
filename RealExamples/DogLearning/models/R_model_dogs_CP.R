# Complete pooling version of the Dog Learning model.
# This model is very similar to the original BUGS example (Example Vol I): 
# Dogs: loglinear model for binary data. But, slight modifications to priors.
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
      log(p[i,t]) <- alpha*xa[i,t] + beta*xs[i,t]
      # Bernoulli likelihood for failure (got shocked, y = 1)
      y[i,t] ~ dbern(p[i,t])
      # log pointwise density
      lpd[i,t] <- log(p[i,t])*y[i,t] + log(1-p[i,t])*(1-y[i,t])      
    }
  }
  
  # Vague priors, log scale
  log.alpha ~ dnorm(0, 0.0001)
  log.beta ~ dnorm(0, 0.0001)
  alpha <- -exp(log.alpha)
  beta <- -exp(log.beta)
}