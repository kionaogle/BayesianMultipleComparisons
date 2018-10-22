# Complete pooing model for orange-tree growth model.

model{
  for(i in 1:N){
    # Normal likelihood for circumference data (Y)
    Y[i] ~ dnorm(mu[i], tau)
    # Pointwise predictive density
    lpd[i] <- logdensity.norm(Y[i], mu[i], tau)
    # let b[1,i] = Ymax[i] (maximum circumference)
    # let b[2,i] = Ymin[i] (circumference at time of planting)
    # let b[3,i] = Ghalf[i] (growth rate at half max circumference)
    # Non-linear, logistic growth curve for predicted (mean) circumference:
    mu[i] <- b[1,Tree[i]]/(1+exp(-4*b[3,Tree[i]]/b[1,Tree[i]]*Time[i])*
                             (b[1,Tree[i]]/b[2,Tree[i]]-1))
  }
  # Set tree-level parameters equal to population-level parameters:
  for(p in 1: Nparms){
    for(t in 1:Ntree){
      # Log scale
      log.b[p,t] <- log.mu.b[p]
      # Regular scale
      b[p,t] <- mu.b[p]
    }
  }
  # Relatively non-informative priors for overall, pop'n level parameters:
  for(i in 1:Nparms){
    log.mu.b[i] ~ dnorm(0, 0.0001)
    mu.b[i] <- exp(log.mu.b[i])
  }
  
  # Conjugate gamma prior for tau (precision), calculate sig (standard deviation)
  tau ~ dgamma(0.01, 0.01)
  sig <- pow(tau, -0.5)
}