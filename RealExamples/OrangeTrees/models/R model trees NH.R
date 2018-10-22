# Non-hierarhical model for orange-tree growth model.

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
  # Independent, non-hierachical, relatively non-informative priors
  # for tree-level growth parameters:
  for(p in 1: Nparms){
    for(t in 1:Ntree){
      # Log scale
      log.b[p,t] ~ dnorm(0, 0.0001)
      # Regular scale
      b[p,t] <- exp(log.b[p,t])
    }
  }
  # Conjugate gamma prior for tau (precision), calculate sig (standard deviation)
  tau ~ dgamma(0.01, 0.01)
  sig <- pow(tau, -0.5)
}