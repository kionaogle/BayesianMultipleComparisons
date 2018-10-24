# The forcal, Hierarchical Bayesian (HB) version of the "allometric scaling" model.

data{
  # Create data matrix for multivariate likelihood / data model. The data include the variables 
  # LogLength = log(l), LogMass = log(M). 
  # Put all data in the data vector Y, which is equivalent to the observation
  # vector for [log(l) log(M)] in eqn (2) in the main text of Price et al. (2009).
  for (i in 1:N) {
    Y[i, 1] <- LogLength[i]
    Y[i, 2] <- LogMass[i]
  }
}

model
{
  # Loop through each observation i in the dataset:
  for(i in 1:N){
    # The likelihood (sampling distribution for loglength and logmass) is a
    # multivariate normal distribution with mean mu and precision matrix Omega.
    # The original paper ran several different models, with different assumptions about
    # the numerical value of the scaling coefficiencts. Here, we implement the SPAM model,
    # or the "specialized" model, which reflects a "standard" hierarchical model
    # for species-specific parameters.
    Y[i,1:2] ~ dmnorm(mu[i,1:2], Omega[1:2,1:2])
    
    # Define the mean vector (i.e., scaling model that relates the true or latent variables).
    # alpha is the species-specific normalizing constant, and beta is the species-specific
    # scaling exponent.
    for(k in 1:2){
      mu[i,k] <- alpha[SP[i],k] + beta[SP[i],k]*Lrho[i]
    }
    
    # Berkson model for "true" or latent log diameter (Lrho). 
    # Lrho varies about measured LogDiameter, and tauD is the precision (1/variance) that 
    # describes measurement errror. The Berkson model is for "LatentRho".
    LatentRho[i] ~ dnorm(LogDiameter[i], tauD)
    Lrho[i] <- LatentRho[i]
    
    for(k in 1:2){
      # difference between Y and mu
      Ydiff[i,k]<-Y[i,k]-mu[i,k]
    }
    # log pointwise density
    lpd[i]<-0.5*logdet(Omega[1:2,1:2])-nVars/2*log(2*pi)-t(Ydiff[i,1:2])%*%Omega[1:2,1:2]%*%(Ydiff[i,1:2])/2 
    #nVars=2, set as data
    
    # lpd could also be computed as (log_lik):
    #first using the logdensity.mnorm function in JAGS 4.0 and higher
    #log_lik[i]<-logdensity.mnorm(Y[i,], mu[i,], Omega)
    #nVars=2, set as data
  }
  
  # Hierarchical priors for the species-specific parameters.
  for(j in 1:Nsp){ # species loop
    for(k in 1:2){ # dependent variable loop
      # Scaling exponent:
      beta[j,k] ~ dnorm(mu.beta[k],tau.beta[k])
      # Normalizing constant (intercept):
      alpha[j,k] ~ dnorm(mu.alpha[k], tau.alpha[k])
    }
  }
  
  # Hyperpriors for global parameters and priors for variance/precision terms:
  Omega[1:2,1:2] ~ dwish(R[1:2,1:2], 2)
  Sigma[1:2,1:2] <- inverse(Omega[1:2,1:2])
  for(k in 1:2){
    mu.alpha[k] ~ dnorm(0,0.00001)
    tau.alpha[k] ~ dgamma(0.01,0.001)
    sig.alpha[k] <- sqrt(1/tau.alpha[k])
    mu.beta[k] ~ dnorm(0,0.00001)
    tau.beta[k]  ~ dgamma(0.01,0.001) 
    sig.beta[k] <- sqrt(1/tau.beta[k])
  }
  
  # Semi-infomative prior for precision of logDiameter measurement error.
  sigD ~ dlnorm(-4.135, 2)
  tauD <- pow(sigD, -2)		
  
  # Pairwise comparisons of species:
  # Compute for each species pair (n = Ncomp), for each of 2 (k loop) dependent variables:
  for(j in 1:Ncomp){
    for(k in 1:2){
      alpha.diff[j,k] <- alpha[pairID[j,1],k] - alpha[pairID[j,2],k]
      beta.diff[j,k] <- beta[pairID[j,1],k] - beta[pairID[j,2],k]
      # Posterior means of p.alpha and p.beta give Bayesian p-values for comparisons
      p.alpha[j,k] <- step(alpha.diff[j,k])
      p.beta[j,k] <- step(beta.diff[j,k])
    }
  }
}