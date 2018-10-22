data{
  # Create data matrix for multivariate likelihood / data model.The cata include the variables 
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
    # The likelihood (sampling distribution for loglength, logmass, logarea) is a
    # multivariate normal distribution with mean mu and precision matrix Omega.
    # The original paper ran several different models, with different assumptions about
    # the numerical value of the scaling coefficiencts. Here, we implement the SPAM model,
    # or the "specialized" model which reflects a "standard" hierarhical model
    # for species-specific parameters.
    Y[i,1:2] ~ dmnorm(mu[i,1:2], Omega[1:2,1:2])
    
    # Define the mean vector (i.e., scaling model that relates the true or latent variables).
    # alpha is the species-specific normalizing constant, and beta is the species-specific
    # scaling exponent.
    for(k in 1:2){
      mu[i,k] <- alpha[SP[i],k] + beta[SP[i],k]*Lrho[i]
    }

    # Berkson model for "true" or latent log diameter (Lrho). Note, the data are in terms of 
    # diameter, and the use of diameter vs radius will not affect the scaling exponents.
    # Lrho varies about measured LogDiameter, and tauD is the precision (1/variance) that 
    # describes measurement errror. The Berkson model is for "LatentRho".
    LatentRho[i] ~ dnorm(LogDiameter[i], tauD)
    Lrho[i] <- LatentRho[i]
  }
  
  # Hierarchical priors for the species-specific parameters.
  for(j in 1:Nsp){
    for(k in 1:2){
      # Set scaling exponent:
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
  
  # Semi-infomative prior for precision of logRadius measurement error.
  sigD ~ dlnorm(-4.135, 2)
  tauD <- pow(sigD, -2)		
  
  # Pairwise comparisons of species:
  # Compute for each species pair (n = Ncomp), for each of 2 (k loop) variables:
  for(j in 1:Ncomp){
    for(k in 1:2){
      alpha.diff[j,k] <- alpha[pairID[j,1],k] - alpha[pairID[j,2],k]
      beta.diff[j,k] <- beta[pairID[j,1],k] - beta[pairID[j,2],k]
      # P-value code:
      p.alpha[j,k] <- step(alpha.diff[j,k])
      p.beta[j,k] <- step(beta.diff[j,k])
    }
  }
}