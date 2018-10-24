# The Non-hierarchical (NH) version of the "allometric scaling" model.

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

# Non-hierarhical model for species-level parameters:

model
{
  # Loop through each observation i in the dataset:
  for(i in 1:N){
    # The likelihood (sampling distribution for loglength and logmass) is a
    # multivariate normal distribution with mean mu and precision matrix Omega.
    # The original paper ran several different models, with different assumptions about
    # the numerical value of the scaling coefficiencts. Here, we implement the SPAM model,
    # or the "specialized" model, which assumes species-specific parameters, but these
    # are not modeled hierarchically.
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
    #nVars=2, set as data
    lpd[i]<-0.5*logdet(Omega[1:2,1:2])-nVars/2*log(2*pi)-t(Ydiff[i,1:2])%*%Omega[1:2,1:2]%*%(Ydiff[i,1:2])/2 
    # lpd could also be computed as (log_lik):
    #first using the logdensity.mnorm function in JAGS 4.0 and higher
    #log_lik[i]<-logdensity.mnorm(Y[i,], mu[i,], Omega)
  }
  
  # NON-HIERARCHICAL priors for the species-specific parameters.
  for(j in 1:Nsp){ # species loop
    for(k in 1:2){ # dependent variable loop
      # Scaling exponent:
      beta[j,k] ~ dnorm(0,0.00001)
      # Normalizing constant (intercept):
      alpha[j,k] ~ dnorm(0,0.00001)
    }
  }
  
  # Hyperpriors for global parameters and priors for variance/precision terms:
  Omega[1:2,1:2] ~ dwish(R[1:2,1:2], 2)
  Sigma[1:2,1:2] <- inverse(Omega[1:2,1:2])
  
  # Semi-infomative prior for precision of logDiameter measurement error.
  sigD ~ dlnorm(-4.135, 2)
  tauD <- pow(sigD, -2)		
}