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
  
  # provide values for alpha and beta, tauD, Omega
}