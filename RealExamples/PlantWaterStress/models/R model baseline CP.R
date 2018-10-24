# Bayesian model that fits a SAM model to baseline water potential data:
# Complete pooling version

model{
  for(i in 1:N){
    # normal likelihood for water potential data
    wp[i] ~ dnorm(mu[i], tau)
    # Log pointwise density, monitor for WAIC pD:
    lpd[i] <- logdensity.norm(wp[i], mu[i], tau)
    
    # linear model; the shrub-specific coefficients will be set equal to
    # global-level parameters to achieve complete pooling (see below)
    mu[i] <- b[1,shrubID[i]] + 
      b[2,shrubID[i]]*(maxDant[i]-meanmaxD) + 
      b[3,shrubID[i]]*(swc30ant[i]-meanswc30) + 
      b[4,shrubID[i]]*(swc60ant[i]-meanswc60) + 
      b[5,shrubID[i]]*(maxDant[i]-meanmaxD)*(swc30ant[i]-meanswc30) + 
      b[6,shrubID[i]]*(maxDant[i]-meanmaxD)*(swc60ant[i]-meanswc60) + 
      b[7,shrubID[i]]*(swc30ant[i]-meanswc30)*(swc60ant[i]-meanswc60)
    
    # Calculate antecedent variables
    maxDant[i] <- sum(maxDTemp[i,])
    swc30ant[i] <- sum(swc30Temp[i,])
    swc60ant[i] <- sum(swc60Temp[i,])
    for(t in 1:7){#for each past time step 
      maxDTemp[i,t] <- maxD[i,t]*wA[t]
      swc30Temp[i,t] <- swc30[i,t]*wB[t]
      swc60Temp[i,t] <- swc60[i,t]*wC[t]
    }
  }
  
  # Sum of the deltas for each covariate
  sumA <- sum(deltaA[])
  sumB <- sum(deltaB[])
  sumC <- sum(deltaC[])
  
  # Employ "delta trick" to give vector of weights dirichlet priors
  for(t in 1:7){#for each time step 
    wA[t]<-deltaA[t]/sumA
    deltaA[t]~dgamma(alphaA[t],1)
    wB[t]<-deltaB[t]/sumB
    deltaB[t]~dgamma(alphaB[t],1)
    wC[t]<-deltaC[t]/sumC
    deltaC[t]~dgamma(alphaC[t],1)
  }
  
  # Complete pooling of shrub level parameters; i.e., set shrub-level
  # parameters equal to global parameters.
  for(k in 1:7){#for each parameter
    for(j in 1:8){#for each shrub
      b[k,j] <- mu.b[k]
    }
  }
  
  #root node priors for mu.b; pop'n level parameters
  for(k in 1:7){#for each parameter 
    #diffuse normal priors for mu.b's
    mu.b[k] ~ dnorm(0, 0.00001)
  }
  
  # Vague gamma prior for precision
  tau ~ dgamma(0.01, 0.01)
  sig <- 1/sqrt(tau)
  
}

