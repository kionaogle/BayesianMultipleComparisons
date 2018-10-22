#SAM model for baseline water potential
model{
  for(i in 1:N){
    # normal likelihood
    wp[i] ~ dnorm(mu[i], tau)
    # Log pointwise density, monitor for WAIC pD:
    lpd[i]<-logdensity.norm(wp[i], mu[i], tau)
    
    #linear model
    mu[i] <- b[1,shrubID[i]] + 
      b[2,shrubID[i]]*(maxDant[i]-meanmaxD) + 
      b[3,shrubID[i]]*(swc30ant[i]-meanswc30) + 
      b[4,shrubID[i]]*(swc60ant[i]-meanswc60) + 
      b[5,shrubID[i]]*(maxDant[i]-meanmaxD)*(swc30ant[i]-meanswc30) + 
      b[6,shrubID[i]]*(maxDant[i]-meanmaxD)*(swc60ant[i]-meanswc60) + 
      b[7,shrubID[i]]*(swc30ant[i]-meanswc30)*(swc60ant[i]-meanswc60)
    
    #calculating antecedent variables
    maxDant[i] <- sum(maxDTemp[i,])
    swc30ant[i] <- sum(swc30Temp[i,])
    swc60ant[i] <- sum(swc60Temp[i,])
    for(t in 1:7){#for each time step 
      maxDTemp[i,t] <- maxD[i,t]*wA[t]
      swc30Temp[i,t] <- swc30[i,t]*wB[t]
      swc60Temp[i,t] <- swc60[i,t]*wC[t]
    }
  }
  
  #sum of the deltas for each covariate
  sumA <- sum(deltaA[])
  sumB <- sum(deltaB[])
  sumC <- sum(deltaC[])
  
  #Employing "delta trick" to give vector of weights dirichlet priors
  for(t in 1:7){#for each time step 
    wA[t]<-deltaA[t]/sumA
    deltaA[t]~dgamma(alphaA[t],1)
    wB[t]<-deltaB[t]/sumB
    deltaB[t]~dgamma(alphaB[t],1)
    wC[t]<-deltaC[t]/sumC
    deltaC[t]~dgamma(alphaC[t],1)
  }
  
  #hierachical normal priors for predawn parameters; species level parameters
  for(k in 1:7){#for each parameter
    for(j in 1:8){#for each shrub
      b[k,j] ~ dnorm(mu.b[k], tau.b[k])
    }
  }
  
  #root node priors for mu.b and tau.b; pop'n level parameters
  for(k in 1:7){#for each parameter 
    #diffuse normal priors for mu.b's
    mu.b[k] ~ dnorm(0, 0.00001)
    #sig.b's parameterized with a folded t distribution, 2 degrees of freedom
    tau.eps.b[k] ~  dt(0,Bb[k],2)
    sig.b[k] <- abs(tau.eps.b[k])
    tau.b[k] <- pow(sig.b[k], -2)
    #params for folded t; set as Ab data
    Bb[k] <- 1/(Ab[k]*Ab[k])
  }
  
  #diffuse gamma priors for precision
  tau ~ dgamma(0.01, 0.01)
  sig <- 1/sqrt(tau)
  
  # Relevant to multiple comparisons: compute pairwise differences in
  # all regression parameters among treatment groups:
  # Control shrubs are: 3,4,5,7
  # Treatment shrubs are: 1,2,6,8
  for(p in 1:7){
    b.control[p] <- (sum(b[p,3:5])+sum(b[p,7]))/4
    b.watered[p] <- (sum(b[p,1:2])+sum(b[p,6])+sum(b[p,8]))/4
    b.diff[p] <- b.watered[p] - b.control[p]
    p.b[p] <- step(b.diff[p])
  }  
}

