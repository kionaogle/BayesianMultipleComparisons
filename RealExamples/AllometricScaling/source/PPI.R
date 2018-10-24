################################################################################
# Function to calculate PPI from HB, NH, and CP coda objects using pWAIC2
# Originally calculations written January 20, 2016 by Michael Fell, Jessica Guo and Drew Peltier.  
# FUnction written September 27, 2018 by Jessica Guo and Mike Fell
# Email: jsg88@nau.edu
#
################################################################################
PPI <- function(codahb, codanh, codacp, name="lpd", pwAIC2=F){
  #codahb is coda object of hierarchical Bayesian model version
  #codanh is coda object of non hierarchical Bayesian model version
  #codacp is coda object of complete pooling Bayesian model version
  #name is the name of monitored log predictive density, default = "lpd"
  #pwAIC2 is whether to return all 3 model versions' pWAIC2, default = F

  if(!(length(codahb)==length(codanh) & length(codahb)==length(codacp))){
    stop("length (number of chains) of coda objects must be equal")
  }

  n.chains <- length(codahb)
  
  #index, a function  that takes the effective number of parameters 
  #(no matter how you calculate it), and converts to an index between 0 and 1
  #0 is no pooling, 1 is complete pooling
  index <- function(PD, PD.nh, PD.cp){
    PDind <- 1-(PD-PD.cp)/(PD.nh-PD.cp)
    PDind2 <- (PD.nh-PD)/(PD.nh-PD.cp)
    PDind
  }

  #
  ind.hb<-grep(pattern = paste0("^", name, "\\["), x=colnames(codahb[[1]]), perl=TRUE)
  ind.nh<-grep(pattern = paste0("^", name, "\\["), x=colnames(codanh[[1]]), perl=TRUE)
  ind.cp<-grep(pattern = paste0("^", name, "\\["), x=colnames(codacp[[1]]), perl=TRUE)
  
 
#HB
#pWAIC2 = variance of the individual terms in the log predictive density summed over n data points
#paste together the chains
  HB<-numeric(0)
  for(i in 1:n.chains){
    HB<-rbind(HB, codahb[[i]][,ind.hb])
  }
#take the variance of each lpd column, sum for pWAIC2
  pwaicHB<-sum(apply(HB, 2, FUN=var))


#NH
  #pWAIC2 = variance of the individual terms in the log predictive density summed over n data points
  #paste together the chains
  NH<-numeric(0)
  for(i in 1:n.chains){
    NH<-rbind(NH, codanh[[i]][,ind.nh])
  }
  #take the variance of each lpd column, sum for pWAIC2
  pwaicNH<-sum(apply(NH, 2, FUN=var))
  
#CP
  #pWAIC2 = variance of the individual terms in the log predictive density summed over n data points
  #paste together the chains
  CP<-numeric(0)
  for(i in 1:n.chains){
    CP<-rbind(CP, codacp[[i]][,ind.cp])
  }
  #take the variance of each lpd column, sum for pWAIC2
  pwaicCP<-sum(apply(CP, 2, FUN=var))

#use the index function to calculate PPI from the pWAIC2 of each model version
  if(pwAIC2==F){print(index(pwaicHB, pwaicNH, pwaicCP))} else
  {print(c(pwaicHB, pwaicNH, pwaicCP, index(pwaicHB, pwaicNH, pwaicCP)))}
}
