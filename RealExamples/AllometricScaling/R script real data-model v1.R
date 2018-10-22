# setwd("C:\\Users\\mkf58\\Google Drive\\Oglelab_NAU_JOB\\NAU_JOB\\MultiComp\\Review_Example\\MulticompExample\\Price_scaling")

# A horrible recursive function to find groups
groupcheck <- function(datain, testvars){
    dataout <- integer(0)
    
    test <- t(apply(X=datain[,4:5], MARGIN=1, FUN=match, table=testvars))
    test <- apply(X=!is.na(test), MARGIN=1, FUN=sum)
    test <- datain[which(test == 1), 4:5]

    if(!is.null(dim(test)) && nrow(test) != 0){
        for(i in 1:1){
            newvar <- test[i,which(test[i,] %in% testvars == FALSE)]
            test2 <- test
            
            test2[,1] <- test2[,1] %in% newvar  
            test2[,2] <- test2[,2] %in% newvar  
            test2 <- test[which(rowSums(test2)==1),][!test2[which(rowSums(test2)==1),]]
            # make sure there are no groupings of 4.
            if(sum(test2 %in% testvars) == length(testvars)){
                testvars <- c(testvars, newvar)
                dataout <- unlist(groupcheck(datain=datain, testvars=testvars))
            }else{
                dataout <- unlist(testvars)
            }
        }
    }else{
        dataout <- unlist(testvars)
    }
    
    return(dataout)
}

if("rjags" %in% installed.packages()[,1] == FALSE){
	install.packages("rjags")
	library(rjags)
} else{
	library(rjags)
}
if("mcmcplots" %in% installed.packages()[,1] == FALSE){
	install.packages("mcmcplots")
	library(mcmcplots)
} else{
	library(mcmcplots)
}

install.packages("devtools")
library(devtools)
devtools::install_github("fellmk/PostJAGS/postjags")
library(postjags)

chuck = read.csv("ChuckDataMatrix.csv")

Nsp = 49
Ncomp = Nsp*(Nsp-1)/2

# Indices for species pairwise comparisons; number of comparisons = (Nsp)*(Nsp-1)/2
pairID = matrix(data =NA, nrow = Ncomp, ncol=2)
ii = 0
for(s1 in 1:(Nsp-1)){
  for(s2 in (s1+1):Nsp){
    ii = ii +1
    pairID[ii,1] = s1
    pairID[ii,2] = s2
  }
}

data = list(N=1162, Nsp=Nsp, R=structure(.Data=c(1,0,0,1),.Dim=c(2,2)), pairID = pairID, Ncomp = Ncomp,
            SP = chuck$SP, LogMass = chuck$LogMass, LogLength = chuck$LogLength, LogDiameter = chuck$LogDiameter,
            nVars=2, pi=pi)

source("R inits chuck v1.R")


n.chains = 3

# Hierarchical model:
jm.hb <- jags.model(file = "models/R model HB_lpd.R", data = data, inits = chuck_inits, n.chains = n.chains, n.adapt = 1000)

# Complete pooling model:
jm.cp <- jags.model(file = "models/R model CP_lpd.R", data = data, inits = chuck_inits, n.chains = n.chains, n.adapt = 1000)

# Non-hierarchical model:
jm.nh <- jags.model(file = "models/R model NH_lpd.R", data = data, inits = chuck_inits, n.chains = n.chains, n.adapt = 1000)

################################## Added by MKF ###############################
# 8/31/2018
# Burn In
jm.hb.names <- c(variable.names(jm.hb), "deviance")
jm.hb.names <- c(jm.hb.names, "alpha.diff", "beta.diff", "p.alpha", "p.beta", "mu.beta", "tau.beta",
                 "mu.alpha", "tau.alpha", "sig.alpha", "sig.beta")
jm.cp.names <- c(variable.names(jm.cp), "deviance")
jm.nh.names <- c(variable.names(jm.nh), "deviance")

load.module("dic")

##### Run Hierarchical model #####
runtime.hb <- system.time(
{
update(jm.hb, n.iter=10000)
codahb =coda.samples(jm.hb,variable.names=jm.hb.names,
                  n.iter=2000000, thin=1000)
}
)
print(runtime.hb)
#save(codahb, file = "Output/jm_hb4.Rdata")
				  
##### Run Complete pooling model #####
runtime.cp <- system.time(
{
update(jm.cp, n.iter=10000)
codacp =coda.samples(jm.cp,variable.names=jm.cp.names,
                  n.iter=2000000, thin=1000)
}
)		
#save(codacp, file = "Output/jm_cp3.Rdata")	
print(runtime.cp)

##### Run Non-hierarchical model #####
runtime.nh <- system.time(
{
update(jm.nh, n.iter=10000)
codanh =coda.samples(jm.nh,variable.names=jm.nh.names,
                  n.iter=2000000, thin=1000)
}
)	
#save(codanh, file = "Output/jm_nh3.Rdata")
print(runtime.nh)			  

##### Save model objects to output folder #####				  
# Moved so files are saved after the model runs
#save(codahb, file = "Output/jm_hb4.Rdata")
#save(codacp, file = "Output/jm_cp3.Rdata")
#save(codanh, file = "Output/jm_nh3.Rdata")

##### Save model objects to output folder #####				  
load(file = "Output/jm_hb4.Rdata")
load(file = "Output/jm_cp3.Rdata")
load(file = "Output/jm_nh3.Rdata")

codahb[[1]][1:5,1:5]

##### Jess Code #####
sumhb<-coda.fast(coda=codahb, chains=3)
sumcp<-coda.fast(coda=codacp, chains=3)
sumnh<-coda.fast(coda=codanh, chains=3)

sumhb[which(rownames(sumhb)=="deviance"),]
-2*sum(sumhb[which(rownames(sumhb)=="log_lik[1]"):which(rownames(sumhb)=="log_lik[1162]"),1])
-2*sum(sumhb[which(rownames(sumhb)=="lpd[1]"):which(rownames(sumhb)=="lpd[1162]"),1])

sumcp[which(rownames(sumcp)=="deviance"),]
-2*sum(sumcp[which(rownames(sumcp)=="log_lik[1]"):which(rownames(sumcp)=="log_lik[1162]"),1])
-2*sum(sumcp[which(rownames(sumcp)=="lpd[1]"):which(rownames(sumcp)=="lpd[1162]"),1])

sumnh[which(rownames(sumnh)=="deviance"),]
-2*sum(sumnh[which(rownames(sumnh)=="log_lik[1]"):which(rownames(sumnh)=="log_lik[1162]"),1])
-2*sum(sumnh[which(rownames(sumnh)=="lpd[1]"):which(rownames(sumnh)=="lpd[1162]"),1])


############################# Get parameters###################################
#
# This section of code can likely be generalized to a function for creating
# inputs for running a new model. 9/4/2018
#
###############################################################################

#sumtime <- system.time(             
#    sumhbold <- summary(codahb)
#)

fasttime <- system.time(
    sumhb <- coda.fast(chains=3, thin=1, coda=codahb)
)
sumhbtemp <- sumhb

sumhb <- as.matrix(sumhb)

varlist <- c("N", "Omega", "alpha", "beta", "SP", "LogDiameter", 
             "tauD")# A list of variables needed

##parms <- initfind(codahb, iteration\=0, OpenBUGS=FALSE)

#summary.names <- row.names(sumhb$statistics)
summary.names <- row.names(sumhb)
out.names <- character(0)
for(snames in summary.names){
	snames <- unlist(strsplit(x=snames, split="\\[", perl=TRUE))[1]
	out.names <- c(out.names, snames)
}

# This will provide a warning if the output is erroneous.
if(length(out.names) != length(summary.names)){
	stop("The length of the names vectors for the posterior and the selection 
	list are not equal.")
}

#parms <- sumhb$statistics[which(out.names %in% varlist),]
parms <- sumhb[which(out.names %in% varlist),]
#parms <- parms[,which(colnames(parms) == "Mean")]
parms <- parms[,which(colnames(parms) == "mean")]

# Trick initfind into accepting summary values as inputs 
parms <- as.matrix(parms)
parms <- t(parms)
parms <- mcmc(data=list(parms), thin=1)
parms <- initfind(mcmcin=parms)

# If parameter list is already restricted then use parms as is without 
# running removevars.
parms <- parms[[2]]
parms <- parms[[1]]

################################ Group Code ###################################
#
# This section of code is for creating groups that are equal. 9/5/2018
#
###############################################################################
sumhb <- sumhbtemp
#source("src/MCMCrestart_v2_3.R")
#source("src/coda_functions/coda_fast.R")
varnames <- c("^p.alpha\\[", "^p.beta\\[")
parmnames <- c("^alpha\\[", "^beta\\[")
probcutoff <- 0.8

for(j in 1:length(varnames)){
    varname <- varnames[j]
    parmname <- parmnames[j]

    # Get the data for the parameter
    parmvals <- rownames(sumhb)
    parmvals <- grep(parmname, parmvals, perl=TRUE)
    parmvals <- sumhb[parmvals,]

    # Get the p-values for the parameter
    pvalue <- rownames(sumhb)
    pvalue <- grep(varname, pvalue, perl=TRUE)
    pvalue <- sumhb[pvalue,]

    # apply pvalue correction for less than 0.5
    pvalue[pvalue$mean > 0.5,]$mean <- 1 - pvalue[pvalue$mean > 0.5,]$mean

    # Get index values from names
    indexvals <- integer(0)
    for(name in rownames(pvalue)){
        name <- unlist(strsplit(x=name, split="\\["))[2]
        name <- unlist(strsplit(x=name, split="\\]"))[1]
        name <- unlist(strsplit(x=name, split=","))
        name <- as.integer(name)
        indexvals <- rbind(indexvals, name)
    }
    # Remove rownames and add species ids to columns 3 and 4
    indexvals <- cbind(indexvals, pvalue$mean, rbind(pairID, pairID))
    rownames(indexvals) <- NULL
    colnames(indexvals) <- c("comp", "k", "mean.p", "sp1", "sp2")
    indexvals <- data.frame(indexvals)

    for(k in 1:2){ #k indexes alpha and beta
        # indexvals[indexvals$k==k,]
        cutoff <- quantile(x=indexvals[indexvals$k==k,]$mean.p, probs=c(probcutoff))
        cutoff <- 0.25
        datafind <- which(indexvals[indexvals$k==k,]$mean.p>cutoff)
        datafind <- indexvals[indexvals$k==k,][datafind,]
        datafind <- datafind[order(datafind$mean, decreasing=TRUE),]


        # Test group finding function
        #head(datafind)
        #testoutput <- groupcheck(datafind, testvars = datafind[1, 4:5])

        infunc <- function(vals, compvals){
            return(vals %in% compvals)
        }

        spfind <- integer(0)
        spmax <- 49
        splist <- 1:data$Nsp
        spout <- list()
        i <- 1
        # Find the groupings 
        while(length(splist) != 0 & i < spmax & nrow(datafind) > 0){
            sptest <- datafind[1, 4:5]
            sptest <- groupcheck(datafind, testvars = datafind[1, 4:5])
            splist <- splist[-which(splist %in% sptest)]
            spfind <- c(spfind, sptest)
            # remove comparisons from datafind for species selected above
            if(length(sptest) > 0 & sum(is.na(sptest)==0)){
                datafind <- datafind[-which(datafind$sp1 %in% spfind),]
                # if the which results in 0 it is bad.
                if(length(which(datafind$sp2 %in% spfind)) > 0){
                    datafind <- datafind[-which(datafind$sp2 %in% spfind),]
                }
            }
            
            if(sum(is.na(sptest))>0){
                break
            }
            
            if(length(sptest) > 0){
                spout[[i]] <- sptest
            }
            i <- i + 1
        }

        # Create new parametervalues
        
        parmname <- unlist(strsplit(x=parmname, split="\\\\", perl=TRUE))[1]
        parmindex <- grep(pattern=parmname, x=names(parms))
        parms[[parmindex]]
        for(group in spout){
            print(group)
            parms[[parmindex]][group,k] <- mean(parms[[parmindex]][group,k])
        }

    } # End of k loop
} # End of j loop

parms.out <- parms
save(parms.out, file="price_scaling_simulation_parameters.Rdata")

############################# Simulate Data ###################################
#
# This section of code can likely be generalized to a function for creating
# inputs for running a new model. 9/4/2018
#
##### Try model run! #####
# Hierarchical model:
###############################################################################
n.chains <- 1
n.iters <- 100
base.file.name <- "Output/datasets/Example1"
parmsave <- c("Y")
sim.model <- "models/R_model_HB_simulate_data_mf.R"

jm.hb.gen <- jags.model(file = sim.model, data = parms, n.chains = n.chains, n.adapt = 1000)

jm.hb.names <- c(variable.names(jm.hb.gen))
# Run the model to get new data!
runtime.hb.gen <- system.time(
{
update(jm.hb.gen, n.iter=1000)
codahbnew =coda.samples(jm.hb.gen,variable.names=jm.hb.names,
                  n.iter=n.iters, thin=1)
}
)

#save(codahbnew, file="newsamples.Rdata")
load("newsamples.Rdata")

for(iter in 1:n.iters){
    if(length(codahbnew) > 1){
        stop("CODA for simulated data should have one chain. If you chose to use 
multiple chains please break up the coda object into three objects and process 
them serially.")
    }
    
    # Extract data for each iteration and save it out.
    parm.sets <- initfind(mcmcin=codahbnew, iteration=iter)
    parm.sets <- removevars(initsin=parm.sets, variables=which(!(parm.sets[[1]] %in% parmsave)))
    parm.sets <- parm.sets[[2]][[1]]
    # Save the data
    save(parm.sets, file=paste0(base.file.name, "_", iter, ".Rdata"))
}









