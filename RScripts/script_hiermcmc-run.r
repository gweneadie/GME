# once a covariance matrix is found and the chain has converged, use this script to take a lot of samples
# call it from a regular terminal, not R-Studio

# file IDs for the covariance matrices and the burnin
adaptiveDate = "2019-04-29"
burninDate = "2019-04-29"

# load 'snow' package, random number generator called 'rlecuyer'
library( snow )
library( rlecuyer )
library( RColorBrewer )
library( emdbook )
library( moments )
library( pracma )
library( BMS )


### load other code ############################
# source model functions, priors, etc
source("script_set-up.r")

####################################################################

# -- Set up a cluster on local machine --
# number of local processors
nProc <- 3

# Set up a local cluster
cl <- makeSOCKcluster( rep.int( "localhost", nProc ) )

# Set the random number sequence on the cluster
clusterSetupRNGstream( cl, seed = c(11,102,33,410,35,671) )


# load the functions and libraries needed to run mcmcFun in the nodes of the cluster
invisible( clusterEvalQ(cl, library(MASS, coda) ) )

# source files needed in the nodes of the cluster
invisible( clusterEvalQ(cl, source("constants_TmatrixSun.r")) )
invisible( clusterEvalQ(cl, source("function_energy.r")) )
invisible( clusterEvalQ(cl, source("functions_conversions.r")) )
invisible( clusterEvalQ(cl, source("functions_deason.r")) )
invisible( clusterEvalQ(cl, source("functions_observe-priors.r")) )
invisible( clusterEvalQ(cl, source("functions_parameter-transformations.r")) )
invisible( clusterEvalQ(cl, source("functions_priors.r")) )

# source hierarchical Bayesian sampling code
invisible( clusterEvalQ(cl, source("GME-hierarchical.r")) )

# script that loads data, set prior parameter values, etc.

# function for proposal distribution for r, vlos, etc.
mypropDF.RPMVlos = function(sig,...)  mvrnorm(Sigma = sig, mu=c(0,0,0,0))

# proposal distribution for model parameters
mypropDF <- function(covmat, ...) mvrnorm(Sigma=covmat, mu=c(0,0,0,0) )


iterations=5e3

covmatrices = readRDS(file = paste("../Results/adaptive-burnin-", filename, "_", adaptiveDate, sep="") )

nuiscovmat = covmatrices[[1]]
propDFcov = covmatrices[[2]]

# load initial burn-in of all chains
sim.list = readRDS(file = paste("../Results/converged", filename, adaptiveDate, iterations, "iterations", nProc, "chains", burninDate, sep="-") )



# run 10,000 samples 3 times and save as three files, for each rmin
  
# set number of samples wanted
iterations = 1e4
  
for( j in 1:3 ){
    parallel.chains = lapply( X=sim.list , FUN=function(x) as.mcmc(x$chain) )
    
    # get the last parameter values for each chain
    # initial parameters needs to be a list
    # change the initial values, using the last position of the previous chain
    newstart.list = init.list
    n = nrow(parallel.chains[[1]])
    for( i in 1:nProc ){
      newstart.list[[i]][[1]] <- as.numeric( parallel.chains[[i]][n, 1:npars] )
      newstart.list[[i]][[2]] <- sim.list[[i]]$last.Data
    }
    
    # remove old results
    rm(sim.list)
    rm(parallel.chains)

    # run chains
    sim.list <- parLapply(cl=cl, N=iterations, dat=mydata, x=newstart.list, fun=HierGalacticMCMC, n.pars=npars,
                          logPDF=logDF, DF=modelDF, pot=modelpot, propDF=mypropDF,
                          priors=prior.wrapper, transform.pars=transform.func,
                          pm.par=yesnopm,
                          PMra.priors=priorfunc.PMra, PMdec.priors=priorfunc.PMdec,
                          Vlos.priors=priorfunc.Vlos, Rgc.priors=priorfunc.Rgc,
                          PMbounds = myPMbounds, Vlosbounds = myVlosbounds, Rgcbounds = myRgcbounds,
                          propDF.RPMVlos = mypropDF.RPMVlos, propsd.RPMVlos.empir=nuiscovmat,
                          Tmatrix=myTmatrix,
                          priorfuncs = list(singleunif.prior, gaus.prior, alphaprior, singleunif.prior),
                          ppars = list(phiobounds, c(gammamean, gammasd), apriorpars, bbounds),
                          countSun=TRUE, solarmotion=sunmotion, velocityLSR=vLSR, savedatapars = TRUE,
                          covmat=propDFcov)
    

    # save object
        saveRDS(object = sim.list, file = paste("../Results/final", filename, adaptiveDate, burninDate, iterations, "iterations", "thinning20", nProc, "chains", j, Sys.Date(), sep="-"))
    
}

stopCluster(cl)


