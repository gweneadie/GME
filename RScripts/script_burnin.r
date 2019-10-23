## SETUP ##
source("script_set-up.r")

# provide the date that the adaptive chain was created
adaptiveDate = "2019-04-29"

# number of iterations for burn-in
iterations = 5e2

# function for proposal distribution for r, vlos, etc.
mypropDF.RPMVlos = function(sig,...){  mvrnorm(Sigma = sig, mu=c(0,0,0,0)) }

# proposal distribution for model parameters
mypropDF <- function(covmat, ...){ mvrnorm(mu = c(0,0,0,0), Sigma = covmat) }


# load(file = paste("../../PaperIII-chains/covmatrices-", filename, sep=""))
burnin = readRDS(file = paste("../Results/adaptive-burnin-", filename, "_", adaptiveDate, sep=""))

nuiscovmat = burnin[[1]]
propDFcov = burnin[[2]]
rm(burnin)


###################################################################
library(snow)
# -- Set up a cluster on local machine --
# number of local processors
nProc <- 3

# Set up a local cluster with nProc processors
cl <- makeSOCKcluster( rep.int( "localhost", nProc ) )

# Set a random number sequence on the cluster
clusterSetupRNGstream( cl, seed = c(111,12,313) )

# load the libraries, constants, and functions needed to run GME in the nodes of the cluster (see ?list.files to understand the pattern argument)
invisible( clusterEvalQ(cl, library(MASS, coda) ) )

# get a vector of the files that need to be sourced
# filestosource = list.files(pattern = "^[cf]")

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


# Make the markov chains!
sim.list <- parLapply(cl=cl, N=iterations, dat=mydata, x=init.list, fun=HierGalacticMCMC, n.pars=npars,
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

# sim.list is a list with the chains AND the acceptance rates.
# put the chains and acceptance rates into their own objects.
parallel.chains <- lapply( X=sim.list , FUN=function(x) as.mcmc(x$chain) )
acceptance.rates <- lapply(X=sim.list, FUN=function(x) x$acceptance.rate)

# print the acceptance rates to check them
print(acceptance.rates)

# Calculate the Gelman-Rubin statistic, set minimum Rhat acceptable, plot
# require min 1.1 for model pars' point and conf. int.
# require min 1.1 for nuis pars' point
# require min 1.2 for nuis pars' conf. int
Rhats <- gelman.diag(parallel.chains)[[1]]
min.Rhat = 1.1
min.nuisRhat = 1.2
good.parsRhats <- Rhats[ 1:npars, ] < min.Rhat
good.nuisRhats <- cbind((Rhats[ (npars+1):(npars + 4*nrow(mydata)), 1] < min.Rhat), (Rhats[ (npars+1):(npars + 4*nrow(mydata)), 2] < min.nuisRhat))
good.Rhats = rbind(good.parsRhats, good.nuisRhats)
# check for convergence
converged.tf <- all( colSums(good.Rhats) == (npars + 4*nrow(mydata)) )

# if chains have already converged, then run for 5000 steps with 30 thinning (150 000 steps)
if( converged.tf ){
  converged.tf=FALSE
  stopcount = 2
}else{
  stopcount=10
}

count=1
iterations = 5000

# if the chains do not have good Rhat values then set a counter...
count=1

# ... and start running a while loop until things have converged, or until a maximum of 10 loops
while( !converged.tf & count!=10 ){
  
  # count how many iterations it takes to converge
  count = count + 1
  
  # get the last parameter values for each chain
  # initial parameters needs to be a list
  # change the initial values, using the last position of the previous chain
  newstart.list = init.list
  
  n = nrow(parallel.chains[[1]])
  
  for( i in 1:nProc ){
    newstart.list[[i]][[1]] <- as.numeric( parallel.chains[[i]][n, 1:npars] )
    newstart.list[[i]][[2]] <- sim.list[[i]]$last.Data
  }
  
  rm(parallel.chains)
  rm(sim.list)
  
  # run chains again
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
  
  # put the chains and acceptance rates into their own objects.
  parallel.chains = lapply( X=sim.list , FUN=function(x) as.mcmc(x$chain) )
  acceptance.rates = lapply(X=sim.list, FUN=function(x) x$acceptance.rate)
  acceptnuis.rates = lapply(X=sim.list, FUN=function(x) x$acceptRPMVlos.rate)
  
  # calculate Rhat for the chains, at different intervals
  Rhats = gelman.diag(parallel.chains)[[1]]
  
  # find out if all chains have Rhats< min Rhat for all intervals
  good.parsRhats <- Rhats[ 1:npars, ] < min.Rhat
  good.nuisRhats <- cbind((Rhats[ (npars+1):(npars + 4*nrow(mydata)), 1] < min.Rhat), (Rhats[ (npars+1):(npars + 4*nrow(mydata)), 2] < min.nuisRhat))
  good.Rhats = rbind(good.parsRhats, good.nuisRhats)
  
  converged.tf <- all( colSums(good.Rhats) == (npars + 4*nrow(mydata)) )
  
} # end while loop

stopCluster(cl)

saveRDS(sim.list, file = paste("../Results/converged", filename, adaptiveDate, iterations, "iterations", nProc, "chains", Sys.Date(), sep="-") )



