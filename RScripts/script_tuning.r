# Adaptive mcmc run for tuning chains

# source file that sets-up parameters like file name, etc.
source("script_set-up.r")

# from here on 'mydata' object used for analysis

# initial proposal distribution for model parameters (zeros in off-diagonals)
mypropDF <- function(x, ...) rnorm( 4, mean=c(0,0,0,0) , sd=c(0.07, 0.05, 0.01, 0.08))
# function for proposal distribution for r, vlos, etc.
mypropDF.RPMVlos = function(sig,...){  mvrnorm(Sigma = sig, mu=c(0,0,0,0)) }

# set-up covariance matrices for jumping distribution
y = 1e-6
mypropsd.RPMVlos = lapply( X=1:nrow(mydata), FUN=function(x) matrix(c(y,0,0,0, 
                                                                      0,y,0,0,
                                                                      0,0,y,0,
                                                                      0,0,0,y), nrow=4) )

test = adjust.propsd(init=init.list[[1]], dat=mydata, DF=modelDF, pot=modelpot, n.pars=npars, 
                     priors=prior.wrapper, propDF=mypropDF, N=1e3, logPDF=logDF, yourpatience=5,
                     transform.pars=transform.func,
                     PMra.priors= priorfunc.PMra, PMdec.priors=priorfunc.PMdec, Vlos.priors=priorfunc.Vlos, Rgc.priors=priorfunc.Rgc,
                     PMbounds = myPMbounds, Vlosbounds = myVlosbounds, Rgcbounds = myRgcbounds,
                     propDF.RPMVlos = mypropDF.RPMVlos, propsd.RPMVlos.empir=mypropsd.RPMVlos,
                     pm.pars=yesnopm,
                     priorfuncs = list(singleunif.prior, gaus.prior, alphaprior, singleunif.prior),
                     ppars = list(phiobounds, c(gammamean, gammasd), apriorpars, bbounds),
                     countSun=TRUE, solarmotion=sunmotion, velocityLSR=vLSR, savedatapars=TRUE
) 


# get covariance matrices for rgc, vlos, and proper motions
newcovmat = test[[1]]

# replace the 0 covariance matrices with much smaller step size
thebadcovmats = unlist(lapply(X = newcovmat, FUN = function(x) all(x==matrix(rep(0,16), nrow=4)) ))

if(any(thebadcovmats)){
  y = 1e-9
  newcovmat[thebadcovmats] = lapply(1:sum(thebadcovmats), 
                                    FUN = function(x) matrix(c(y,0,0,0,
                                                               0,y,0,0,
                                                               0,0,y,0,
                                                               0,0,0,y), nrow=4))
}

# get covariance matrix for model parameters
propDFcov = cov(test[[4]][,1:4])
dimnames(propDFcov)=NULL

init.list[[1]] = test[[2]]

# new proposal distribution FUNCTION that uses newcovmat
mypropDF <- function(covmat, ...){ mvrnorm(mu = c(0,0,0,0), Sigma = covmat) }


# run again but with the adaptiveBurnin function
test2 = adaptiveBurnin(init=init.list[[1]], dat=mydata, DF=modelDF, pot=modelpot, n.pars=npars, 
                       priors=prior.wrapper, propDF=mypropDF, N=1e3, logPDF=logDF, yourpatience=3,
                       transform.pars=transform.func,
                       PMra.priors= priorfunc.PMra, PMdec.priors=priorfunc.PMdec, Vlos.priors=priorfunc.Vlos, Rgc.priors=priorfunc.Rgc,
                       PMbounds = myPMbounds, Vlosbounds = myVlosbounds, Rgcbounds = myRgcbounds,
                       propDF.RPMVlos = mypropDF.RPMVlos, propsd.RPMVlos.empir=newcovmat,
                       pm.par=yesnopm,
                       priorfuncs = list(singleunif.prior, gaus.prior, alphaprior, singleunif.prior),
                       ppars = list(phiobounds, c(gammamean, gammasd), apriorpars, bbounds),
                       countSun=TRUE, solarmotion=sunmotion, velocityLSR=vLSR,
                       savedatapars=TRUE, parcovmat=propDFcov)

# get covariance matrices for rgc, vlos, and proper motions
newcovmat = test2[[1]]

# replace the 0 covariance matrices with much smaller step size
thebadcovmats = unlist(lapply(X = newcovmat, FUN = function(x) all(x==matrix(rep(0,16), nrow=4)) ))

if(any(thebadcovmats)){
  y = 1e-10
  newcovmat[thebadcovmats] = lapply(1:sum(thebadcovmats), FUN = function(x) matrix(c(y,0,0,0,
                                                                                     0,y,0,0,
                                                                                     0,0,y,0,
                                                                                     0,0,0,y), nrow=4))
  
}

# get covariance matrix for model parameters and initial values for next run
propDFcov = test2[[2]]
dimnames(propDFcov)=NULL
init.list[[1]] = test2[[3]]

# plot to check chains
plot(as.mcmc(test2[[5]]))

# run AGAIN if necessary 
test3 = adaptiveBurnin(init=init.list[[1]], dat=mydata, DF=modelDF, pot=modelpot, n.pars=npars, 
                       priors=prior.wrapper, propDF=mypropDF, N=1e3, logPDF=logDF, yourpatience=10,
                       transform.pars=transform.func,
                       PMra.priors= priorfunc.PMra, PMdec.priors=priorfunc.PMdec, Vlos.priors=priorfunc.Vlos, Rgc.priors=priorfunc.Rgc,
                       PMbounds = myPMbounds, Vlosbounds = myVlosbounds, Rgcbounds = myRgcbounds,
                       propDF.RPMVlos = mypropDF.RPMVlos, propsd.RPMVlos.empir=newcovmat,
                       pm.par=yesnopm,
                       priorfuncs = list(singleunif.prior, gaus.prior, alphaprior, singleunif.prior),
                       ppars = list(phiobounds, c(gammamean, gammasd), apriorpars, bbounds),
                       countSun=TRUE, solarmotion=sunmotion, velocityLSR=vLSR,
                       savedatapars=TRUE, parcovmat=propDFcov)


# get covariance matrices for rgc, vlos, and proper motions
newcovmat = test3[[1]]

# replace the 0 covariance matrices with much smaller step size
thebadcovmats = unlist(lapply(X = newcovmat, FUN = function(x) all(x==matrix(rep(0,16), nrow=4)) ))

if(any(thebadcovmats)){
  y = 1e-9
  newcovmat[thebadcovmats] = lapply(1:sum(thebadcovmats), FUN = function(x) matrix(c(y,0,0,0,
                                                                                     0,y,0,0,
                                                                                     0,0,y,0,
                                                                                     0,0,0,y), nrow=4))
  
}


# get covariance matrix for model parameters using Roberts & Rosenthal method
propDFcov = test3[[2]]
dimnames(propDFcov)=NULL
init.list[[1]] = test3[[3]]

# plot to check chains
plot(as.mcmc(test3[[5]]))




filename

saveRDS(object = test3, file = paste("../Results/adaptive-burnin-", filename, "_", Sys.Date(), sep=""))



