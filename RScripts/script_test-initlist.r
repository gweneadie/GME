# initial values for parameters: phi_o, gamma, alpha, beta
start.list = list( c(45, 0.4, 3.2, 0.5), c(50, 0.35, 3.1, 0), c(48, 0.42, 3.6, 0.4))
# transform to space used for sampling
par.init <- lapply( start.list, FUN=inv.transform )


# initial guesses for pm's (take fraction of measurements)
initdata = mydata
initdata[miss.pm, c("PMra", "PMdec")] = mean( na.omit(mydata$PMra) )/10

# all initial parameters needs to be a list together
init.list = lapply(X = par.init, FUN = function(x, idata) list( x, idata ), idata=initdata )

# initial proposal distribution for model parameters (zeros in off-diagonals)
mypropDF <- function(x, ...) rnorm( 4, mean=c(0,0,0,0) , sd=c(0.07, 0.05, 0.01, 0.08))

# function for proposal distribution for r, vlos, etc.
mypropDF.RPMVlos = function(sig,...){  mvrnorm(Sigma = sig, mu=c(0,0,0,0)) }

y = 1e-6
mypropsd.RPMVlos = lapply( X=1:nrow(mydata), FUN=function(x) matrix(c(y,0,0,0, 
                                                                      0,y,0,0,
                                                                      0,0,y,0,
                                                                      0,0,0,y), nrow=4) )

# start running a single Markov Chain to test the initial parameter values
# it will return an error if the initial values make the DF negative or Inf

# number of iterations in test chain
iterations = 10


test <- HierGalacticMCMC(dat=mydata, init=init.list[[3]], N=iterations, 
                         n.pars=npars,
                         logPDF=logDF, DF=modelDF, pot=modelpot, 
                         propDF=mypropDF,
                         priors=prior.wrapper, transform.pars=transform.func,
                         pm.par=yesnopm,
                         PMra.priors=priorfunc.PMra, PMdec.priors=priorfunc.PMdec,
                         Vlos.priors=priorfunc.Vlos, Rgc.priors=priorfunc.Rgc,
                         PMbounds = myPMbounds, Vlosbounds = myVlosbounds, Rgcbounds = myRgcbounds,
                         propDF.RPMVlos = mypropDF.RPMVlos, propsd.RPMVlos.empir=mypropsd.RPMVlos,
                         Tmatrix=myTmatrix,
                         priorfuncs = list(singleunif.prior, gaus.prior, alphaprior, singleunif.prior),
                         ppars = list(phiobounds, c(gammamean, gammasd), apriorpars, bbounds),
                         countSun=TRUE, solarmotion=sunmotion, velocityLSR=vLSR, savedatapars = TRUE,
                         covmat=propDFcov)
