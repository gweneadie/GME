# function to adjust vt proposal distributions

# beta.AM    small positive constant for adaptive metropolis part, default 0.05
#
adaptiveBurnin <- function( acceptrange = c(0.25, 0.35), acceptnuisrange=c(0.25, 0.40), changestep, yourpatience=25,
                           init, dat, DF, pot, n.pars, priors, propDF, N, logPDF=logDF, transform.pars,
                           pm.par=TRUE,
                           PMra.priors, PMdec.priors, Vlos.priors, Rgc.priors,
                           propDF.RPMVlos, propsd.RPMVlos.empir,
                           Tmatrix=myTmatrix, parcovmat, ...){

  # start the count
  count = 1
  # get range for acceptance rates
  accept.lower <- acceptnuisrange[1]
  accept.upper <- acceptnuisrange[2]
  acceptpar.lower <- acceptrange[1]
  acceptpar.upper <- acceptrange[2]
  
  # make sequence of numbers for tracers and number of pars
  seq.tracers = 1:nrow(dat)
  seq.pars = 1:n.pars
  
  init.try <- init
  
  # run chain for 100 steps
  chaintry <- HierGalacticMCMC( init=init.try, dat, DF, pot, n.pars, priors, propDF, N=N, logPDF=logDF, transform.pars=transform.pars,
                                pm.par=TRUE,
                                PMra.priors, PMdec.priors, Vlos.priors, Rgc.priors,
                                propDF.RPMVlos, propsd.RPMVlos.empir,
                                Tmatrix=myTmatrix, covmat = parcovmat, ...)
  
  # get acceptance rates for nuisance pars and model pars
  acceptRPMVlos.rates <- chaintry$acceptRPMVlos.rate
  accept.rates <- chaintry$acceptance.rate
  
  # see if they are all good
  RPMVlosaccept.good <- all( (acceptRPMVlos.rates > accept.lower) & (acceptRPMVlos.rates < accept.upper) )
  parsaccept.good <- all( (accept.rates > acceptpar.lower) & (accept.rates < acceptpar.upper) )
  
  # change the initial values, using the last position of the previous chain
  init.try[[1]] <- as.numeric( chaintry$chain[N, 1:n.pars] )
  init.try[[2]][, c("Rgc", "PMra", "PMdec", "Vlos")] <- matrix( chaintry$chain[N, ( (n.pars+1):(4*nrow(dat)+n.pars))], nrow=nrow(dat))
  
  # make object to update covariance matrices that need changing
  new.covmatrices = propsd.RPMVlos.empir
  new.parcovmat = parcovmat

  # make combineghains object in case while loop isn't used
  combinechains=NULL
  
  # while the acceptances are not good
  while( ( (!RPMVlosaccept.good ) | (!parsaccept.good) ) & count<yourpatience ){

    # find problem rates for r, proper motions and vlos
    smallRPMVlos <- (acceptRPMVlos.rates < accept.lower)
    bigRPMVlos <- (acceptRPMVlos.rates > accept.upper)

    # vector of covariance matrices to change
    changethese = seq.tracers[smallRPMVlos | bigRPMVlos]

    if( count==1 ){
      
      # calculate the new covariance matrix for the PMVlos chains, using all the previous steps
      for( m in changethese){
        new.covmatrices[[m]] = cov( chaintry$chain[, paste(c("Rgc", "PMra", "PMdec", "Vlos"), m, sep="")]  )
      }
      
      # if the model parameter acceptance is not good, then calculate new covariance matrix
      if( !parsaccept.good ){
        new.parcovmat = cov( chaintry$chain[, 1:npars] )
        dimnames(new.parcovmat) = NULL
      }
      
    }else{
      # calculate the new covariance matrix for the PMVlos chains, using all the previous steps
      for( p in changethese){
        new.covmatrices[[p]] = cov( combinechains[, paste(c("Rgc", "PMra", "PMdec", "Vlos"), p, sep="")]  )
        
        # if there are some bad ones that have gone to zero, fix them
        # replace the 0 covariance matrices with much smaller step size
        thebadcovmats = unlist(lapply(X = new.covmatrices, FUN = function(x) all(x==matrix(rep(0,16), nrow=4)) ))
        
        if(any(thebadcovmats)){
          y = 1e-10
          new.covmatrices[thebadcovmats] = lapply(1:sum(thebadcovmats), FUN = function(x) matrix(c(y,0,0,0,
                                                                                             0,y,0,0,
                                                                                             0,0,y,0,
                                                                                             0,0,0,y), nrow=4))
          
        }
        
      }
      
      # if the model parameter acceptance is not good, then calculate new covariance matrix
      if( !parsaccept.good ){
        new.parcovmat = cov( combinechains[, 1:npars] )
        dimnames(new.parcovmat) = NULL
      }
      
    }

    # retry chain with new Rgc step size and new covariance matrix
    chaintryagain <- HierGalacticMCMC(  init=init.try, dat=dat, DF, pot, n.pars, priors, propDF, N, logPDF=logDF, transform.pars,
                                   pm.par=TRUE,
                                   PMra.priors, PMdec.priors, Vlos.priors, Rgc.priors,
                                   propDF.RPMVlos, propsd.RPMVlos.empir=new.covmatrices,
                                   Tmatrix=myTmatrix, covmat  = new.parcovmat, ...)
 
    # get acceptance rates
    acceptRPMVlos.rates <- chaintryagain$acceptRPMVlos.rate
    accept.rates <- chaintryagain$acceptance.rate
    
    # see if they are all good
    RPMVlosaccept.good <- all( (acceptRPMVlos.rates > accept.lower) & (acceptRPMVlos.rates < accept.upper) )
    parsaccept.good <- all( (accept.rates > acceptpar.lower) & (accept.rates < acceptpar.upper) )
    
    # change the initial values, using the last position of the previous chain
    init.try[[1]] <- as.numeric( chaintryagain$chain[N, 1:n.pars] )
    init.try[[2]][, c("Rgc", "PMra", "PMdec", "Vlos")] <- matrix( chaintryagain$chain[N, ( (n.pars+1):(4*nrow(dat)+n.pars))], nrow=nrow(dat))
  
    # add new chain to old chain
    if( count==1 ){
      combinechains = rbind(chaintry$chain, chaintryagain$chain )
    }else{
      combinechains = rbind(combinechains, chaintryagain$chain)
    }
    # Patience count
    count = count + 1
    
  }
  
  if(yourpatience==1){
    chaintryagain = chaintry
    combinechains = chaintry$chain
  }
  

  out = list(new.covmatrices, new.parcovmat, init.try, count, combinechains, chaintryagain$acceptRPMVlos.rate, chaintryagain$acceptance.rate)
  
  out
}

  
  

  
