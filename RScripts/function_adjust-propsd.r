# function to adjust vt proposal distributions

# beta.AM    small positive constant for adaptive metropolis part, default 0.05
#
adjust.propsd <- function( acceptrange=c(0.20, 0.30), changestep, yourpatience=25,
                           init, dat, DF, pot, n.pars, priors, propDF, N, logPDF=logDF, transform.pars,
                           pm.pars,
                           PMra.priors, PMdec.priors, Vlos.priors, Rgc.priors,
                           propDF.RPMVlos, propsd.RPMVlos.empir,
                           Tmatrix=myTmatrix, ...){

  # start the count
  count = 1
  # get range for acceptance rates
  accept.lower <- acceptrange[1]
  accept.upper <- acceptrange[2]
  
  # make sequence of numbers for tracers
  seq.tracers = 1:nrow(dat)
  
  init.try <- init
  
  # run chain for 100 steps
  chaintry <- HierGalacticMCMC( init=init.try, dat, DF, pot, n.pars, priors, propDF, N=N, logPDF=logDF, transform.pars=transform.pars,
                                pm.par=pm.pars,
                                PMra.priors, PMdec.priors, Vlos.priors, Rgc.priors,
                                propDF.RPMVlos, propsd.RPMVlos.empir,
                                Tmatrix=myTmatrix, ...)
  
  # get acceptance rates
  acceptRPMVlos.rates <- chaintry$acceptRPMVlos.rate
  
  # see if they are all good
  RPMVlosaccept.good <- all( (acceptRPMVlos.rates > accept.lower) & (acceptRPMVlos.rates < accept.upper) )
  
  # change the initial values, using the last position of the previous chain
  init.try[[1]] <- as.numeric( chaintry$chain[N, 1:n.pars] )
  init.try[[2]][, c("Rgc", "PMra", "PMdec", "Vlos")] <- matrix( chaintry$chain[N, ( (n.pars+1):(4*nrow(dat)+n.pars))], nrow=nrow(dat))
  
  # make object to update covariance matrices that need changing
  new.covmatrices = propsd.RPMVlos.empir

  # make combineghains object in case while loop isn't used
  combinechains=NULL
  
  # while the acceptances are not good
  while( (!RPMVlosaccept.good ) & count<yourpatience ){

     # find problem rates for proper motions and vlos
    smallRPMVlos <- (acceptRPMVlos.rates < accept.lower)
    bigRPMVlos <- (acceptRPMVlos.rates > accept.upper)
     
    # vector of covariance matrices to change
    changethese = seq.tracers[smallRPMVlos | bigRPMVlos]
      
    if( count==1 ){
      # calculate the new covariance matrix for the PMVlos chains, using all the previous steps
      for( m in changethese){
        new.covmatrices[[m]] = cov( chaintry$chain[, paste(c("Rgc", "PMra", "PMdec", "Vlos"), m, sep="")]  )
      }
    }else{
      # calculate the new covariance matrix for the PMVlos chains, using all the previous steps
      for( p in changethese){
        new.covmatrices[[p]] = cov( combinechains[, paste(c("Rgc", "PMra", "PMdec", "Vlos"), p, sep="")]  )
      }
    }

    # retry chain with new Rgc step size and new covariance matrix
    chaintryagain <- HierGalacticMCMC(  init=init.try, dat=dat, DF, pot, n.pars, priors, propDF, N, logPDF=logDF, transform.pars,
                                   pm.par=pm.pars,
                                   PMra.priors, PMdec.priors, Vlos.priors, Rgc.priors,
                                   propDF.RPMVlos, propsd.RPMVlos.empir=new.covmatrices,
                                   Tmatrix=myTmatrix, ...)
 
    # get acceptance rates
    acceptRPMVlos.rates <- chaintryagain$acceptRPMVlos.rate
        
    # see if they are all good
    RPMVlosaccept.good <- all( (acceptRPMVlos.rates > accept.lower) & (acceptRPMVlos.rates < accept.upper) )
  

    # change the initial values, using the last position of the previous chain
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
  

  out = list(new.covmatrices, init.try, count, combinechains, chaintryagain$acceptRPMVlos.rate)
  
  out
}

  
  

  
