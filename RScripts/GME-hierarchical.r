# Gibbs Hybrid MCMC script

# logPDF   = function, calculates log of the distribution function DF
# DF       = the distribution function of choice
# priors   = the prior function(s)
# propDF   = proposal distribution for parameters
# dat      = matrix of Data 
# n.pars   = number of parameters in model
# N        = length of chain
# init     = LIST of [[1]] inital parameter values and [[2]] initial data parameter values
# vt.par   = TRUE/FALSE indicating if tangential velocities need to be treated as parameters
#          = TRUE, dat = (radius, line of sight velocity); 
#          = FALSE, dat = (radius, speed) OR (radius, line of sight velocity, transerve velocity)
# vt.propsd= standard deviation for vt parameter proposal distribution 
# vt.priors= the prior function(s) for the vt parameters
# n.vt     = number of vt parameters
# vt.notmeasured = TRUE/FALSE indicating if some vt's were not measured
# vt.init  = initial vt parameter values

# default proposal distribution for vt is given by
#  rnorm( n.vt, mean=0, sd = vt.propsd )

# default prior on vt (uniform on vt^2)
# unifvt2.prior <- function( vt, ... ){
#   
#   if( vt<0 ){ return( 0 ) }else{  return( vt ) }
#   
# }

# dat must be in a dataframe and in the following form, and NOT adjust for solar motion:

# RA    (decimal degrees)
# DEC   (decimal degrees)
# PMra  (arcsec/yr)
# PMdec (arcsec/yr)
# Vlos  (km/s)
# plx   (arcsec)
# Rgc   (kpc)
# cosphi
# sinphi
# ePMra (arsec/yr)
# ePMdec(arcsec/yr)
# eVlos (km/s)
# eRgc  (kpc)


# propsd.PMVlos   list of covariance matrices for proper motion and line-of-sight velocity jumping distribution

# transformation matrix (up-to-date transformation matrix using J2000 epoch)
# RIGHT-HANDED UVW

# ---------------- GME ------------------ #
HierGalacticMCMC = function( init, dat, DF, pot, n.pars, priors, propDF, N, logPDF=logDF, 
                             transform.pars,
                  pm.par=FALSE,
                  PMra.priors, PMdec.priors, Vlos.priors, Rgc.priors,
                  propDF.RPMVlos, propsd.RPMVlos.empir,
                  Tmatrix=myTmatrix, countSun, solarmotion, velocityLSR,
                  savedatapars, progressBar=TRUE, thinning = 1, ...){    
  # ... arguments could be e.g. PMbounds, Vlosbounds, Rgcbounds, propsd.PMVlos,  ... etc.
  
  
  # get the number of data points
  n.dat = nrow(dat)
  # replicate original data --> Data, so dat doesn't get overwritten
  # RA, DEC, PMra, PMdec, Vlos, plx, Rgc, cosphi, sinphi,
  # ePMra, ePMdec, eVlos, eRgc
  Data = dat
  
  Data.init = init[[2]]
  init = init[[1]]

  # calculate A, B-matrix and inverse-B matrix, and transformation matrix for velocities
  # (for going from PMs, to Cylindrical i.e. U,V,W to Spherical i.e. vr, vtheta, vphi)
  Bmatrix.list = vector("list",length=n.dat)
  TMatrixVels.list = vector("list", length=n.dat)

  for(i in 1:n.dat){
    Amatrix.temp = Amatrix(RA = dat$RA[i], DEC = dat$DEC[i])
    Bmatrix.list[[i]] = Tmatrix %*% Amatrix.temp
    TMatrixVels.list[[i]] = Transform.Matrix(dat = dat[i, ])
  }

  # object to keep r, vr, vt matrix. This will change at each iteration
  dat.current = matrix(nrow=n.dat, ncol=3)

  # Transform Data.init line-of-sight and proper motion measurements
  # to spherical Galactocentric measurements
  # Data.init has NA's filled in with initial guesses for missing values
  
  for(j in 1:n.dat){
    UVW.temp = PMtoUVW(Bmatrix=Bmatrix.list[[j]], dat.=Data.init[j, ], countSun=countSun, solarmotion =solarmotion, velocityLSR=velocityLSR )  
    # adjust for rotation of Galaxy
    
    PiThetaZ.temp = UVWtoPTZ( UVW=UVW.temp, dat.=Data.init[j, ])
    # transform velocities in right-handed X, Y, Z  to right-handed spherical galactocentric
    vr.vtheta.vphi = ( TMatrixVels.list[[j]] %*% PiThetaZ.temp )/100 # put into units of 100km/s for DF
    # save as r, vr, vt
    dat.current[j, ] = c(Data.init$Rgc[j], vr.vtheta.vphi[1, ], sqrt(vr.vtheta.vphi[2, ]^2 + vr.vtheta.vphi[3, ]^2) )
  }
  
  ### check the initial model and nuisance parameter values
  
  # get logDF values for initial model parameters
  test.init = logPDF( pars=init, dat=dat.current, DF=DF, pot=pot, transform.pars=transform.pars )
  
  # get log prior values for initial model parameters
  logpriors = sum(log( priors( pars = transform.pars(init), ... ) ) )
  
  # need to check that initial values are OK. 
  # if missing PMs then calculate log of nuisance priors first way, if not missing any then calculate log of nuisance priors the second way
  if( pm.par ){

    # find the unknown pm's (note if PMra not known, then PMdec not known either)
    pm.unknown <- is.na(dat[, "PMra"])
    n.pm <- sum(pm.unknown)
    pm.unknown <- which(pm.unknown[])

    # get prior values for initial parameters for Rgc, PMs, and Vlos
    lognuispriors = ( ( sum( log( dnorm(x = dat$Rgc, mean=Data.init$Rgc, sd=dat$eRgc) ) )
        + sum( log( dnorm(x = dat$Vlos, mean=Data.init$Vlos, sd=dat$eVlos) ) )
        + sum( log( dnorm(x = dat$PMra[-pm.unknown], mean=Data.init$PMra[-pm.unknown], sd=dat$ePMra[-pm.unknown]) ) )
        + sum( log( dnorm(x = dat$PMdec[-pm.unknown], mean=Data.init$PMdec[-pm.unknown], sd=dat$ePMdec[-pm.unknown]) ) )
      + sum( log( Rgc.priors(Rgc=Data.init$Rgc, ...) ) )
      + sum( log( PMra.priors(PMra=Data.init$PMra[-pm.unknown], ... ) ) )
      + sum( log( PMdec.priors(PMdec=Data.init$PMdec[-pm.unknown], ... ) ) )
      + sum( log( Vlos.priors(Vlos=Data.init$Vlos, ... ) ) )
      ) )

  
    }else{ # otherwise, if no PMs missing, then...

    lognuispriors = ( sum( log( dnorm(x = dat$Rgc, mean=Data.init$Rgc, sd=dat$eRgc) ) )
                      + sum( log( dnorm(x = dat$Vlos, mean=Data.init$Vlos, sd=dat$eVlos) ) )
                      + sum( log( dnorm(x = dat$PMra, mean=Data.init$PMra, sd=dat$ePMra) ) )
                      + sum( log( dnorm(x = dat$PMdec, mean=Data.init$PMdec, sd=dat$ePMdec) ) )
                      + sum( log( Rgc.priors(Rgc=Data.init$Rgc, ...) ) )
                      + sum( log( PMra.priors(PMra=Data.init$PMra, ... ) ) )
                      + sum( log( PMdec.priors(PMdec=Data.init$PMdec, ... ) ) )
                      + sum( log( Vlos.priors(Vlos=Data.init$Vlos, ... ) ) )
                    )
  }
    
  # now check that the initial values are OK
  if( any( !is.finite(logpriors) ) ){ stop("bad initial model parameters for priors") }

    if( any( !is.finite(lognuispriors) ) ){ stop("bad initial Rgc, PMs, and/or Vlos parameters") }

  if( any( !is.finite(test.init) ) ){  stop("bad initial model parameters")   } 

  # Assuming all the initial values are OK,
  # set up the parameter chain
  chain = matrix( ncol = n.pars , nrow = N )
  chain[1, ] = init
  accept = 0

  if( savedatapars ){ # set up the Rgc, PMra, PMdec, and Vlos parameter chains
   
     r.chain = matrix( ncol=n.dat, nrow=N )
    PMra.chain = matrix( ncol = n.dat, nrow = N )
    PMdec.chain = matrix( ncol = n.dat, nrow = N )
    Vlos.chain = matrix( ncol = n.dat, nrow=N)
  
    r.chain[1, ] = Data.init$Rgc
    PMra.chain[1, ] = Data.init$PMra
    PMdec.chain[1, ] = Data.init$PMdec
    Vlos.chain[1, ] = Data.init$Vlos

    # keep track of the proper motion and V_los parameter acceptance rates
    RPMVlos.accept = rep(0, length.out=n.dat)
    
  }

  if( progressBar ){
    # make a progress bar
    pb <- txtProgressBar(min = 0, max = N*thinning, style = 3)
  }

  # make the Markov Chain
  for(i in 2:N*thinning){
    
    # draw new values for model parameters, by taking a step from the previous pars
    xtry = init + propDF(...)
   
    # calculate the priors for parameters
    logpriors.init = sum( log( priors( pars = transform.pars( init ), ... ) ) )
    logpriors.try = sum( log( priors( pars = transform.pars( xtry ), ... ) ) )

    # accept/reject the (M,a) pair
    difflog = ( logPDF( pars=xtry, dat=dat.current, DF=DF,  pot=pot, transform.pars ) + logpriors.try
                  - ( logPDF( pars=init, dat=dat.current, DF=DF, pot=pot, transform.pars ) + logpriors.init ) )

    ###    
 if( !is.numeric(difflog)){browser()}
    ###
    
    if( difflog > 0 | ( exp(difflog) > runif(1) ) ){
        
      if(is.whole(i/thinning)){  chain[i/thinning, ] = xtry  }
      
      init = xtry
      accept = accept + 1
        
    }else{ 
      
      if(is.whole(i/thinning)){  chain[i/thinning, ] = init  }
      
      }
  
# ---------------------------------------
    # accept/reject r parameters based on current M, a, PMra, PMdec, parameters
    # priors for the r's not jth tracer will cancel, b/c r's the same in both cases
#     r.init = dat.current[, 1]
#     logr.try = log(r.init) + rnorm(n.dat, mean=0, sd=propsd.Rgc)
#     r.try = exp(logr.try)

    # set up object to update
    dat.try = dat.current
    Data.try = Data.init

    for( k in 1:n.dat){
      
      
      # draw new proper motions and line of sight velocities
      RPMVlos.try = Data.init[k, c("Rgc", "PMra", "PMdec", "Vlos")] + propDF.RPMVlos(sig=propsd.RPMVlos.empir[[k]], ...)
      RPMVlos.try = as.numeric(RPMVlos.try)
      Data.try[k, c("Rgc", "PMra", "PMdec", "Vlos")] = RPMVlos.try
            
      # -----------------------
      # calculate new vr and vt
        UVW.temp = PMtoUVW(Bmatrix=Bmatrix.list[[k]], dat.=Data.try[k, ], countSun=countSun, solarmotion=solarmotion, velocityLSR=velocityLSR)
        #UVW.temp[2, ] = UVW.temp[2, ] + 220 
        PiThetaZ.temp = UVWtoPTZ( UVW = UVW.temp, dat.=Data.try[k, ])
        v.temp = ( TMatrixVels.list[[k]] %*% PiThetaZ.temp )/100 # put into units of 100km/s for DF
        # save as r, vr, vt
        dat.try[k, ] = c(Data.try$Rgc[k], v.temp[1, ], sqrt(v.temp[2, ]^2 + v.temp[3, ]^2) )
      # -----------------------      
      
    # calculate the logPDF at the new place in parameter space and at the old place in parameter space
    if(pm.par){
          
      
      if( any(k==pm.unknown) ){
        
        # calculate logPDF for the new value of vt
        newlogPDF = ( log( dnorm( x = dat$Vlos[k], mean=Data.try$Vlos[k], sd=dat$eVlos[k] ) )
                      + log( dnorm( x = dat$Rgc[k], mean=Data.try$Rgc[k], sd=dat$eRgc[k]) )
                      
                      + logPDF( pars=init, dat=dat.try, DF=DF, pot=pot, transform.pars=transform.pars ) 
                      
                      + log( Rgc.priors(Rgc=Data.try$Rgc[k], ...) )
                      + log( PMra.priors(PMra=Data.try$PMra[k], ... ) )
                      + log( PMdec.priors(PMdec=Data.try$PMdec[k], ... ) )
                      + log( Vlos.priors(Vlos=Data.try$Vlos[k], ... ) )
                      
        )
        
        # calculate the current logPDF value for jth particle
        oldlogPDF = ( log( dnorm( x = dat$Vlos[k], mean=Data.init$Vlos[k], sd=dat$eVlos[k] ) )
                      + log( dnorm( x = dat$Rgc[k], mean=Data.init$Rgc[k], sd=dat$eRgc[k]) )
                      
                      + logPDF( pars=init, dat=dat.current, DF=DF, pot=pot, transform.pars=transform.pars)
                      
                      + log( Rgc.priors(Rgc=Data.init$Rgc[k], ...) )
                      + log( PMra.priors(PMra=Data.init$PMra[k], ... ) )
                      + log( PMdec.priors(PMdec=Data.init$PMdec[k], ... ) ) 
                      + log( Vlos.priors(Vlos=Data.init$Vlos[k], ... ) )
        )
        
      
        }else{
        newlogPDF = ( log( dnorm(x = dat$PMra[k], mean=Data.try$PMra[k], sd=dat$ePMra[k]) )
                      + log( dnorm(x = dat$PMdec[k], mean=Data.try$PMdec[k], sd=dat$ePMdec[k]) )
                      + log( dnorm(x = dat$Vlos[k], mean=Data.try$Vlos[k], sd=dat$eVlos[k]) )
                      + log( dnorm( x = dat$Rgc[k], mean=Data.try$Rgc[k], sd=dat$eRgc[k]))
                      
                      + logPDF( pars=init, dat=dat.try, DF=DF, pot=pot, transform.pars=transform.pars ) 
                      
                      + log( Rgc.priors(Rgc=Data.try$Rgc[k], ...) )
                      + log( PMra.priors(PMra=Data.try$PMra[k], ... ) )
                      + log( PMdec.priors(PMdec=Data.try$PMdec[k], ... ) )
                      + log( Vlos.priors(Vlos=Data.try$Vlos[k], ... ) )
        )
        
        # calculate the current logPDF value for jth particle
        oldlogPDF = ( log( dnorm(x = dat$PMra[k], mean=Data.init$PMra[k], sd=dat$ePMra[k]) )
                      + log( dnorm(x = dat$PMdec[k], mean=Data.init$PMdec[k], sd=dat$ePMdec[k]) )
                      + log( dnorm(x = dat$Vlos[k], mean=Data.init$Vlos[k], sd=dat$eVlos[k]) )
                      + log( dnorm( x = dat$Rgc[k], mean=Data.init$Rgc[k], sd=dat$eRgc[k]) )
                      
                      + logPDF( pars=init, dat=dat.current, DF=DF, pot=pot, transform.pars = transform.pars) 
                      
                      + log( Rgc.priors(Rgc=Data.init$Rgc[k], ...) )
                      + log( PMra.priors(PMra=Data.init$PMra[k], ... ) )
                      + log( PMdec.priors(PMdec=Data.init$PMdec[k], ... ) ) 
                      + log( Vlos.priors(Vlos=Data.init$Vlos[k], ... ) )
        )
        }
      
      }else{
        
        newlogPDF = ( log( dnorm(x = dat$PMra[k], mean=Data.try$PMra[k], sd=dat$ePMra[k]) )
                      + log( dnorm(x = dat$PMdec[k], mean=Data.try$PMdec[k], sd=dat$ePMdec[k]) )
                      + log( dnorm(x = dat$Vlos[k], mean=Data.try$Vlos[k], sd=dat$eVlos[k]) )
                      + log( dnorm( x = dat$Rgc[k], mean=Data.try$Rgc[k], sd=dat$eRgc[k]))
                      
                      + logPDF( pars=init, dat=dat.try, DF=DF, pot=pot, transform.pars=transform.pars ) 
                      
                      + log( Rgc.priors(Rgc=Data.try$Rgc[k], ...) )
                      + log( PMra.priors(PMra=Data.try$PMra[k], ... ) )
                      + log( PMdec.priors(PMdec=Data.try$PMdec[k], ... ) )
                      + log( Vlos.priors(Vlos=Data.try$Vlos[k], ... ) )
        )
        
        # calculate the current logPDF value for jth particle
        oldlogPDF = ( log( dnorm(x = dat$PMra[k], mean=Data.init$PMra[k], sd=dat$ePMra[k]) )
                      + log( dnorm(x = dat$PMdec[k], mean=Data.init$PMdec[k], sd=dat$ePMdec[k]) )
                      + log( dnorm(x = dat$Vlos[k], mean=Data.init$Vlos[k], sd=dat$eVlos[k]) )
                      + log( dnorm( x = dat$Rgc[k], mean=Data.init$Rgc[k], sd=dat$eRgc[k]) )
                      
                      + logPDF( pars=init, dat=dat.current, DF=DF, pot=pot, transform.pars = transform.pars) 
                      
                      + log( Rgc.priors(Rgc=Data.init$Rgc[k], ...) )
                      + log( PMra.priors(PMra=Data.init$PMra[k], ... ) )
                      + log( PMdec.priors(PMdec=Data.init$PMdec[k], ... ) ) 
                      + log( Vlos.priors(Vlos=Data.init$Vlos[k], ... ) )
        )
        
      }
      
        
      # find difference in logPDF values
      difflog = newlogPDF - oldlogPDF

      # if difference in logs is positive, or if the ratio of the posteriors is greater than an random number generated from a uniform distribution between 0 and 1, then accept the new parameters
      ###
      if( !is.numeric(difflog) | is.na(difflog) ){browser()}
      ###
      if( difflog > 0 | (exp(difflog) > runif(1)) ){
    
        if( savedatapars ){
      
          if(is.whole(i/thinning)){ 
            
            r.chain[i/thinning, k] = Data.try$Rgc[k]
            PMra.chain[i/thinning, k] = Data.try$PMra[k]
            PMdec.chain[i/thinning, k] = Data.try$PMdec[k]
            Vlos.chain[i/thinning, k] = Data.try$Vlos[k]
            
          }
          
          RPMVlos.accept[k] = RPMVlos.accept[k] + 1
          
        }    
        
        dat.current[k, ] <- dat.try[k, ]
        Data.init[k, ] <- Data.try[k, ]
                
      }else{ # if not, then save the old values
        
        if( savedatapars ){
          
          if(is.whole(i/thinning)){
           
            r.chain[i/thinning, k] = Data.init$Rgc[k]
            PMra.chain[i/thinning, k] = Data.init$PMra[k]
            PMdec.chain[i/thinning, k] = Data.init$PMdec[k]
            Vlos.chain[i/thinning, k] = Data.init$Vlos[k]
          }
        }
        
        dat.try[k, ] = dat.current[k, ]
        
      }
    }

# ---------------------------------------
    
    if( progressBar ){ setTxtProgressBar(pb, i/thinning) }
    
  } # close for loop (N-1 iterations --- once this is done the chain is done)


    if( savedatapars ){
      chain.colnames = c( paste("par", 1:n.pars, sep=""), paste("Rgc", 1:n.dat, sep=""), paste("PMra", 1:n.dat, sep=""),
                          paste("PMdec", 1:n.dat, sep=""), paste("Vlos", 1:n.dat, sep="") )
      
      allchains = cbind(chain, r.chain, PMra.chain, PMdec.chain, Vlos.chain)
      
    
      }else{
      chain.colnames = paste("par", 1:n.pars, sep="")
      allchains = chain
    }
    
    colnames(allchains) = chain.colnames
    
    if( savedatapars ){
      out = list( chain=allchains, acceptance.rate=accept/N,
                  acceptRPMVlos.rate = RPMVlos.accept/N, last.Data=Data.init )
      
      print( paste( "R-PM-Vlos acceptance rate = " , out$acceptRPMVlos.rate , sep="" ) )
    
      }else{
      out = list( chain=allchains, acceptance.rate=accept/N, last.Data=Data.init )
    }
    
  
  print( paste( "acceptance rate = ",out$acceptance.rate, sep="" ) )
  out
  
} # close function
