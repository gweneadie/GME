# functions for possible prior distribtions on model parameters

######## wrapper for prior when different prior distributions are being used for model parameters (e.g., phi_o, gamma, alpha, and beta)

prior.wrapper = function(pars, priorfuncs, ppars, ...){
  
  if( length(priorfuncs) != length(pars) ){ stop("Error: number of prior functions doesn't match number of parameters")}
  
  # make object for output
  priorvalues = rep(NA, length=length(pars))
  
  for( i in 1:length(pars) ){
    priorvalues[i] = priorfuncs[[i]](pars[i], ppars[[i]])
  }
  
  priorvalues
  
}


# prior on vt (uniform on vt^2)
unifvt2.prior <- function( vt, ... ){
  
  if( any(vt<0) ){ return( 0 ) }else{  return( vt ) }
  
}


######### Uniform (truncated) prior function
singleunif.prior <- function( pars, ppars ){
  
  # bounds of truncated prior
  par.min = ppars[1]
  par.max = ppars[2]
  
  # if parameter is held fixed...
  if( par.min==par.max ){
    priorvalue = 1
  }else{
    # calculated uniform prior value
    priorvalue = 1/(par.max - par.min) * ( par.min <= pars  &  pars <= par.max )
  }
  
  priorvalue
  
}

########### Gaussian prior
gaus.prior <- function( pars, ppars ){
  
  dnorm(x = pars, mean = ppars[1], sd = ppars[2])
  
}

########### Gamma prior (in regular space) on alpha

alphaprior <- function( pars, ppars ){
  
  rmin = ppars[1]
  myc = ppars[2]
  p = ppars[3]
  n = ppars[4]
  nco = ppars[5]
  
  eta = pars - 3
  
  dgamma(x = eta, shape = (n+myc), scale=1/(nco+p))
  
}

########## alpha prior, using other data not used in the analysis
# Pareto cumulative distribution function
r.inverseParetocdf = function(x, A=2.5, xmin){  xmin/((1-x)^(1/A)) }

# Pareto distribution function
Paretopdf = function(x, A=2.5, xmin){ A*xmin^A/x^(A+1) }


############ General gamma prior

gamma.prior <- function( pars, ppars ){
  
  eta = pars - 3
  
  dgamma(x = eta, shape = ppars[1], rate = ppars[2] )
  
}

