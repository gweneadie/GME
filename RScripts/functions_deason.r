# Deason, Evans model

# total gravitational potential, the first column of dat must be the Galactocentric distance r
pot.Deason = function(dat, pars){
  p = pars[1] # potential constant
  g = pars[2] # potential slope
  a = pars[3] # density slope
  
  # if dat is a vector, turn it into a matrix
  if( is.vector(dat) ){ dat = matrix( dat , ncol=length(dat) ) }
    
  p/(dat[, 1]^g)
  
}


# tracer number density
n.r = function(r, pars){  1/(r^pars[3])  }
  
# cumulative mass profile
M.r = function(r, pars){
  
  if( !is.vector(pars) ){
    p = pars[,1] # potential constant
    g = pars[,2] # potential slope
  }else{
    p = pars[1]
    g = pars[2]  }

  g*p*r^(1-g)
  
}

# dark matter
r.M = function(M, pars){
  p = pars[1] # potential constant parameter
  g = pars[2] # potential slope parameter
  
  (M/(g*p))^(1/(1-g))
  
}

# log of the distribution function
logDF.Deason = function(pars, dat, pot=pot.Deason, DF=NULL, transform.pars=default.transform){
  
  pars = transform.pars( pars )
  
  p = pars[1] # potential constant
  g = pars[2] # potential slope
  a = pars[3] # density slope
  b = pars[4] # anisotropy parameter
  
  if( ( a - 2*b ) <= ( g*(-b + 0.5) ) ){ return(-Inf) }
  
  if( b >=1 ){ return(-Inf) }
  
  if( is.vector(dat) ){ dat = matrix( dat, ncol=3 ) }
  
  Es = relE(dat = dat, pars = pars, pot = pot)

  if( any( dat[, 3] < 0 ) ){ return(-Inf) }
  if( any( (Es < 0) | !is.finite(Es) | !is.finite(log(dat[, 3])) ) ){ return(-Inf) }
 
  logeta = (2*b/g - a/g)*log(p) - log( pi*sqrt(pi)*2^( -b+3/2 ) ) - lgamma(x = (1-b)) +
    lgamma( x = (a/g - 2*b/g + 1) ) - lgamma( x = ( b*(g-2)/g + a/g - 1/2 ) )
    
  output = logeta - 2*b*log( dat[, 1]*dat[, 3] ) + ( b*(g-2)/g + a/g - 3/2 )*log(Es)
  
  sum(output)
  
}
##################

# function to check if a is a whole number. Used in mymcmc for thinning chains
is.whole <- function(a) { 
  (is.numeric(a) && floor(a)==a) ||
    (is.complex(a) && floor(Re(a)) == Re(a) && floor(Im(a)) == Im(a))
}

# plotting DF function
plotlogDF.Deason = function(pars, dat, pot=pot.Deason, DF=NULL, transform.pars=default.transform){
  
  pars = transform.pars( pars )
  
  p = pars[1] # potential constant
  g = pars[2] # potential slope
  a = pars[3] # density slope
  b = pars[4] # anisotropy parameter
  
  if( ( a - 2*b ) <= ( g*(-b + 0.5) ) ){ return(-Inf) }
  
  if( b >=1 ){ return(-Inf) }
  
  if( is.vector(dat) ){ dat = matrix( dat, ncol=3 ) }
  
  Es = relE(dat = dat, pars = pars, pot = pot)
  
  if( any( dat[, 3] < 0 ) ){ return(-Inf) }
  if( any( (Es < 0) | !is.finite(Es) | !is.finite(log(dat[, 3])) ) ){ return(-Inf) }
  
  logeta = (2*b/g - a/g)*log(p) - log( pi*sqrt(pi)*2^( -b+3/2 ) ) - lgamma(x = (1-b)) +
    lgamma( x = (a/g - 2*b/g + 1) ) - lgamma( x = ( b*(g-2)/g + a/g - 1/2 ) )
  
  output = logeta - 2*b*log( dat[, 1]*dat[, 3] ) + ( b*(g-2)/g + a/g - 3/2 )*log(Es)
  
 output 
}

# plotting isotopric DF function using energy as variable
logDFenergy.Deason = function(Es, pars, pot=pot.Deason){

  p = pars[,1] # potential constant
  g = pars[,2] # potential slope
  a = pars[,3] # density slope
  b = pars[,4] # anisotropy parameter
  
  if( any( ( a - 2*b ) <= ( g*(-b + 0.5) ) ) ){ return(-Inf) }
  
  if( any(b >=1) ){ return(-Inf) }
  
  logeta = (2*b/g - a/g)*log(p) - log( pi*sqrt(pi)*2^( -b+3/2 ) ) - lgamma(x = (1-b)) +
    lgamma( x = (a/g - 2*b/g + 1) ) - lgamma( x = ( b*(g-2)/g + a/g - 1/2 ) )
  
  output = logeta + ( b*(g-2)/g + a/g - 3/2 )*log(Es)
  
  output 
}

Watkins.Mest = function(dat, Watpars){
  
  a=Watpars[1] ; g = Watpars[2]; b=Watpars[3]
  rout = max(dat[, 1])
  Const = (a+g-2*b)*rout^(1-a)/(3-2*b)
  
  2.325e9 * Const * mean( ( dat[, 2]^2 + dat[, 3]^2 )*dat[, 1]^a )
  
}

# virial radius and virial mass
rvir = function( chain, H=0.678e-3 ){
  
  phio = chain[, 1]
  gam = chain[, 2]
  
  ( gam*phio/(100*H^2) )^( 1/(2+gam) )

}

# function for plotting Deason potenial
plotpot.deason = function(r, pars){
  
  if( !is.vector(pars) ){
    p = pars[, 1]
    g = pars[, 2]
  }else{
    p = pars[1] # potential constant
    g = pars[2] # potential slope
  }  
  
  p/(r^g)
  
}


