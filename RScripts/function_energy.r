# function to calculate the specific energy of a tracer in Galactocentric coordinates, given a model for the gravitational potential

# particle's relative specific energy. dat must be in the form r, vr, and vt
# where r is the Galactocentric distance, vr is the radial velocity, and vt = sqrt(v_theta^2 + v_phi^2)

relE <- function( dat, pars, pot){
  
  # if dat is a vector, turn it into a matrix
  if( is.vector(dat) ){ dat = matrix( dat , ncol=length(dat) ) }
  
  vr = dat[,2]
  vt = dat[,3]
  
  -( vr^2 + vt^2 )/2  + pot(dat=dat, pars=pars)
  
}
