### PRIORS ###
# bounds for PM and Vlos
priorfunc.PMra = function( PMra, PMbounds, ... ){
  PMmin = PMbounds[1]
  PMmax = PMbounds[2]

  1/( PMmax - PMmin ) * ( PMmin <= PMra  &  PMra <= PMmax )

}

priorfunc.PMdec = function( PMdec, PMbounds, ...  ){
  PMmin = PMbounds[1]
  PMmax = PMbounds[2]

  1/( PMmax - PMmin ) * ( PMmin <= PMdec  &  PMdec <= PMmax )

}

priorfunc.Vlos = function( Vlos, Vlosbounds, ...  ){
  Vlosmin = Vlosbounds[1]
  Vlosmax = Vlosbounds[2]

  1/( Vlosmax - Vlosmin ) * ( Vlosmin <= Vlos & Vlos <= Vlosmax )
}

priorfunc.Rgc = function( Rgc, Rgcbounds, ...  ){
  Rgcmin = Rgcbounds[1]
  Rgcmax = Rgcbounds[2]

  1/( Rgcmax - Rgcmin ) * ( Rgcmin <= Rgc & Rgc <= Rgcmax )
}





# Uniform Jeffery's Prior (in regular space).
# Bounds must be a vector of length two (lower, upper)
unif.Jeffprior <- function( pars, bounds, ... ){
  
  par.min = bounds[,1]
  par.max = bounds[,2]
  
  pars/( par.max - par.min ) * ( par.min <= pars  &  pars <= par.max )
  
}


# # uniform prior for r , vr
# r.prior <- function( pars, rbounds, ... ){
#   
#   par.min = rbounds[,1]
#   par.max = rbounds[,2]
#   
#   1/( par.max - par.min ) * ( par.min <= pars  &  pars <= par.max )
#   
# }
# 
# vr.prior <- function( pars, vrbounds, ...){
#   par.min = vrbounds[,1]
#   par.max = vrbounds[,2]
#   
#   1/( par.max - par.min ) * ( par.min <= pars  &  pars <= par.max )
#   
# }
