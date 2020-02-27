massfunc = function(chain, r, credregion){
  
  # get the mass at that radius
  Ms = M.r(r, pars = chain)
  
  # calculate the credible region for th emass at that radius, return in 1e12 solar mass unite
  calc.CredInt(chain = as.matrix(Ms*2.325e9), colnum  = 1, cred.reg = credregion)/1e12

}
