# function with transformation functions for sampling

# for transforming parameters back to regular space
transform.func = function(x){
  c( exp(x[1:3]), (1-exp(x[4])) )
}

# parametrization for sampling
inv.transform = function(x){
  c( log(x[1:3]), log(1 - x[4]) )
}

# defualt if no parameter transformation
default.transform = function(x){
  x
}

# for transforming chains back to regular space
transform.chains = function(x, vtpars=FALSE, npars=NULL, nvtpars=NULL){
  if(vtpars){
    cbind( exp(x[, 1:3]), (1-exp(x[, 4])), x[, 5:(npars+nvtpars)] )
  }else{
    cbind( exp(x[, 1:3]), (1-exp(x[, 4])))
  }
}
