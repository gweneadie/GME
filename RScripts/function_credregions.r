library(emdbook)

# function to calculate M(R<r) credible regions

# Find the 50%, 75%, and 95% credible regions (default) at a specific radius by setting rmin and rmax as the same value, and n.r=1

# Alternatively, get a set of credible regions for a bunch of radii by setting rmin<rmax and n.r = number of radii you want

calc.MrCredInt = function(chain, rmax, regions=c(0.5, 0.75, 0.95), modelMr=M.r, npars=4, rmin=1e-3,  n.r=1e2){
  
  # r.values in a log-sequence
  r.values = lseq(rmin, rmax, length.out=n.r)
  
  # credible region bounds
  drop.perc = (1 - regions)/2
  
  # columns = r , rows = M(r) in units of 10^12 solar masses
  Mr.values = 2.325e-3*sapply(X=r.values, FUN=modelMr, pars= chain[ ,1:npars] )
  colnames(Mr.values) = paste( "r", seq(1,length(r.values)), sep="")

  # sort the M(r) values
  sorted.Mr = apply(X=Mr.values, MARGIN=2, FUN=sort)
  
  if( n.r==1 ){
    lower.creds = matrix( sapply( X=drop.perc, FUN=function(x) ( sorted.Mr[ -( 1:(x*nrow(sorted.Mr)) ) , ] )[1] ),
                            nrow=1, ncol=length(regions))
    colnames(lower.creds) = as.character(100*regions)
    
    upper.creds = matrix(apply( X=cbind( drop.perc, regions ), MARGIN=1,
                         FUN=function(x) ( sorted.Mr[ -( 1:( sum(x) * nrow( sorted.Mr ) ) ), ] )[1]),
                         nrow=1, ncol=length(regions))
    colnames(upper.creds) = colnames(lower.creds)
    
  }else{
    lower.creds = sapply( X=drop.perc, FUN=function(x) ( sorted.Mr[ -( 1:(x*nrow(sorted.Mr)) ) , ] )[1, ] )
    colnames(lower.creds) = as.character(100*regions)

    upper.creds = apply( X=cbind( drop.perc, regions ), MARGIN=1,
                         FUN=function(x) ( sorted.Mr[ -( 1:( sum(x) * nrow( sorted.Mr ) ) ), ] )[1, ] )
    colnames(upper.creds) = colnames(lower.creds)
   
  }
  

  meanMs = colMeans(Mr.values)
  medianMs= apply(X = Mr.values, MARGIN = 2, FUN = median)
  
  # output the credible regions
  output = list(lower.creds, upper.creds, meanMs, medianMs)

  # name the lists
  names(output) = c("lower", "upper", "meanM", "medianM")

  output
  
}

# function to calcualte credible region for any parameter from the Markov chain. Assume the Markov chain is a matrix with each column being a different parameter

calc.CredInt = function( chain, cred.reg, colnum ){
  # length of chain
  n = nrow(chain)
  
  # calc. how much to drop to get credible interval
  drop.percent = (1 - cred.reg)/2
  drop.n = drop.percent*n
  
  # sort the parameter values
  sorted = sort(chain[, colnum])
  
  # get credible regions
  credreg = sorted[-((n - drop.n):n)]
  credreg = credreg[-(1:drop.n)]
  
  # return bounds of credible interval, median, mean, sd
  output = c(range(credreg), mean(chain[, colnum]), median(chain[, colnum]), sd(chain[, colnum]))
  names(output)= c("credrangelower", "credrangeupper", "mean", "median", "sd")
  output
  
}


########## calculate credible regions for isotropic DF

calc.logDFregions = function(chain, Emin=1e-5, Emax=40, regions=c(0.5, 0.75, 0.95), modelDF=logDFenergy.Deason, npars=4, relE=relE, pot=pot.Deason, n.E=1e2){
  
  # energy values
  E.values = lseq(Emin, Emax, length.out=n.E)
  
  # credible region bounds
  drop.perc = (1 - regions)/2

  # columns = r , rows = M(r) in units of 10^12 solar masses
  logDF.values = sapply(X=E.values, FUN=modelDF, pars=chain[ ,1:npars], pot=pot)
  colnames(logDF.values) = paste( "E", seq(1,n.E), sep="")
  
  # sort the log(f(E)) values
  sorted.logDF = apply(X=logDF.values, MARGIN=2, FUN=sort)
  
  if( n.E==1 ){
    lower.creds = matrix( sapply( X=drop.perc, FUN=function(x) ( sorted.logDF[ -( 1:(x*nrow(sorted.logDF)) ) , ] )[1] ),
                          nrow=1, ncol=length(regions))
    colnames(lower.creds) = as.character(100*regions)
    
    upper.creds = matrix(apply( X=cbind( drop.perc, regions ), MARGIN=1,
                                FUN=function(x) ( sorted.logDF[ -( 1:( sum(x) * nrow( sorted.logDF ) ) ), ] )[1]),
                         nrow=1, ncol=length(regions))
    colnames(upper.creds) = colnames(lower.creds)
    
  }else{
    lower.creds = sapply( X=drop.perc, FUN=function(x) ( sorted.logDF[ -( 1:(x*nrow(sorted.logDF)) ) , ] )[1, ] )
    colnames(lower.creds) = as.character(100*regions)
    
    upper.creds = apply( X=cbind( drop.perc, regions ), MARGIN=1,
                         FUN=function(x) ( sorted.logDF[ -( 1:( sum(x) * nrow( sorted.logDF ) ) ), ] )[1, ] )
    colnames(upper.creds) = colnames(lower.creds)
    
  }
  
  
  meanlogDF = colMeans(logDF.values)
  medianlogDF= apply(X = logDF.values, MARGIN = 2, FUN = median)
  
  # output the credible regions
  output = list(lower.creds, upper.creds, meanlogDF, medianlogDF)
  
  # name the lists
  names(output) = c("lower", "upper", "mean", "median")
  
  output
  
}

