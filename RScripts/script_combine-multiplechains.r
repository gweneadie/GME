# script to combine chains 

source("functions_deason.r")

source("script_set-up.r")


adaptiveDate = "2019-04-29"
burninDate = "2019-04-29"
chaindate = "2019-04-29"
hiermcmcDate = "2019-04-29"

nfiles = 3

rmins = 20

savenuispars = FALSE
n.pars = 4


iterations = 1e4
nProc = 3

# for each rmin value
for( j in 1:length(rmins) ){

    for( i in 1:nfiles ){
    file.name = paste("../Results/final", filename, adaptiveDate, burninDate, iterations, "iterations", "thinning20", nProc, "chains", j, hiermcmcDate, sep="-")

      
    temp = readRDS(file.name)
    nnuis = ncol(temp[[1]]$chain) - n.pars
    
    if(savenuispars){
      temp = lapply(X = temp, FUN = function(x) x$chain[, (n.pars+1):(nnuis+n.pars)] )  
    }else{
      temp = lapply(X = temp, FUN = function(x) x$chain[, 1:n.pars] )
    }
    temp = do.call(rbind, args=temp)
 
    if( i ==1 ){ chain = temp }else{
      chain = rbind(chain, temp)
    }
      
    }
  
  rm(temp)

  if(savenuispars){
    # begin filename with
    beginname = "nuis"    
  }else{
    # save only model parameters (transform model parameters)
     chain = transform.chains(x = chain, vtpars = FALSE, npars=n.pars, nvtpars = nnuis )    
     # begin filename with
     beginname = "modelpars"
  }
  
  # save the combined chain
    saveRDS(object = chain, file = paste("../Results/", beginname, "-", nfiles, "-chainscombined", "-", filename, "-", adaptiveDate, "_", burninDate, "_", chaindate, "_", hiermcmcDate, "_", Sys.Date(), sep=""))
    
  rm(chain)

}

