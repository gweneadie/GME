# file that sets prior distributions, data to be used, initial parameter values, proposal distributions, etc.

# This file is used by many other files and should not be changed mid-way through an analysis

# libraries and other scripts needed
library(MASS)
library(coda)

# source r files/scripts with constants, functions, priors, GME code, etc
source("constants_TmatrixSun.r")
source("function_adaptive-burnin.r")
source("function_adjust-propsd.r")
source("functions_conversions.r")
source("functions_deason.r")
source("function_energy.r")
source("functions_observe-priors.r")
source("functions_parameter-transformations.r")
source("functions_priors.r")
source("GME-hierarchical.r")

# minimum distance tracer objects must have from the center of the Galaxy to be included in the analysis
rmin = 10
rminname = "10"

# name of the data file that will be used (this is used later in naming conventions)
file.name = "GCdata_noSgrGCsnoCraternoFSR_GaiaHSTPROMO_2019-02-27_2019-02-27"
# file.name = "GCdata_noSgrGCsnoCraternoFSR_Vasiliev_2019-02-27_2019-02-27"


# load data that is already in RDS format
tracerdata = readRDS( paste("../DataRFormat/", file.name, sep="") )

tracerdata = readRDS( paste("../../Data-tracers/TracerData_Rformat/", file.name, sep=""))

# Define the model you want to use and provide the names of the functions that will get called
model.name = "Deason"  
modelDF = NULL
logDF = logDF.Deason
modelpot = pot.Deason
modelMr = M.r
# number of model parameters
npars = 4

# bounds for single uniform prior function
phiobounds = c(1,200)    # phi_o
bbounds = c(-1,1)      # beta

# gamma parameter mean and variance
gammasd = 0.06
gammamean = 0.5

## (alpha will be set later using unused data)

# bounds for proper motions, line of sight velocities, distances
myPMbounds = c(-0.1, 0.1)   # arcsec/yr
myVlosbounds = c(-1e3, 1e3) # km/s
myRgcbounds = c(0, 3e2)     # kpc

# initial values for parameters: Phi_o, gamma, alpha, beta
# This is currently set up to run 3 independent Markov Chains---- that's why there are 3 sets of parameter values, the first set being Phi_o=150, gamma = 0.4, alpha = 3.2, and beta = 0.5.
start.list = list( c(100, 0.4, 3.2, 0.5), c(50, 0.35, 3.1, 0), c(80, 0.42, 3.6, 0.4))

# these were the values used for Gaia and HSTPROMO

# transform these values to space used for sampling (defined by inv.transform function which is set in )
par.init <- lapply( start.list, FUN=inv.transform )


# remove the ones that don't have a line-of-sight velocity
mydata = tracerdata[!is.na(tracerdata$Vlos), ]

# take only GCs beyond rmin value
mydata = mydata[mydata[,"Rgc"]>rmin, ]
# remove rownames because they're annoying later
rownames(mydata)=NULL

# how many are missing proper motions?
no.pm = sum(is.na(mydata$PMra))
# which ones are missing proper motions?
miss.pm = which(is.na(mydata$PMra))

# yes/no vector determining if any proper motions are missing at all
yesnopm = any(is.na(mydata$PMra))
# 

# create a copy of the data that will be re-worked to provide the initial values of the nuisance parameters
initdata = mydata

# initial guesses for pm's (take fraction of measurements to avoid making things unbound, which the DF doesn't allow)
initdata[miss.pm, c("PMra", "PMdec")] = mean( na.omit(mydata$PMra) )/10

# all initial parameters (both nuisance and model) need to be in a list
init.list = lapply(X = par.init, FUN = function(x, idata) list( x, idata ), idata=initdata )

# parameters for prior on alpha, using unused data
TracerNotUsed = tracerdata[ is.na(tracerdata$Vlos) & tracerdata$Rgc>rmin, ]
rvals = TracerNotUsed$Rgc[TracerNotUsed$Rgc>rmin]
nextra = length(rvals)
nco = sum( log(rvals/rmin))
rseq = seq(0,4, length.out=1e3)
myc=1e-3
p=1e-3
mu=(myc+nextra)/(nco+p)
apriorpars = c(rmin, myc, p, nextra, nco)

apriorpars = c(10, 0.001, 0.001, 3.0, 4.5)

# part of identifying filename that will be used in outputs
filename = paste(file.name, model.name, "rmin", rminname, sep="-")

