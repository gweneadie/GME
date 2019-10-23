# quick script to calculate mass

source("function_credregions.r")
source("functions_deason.r")

# load the posterior distribution parameter you'd like to investigate
final.chain = readRDS(file = "../Results/modelpars-3-chainscombined-GCdata_noSgrGCsnoCraternoFSR_Vasiliev_2019-02-27_2019-02-27-Deason-rmin-15-2019-04-26_2019-04-26_2019-04-26_2019-04-26_2019-04-26")


# get the 50%, 75%, and 95% credible regions at a specific radius by setting rmin and rmax as the same value, and n.r=1

calc.MrCredInt(chain = final.chain, rmax = 50, rmin=50, n.r=1)

# calculate the virial radius from every set of parameter values, given H_0 = 73
rvirials = rvir(final.chain, H = 0.73e-3)

# get the credible regions for the virial radius
calc.CredInt(chain = as.matrix(rvirials), colnum = 1, cred.reg = 0.95)

# calculate the posterior distribution for M200
Mvirials = M.r(r = rvirials, pars = final.chain)

# then get the credibl regions for M200
calc.CredInt(chain = as.matrix(Mvirials*2.325e9), colnum = 1, cred.reg = 0.950)

