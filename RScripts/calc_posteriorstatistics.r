# script to calculate stats from posterior distributions


library(coda)
library(xtable)

source("function_credregions.r")
source("functions_deason.r")


# load final chains

#########################

# Gaia + HSTPROMO, r>15kpc
final.chainGaiaHSTPROMO = readRDS(file = "../Results/modelpars-3-chainscombined-GCdata_noSgrGCsnoCrater_GaiaHSTPROMO_2018-10-05_2018-10-05-Deason-rmin-15-2018-10-05_2018-10-05_2018-10-05_2018-10-05_2018-10-06")

# Gaia + HSTPROMO, r>20kpc
final.chainGaia20 = readRDS(file = "../Results/modelpars-3-chainscombined-GCdata_noSgrGCsnoCrater_GaiaHSTPROMO_2018-10-05_2018-10-05-Deason-rmin-20-2018-10-08_2018-10-08_2018-10-08_2018-10-08_2018-10-09")



# label the columns
colnames(final.chainGaiaHSTPROMO) = c("Phi", "gamma", "alpha", "beta")
colnames(final.chainVasiliev) = c("Phi", "gamma", "alpha", "beta")

# rvirial estimate with Ho = 0.73e-3 (G=1 units)
r200GaiaHST = rvir(final.chainGaiaHSTPROMO, H = 0.73e-3)
r200Vas = rvir(final.chainVasiliev, H=0.73e-3)

# mass estimate
M200GaiaHST = M.r(r = r200GaiaHST, pars = final.chainGaiaHSTPROMO)
M200Vas = M.r(r = r200Vas, pars = final.chainVasiliev)

all.M200GaiaHST = M200GaiaHST*2.325e-3
all.M200Vas = M200Vas*2.325e-3

# calculate effective sample size
effectiveSize(as.mcmc(all.M200GaiaHST))
effectiveSize(as.mcmc(all.M200Vas))

# calculate summary statistics
sM200Gaia = summary(as.mcmc(all.M200GaiaHST))
sM200Vas = summary(as.mcmc(all.M200Vas))

# calculate credible regions for M200
calc.CredInt(chain = as.matrix(all.M200GaiaHST), cred.reg = c(0.95), colnum = 1)

calc.CredInt(chain = as.matrix(all.M200Vas), cred.reg = c(0.95), colnum = 1)

