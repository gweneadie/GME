# Script to plot Galactocentric velocity components

# load constants and conversions needed
source("constants_TmatrixSun.r")
source("functions_conversions.r")

GaiaHST = readRDS(file = "../DataRFormat/GCdata_noSgrGCsnoCraternoFSR_Vasiliev_2019-02-27_2019-02-27")

# remove rows without proper motions
GaiaHSTcomplete = GaiaHST[ which(!is.na(GaiaHST$PMra)), ]
GaiaHSTcomplete = GaiaHSTcomplete[which(!is.na(GaiaHSTcomplete$Vlos)), ]

# get the number of data points
n.dat = nrow(GaiaHSTcomplete)

# calculate A, B-matrix and inverse-B matrix, and transformation matrix for velocities
# (for going from PMs, to Cylindrical i.e. U,V,W to Spherical i.e. vr, vtheta, vphi)
Bmatrix.list = vector("list",length=n.dat)
TMatrixVels.list = vector("list", length=n.dat)

# pre-make empty columns
GaiaHSTcomplete[, c("myU", "myV", "myW", "myPI", "myTHETA", "myZ", "vr", "vtheta", "vphi")] = NA

for(i in 1:n.dat){
  
  Amatrix.temp = Amatrix(RA = GaiaHSTcomplete$RA[i], DEC = GaiaHSTcomplete$DEC[i])
  
  Bmatrix.list[[i]] = myTmatrix %*% Amatrix.temp
  
  TMatrixVels.list[[i]] = Transform.Matrix(dat= GaiaHSTcomplete[i, ])
  
  GaiaHSTcomplete[i, c("myU", "myV", "myW")] = PMtoUVW(Bmatrix=Bmatrix.list[[i]], dat=GaiaHSTcomplete[i, ],  countSun = TRUE, solarmotion = sunmotion, velocityLSR=vLSR)
  
  UVW.temp = GaiaHSTcomplete[i, c("myU", "myV", "myW")] + vLSR
  
  GaiaHSTcomplete[i, c("myPI", "myTHETA", "myZ")] = UVWtoPTZ( UVW = t(UVW.temp), dat.=GaiaHSTcomplete[i, ])

  # transform to Pi, Theta, Z components
  PiThetaZ.temp = UVWtoPTZ( UVW=t(UVW.temp), dat.=GaiaHSTcomplete[i, ])
  
  # transform to spherical components
  vr.vtheta.vphi = ( TMatrixVels.list[[i]] %*% PiThetaZ.temp )
  
  # save as r, vr, vt
  GaiaHSTcomplete[i, c("vr", "vtheta", "vphi") ] = vr.vtheta.vphi
  
}


plot(GaiaHSTcomplete$Rgc, with(GaiaHSTcomplete, expr = sqrt(vr^2 + vtheta^2 + vphi^2)), ylab="|v| (km/s)", xlab = expression(r[galactic]), ylim=c(0, 750), xlim=c(0, 110))

identify(x = GaiaHSTcomplete$Rgc, with(GaiaHSTcomplete, expr = sqrt(vr^2 + vtheta^2 + vphi^2)), labels = GaiaHSTcomplete$FullName)
