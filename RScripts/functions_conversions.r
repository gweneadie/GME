# conversions
library(astrolibR)

# wrapper for ten function, which will be needed to transform RA and Dec into decimals
wrap.ten = function(RAdec, RA=TRUE){
  if(RA){
    out = ten(dd=RAdec[1], mm=RAdec[2], ss=RAdec[3])*15
  }else{
    out = ten(dd=RAdec[1], mm=RAdec[2], ss=RAdec[3])
    if(RAdec[1]<0){ out = -out}
  }
  out
}

# function to get angles for spherical coordinates and for transforming from cylindrical to spherical
  # phi measured down from the positive z-axis
  # theta measured from the positive x-axis in x-y plane
getangle = function( datfram){
  
  # projected distance onto the x-y plane
  projR = sqrt(datfram$Xgc^2 + datfram$Ygc^2)
  
  # 3D distance from origin to object
  fullR = sqrt(datfram$Xgc^2 + datfram$Ygc^2 + datfram$Zgc^2)
  
  # phi is measured down from the positive z axis
  sinphi = projR/fullR
  cosphi = datfram$Z/fullR
  
  # theta is measured from the positive x axis in the x-y plane
  costheta = datfram$Xgc/projR
  sintheta = datfram$Ygc/projR
  
  cbind(cosphi, sinphi, costheta, sintheta)
}


# function to convert Sun-centered X, Y, Z, to Galactic-centered X, Y, Z
  # units: kpc
  # sun X value 8.27 kpc to change to a Galactocentric coordinate system
  # sun Z value 0.0205 kpc from Humphrey's & Larsen (1995)
XYZtoGalXYZ = function( XYZ, adjustSunXYZ){
  XYZ - matrix(data = adjustSunXYZ, nrow=nrow(XYZ), ncol=3, byrow = TRUE)
}


# ---------- My version of the gal_uvw code
# This is more efficient and IMO, clearer.


# Coordinate matrix
Amatrix = function( RA, DEC ){
  
  degtorad = pi/180
  cosa = cos(RA*degtorad)
  sina = sin(RA*degtorad)
  cosd = cos(DEC*degtorad)
  sind = sin(DEC*degtorad)
  
  matrix( data = c( cosa*cosd, -sina, -cosa*sind,
                    sina*cosd,  cosa, -sina*sind,
                    sind,       0   ,  cosd),
          ncol=3, byrow=TRUE)
  
}

# function to transform proper motions into UVW, accounting for solar motion
# Solar motion here is in RIGHT-handed system
# Bmatrix   T * A 
# dat       data frame with (column names)
# plx     parallax                                          arcsec
# Vlos    radial velocity                                   km/s
# PMra    proper motion in ra, corrected for declination    arcsec/yr
# PMdec   proper motion in dec                              arcsec/yr
# RA      right ascension                                   decimal degrees
# DEC     declination                                       decimal degrees
# k         the equivalent of 1 AU/yr in units of km/s


PMtoUVW = function(Bmatrix, dat., k=4.74057, countSun=FALSE, solarmotion=NULL, velocityLSR=NULL){
  
  uvw = Bmatrix %*% rbind(dat.$Vlos, k*dat.$PMra/dat.$plx, k*dat.$PMdec/dat.$plx)
  
  # make uvw into a matrix, where each row is u, v, w
  uvw = matrix(data = uvw, nrow = 3, byrow = FALSE )
 
  if(countSun){
    output = uvw + solarmotion + velocityLSR
  }else{
    output = uvw
  }
  output
  
}

# Transform U, V, W space velocities (already accounted for rotation of Galaxy)
# to Pi, Theta, Z
UVWtoPTZ = function( UVW, dat. ){
  
  projR = sqrt( dat.$Xgc^2 + dat.$Ygc^2 )
  
  Pi = UVW[1, ]*dat.$Xgc/projR + UVW[2, ]*dat.$Ygc/projR
  Theta = - UVW[1, ]*dat.$Ygc/projR + UVW[2, ]*dat.$Xgc/projR
  Z = UVW[3, ]
  
  c(Pi, Theta, Z)
  
}

# function to calculate transformation matrix for going from cylindrical Galactocentric coordinates (CGC) to spherical Galactocentric coordinates (SGC) takes in UVW RIGHT-handed coordinates, converts to RIGHT-handed SGC coordinates

Transform.Matrix = function(dat){
  
  # phi is measured down from the positive z-axis
  matrix( c(dat$sinphi, 0 , dat$cosphi,
            0, 1, 0,
            dat$cosphi, 0, -dat$sinphi), byrow=TRUE, ncol=3)
}


