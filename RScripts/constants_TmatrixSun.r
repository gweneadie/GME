# constants for transforming from Galactocentric to heliocentric

# set the transformation matrix (up-to-date transformation matrix using J2000 epoch)
# this will give RIGHT-HANDED UVW, which is what we want here
myTmatrix = matrix(data = c(-0.0548755604, -0.8734370902 , -0.4838350155,
                            0.4941094279, -0.44482963, 0.7469822445,
                            -0.867666149, -0.1980763734, 0.4559837762),
                   nrow = 3, ncol = 3, byrow = TRUE)

# current solar motion from Schonrich, Binney, and Dehnen (2010), MNRAS 403 (4)
sunmotion = rbind(11.1, 12.24, 7.25)
# velocity of the LSR
vLSR = rbind(0,220,0)
