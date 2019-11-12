# set prior distributions for model parameters

# calculate the hyper-parameters for prior on alpha, using data from *other* studies (e.g. simulations) not used in the present analysis.

# uncomment and put in appropriate filename to load the data
# TracerNotUsed = readRDS("../Data/filename")

# get the r values of these unused data and only use ones beyond rmin
rvals = TracerNotUsed$Rgc[TracerNotUsed$Rgc>rmin]

# the number of tracers used to determine the prior
nextra = length(rvals)

# calculate and set the hyperparameter values for the hyperprior distribution on alpha (which is a gamma distribution with shape = (n+myc), scale=1/(nco+p))
nco = sum( log(rvals/rmin))
myc=1e-3
p=1e-3

# assign the hyperparameter values to a single object
apriorpars = c(rmin, myc, p, nextra, nco)
