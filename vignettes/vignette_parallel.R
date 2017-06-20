devtools::load_all()


# a vignette showing how to use parallel functions
#library(aldodevel)

gsize <- 100
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
locs <- simulateGrid(nvec,jittersize=0)

# covariance function and parameters
covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, sig2noise = 1)

# simulate some data
#y <- simulateData(locs,covparms,covfun)
y <- rnorm(n)

# generate an ordering
ord <- orderMaxMinLocal(locs)
#ord <- 1:n

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]

# find the ordered m nearest neighbors
m <- 30
NNarray <- findOrderedNNfast(locsord,m)
system.time( sv0 <- orderedCompLik(covparms,covfun,yord,locsord,NNarray) )
system.time( sv1 <- OrderedCompLik(covparms,yord,locsord,NNarray) )
sv0-sv1




# try some parallel stuff




















