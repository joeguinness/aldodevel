

# a vignette showing how to use parallel functions
library(aldodevel)

gsize <- 40
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
locs <- simulateGrid(nvec,jittersize=0)

# covariance function and parameters
covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 1/2, sig2noise = 1)

# simulate some data
#y <- simulateData(locs,covparms,covfun)
y <- rnorm(n)
#y <- rep(0,n)

# generate an ordering
#ord <- orderMaxMinLocal(locs)
ord <- 1:n

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]


# find the ordered m nearest neighbors
m <- 10
NNarray <- findOrderedNNfast(locsord,m)
system.time( sv0 <- orderedCompLik(covparms,covfun,yord,locsord,NNarray) )
system.time( sv1 <- OrderedCompLik(covparms,yord,locsord,NNarray) )
sv0-sv1

#j <- 14
#sv2 <- t(solve(chol(covfun( locsord[rev(NNarray[j,]),],covparms ))))
round(sv1-sv2,4)

round(sv1,3)
round(sv2,3)
sv1
locsord[rev(NNarray[j,]),]
locsord[rev(NNarray[j,]),] - sv1

sv3 <- t(chol(covfun( locsord[rev(NNarray[j,]),],covparms )))
sv3





