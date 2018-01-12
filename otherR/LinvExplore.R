
# testing for Linv functions
devtools::load_all(".")

# grid size for data locations
gsize <- 30
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# covariance function and parameters
# covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 1/2, nugget = 0)

# generate an ordering
ord <- order_maxmin(locs)

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]

# find the ordered m nearest neighbors
m <- 30
NNarray <- findOrderedNN_kdtree(locsord,m)


#
system.time( Linv1 <- getLinvEntries(covparms,"maternIsotropic", locsord,NNarray) )
system.time( Linv2 <- vecchiaLinverse(covparms,"maternIsotropic",locsord,NNarray) )
