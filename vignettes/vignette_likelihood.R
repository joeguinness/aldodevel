

#install.packages("/Users/guinness/Dropbox/research/aldodevel_0.1.0.tar.gz",
#                 repos = NULL, type = "source" )


# a short vignette demonstrating how to use the functions
library("aldodevel")

# grid size for data locations
gsize <- 100
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# covariance function and parameters
covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 1/2, sig2noise = 1)

# simulateData does full covariance calculation. beware!
y <- rnorm(n)

# generate an ordering
ord <- sample(n)

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]

# find the ordered m nearest neighbors
m <- 20
NNarray <- findOrderedNN_kdtree(locsord,m)

# automatically group the observations
NNlist <- groupNN3(NNarray)
# NNlist[1:10]  # list elements are subsets of rows of NNarray

# compute the ungrouped and grouped likelihoods
system.time(  ll_ungrouped <-
    OrderedCompLik(covparms,"maternIsotropic",yord,locsord,NNarray)      )
system.time(  ll_grouped   <-
    OrderedGroupCompLik(covparms,yord,locsord,NNlist)  )

# double check answers with slower R implementation
ll_ungrouped - orderedCompLik(covparms, covfun,yord,locsord,NNarray )
ll_grouped - orderedGroupCompLik(covparms, covfun,yord,locsord,NNlist )


