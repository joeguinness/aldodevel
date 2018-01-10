

#install.packages("/Users/guinness/Dropbox/research/aldodevel_0.1.0.tar.gz",repos = NULL, type = "source" )


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
# covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 1/2, nugget = 0)

# simulateData does full covariance calculation. beware!
y <- rnorm(n)

# generate an ordering
ord <- order_maxmin(locs)

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]

# find the ordered m nearest neighbors
m <- 30
NNarray <- findOrderedNN_kdtree(locsord,m)

# automatically group the observations
NNlist <- groupNN(NNarray)
# NNlist[1:10]  # list elements are subsets of rows of NNarray

# compute the ungrouped and grouped likelihoods
t1 <- proc.time()[3]
ll_ungrouped <- vecchiaLik(covparms,"maternIsotropic",yord,locsord,NNarray)
t2 <- proc.time()[3]
ll_ungrouped_chol <- vecchiaLik_function(covparms,"maternIsotropic",yord,locsord,NNarray)
t3 <- proc.time()[3]
ll_ungrouped - ll_ungrouped_chol
t2-t1
t3-t2
# > ll_ungrouped
# [1] -21133.47


# figure out why _function is slower
t4 <- proc.time()[3]
ll_grouped <- vecchiaLik_grouped_function(covparms,"maternIsotropic",yord,locsord,NNlist)
t5 <- proc.time()[3]
ll_grouped_chol <- vecchiaLik_grouped(covparms,yord,locsord,NNlist)
t6 <- proc.time()[3]
ll_grouped - ll_grouped_chol
t5-t4
t6-t5




