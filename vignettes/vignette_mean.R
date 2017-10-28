

# a short vignette demonstrating how to use the functions
#library(aldodevel)
rm(list=ls())
# load all package functions
devtools::load_all(".")


# grid size for data locations
gsize <- 230
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
locs <- simulateGrid(nvec,jittersize=0)
#plot(locs[,1],locs[,2])

# covariance function and parameters
covfun <- maternIsotropic
covparms <- c(variance = 4, range = 0.1, smoothness = 1/2, sig2noise = 1)

# simulate some data and plot a map of it
# uses Cholesky method so don't try this for large n
X <- cbind( rep(1,n), locs[,1], locs[,2] )
beta <- c(1,6,0)
# simulateData does full covariance calculation. beware!
#y <- X %*% beta + simulateData(locs,covparms,covfun)
y <- rnorm(n)
#image( matrix(y,nvec) )

# generate an ordering and plot the first n/8
ord <- orderMaxMinLocal(locs)
#ord <- 1:n
n0 <- round(n/8)
#plot( locs[ord[1:n0],1],locs[ord[1:n0],2] )

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]
Xord <- X[ord,]

# find the ordered m nearest neighbors
m <- 15
NNarray <- findOrderedNN_kdtree(locsord,m)

# automatically group the observations
NNlist <- groupNN3(NNarray)
#NNlist[1:10]  # list elements are subsets of rows of NNarray


# compute ungrouped ("apply" and "Rcpp" implementations) and grouped ordered composite logliks
# system.time(  ll0 <- mvnMargLik(covparms,covfun,y,locs) # only do this if n is small
# "apply" implementation (no grouping)
system.time(  ll1 <- orderedCompLik(covparms,covfun,yord,locsord,NNarray)      )
# "apply" implementation (grouping)
system.time(  ll2 <- orderedGroupCompLik(covparms,covfun,yord,locsord,NNlist)  )
# Rcpp implementation (no grouping)
system.time(  ll3 <- OrderedCompLik(covparms,yord,locsord,NNarray)             )
ll1 - ll3 # should be zero
# Rcpp implementation (grouping)
system.time(  ll4 <- OrderedGroupCompLik(covparms,yord,locsord,NNlist)  )
ll2-ll4




# stuff below doesn't work quite yet



# an attempt to write a wrapper function to do all of this stuff:
# ordering, finding neighbors, maximizing parameters
# only covariance function implemented is isotropic matern
# interesting thing: block independent likelihood is
# good at finding parameter estimates that nearly maximize
# vecchia's likelihood approximation
# This function uses block independent likelihood to quickly get starting values
# for an optimization of Vecchia's approximation
system.time( result <- fitmodel(y,X,locs,maternIsotropic,numneighbors=30,fixedparameters=c(1,NA,NA,1))  )
#system.time( result2 <- fitmodelpp(y,X,locs,maternIsotropic,numneighbors=30,fixedparameters=c(1,NA,NA,1))  )


# see how long it takes to get estimates with exact likelihood
fixedparameters <- c(NA,NA,NA,1)
notfixedinds <- which(is.na(fixedparameters))  # indices of parms to estimate
linkfun <- list( function(x) log(x), function(x) log(x), function(x) log(x), function(x) log(x)/log(1-x) )
invlinkfun <- list(function(x) exp(x),function(x) exp(x),function(x) exp(x), function(x) exp(x)/(1+exp(x)))
startvals <- rep(0,4)

f2 <- function(x){
    for(j in 1:length(notfixedinds)){
        covparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](x[j])
    }
    ll <- mvnMargLik(covparms,covfun,y,locs)
    return(-ll)
}
system.time({
    result2 <- optim(startvals[notfixedinds],f2,method="Nelder-Mead",control=list(maxit=500,trace=1))
})
outparms <- fixedparameters
# transform back to original parameter domain with inverse link
for(j in 1:length(notfixedinds)){
    outparms[notfixedinds[j]] <- invlinkfun[[notfixedinds[j]]](result2$par[j])
}
result2$par <- outparms


# compare the estimates
# approximate methods
result$par
# exact likelihood
result2$par











