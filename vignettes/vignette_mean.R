

#devtools::install_github("joeguinness/aldodevel")

# a short vignette demonstrating how to use the functions
#library("aldodevel")


# grid size for data locations
gsize <- 40
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# covariance function and parameters
covfun <- maternIsotropic
covparms <- c(variance = 1, range = 0.05, smoothness = 1/2, nugget = 0)

# simulateData does full covariance calculation. beware!
X <- as.matrix( cbind(rep(1,n), x1) )
truebeta <- c(1,1)

y <- simulateData(locs,covparms,covfun) + X %*% truebeta
#y <- rnorm(n)
fields::image.plot(array(y,nvec))

# generate an ordering
ord <- sample(n)

# define ordered locations and observations
locsord <- locs[ord,]
yord <- y[ord]
Xord <- X[ord,,drop=FALSE]

# find the ordered m nearest neighbors
m <- 25
NNarray  <- findOrderedNN_kdtree(locsord,m)

# define function to plug into optim
funtomin <- function( logparms, m ){
    parms <- exp(logparms)
    ll <- proflik(parms,"maternIsotropic",yord,Xord,locsord,NNarray[,1:(m+1)])
    return(-ll)
}

# minimize the negative profile loglikelihood
# starting values
startparms <- c(0.2,1/2,1/2)
# move towards maximum with rough approximation (small m)
result0 <- optim(par = log(startparms), fn = funtomin, control = list(trace=0), m = 10 )
# sharper approximation (larger m)
result <- optim(par = result0$par, fn = funtomin, control = list(trace=0), m = 25 )

# look at the results
subparms <- exp(result$par)
fitinfo <- proflik(subparms,"maternIsotropic",yord,Xord,locsord,NNarray,returnparms = TRUE)
fitinfo
