

#devtools::install_github("joeguinness/aldodevel")

# a short vignette demonstrating how to use the functions
library("aldodevel")

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
covparms <- c(variance = 1, range = 0.1, smoothness = 1/2, nugget = 0)

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
m <- 15
NNarray <- findOrderedNN_kdtree(locsord,m)

# profile out mean and variance parameters
proflik <- function(yord,Xord,subparms,covfun,locsord,NNarray,returnparms = FALSE){

    n <- length(yord)
    covparms1 <- c(1,subparms)
    LinvEntries <- getLinvEntries(covparms1,"maternIsotropic",locsord,NNarray)
    B <- array(NA, dim(Xord))
    for(j in 1:ncol(X)){
        B[,j] <- LinvMultFromEntries(LinvEntries,Xord[,j],NNarray)
    }
    infomat <- crossprod(B)
    z <- LinvMultFromEntries(LinvEntries,yord,NNarray)
    beta <- solve( infomat, crossprod(B,z) )
    resids <- yord - Xord %*% beta
    z_resids <- LinvMultFromEntries(LinvEntries,resids,NNarray)
    sigmasq <- c( crossprod(z_resids)/n )

    logdet <- -2*sum(log(LinvEntries[,1])) + n*log(sigmasq)
    quadform <- n
    profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    if( !returnparms ){
        return(profloglik)
    }
    if( returnparms ){
        betacovmat <- sigmasq*solve(infomat)
        return(list(loglik = profloglik, covparms = c(sigmasq,subparms),
                    beta = beta, betacovmat = betacovmat))
    }
}

# define function to plug into optim
funtomin <- function( logparms ){

    parms <- exp(logparms)
    ll <- proflik(yord,Xord,parms,covfun,locsord,NNarray)
    return(-ll)

}

# minimize the negative profile loglikelihood
startparms <- c(0.2,1/2,1/2)
result <- optim(par = log(startparms), fn = funtomin, control = list(trace=2) )

# look at the results
subparms <- exp(result$par)
fitinfo <- proflik(yord,Xord,subparms,covfun,locsord,NNarray,returnparms = TRUE)
fitinfo
