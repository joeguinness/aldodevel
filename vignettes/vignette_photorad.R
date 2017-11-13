

# fit covariance parameters to the
# photosynthetically available (active) radiation data
library("aldodevel")

# load in the data
load("datasets/pardata.RData")

# overwrite image.plot to get axes right
implot <- function(mat, ... ){
    fields::image.plot( t( mat[ nrow(mat):1 , ] ) , ... )
}

# plot the data and get dimensions
implot(pardata)
nvec <- dim(pardata)
n1 <- nvec[1]
n2 <- nvec[2]

# reshape the data into vectors
observed <- !is.na(pardata)
locs <- as.matrix( expand.grid( 1:n1, 1:n2 ) )[observed,]
parvec <- pardata[ observed ]
n <- length(parvec)
n # possible to do exact max lik, but very slow

# get an ordering
ord <- spacefill_kdtree(locs)
#ord <- 1:n

# reorder data
locsord <- locs[ord,]
parord <- parvec[ord]

# get ordered nearest neighbors
m <- 10
NNarray <- findOrderedNN_kdtree( locsord, m )
NNlist <- groupNN3(NNarray)

# plot examples of the neighbors
par(mfrow=c(1,3),mar=c(1,1,3,1))
for( k in c(100,500,5000) ){
    plot( locsord[,2], -locsord[,1], pch=16, cex = 0.3, col = "gray", axes = FALSE )
    points( locsord[1:k,2], -locsord[1:k,1], pch=16, cex = 0.8, col = "black" )
    points( locsord[NNarray[k,2:(m+1)],2], -locsord[NNarray[k,2:(m+1)],1], pch = 16, cex=0.9, col = "blue")
    points( locsord[k,2], -locsord[k,1], pch = 16, cex = 0.9, col = "magenta" )
    mtext(paste("i =",k), side = 3, line = 1 )
}


# maximize a likelihood

# function to minimize
funtomin <- function( logparms ){

    parms <- exp(logparms)
    parms[4] <- exp(logparms[4])/( 1 + exp(logparms[4]) )
    loglik <- OrderedCompLik(parms, "maternIsotropic", parord - mean(parord), locsord, NNarray )
    #loglik <- OrderedGroupCompLik(parms, parord - mean(parord), locsord, NNlist )
    return(-loglik)

}

# use optim to get minimim
# control = list(trace=5) tells us to print out info
startvals <- rep(0,4)
# startvals <- result$par
t1 <- proc.time()[3]
result <- optim( startvals, funtomin, control = list(trace=5, maxit=1000) )
totaltime <- proc.time()[3] - t1

result

# get estimated parameters
logparms <- result$par
parms <- exp(logparms)
parms[4] <- exp(logparms[4])/( 1 + exp(logparms[4]) )
parms

# how long did it take
totaltime

# how long per function evaluation
totaltime/result$counts[1]









# write a routine to start with small m, then increase
# and monitor parameter values for convergence
fitmodel <- function( yord, locsord, startvals, mvec = seq(5,30,by=5) ){

    number_m <- length(mvec)
    allparms <- array( NA, c(number_m, length(startvals) ) )

    funtomin <- function( logparms ){
        parms <- exp(logparms)
        parms[4] <- exp(logparms[4])/( 1 + exp(logparms[4]) )
        loglik <- OrderedCompLik(covparms = parms, covfun_name = "maternIsotropic",
                                 y = parord - mean(parord), locs = locsord,
                                 NNarray = NNarray )
        return(-loglik)
    }

    for(j in 1:number_m){
        if(j > 1) startvals <- result$par
        NNarray <- findOrderedNN_kdtree( locsord, mvec[j] )

        t1 <- proc.time()[3]
        result <- optim( par = startvals, fn = funtomin )
        totaltime <- proc.time()[3] - t1
        print(totaltime)

        logparms <- result$par
        parms <- exp(logparms)
        parms[4] <- exp(logparms[4])/( 1 + exp(logparms[4]) )
        allparms[j,] <- parms
    }
    return(list(parms=allparms,result=result))
}


fit <- fitmodel(yord,locsord, startvals = rep(0,4) )

fit$parms


