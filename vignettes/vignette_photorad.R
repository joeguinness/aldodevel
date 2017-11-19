
#devtools::install_github("joeguinness/aldodevel")


# fit covariance and mean parameters to the
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
X <- cbind( rep(1,n), locs )



# get an ordering
ord <- spacefill_kdtree(locs)
#ord <- 1:n

# reorder data
locsord <- locs[ord,]
parord <- parvec[ord]
Xord <- X[ord,,drop=FALSE]

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




# maximize a profile likelihood
funtomin <- function( logparms ){
    parms <- exp(logparms)
    pl <- proflik(parord,Xord,parms,covfun,locsord,NNarray)
    return(-pl)
}

# use optim to get minimim
# control = list(trace=5) tells us to print out info
startvals <- rep(0,3)
# startvals <- result$par
t1 <- proc.time()[3]
result <- optim( par = startvals, fn = funtomin, control = list(trace=5) )
totaltime <- proc.time()[3] - t1


# another call to proflike with returnparms = TRUE
# in order to get betahat and sigmasq
subparms <- exp(result$par)
fitinfo <- proflik(parord,Xord,subparms,covfun,locsord,NNarray,returnparms = TRUE)
fitinfo

# how long did it take
totaltime

# how long per function evaluation
totaltime/result$counts[1]


# function to do all of this stuff
t1 <- proc.time()
fitinfo <- fitmodel(parvec,locs,X)
t2 <- proc.time()
t2-t1
fitinfo






