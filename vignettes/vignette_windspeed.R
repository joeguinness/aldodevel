
# analyze the averaged jason3 wind speed data
#library("aldodevel")
devtools::load_all(".")


# loads in object
load("datasets/jason3_windspeed_avg.RData")
attach(dframe)
n <- length(windspeed)

# quick plot
fields::quilt.plot(lon,lat,windspeed)

# higher resolution to see satellite path
fields::quilt.plot(lon,lat,windspeed,nx=400,ny=200)

# get xyz coordinates
lonrad <- lon*2*pi/360
latrad <- (lat+90)*2*pi/360
x <- sin(latrad)*cos(lonrad)
y <- sin(latrad)*sin(lonrad)
z <- cos(latrad)
xyz <- cbind(x,y,z)
coords <- cbind(xyz,time)

# get ordering
ord <- order_maxmin(xyz)
plot(y[ord[1:1000]],z[ord[1:1000]])
xyzord <- xyz[ord,]
coordsord <- coords[ord,]
windspeedord <- windspeed[ord]

# get ordering from spatial locations
NNarray <- findOrderedNN_kdtree(xyzord,m=20)
X <- as.matrix( rep(1,length(windspeed)) )
Xord <- X[ord,,drop=FALSE] # don't actually need this since all ones


## first analysis: ignore temporal dimension

# using maternSphere
funtomax1 <- function( logparms ){
    parms <- exp(logparms)
    parms[5] <- logparms[5]
    loglik <- vecchiaLik_function(parms[1:4],"maternSphere",windspeedord-parms[5],
                                  cbind(lon[ord],lat[ord]),NNarray)
    return(-loglik)
}

startparms <- c(10,0.1,0.8,0.001,12)
fit1 <- optim(c(log(startparms[1:4]),startparms[5]),funtomax1,control=list(trace=5,maxit=500))


## second analysis: space-time covariance

# using maternSphereTime
funtomax2 <- function( logparms ){
    parms <- exp(logparms)
    parms[6] <- logparms[6]
    loglik <- vecchiaLik_function(parms[1:5],"maternSphereTime",windspeedord-parms[6],
                                  cbind(lon[ord],lat[ord],time[ord]),NNarray)
    return(-loglik)
}
startparms <- c(10,0.1,6e4,0.8,0.001,12)
fit2 <- optim(c(log(startparms[1:5]),startparms[6]),funtomax2,control=list(trace=5,maxit=500))








parms <- exp(fit1$par)
funtomax1( log(startparms) )
funtomax2( log(startparms) )



funtomax <- function( logparms ){

    # logparms are logarithm of parms
    parms <- exp(logparms)
    # parms[1] = spatial range
    # parms[2] = temporal range
    # parms[3] = smoothness
    # parms[4] = nugget
    scaledcoords <- coordsord
    scaledcoords[,1:3] <- coordsord[,1:3]/parms[1]
    scaledcoords[,4] <- coordsord[,4]/parms[2]

    maternparms <- c( 1, parms[3], parms[4] )
    loglik <- proflik(windspeedord, Xord, maternparms, "maternIsotropic", scaledcoords, NNarray )
    #loglik <- OrderedCompLik(maternparms, "maternIsotropic", windspeed0ord, scaledcoords, NNarray )
    return(-loglik)
}

startparms <- c(0.1,6e4,0.8,0.001)
result <- optim(log(startparms),funtomax,control=list(trace=5,maxit=100))

parms <- exp(result$par)
scaledcoords <- coordsord
scaledcoords[,1:3] <- coordsord[,1:3]/parms[1]
scaledcoords[,4] <- coordsord[,4]/parms[2]
maternparms <- c(1,parms[3],parms[4])
t1 <- proc.time()
fitinfo <- proflik(windspeedord,Xord,maternparms,"maternIsotropic",scaledcoords,NNarray,returnparms=TRUE)
proc.time()-t1
fitinfo


#               #
#  Predictions  #
#               #

# get prediction locations
mediantime <- median(time) + 1000
latgrid <- seq( min(lat), max(lat), length.out = 60 )
longrid <- seq( 0, 360, length.out = 121)[1:120] # so no locations repeated
lonradgrid <- longrid*2*pi/360
latradgrid <- (latgrid+90)*2*pi/360

lonlatgrid <- as.matrix( expand.grid(lonradgrid,latradgrid) )
xyzgrid <- cbind( sin(lonlatgrid[,2])*cos(lonlatgrid[,1]),
                  sin(lonlatgrid[,2])*sin(lonlatgrid[,1]),
                  cos(lonlatgrid[,2])   )
n_pred <- nrow(coords_pred)
coords_pred <- cbind( xyzgrid, rep(mediantime, n_pred) )

# reorder prediction locations
ord_pred <- spacefill_kdtree( xyzgrid )
coordsord_pred <- coords_pred[ord_pred,]

# put all coordinates together
coordsord_all <- rbind( coordsord, coordsord_pred )

# get nearest neighbor array (in space only)
NNarray_all <- findOrderedNN_kdtree(coordsord_all[,1:3],m=60)

# scale them by fitted parameters
scaledcoords_all <- coordsord_all
scaledcoords_all[,1:3] <- coordsord_all[,1:3]/parms[1]
scaledcoords_all[,4] <- coordsord_all[,4]/parms[2]

windspeedord_withzeros <- c( windspeedord - c(fitinfo$beta), rep(0,n_pred) )

# get entries of Linv for obs locations and pred locations
LinvEntries_all <- getLinvEntries(fitinfo$covparms,"maternIsotropic",scaledcoords_all,NNarray_all)

# compute conditional expectation
inds1 <- (1:n)
inds2 <- (n+1):(n+n_pred)
v1 <- LinvMultFromEntries(LinvEntries_all, windspeedord_withzeros, NNarray_all )
v1[inds1] <- 0
v2 <- -LMultFromEntries(LinvEntries_all,v1,NNarray_all)

condexp <- v2[inds2] + c(fitinfo$beta)
condexp[ord_pred] <- v2[inds2]

condexp_array <- array( condexp, c(length(longrid),length(latgrid)) )

fields::image.plot(condexp_array)




