
# analyze the averaged jason3 wind speed data
library("aldodevel")

# loads in object
load("datasets/jason3_windspeed_avg.RData")
attach(dframe)

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

# subtract off mean
windspeed0 <- windspeed - mean(windspeed)

# get ordering
ord <- spacefill_kdtree(xyz)
plot(y[ord[1:1000]],z[ord[1:1000]])

coordsord <- coords[ord,]
windspeed0ord <- windspeed0[ord]

# get ordering from spatial locations
NNarray <- findOrderedNN_kdtree(xyz[ord,],m=20)

funtomax <- function( logparms ){

    parms <- exp(logparms)
    parms[5] <- exp(logparms[5])/(1+exp(logparms[5]))
    scaledcoords <- coordsord
    scaledcoords[,1:3] <- coordsord[,1:3]/parms[2]
    scaledcoords[,4] <- coordsord[,4]/parms[3]

    maternparms <- c( parms[1], 1, parms[4], parms[5] )
    loglik <- OrderedCompLik(maternparms, "maternIsotropic", windspeed0ord, scaledcoords, NNarray )
    return(-loglik)
}


result <- optim(rep(0,5),funtomax,control=list(trace=5,maxit=1000))


