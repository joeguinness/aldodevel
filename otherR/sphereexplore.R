
# analyze the averaged jason3 wind speed data
#library("aldodevel")
devtools::load_all(".")

# loads in object
load("datasets/jason3_windspeed_avg.RData")
attach(dframe)
n <- length(windspeed)

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
plot(y[ord[1:200]],z[ord[1:200]])
xyzord <- xyz[ord,]

coordsord <- coords[ord,]
windspeedord <- windspeed[ord]

# get ordering from spatial locations
NNarray <- findOrderedNN_kdtree(xyz[ord,],m=20)



startparms <- c(10,0.1,0.8,0.001)

t1 <- proc.time()

loglik1 <- vecchiaLik_function(startparms,"maternSphere",windspeedord,cbind(lon[ord],lat[ord]),NNarray)
loglik1

proc.time() - t1

t1 <- proc.time()

loglik2 <- vecchiaLik_function(startparms,"maternIsotropic",windspeedord,xyzord,NNarray)
loglik2

proc.time() - t1


