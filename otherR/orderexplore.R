
# explore ordering timings
gsize <- 200
nvec <- c(gsize,gsize)
n <- prod(nvec)

# generate data locations and plot them
x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# generate an ordering
system.time( ord1 <- orderMaxMinLocal(locs) )
system.time( ord2 <- order_maxmin(locs) )
system.time( ord3 <- order_maxmin2(locs) )

par(mfrow=c(1,3))
ntoplot <- 421

ord <- ord1[1:ntoplot]
plot( locs[ord,1], locs[ord,2] )
ord <- ord2[1:ntoplot]
plot( locs[ord,1], locs[ord,2] )
ord <- ord3[1:ntoplot]
plot( locs[ord,1], locs[ord,2] )







