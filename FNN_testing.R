system.time( unordNNarray <- FNN::get.knnx( xyzord, xyzord[floor(3*n/4+1):n,], k = 2*m )$nn.index )
sv <- sapply( 1:nrow(unordNNarray), function(j) sum( unordNNarray[j,] < j+floor(3*n/4+1) ) < m )
sum(sv)


sv <- matrix(NA,n,m)
times <- rep(NA,10)
for(k in 1:10){
    times[k] <- system.time({for(j in (m+1):(m+k*1000)){
        sv[j,] <- FNN::get.knnx( latlonord[1:(j-1),], latlonord[j,,drop=FALSE], k = m, algorithm = "kd_tree" )$nn.index
    }
    })[[3]]}

geteach <- function( j ){
    out <- FNN::get.knnx( latlonord[1:(j-1),,drop=FALSE], latlonord[j,,drop=FALSE], k = m )$nn.index
}

sv <- lapply( (m+1):n, geteach )
