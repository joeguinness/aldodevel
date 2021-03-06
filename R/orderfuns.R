

# This is the O(n^2) algorithm for AMMD ordering. Start with a point
# in the middle, then propose a random set of the remaining points
# (of size 'numpropose') and choose the one that has maximum minimum
# distance to the already selected points. set 'numpropose' to n
# to get the exact maximum minimum distance ordering
#' @importFrom matrixStats colMins
#' @importFrom matrixStats rowMins
orderMaxMinFast <- function( locs, numpropose ){

    n <- nrow(locs)
    d <- ncol(locs)
    remaininginds <- 1:n
    orderinds <- rep(0L,n)
    # pick a center point
    mp <- matrix(colMeans(locs),1,d)
    distmp <- fields::rdist(locs,mp)
    ordermp <- order(distmp)
    orderinds[1] = ordermp[1]
    remaininginds <- remaininginds[remaininginds!=orderinds[1]]
    for( j in 2:(n-1) ){
        randinds <- sample(remaininginds,min(numpropose,length(remaininginds)))
        distarray <-  fields::rdist(locs[orderinds[1:j-1],,drop=FALSE],locs[randinds,,drop=FALSE])
        bestind <- which(matrixStats::colMins(distarray) ==  max( matrixStats::colMins( distarray ) ))
        orderinds[j] <- randinds[bestind[1]]
        remaininginds <- remaininginds[remaininginds!=orderinds[j]]
    }
    orderinds[n] <- remaininginds
    orderinds
}


# ordered by distance to some point loc0
orderDist2Point <- function( locs, loc0 ){
    d <- ncol(locs)
    loc0 <- matrix(c(loc0),1,d)
    distvec <- rdist(locs,loc0)
    orderinds <- order(distvec)
}

# ordered by distance to the center
orderMiddleOut <- function( locs ){
    d <- ncol(locs)
    loc0 <- matrix(colMeans(locs),1,d)
    orderDist2Point(locs,loc0)
}


orderByCoordinate <- function( locs, coordinate ){
    # coordinate can be a single coordinate in {1,2,..,d}
    # in this case, the points are ordered increasingly
    # in that coordinate. if coordinate is a vector
    # of coordinates, the points are ordered increasingly
    # in the sum of the coordinates dicated by the vector
    order(rowSums(locs[,coordinate,drop=FALSE]))
}



# this is the O(n log n) algorithm for finding an AMMD ordering
# break domain into grid boxes, then order the grid boxes with
# an MMD or AMMD ordering. The loop over the grid boxes in order,
# each time selecting the one point within each grid box that has
# maximum minimum distance to all points in the grid box and
# neighboring grid boxes.
#' @export
orderMaxMinLocal <- function(locs){

    n <- dim(locs)[1]
    # number of grid boxes in each dimension.
    # could be changed but probably should remain proportional to sqrt(n)
    nside <- ceiling( 1/8*sqrt(n) )

    # like in NN search, indcube holds the indices of points within
    # each gridbox
    eps <- sqrt(.Machine$double.eps)
    indcube <- array(0,c(nside,nside,ceiling(3*n/nside^2)))  # perhaps change to 5
    lims <- matrix( c(apply(locs,2,min)-eps, apply(locs,2,max)+eps ),2,2 )
    locround <- ceiling( cbind( nside*(locs[,1]-lims[1,1])/(lims[1,2]-lims[1,1]),
                                nside*(locs[,2]-lims[2,1])/(lims[2,2]-lims[2,1]) ) )

    # put the indices in indcube
    for(j in 1:n){
        inds <- locround[j,]
        k <- which( indcube[inds[1],inds[2],] == 0 )[1]
        indcube[inds[1],inds[2],k] <- j
    }


    remainingcube <- indcube
    usedcube <- array(0,dim(indcube))

    # define grid locations and order the grid boxes
    #gridlocs <- simulateGrid(c(nside,nside),0.01)
    gridlocs <- as.matrix(expand.grid(1:nside,1:nside)) +
                            0.01*matrix(2*nside^2,nside^2,2)
    if( nside*nside > 500 ){ # if the number of grid boxes is very large
        # use the function recursively to order them
        gridorder <- orderMaxMinLocal(gridlocs)
    } else {
        # use the quadratic algorithm to order grid boxes
        gridorder <- orderMaxMinFast(gridlocs,20) }


    k <- 1
    totalnumused <- 0
    orderinds <- rep(0,n)
    while( totalnumused < n ){

        # kmod tells us which grid box we are using
        # in this iteration
        kmod <- ((k-1) %% (nside*nside) ) + 1
        gridind <- gridorder[kmod]

        # grid box in (i,j) coordinates
        # used to define and grab neighboring grid boxes
        i1 <- ( (gridind-1) %% nside ) + 1
        i2 <- ceiling( gridind/nside )

        # locations that haven't been selected yet in current grid box
        rem <- remainingcube[i1,i2,remainingcube[i1,i2,] != 0]

        # already selected locations in neighboring boxes
        # we want the next one to be as far as possible from these
        used0 <- c(usedcube[max(1,i1-1):min(nside,i1+1),max(1,i2-1):min(nside,i2+1),])
        used <- used0[used0 != 0]

        if( length(rem) > 0 ){

            if(length(used) > 0){

                # figure out which remaining point (in 'rem') has maximum minimum
                # distance to the already selected points (in 'used')
                distmat <- fields::rdist(locs[rem,,drop=FALSE],locs[used,,drop=FALSE])
                #mindist <- apply(distmat,1,min)  # too slow
                mindist <- matrixStats::rowMins(distmat)
                whichind <- which( mindist == max(mindist) )[1]
                nextind <- rem[whichind]

            } else {

                nextind <- rem[1]

            }

            # update the counter and the ordering vector
            totalnumused <- totalnumused+1
            orderinds[totalnumused] <- nextind

            # update the set of already selected points
            repind <- which(usedcube[i1,i2,] == 0)[1]
            usedcube[i1,i2,repind] <- nextind

            # make sure the selected point is no longer a remaining point
            repind <- which(remainingcube[i1,i2,] == nextind)
            remainingcube[i1,i2,repind] = 0

        }

        # move to the next grid box
        k <- k+1
    }
    return(orderinds)

}



spacefill_kdtree <- function(locs, m = 50){

    # get number of locs
    n <- nrow(locs)
    # m is number of neighbors to search over
    # get the past and future nearest neighbors
    NNall <- FNN::get.knn( locs, k = m )$nn.index
    # pick a random ordering
    ord <- sample(n)
    invord <- order(ord)
    # loop over the first n/4 locations
    # move an index to the end if it is a
    # near neighbor of a previous location
    nmoved <- 0
    for(j in 1:round(n/4)){
        nneigh <- max(1, round( m*(1 - 4*j/n) ) )
        neighbors <- NNall[ord[j],1:nneigh]
        neighbors <- neighbors[ invord[neighbors] > j ]
        nneigh <- length(neighbors)
        newindices <- round( (n-j)*(nneigh:1)/nneigh )

        if(length(neighbors)>0){
            for(k in 1:nneigh){
                curpos <- invord[neighbors[k]]
                newpos <- newindices[k]
                if( newpos > j && newpos > curpos){
                    ord[curpos:(newpos-1)] <- ord[(curpos+1):newpos]
                    ord[newpos] <- neighbors[k]
                    invord[ord[curpos:(newpos-1)]] <- invord[ord[curpos:(newpos-1)]]-1
                    invord[neighbors[k]] <- newpos
                }
            }
        }
    }
    return(ord)
}



order_maxmin <- function(locs){

    # get number of locs
    n <- nrow(locs)
    m <- round(sqrt(n))
    # m is number of neighbors to search over
    # get the past and future nearest neighbors
    NNall <- FNN::get.knn( locs, k = m )$nn.index
    # pick a random ordering
    index_in_position <- c( sample(n), rep(NA,1*n) )
    position_of_index <- order(index_in_position[1:n])
    # loop over the first n/4 locations
    # move an index to the end if it is a
    # near neighbor of a previous location
    #nmoved <- 0
    curlen <- n
    curpos <- 1
    nmoved <- 0
    for(j in 2:(2*n) ){
        nneigh <- round( min(m,n/(j-nmoved+1)) )
        neighbors <- NNall[index_in_position[j],1:nneigh]
        if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
            nmoved <- nmoved+1
            curlen <- curlen + 1
            position_of_index[ index_in_position[j] ] <- curlen
            index_in_position[curlen] <- index_in_position[j]
            index_in_position[j] <- NA
        }
        if( j - nmoved > n/4){cat("breaking \n");  break }
    }
    ord <- index_in_position[ !is.na( index_in_position ) ]
    cat(nmoved)

    return(ord)
}



order_maxmin2 <- function(locs){

    # get number of locs
    n <- nrow(locs)
    m <- 200
    # m is number of neighbors to search over
    # get the past and future nearest neighbors
    NNall <- FNN::get.knn( locs, k = m )$nn.index
    # pick a random ordering
    index_in_position <- c( sample(n), rep(NA,1*n) )
    position_of_index <- order(index_in_position[1:n])
    # loop over the first n/4 locations
    # move an index to the end if it is a
    # near neighbor of a previous location
    #nmoved <- 0
    curlen <- n
    curpos <- 1
    nmoved <- 0
    for(j in 2:(2*n) ){
        nneigh <- round( n/(j-nmoved+1) )
        if(nneigh > m){
            #print(index_in_position[j])
            jj <- index_in_position[j]
            tmp <- FNN::get.knnx( locs[(1:n)[-jj],,drop=FALSE],
                                  locs[jj,,drop=FALSE],
                                  k = nneigh)$nn.index

            neighbors <- c(tmp)
            #print(neighbors)
        } else {
            neighbors <- NNall[index_in_position[j],1:nneigh]
        }
        if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
            nmoved <- nmoved+1
            curlen <- curlen + 1
            position_of_index[ index_in_position[j] ] <- curlen
            index_in_position[curlen] <- index_in_position[j]
            index_in_position[j] <- NA
        }
        if( j - nmoved > n/4){cat("breaking \n");  break }
    }
    ord <- index_in_position[ !is.na( index_in_position ) ]
    cat(nmoved)

    return(ord)
}
