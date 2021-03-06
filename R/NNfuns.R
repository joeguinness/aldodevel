


# functions to find nearest neighbors and do the grouping operations


# naive nearest neighbor finder

findOrderedNN <- function( locs, m ){
     # find the m+1 nearest neighbors to locs[j,] in locs[1:j,]
     # by convention, this includes locs[j,], which is distance 0
     n <- dim(locs)[1]
     NNarray <- matrix(NA,n,m+1)
     for(j in 1:n ){
         distvec <- c(fields::rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
         NNarray[j,1:min(m+1,j)] <- order(distvec)[1:min(m+1,j)]
     }
     return(NNarray)
}





# take in an array of nearest neighbors, and automatically group
# the observations into groups that share neighbors
# this is helpful to speed the computations and improve their accuracy
# #' @export
groupNN <- function(NNarray){
    n <- nrow(NNarray)
    m <- ncol(NNarray)-1

    clust <- vector("list",n)
    for(j in 1:n) clust[[j]] <- j
    for( ell in 2:(m+1) ){  # 2:(m+1)?
        sv <- which( NNarray[,1] - NNarray[,ell] < n )
        for(j in sv){
            k <- NNarray[j,ell]
            if( length(clust[[k]]) > 0){
                nunique <- length(unique(c(NNarray[c(clust[[j]],clust[[k]]),])))

                # this is the rule for deciding whether two groups
                # should be combined
                if( nunique^2 <= length(unique(c(NNarray[clust[[j]],])))^2 + length(unique(c(NNarray[clust[[k]],])))^2 ) {
                    clust[[j]] <- c(clust[[j]],clust[[k]])
                    clust[[k]] <- numeric(0)
                }
            }
        }
    }
    zeroinds <- unlist(lapply(clust,length)==0)
    clust[zeroinds] <- NULL
    NNlist <- lapply(clust,function(inds) NNarray[inds,,drop=FALSE])
    return(NNlist)
}



# take in an array of nearest neighbors, and automatically group
# the observations into groups that share neighbors
# this is helpful to speed the computations and improve their accuracy
# this is the same as groupNN except uses a cubed criterion
#' @export
groupNN3 <- function(NNarray){
    n <- nrow(NNarray)
    m <- ncol(NNarray)-1

    clust <- vector("list",n)
    for(j in 1:n) clust[[j]] <- j
    for( ell in 2:(m+1) ){  # 2:(m+1)?
        sv <- which( NNarray[,1] - NNarray[,ell] < n )
        for(j in sv){
            k <- NNarray[j,ell]
            if( length(clust[[k]]) > 0){
                nunique <- length(unique(c(NNarray[c(clust[[j]],clust[[k]]),])))

                # this is the rule for deciding whether two groups
                # should be combined
                if( nunique^3 <= length(unique(c(NNarray[clust[[j]],])))^3 + length(unique(c(NNarray[clust[[k]],])))^3 ) {
                    clust[[j]] <- c(clust[[j]],clust[[k]])
                    clust[[k]] <- numeric(0)
                }
            }
        }
    }
    zeroinds <- unlist(lapply(clust,length)==0)
    clust[zeroinds] <- NULL
    NNlist <- lapply(clust,function(inds) NNarray[inds,,drop=FALSE])
    return(NNlist)
}


# faster algorithm to find nearest neighbors. This one splits the
# observation domain into grid boxes and searches neighboring
# grid boxes to find neighbors
#' @export
findOrderedNNfast <- function( locs, m ){

    n <- dim(locs)[1]

    # number of grid boxes in each dimension
    nside <- ceiling( 1/10*sqrt(n) ) # perhaps change to 2*sqrt(n)

    eps <- sqrt(.Machine$double.eps)

    # this is a 3D array that will hold the indices of the points
    # within each grid box. This is a bit tricky because we don't
    # know ahead of time how many points are in each grid box. Not sure
    # of the best way to deal with that problem. For now I just pick
    # something big (10*n/nside^2)
    indcube <- array(0,c(nside,nside,ceiling(10*n/nside^2)))

    # rectangular bounds of observation domain
    lims <- matrix( c(apply(locs,2,min)-eps, apply(locs,2,max)+eps ),2,2 )

    # simply round the coordinates to assign them to grid boxes
    locround <- ceiling( cbind( nside*(locs[,1]-lims[1,1])/(lims[1,2]-lims[1,1]),
                                nside*(locs[,2]-lims[2,1])/(lims[2,2]-lims[2,1]) ) )

    # could be faster, but this is not the bottle neck.
    # just puts the indices into indcube, the object that
    # maps grid boxes to observation indices
    for(j in 1:n){
        inds <- locround[j,]
        k <- which( indcube[inds[1],inds[2],] == 0 )[1]
        indcube[inds[1],inds[2],k] <- j
    }

    # initialize
    NNarray <- matrix(NA,n,m+1)

    # loop over all of the grid boxes
    for(i1 in 1:nside){
        for(i2 in 1:nside){

            # number of observations in current grid box
            # maybe better to initialize indcube with NAs
            nk <- sum( indcube[i1,i2,] > 0 )

            # loop over observation indices
            for(k in seq(length.out=nk) ){
                j <- indcube[i1,i2,k]

                # s gives the number of neighboring grid boxes to search
                # s=1 means search current and adjacent 8 grid boxes
                s <- max(1,floor( 1/2*nside^2/j ))
                subind = c()
                nlocal <- 0

                # increase s until we have found "enough" PREVIOUS points
                # in neighboring grid boxes. Must be PREVIOUS in ordering
                while( nlocal < min(j,m*2) ){  # could be m*3?
                    l1 <- max(1,i1-s):min(nside,i1+s)
                    l2 <- max(1,i2-s):min(nside,i2+s)
                    subcube <- indcube[l1,l2,]
                    subind <- subcube[ subcube > 0 & subcube <= j ] # PREVIOUS!
                    nlocal <- length(subind)
                    s <- max(s*2,s+1)
                }
                # pick the closest m of them
                subdistmat <- rdist(locs[j,,drop=FALSE],locs[subind,,drop=FALSE])
                orderind <- order(subdistmat)
                NNarray[j,1:min(j,m+1)] <- as.integer(subind[orderind[1:min(j,m+1)]])
            }
        }
    }
    return(NNarray)
}


# this algorithm is a modification of the kdtree
# algorithm in the FNN package, adapted to the setting
# where the nearest neighbors must come from previous
# in the ordering
#' @export
#' @importFrom FNN get.knnx
findOrderedNN_kdtree <- function(locs,m,mult=2,printsearch=FALSE){

    # number of locations
    n <- nrow(locs)

    # to store the nearest neighbor indices
    NNarray <- matrix(NA,n,m+1)
    # to the first mult*m+1 by brutce force
    NNarray[1:(mult*m+1),] <- findOrderedNN(locs[1:(mult*m+1),],m)

    query_inds <- (mult*m+2):n
    data_inds <- 1:n

    msearch <- m

    while( length(query_inds) > 0 ){

        msearch <- min( max(query_inds), 2*msearch )
        data_inds <- 1:max(query_inds)
        NN <- get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
        less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
        sum_less_than_k <- apply(less_than_k,1,sum)
        ind_less_than_k <- which(sum_less_than_k >= m+1)
        NN_less_than_k <- NN[ind_less_than_k,]

        NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))

        NNarray[ query_inds[ind_less_than_k], ] <- NN_m

        query_inds <- query_inds[-ind_less_than_k]
        if( printsearch ) print(length(query_inds))

    }

    return(NNarray)
}



# this is a minor modification to FNN::get.knnx
# added PACKAGE = "FNN" to the call to .C
# #' @useDynLib FNN
# get.knnx_aldo <- function (data, query, k = 10,
#             algorithm = c("kd_tree", "cover_tree", "CR", "brute"))
# {
#     algorithm <- match.arg(algorithm)
#     if (!is.matrix(data))
#         data <- as.matrix(data)
#     if (!is.numeric(data))
#         stop("Data non-numeric")
#     if (any(is.na(data)))
#         stop("Data include NAs")
#     if (storage.mode(data) == "integer")
#         storage.mode(data) <- "double"
#     if (!is.matrix(query))
#         query <- as.matrix(query)
#     if (!is.numeric(query))
#         stop("Data non-numeric")
#     if (any(is.na(query)))
#         stop("Data include NAs")
#     if (storage.mode(query) == "integer")
#         storage.mode(query) <- "double"
#     n <- nrow(data)
#     m <- nrow(query)
#     d <- ncol(data)
#     p <- ncol(query)
#     if (d != p)
#         stop("Number of columns must be same!.")
#     if (k > n)
#         warning("k should be less than sample size!")
#     Cname <- switch(algorithm, cover_tree = "get_KNNX_cover",
#                     kd_tree = "get_KNNX_kd", CR = "get_KNNX_CR",
#                     brute = "get_KNNX_brute")
#     knnres <- .C(Cname, t(data), t(query), as.integer(k), d,
#                  n, m, nn.index = integer(m * k), nn.dist = double(m *k),
#                  DUP = FALSE, PACKAGE = "FNN")
#     nn.index <- matrix(knnres$nn.index, byrow = T, nrow = m, ncol = k)
#     nn.dist <- matrix(knnres$nn.dist, byrow = T, nrow = m, ncol = k)
#     if (k > n) {
#         nn.index[, (n + 1):k] <- NA
#         nn.dist[, (n + 1):k] <- NA
#     }
#     return(list(nn.index = nn.index, nn.dist = nn.dist))
# }




