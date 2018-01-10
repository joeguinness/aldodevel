

spacefill_kdtree2 <- function(locs){

    # get number of locs
    n <- nrow(locs)

    # pick a random ordering
    set.seed(7)
    ord <- sample(n)
    invord <- order(ord)

    # loop over the ordering
    for(i in 2:round(n/4)){

        # pick out the nneigh nearest neighbors
        # to the ith point in the ordering
        nneigh <- round(n/i)-1
        neighbors <- FNN::get.knnx( locs,
                                    locs[ ord[i],,drop=FALSE ],
                                    k = nneigh )$nn.index

        # disperse the nneigh neighbors of the ith point in
        # the ordering to new indices in the ordering

        # these are the new ordering locations of the neighbors
        newindices <- round( (n-i)*(nneigh:1)/nneigh )

        # loop over nneigh neighbors
        for(k in 2:nneigh){

            # wher ethe kth neighbor currently is
            curpos <- invord[neighbors[k]]
            # its new position
            newpos <- newindices[k]

            # only move it if its new position is later in the ordering
            # and if its new position is later than its current position
            if( newpos > i && newpos > curpos){

                # shift the points in between down by one
                ord[curpos:(newpos-1)] <- ord[(curpos+1):newpos]
                # put the neighbor in the new position
                ord[newpos] <- neighbors[k]
                # update the inverse ordering
                invord[ord[curpos:(newpos-1)]] <- invord[ord[curpos:(newpos-1)]]-1
                invord[neighbors[k]] <- newpos
            }
        }
    }
    return(ord)
}




spacefill_kdtree3 <- function(locs,setlength=10){

    # get number of locs
    n <- nrow(locs)

    # pick a random ordering
    set.seed(7)
    ord <- sample(n)
    invord <- order(ord)

    i <- 2

    # do the first 100 individually
    for(i in 2:100){

        nneigh <- round(n/i)-1
        neighbors <- FNN::knnx.index( locs,
                                      locs[ ord[i],,drop=FALSE ],
                                      k = nneigh )
        newindices <- round( (n-i)*(nneigh:1)/nneigh )
        #print(neighbors)

        for(k in 2:nneigh){
            curpos <- invord[neighbors[k]]
            newpos <- newindices[k]
            if( newpos > i && newpos > curpos){
                ord[curpos:(newpos-1)] <- ord[(curpos+1):newpos]
                ord[newpos] <- neighbors[k]
                invord[ord[curpos:(newpos-1)]] <- invord[ord[curpos:(newpos-1)]]-1
                invord[neighbors[k]] <- newpos
            }
        }
        #print(ord)
    }


    #print(ord)
    for(j in 1:round(n/setlength/4)){
        curinds <- i:(i+setlength-1)
        nneigh <- round(n/i)-1
        #print(nneigh)
        neighbors <- FNN::knnx.index( locs,
                                      locs[ ord[curinds],,drop=FALSE ],
                                      k = nneigh )
        #cat("\n \n")

        keeprows <- c(TRUE,rep(FALSE,setlength-1))
        for(k in 2:setlength){
            if( ! ord[curinds[k]] %in% c(neighbors[keeprows,]) ){
                keeprows[k] <- TRUE
            }
        }

        neighbors <- neighbors[keeprows,,drop=FALSE]
        #print(neighbors[,1:3])
        #print(neighbors)
        #cat("\n\n\n\n")

        for(el in 1:nrow(neighbors)){

            nneigh <- round(n/i)-1
            newindices <- round( (n-i)*(nneigh:1)/nneigh )

            for(k in 2:nneigh){
                curpos <- invord[neighbors[el,k]]
                newpos <- newindices[k]
                if( newpos > i && newpos > curpos){
                    ord[curpos:(newpos-1)] <- ord[(curpos+1):newpos]
                    ord[newpos] <- neighbors[el,k]
                    invord[ord[curpos:(newpos-1)]] <- invord[ord[curpos:(newpos-1)]]-1
                    invord[neighbors[el,k]] <- newpos
                }
            }
            i <- i+1
            #print(ord)
        }
    }
    return(ord)
}

