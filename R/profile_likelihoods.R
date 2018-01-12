# various profile likelioods

# profile out both linear mean parameters and variance parameter
#' @export
proflik_mean_variance <- function(subparms,covfun_name = "maternIsotropic",
                    y,X,locs,NNarray,returnparms = FALSE){

    n <- length(y)
    covparms1 <- c(1,subparms)
    LinvEntries <- vecchiaLinverse(covparms1,covfun_name,locs,NNarray)
    B <- array(NA, dim(X))
    for(j in 1:ncol(X)){
        B[,j] <- LinvMultFromEntries(LinvEntries,X[,j],NNarray)
    }
    infomat <- crossprod(B)
    z <- LinvMultFromEntries(LinvEntries,y,NNarray)
    beta <- solve( infomat, crossprod(B,z) )
    resids <- y - X %*% beta
    z_resids <- LinvMultFromEntries(LinvEntries,resids,NNarray)
    sigmasq <- c( crossprod(z_resids)/n )

    logdet <- -2*sum(log(LinvEntries[,1])) + n*log(sigmasq)
    quadform <- n
    profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    if( !returnparms ){
        return(profloglik)
    }
    if( returnparms ){
        betacovmat <- sigmasq*solve(infomat)
        return(list(loglik = profloglik, covparms = c(sigmasq,subparms),
                    beta = beta, betacovmat = betacovmat))
    }
}


# profile out variance for mean-zero model
#' @export
proflik_variance <- function(subparms,covfun_name = "maternIsotropic",
                                  y,locs,NNarray,returnparms = FALSE){

    n <- length(y)
    covparms1 <- c(1,subparms)
    LinvEntries <- vecchiaLinverse(covparms1,covfun_name,locs,NNarray)
    z <- LinvMultFromEntries(LinvEntries,y,NNarray)
    sigmasq <- c( crossprod(z)/n )

    logdet <- -2*sum(log(LinvEntries[,1])) + n*log(sigmasq)
    quadform <- n
    profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    if( !returnparms ){
        return(profloglik)
    }
    if( returnparms ){
        return(list(loglik = profloglik, covparms = c(sigmasq,subparms) ) )
    }
}



# if for some reason you want to profile out the mean only
#' @export
proflik_mean <- function(parms,covfun_name = "maternIsotropic",
                                  y,X,locs,NNarray,returnparms = FALSE){

    n <- length(y)
    LinvEntries <- vecchiaLinverse(parms,covfun_name,locs,NNarray)
    B <- array(NA, dim(X))
    for(j in 1:ncol(X)){
        B[,j] <- LinvMultFromEntries(LinvEntries,X[,j],NNarray)
    }
    infomat <- crossprod(B)
    z <- LinvMultFromEntries(LinvEntries,y,NNarray)
    beta <- solve( infomat, crossprod(B,z) )
    resids <- y - X %*% beta
    z_resids <- LinvMultFromEntries(LinvEntries,resids,NNarray)
    #sigmasq <- c( crossprod(z_resids)/n )

    logdet <- -2*sum(log(LinvEntries[,1]))
    quadform <- sum( z_resids^2 )
    profloglik <- -1/2*( n*log(2*pi) + logdet + quadform )
    if( !returnparms ){
        return(profloglik)
    }
    if( returnparms ){
        betacovmat <- sigmasq*solve(infomat)
        return(list(loglik = profloglik, covparms = parms,
                    beta = beta, betacovmat = betacovmat))
    }
}
