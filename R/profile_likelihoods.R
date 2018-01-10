# various profile likelioods

# profile out both linear mean parameters and variance parameter
#' @export
proflik <- function(subparms,covfun_name = "maternIsotropic", 
                    y,X,locs,NNarray,returnparms = FALSE){
    
    n <- length(y)
    covparms1 <- c(1,subparms)
    LinvEntries <- getLinvEntries(covparms1,covfun_name,locs,NNarray)
    B <- array(NA, dim(X))
    for(j in 1:ncol(X)){
        B[,j] <- LinvMultFromEntries(LinvEntries,X[,j],NNarray)
    }
    infomat <- crossprod(B)
    z <- LinvMultFromEntries(LinvEntries,y,NNarray)
    beta <- solve( infomat, crossprod(B,z) )
    resids <- yord - X %*% beta
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
