// calculate matern covariances

#include <Rcpp.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/bessel.hpp>

// [[Rcpp::export]]

Rcpp::NumericMatrix MaternFun( Rcpp::NumericMatrix& distmat, Rcpp::NumericVector covparms ){

    int d1 = distmat.nrow();
    int d2 = distmat.ncol();
    int j1;
    int j2;
    Rcpp::NumericMatrix covmat(d1,d2);
    double scaledist;

    double normcon = covparms(0)/(pow(2.0,covparms(2)-1)*tgamma(covparms(2)));

    for (j1 = 0; j1 < d1; j1++){
        for (j2 = 0; j2 < d2; j2++){
            if ( distmat(j1,j2) == 0 ){
                covmat(j1,j2) = covparms(0);
            } else {
                scaledist = distmat(j1,j2)/covparms(1);
                covmat(j1,j2) = normcon*pow( scaledist, covparms(2) )*
                   boost::math::cyl_bessel_k(covparms(2),scaledist);
            }
        }
    }
    return covmat;
}
