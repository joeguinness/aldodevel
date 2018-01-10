#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "covfuns.h"

using namespace std;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
NumericVector vecchiaLik(NumericVector covparms, StringVector covfun_name,
                             NumericVector y,
                             NumericMatrix locs, IntegerMatrix NNarray) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    double cparms[3] = {covparms[0], covparms[1], covparms[2]};
    double nugget = covparms[0]*covparms[3];

    NumericVector ll(1);
    int n = y.length();
    int dim = locs.ncol();

    // number of neighbors + 1
    int m = NNarray.ncol();

    double d;
    double ysub[m];
    std::vector<std::vector<double> > locsub(m, std::vector<double>(dim, 0));
    std::vector<std::vector<double> > L(m, std::vector<double>(m, 0));

    std::vector<double> g (m,0);
    std::vector<double> sig (m,0);


//    vector<double> loc1(dim);
//    vector<double> loc2(dim);

    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            for(k=0;k<dim;k++){ locsub[m-1-j][k] = locs( NNarray(i-1,j)-1, k ); }
            ysub[m-1-j] = y[ NNarray(i-1,j)-1 ];
        }

        // make sure all elements are zero
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ L[k][j] = 0.0; }}

        // get 0,0 entry
        L[0][0] = pow( MaternFunction(0, cparms) + nugget, 0.5 );


        for(j=2; j<m+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<m+1; k++){ g[k-1] = 0.0; }

            // get first j-1 entries of jth row of L
            for(k=0; k<j-1; k++){

                // get distance //
                // i would prefer to do this calculation
                // inside of matern function, i.e. input two
                // locations and return covariance, but compiler
                // seems to have trouble vectorizing that
                d = 0.0;
                for(el=0;el<dim;el++){ d += pow(locsub[k][el]-locsub[j-1][el],2); }
                d = pow( d, 0.5 );

                // compute covariance, would prefer commented versin below
                sig[k] = MaternFunction(d,cparms);
                //sig[k-1] = Matern_from_locs(loc1,loc2,cparms);

                // solve lower triangular system to get L[j-1][k]
                // might be faster for work on L[j-1][k] directly
                // instead of through g
                g[k] = sig[k];
                if(k>0){
                    for(el=0; el<k; el++){
                        g[k] -= L[k][el]*g[el];
                    }
                }
                g[k] = g[k]/L[k][k];
                L[j-1][k] = g[k];
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            L[j-1][j-1] = pow( MaternFunction(0,cparms) + nugget - d, 0.5 );

        }

        // get g = L^{-1}y
        g[0] = ysub[0]/L[0][0];
        for(j=1; j<m; j++){
            g[j] = ysub[j];
            for(k=0; k<j; k++){
                g[j] -= L[j][k]*g[k];
            }
            g[j] = g[j]/L[j][j];
        }

        if(i==m){
            for(j=0; j<m; j++){
                ll(0) += -g[j]*g[j]/2 - log( L[j][j] );
            }
        } else {
            ll(0) += -g[m-1]*g[m-1]/2 - log( L[m-1][m-1] );
        }
    }
    ll(0) += -n*log(2*M_PI)/2;

    return ll;
}
