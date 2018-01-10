#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "covfuns.h"

using namespace std;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
NumericVector vecchiaLik_function(NumericVector covparms, StringVector covfun_name,
                         NumericVector y,
                         NumericMatrix locs, IntegerMatrix NNarray) {

    // utility integers
    int i;
    int j;
    int k;
    int el;

    // covariance parameters
    int nparms = covparms.length();

    double cparms[nparms-1];                                 // number of covariance parameters
    for(k=0; k<nparms-1; k++){ cparms[k] = covparms[k]; }   // assign covariance parameters
    double nugget = covparms[0]*covparms[nparms-1];           // nugget is always last

    NumericVector ll(1);              // loglikelihood to be returned
    int n = y.length();               // length of response
    int dim = locs.ncol();            // dimension of input locations
    int m = NNarray.ncol();           // number of neighbors + 1
    double d;                         // utility double

    // subset of locations and response values
    std::vector<std::vector<double> > locsub(m, std::vector<double>(dim, 0));
    double ysub[m];

    // cholesky factor
    std::vector<std::vector<double> > L(m, std::vector<double>(m, 0));

    // vectors for covariance and cholesky row
    std::vector<double> g (m,0);
    std::vector<double> sig (m,0);

    // the user inputs a covariance function name as a string (covfun_name)
    // then based on the string, we assign a pointer to an allowable
    // c++ covariance function
    double (*p_covfun)(const std::vector<double>* loc1, const std::vector<double>* loc2, double* cparms);
    p_covfun = &MaternIsotropic;

    /*
    if( strncmp(covfun_name,"maternIsotropic") )
    {
        p_covfun = &MaternIsotropic;
    } else {
        cout << "Unrecognized Covariance Function Name \n";
        assert(0);
    }
    */

    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            for(k=0;k<dim;k++){ locsub[m-1-j][k] = locs( NNarray(i-1,j)-1, k ); }
            ysub[m-1-j] = y[ NNarray(i-1,j)-1 ];
        }

        // make sure all elements are zero
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ L[k][j] = 0.0; }}

        // get 0,0 entry
        L[0][0] = pow( (*p_covfun)(&locsub[0], &locsub[0], cparms) + nugget, 0.5 );

        for(j=2; j<m+1; j++){  // j = row of Li

            // initialize g
            for(k=0; k<m; k++){ g[k] = 0.0; }

            // get first j-1 entries of jth row of L
            for(k=0; k<j-1; k++){

                // compute covariance, would prefer commented version below it
                sig[k] = (*p_covfun)(&locsub[k], &locsub[j-1], cparms);

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
            L[j-1][j-1] = pow( (*p_covfun)(&locsub[j-1],&locsub[j-1],cparms) + nugget - d, 0.5 );

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

        // get contribution to likelihood
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
