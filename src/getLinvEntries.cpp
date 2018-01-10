#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "covfuns.h"

using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix getLinvEntries(NumericVector covparms, StringVector covfun_name,
                             NumericMatrix locs, IntegerMatrix NNarray) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    int B;
    double cparms[3] = {covparms[0], covparms[1], covparms[2] };
    double nugget = covparms[0]*covparms[3];

    int n = locs.nrow();
    int dim = locs.ncol();

    // number of neighbors + 1
    int m = NNarray.ncol();

    double d;
    std::vector<std::vector<double> > locsub(m, std::vector<double>(dim, 0));

    double Li[m][m];

    double g[m];
    double sig[m];

    NumericMatrix LinvEntries(n,m);


    // should really do i=m separately. Right now, I need to make
    // sure that the (m+1)th row of NNarray is
    // m+1,m,m-1,...,1
    for(k=0;k<dim;k++){ locsub[0][k] = locs(0,k); }
    LinvEntries(0,0) = pow( MaternIsotropic(&locsub[0], &locsub[0], cparms) + nugget, -0.5 );

    for(i=2; i<n+1; i++){

        // because i could be less than m
        // get maximum value
        B = min(i,m);

        // first, fill in ysub and locsub in reverse order
        for(j=B-1; j>=0; j--){
            for(k=0;k<dim;k++){
                locsub[B-1-j][k] = locs( NNarray(i-1,j)-1, k );
            }
        }
        // initialize Li
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

        // first row has just a single nonzero entry
        Li[1-1][1-1] = pow( MaternIsotropic(&locsub[0], &locsub[0], cparms) + nugget, -0.5 );

        for(j=2; j<B+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<m+1; k++){
                g[k-1] = 0.0;
            }

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){

                // get the distance
                //d = 0.0;
                //for(el=0;el<dim;el++){ d += pow(locsub[k-1][el]-locsub[j-1][el],2); }
                //d = pow( d, 0.5 );

                // get the covariance
                //sig[k-1] = Matern_from_dist(d,cparms);
                sig[k-1] = MaternIsotropic(&locsub[k-1], &locsub[j-1], cparms);


                // solve for g[k-1]
                g[k-1] = 0.0;
                for(el=1; el<k+1; el++){
                    g[k-1] += Li[k-1][el-1]*sig[el-1];
                }

            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            Li[j-1][j-1] = pow( MaternIsotropic(&locsub[j-1],&locsub[j-1],cparms) + nugget - d, -0.5 );

            // now get first j-1 entries jth row of Linverse
            for(k=1; k<j; k++){
                for(el=1; el<j+1; el++){
                    Li[j-1][k-1] += g[el-1]*Li[el-1][k-1];
                }
                Li[j-1][k-1] = -Li[j-1][k-1]*Li[j-1][j-1];
            }

        }

        d = 0.0;
        // fill in ith row of LinvEntries
        for(k=1;k<B+1;k++){
            LinvEntries(i-1,k-1) = Li[B-1][B-k];
        }
    }

    return LinvEntries;
}


// [[Rcpp::export]]
NumericVector LinvMultFromEntries(NumericMatrix LinvEntries, NumericVector z,
                                  IntegerMatrix NNarray) {

    // here, z is the result of Linv * z
    // return x = Linv^T z
    //const double PI = 3.141592653589793238463;
    int i;
    int j;

    int n = z.length();
    NumericVector x(n);
    for(j=0;j<n;j++){ x[j] = 0.0; }

    // number of neighbors + 1
    int m = NNarray.ncol();

    // is it this simple (first m rows)
    for(i=1; i<m; i++){
        for(j=1; j<i+1; j++){
            x( i - 1 ) += z( NNarray(i-1,j-1) - 1 )*LinvEntries(i-1,j-1);
        }
    }

    // all rows after m
    for(i=m; i<n+1; i++){
        for(j=1; j<m+1; j++){
            x( i-1 ) += z( NNarray(i-1,j-1) - 1 )*LinvEntries(i-1,j-1);
        }
    }

    return x;
}



// [[Rcpp::export]]
NumericVector LMultFromEntries(NumericMatrix LinvEntries, NumericVector z,
                                  IntegerMatrix NNarray) {

    // return x = L z
    // by solving (L^{-1})x = z
    int i;
    int j;
    int B;

    int n = z.length();
    NumericVector x(n);

    // number of neighbors + 1
    int m = NNarray.ncol();

    // get entry 0
    x(0) = z(0)/LinvEntries(0,0);

    // get entries 1 through n-1
    for(i=1; i<n; i++){
        B = min(i,m);
        x(i) = z(i);
        for(j=1; j<B; j++){
            x(i) -= LinvEntries(i,j)*x( NNarray(i,j) - 1 );
        }
        x(i) = x(i)/LinvEntries(i,0);
    }

    return x;
}



