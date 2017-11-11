#include <Rcpp.h>
#include <math.h>
#include <iostream>


using namespace std;
using namespace Rcpp;


double covfun3(double d, double *cparms){

    // has special cases for 1/2 and 3/2
    if( d == 0.0 ){
        d = cparms[0];
    } else {
        if( cparms[2] == 0.5 ){
            d = cparms[0]*exp(-d/cparms[1])*cparms[3];
        } else if( cparms[2] == 1.5 ){
            d = cparms[0]*(1+d/cparms[1])*exp(-d/cparms[1])*cparms[3];
        } else {
            double normcon = cparms[0]/(pow(2.0,cparms[2]-1)*Rf_gammafn(cparms[2]));
            d = normcon*pow( d/cparms[1], cparms[2] )*
                Rf_bessel_k(d/cparms[1],cparms[2],1.0);
        }
    }
    return d;
}

// this does the transformation L^{-1}y
// [[Rcpp::export]]
NumericVector whiten3d(NumericVector covparms, NumericVector y,
                               NumericMatrix locs, IntegerMatrix NNarray) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    double cparms[4] = {covparms[0], covparms[1], covparms[2], covparms[3]};

    int n = y.length();
    NumericVector z(n);

    // number of neighbors + 1
    int m = NNarray.ncol();

    double d;
    double ysub[m];
    double locsub[m][3];

    double Li[m][m];

    double g[m];
    double sig[m];



    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            locsub[m-1-j][0] = locs( NNarray(i-1,j)-1, 0 );
            locsub[m-1-j][1] = locs( NNarray(i-1,j)-1, 1 );
            locsub[m-1-j][2] = locs( NNarray(i-1,j)-1, 2 );
            ysub[m-1-j] = y[ NNarray(i-1,j)-1 ];
        }

        // initialize Li
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

        // first row has just a single nonzero entry
        Li[1-1][1-1] = pow( covfun3(0, cparms), -0.5 );

        for(j=2; j<m+1; j++){  // j = row of Li

        // initialize g
        for(k=1; k<m+1; k++){
            g[k-1] = 0.0;
        }

        // get first j-1 entries of jth row of L (not Linverse!)
        for(k=1; k<j; k++){
            d = pow(  pow(locsub[k-1][0] - locsub[j-1][0],2) +
                          pow(locsub[k-1][1] - locsub[j-1][1],2) +
                          pow(locsub[k-1][2] - locsub[j-1][2],2)   , 0.5 );
            sig[k-1] = covfun3(d,cparms);
            g[k-1] = 0.0;
            for(el=1; el<k+1; el++){
                g[k-1] += Li[k-1][el-1]*sig[el-1];
            }
        }

        // get diagonal entry
        d = 0.0;
        for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
        Li[j-1][j-1] = pow( covfun3(0,cparms) - d, -0.5 );

        // now get first j-1 entries jth row of Linverse
        for(k=1; k<j; k++){
            for(el=1; el<j+1; el++){
                Li[j-1][k-1] += g[el-1]*Li[el-1][k-1];
            }
            Li[j-1][k-1] = -Li[j-1][k-1]*Li[j-1][j-1];
        }

        }
        d = 0.0;
        if(i==m){
            for(k=1; k<m+1; k++){
                d = 0.0;
                for(el=1; el<k+1; el++){
                    d += Li[k-1][el-1]*ysub[el-1];
                }
                z(NNarray(m-1,m-k)-1) = d;
            }
        } else {
            d = 0.0;
            for(k=1;k<m+1;k++){
                d += Li[m-1][k-1]*ysub[k-1];
                //printf("%6.3f \n",d);
            }
            z(i-1) = d;
        }
    }

    return z;
}



// this returns the diagonal entries of L^{-1}
// same code as before. this could be simlified to
// just do a cholesky decomposition (instead of inverse cholesky)
// [[Rcpp::export]]
NumericVector getinvcondsd3d(NumericVector covparms, NumericVector y,
                       NumericMatrix locs, IntegerMatrix NNarray) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    double cparms[4] = {covparms[0], covparms[1], covparms[2], covparms[3]};

    int n = y.length();
    NumericVector invcondsd(n);

    // number of neighbors + 1
    int m = NNarray.ncol();

    double d;
    double ysub[m];
    double locsub[m][3];

    double Li[m][m];

    double g[m];
    double sig[m];



    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            locsub[m-1-j][0] = locs( NNarray(i-1,j)-1, 0 );
            locsub[m-1-j][1] = locs( NNarray(i-1,j)-1, 1 );
            locsub[m-1-j][2] = locs( NNarray(i-1,j)-1, 2 );
            ysub[m-1-j] = y[ NNarray(i-1,j)-1 ];
        }

        // initialize Li
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

        // first row has just a single nonzero entry
        Li[1-1][1-1] = pow( covfun3(0, cparms), -0.5 );

        for(j=2; j<m+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<m+1; k++){
                g[k-1] = 0.0;
            }

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){
                d = pow(  pow(locsub[k-1][0] - locsub[j-1][0],2) +
                    pow(locsub[k-1][1] - locsub[j-1][1],2) +
                    pow(locsub[k-1][2] - locsub[j-1][2],2)   , 0.5 );
                sig[k-1] = covfun3(d,cparms);
                g[k-1] = 0.0;
                for(el=1; el<k+1; el++){
                    g[k-1] += Li[k-1][el-1]*sig[el-1];
                }
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            Li[j-1][j-1] = pow( covfun3(0,cparms) - d, -0.5 );

            // now get first j-1 entries jth row of Linverse
            for(k=1; k<j; k++){
                for(el=1; el<j+1; el++){
                    Li[j-1][k-1] += g[el-1]*Li[el-1][k-1];
                }
                Li[j-1][k-1] = -Li[j-1][k-1]*Li[j-1][j-1];
            }

        }

        if(i==m){
            for(k=1; k<m+1; k++){
                invcondsd(k-1) = Li[k-1][k-1];
            }
        } else {
            invcondsd(i-1) = Li[m-1][m-1];
        }
    }

    return invcondsd;
}
