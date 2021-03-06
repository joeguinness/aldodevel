#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "covfuns.h"

using namespace std;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
NumericVector OrderedCompLik(NumericVector covparms, StringVector covfun_name,
                             NumericVector y,
                             NumericMatrix locs, IntegerMatrix NNarray) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    double cparms[4] = {covparms[0], covparms[1], covparms[2], covparms[3]};

    NumericVector ll(1);
    int n = y.length();
    int dim = locs.ncol();

    // number of neighbors + 1
    int m = NNarray.ncol();

    double d;
    double ysub[m];
    //vector< vector <double> > locsub[m][2];
    std::vector<std::vector<double> > locsub(m, std::vector<double>(dim, 0));


    double Li[m][m];

    double g[m];
    double sig[m];


    vector<double> loc1(dim);
    vector<double> loc2(dim);

    //const auto covfun = Matern_from_dist;


    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            for(k=0;k<dim;k++){ locsub[m-1-j][k] = locs( NNarray(i-1,j)-1, k ); }
            ysub[m-1-j] = y[ NNarray(i-1,j)-1 ];
        }


        for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

        loc1 = locsub[0];
        Li[1-1][1-1] = pow( MaternFunction(0, cparms), -0.5 );

        for(j=2; j<m+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<m+1; k++){
                g[k-1] = 0.0;
            }

            //loc1 = locsub[j-1];

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){

                //loc2 = locsub[k-1];
                for(el=0;el<dim;el++){ loc2[el] = locsub[k-1][el]; }

                d = 0.0;
                for(el=0;el<dim;el++){ d += pow(locsub[k-1][el]-locsub[j-1][el],2); }
                d = pow( d, 0.5 );

                sig[k-1] = MaternFunction(d,cparms);
                //sig[k-1] = Matern_from_locs(loc1,loc2,cparms);

                g[k-1] = 0.0;
                for(el=1; el<k+1; el++){
                    g[k-1] += Li[k-1][el-1]*sig[el-1];
                }
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            Li[j-1][j-1] = pow( MaternFunction(0,cparms) - d, -0.5 );

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
                ll(0) += -d*d/2 + log( Li[k-1][k-1] );
            }
        } else {
            d = 0.0;
            for(k=1;k<m+1;k++){
                d += Li[m-1][k-1]*ysub[k-1];
                //printf("%6.3f \n",d);
            }
            ll(0) += -d*d/2 + log( Li[m-1][m-1] );
        }
    }
    ll(0) += -n*log(2*M_PI)/2;

    return ll;
}
