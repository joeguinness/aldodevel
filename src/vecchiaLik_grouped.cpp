#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include "covfuns.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector vecchiaLik_grouped(NumericVector covparms, NumericVector y,
                             NumericMatrix locs, List NNlist ) {

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
    double d;

    int nblock = NNlist.length();

    vector<double> loc1(dim);
    vector<double> loc2(dim);

    // dynamic memory allocation for cholesky
    //double** L;

    //const auto covfun = Matern_from_dist;


    for(i=0; i<nblock; i++){

        // extract the block that we need
        IntegerMatrix NNarray = as<IntegerMatrix>(NNlist[i]);
        int m = NNarray.ncol();
        int r = NNarray.nrow();

        // save the indices of the responses for this block
        vector<int> resp_inds(r,0);
        for(j=0; j<r; j++){ resp_inds[j] = NNarray(j,0); }
        sort( resp_inds.begin(), resp_inds.end() );

        // put all of the indices in a vector
        vector<int> all_inds (r*m,0);
        for(j=0; j<r; j++){
            for(k=0; k<m; k++){
                //cout << NNarray(j,k) << endl;
                all_inds[j*m+k] = NNarray(j,k);
            }
        }

        // get rid of duplicates
        sort(all_inds.begin(),all_inds.end());
        all_inds.erase( unique( all_inds.begin(), all_inds.end() ), all_inds.end() );
        // get rid of first value if negative
        // this is because NAs are converted to negative integers
        if( all_inds[0] < 0 ){
            all_inds.erase( all_inds.begin() );
        }

        int blocksize = all_inds.size();

        //for(j=0; j<r; j++){ cout << resp_inds[j] << " " ; }
        //cout << endl;

        // figure out the indices in all_inds that contain resp_inds
        vector<int> which_resp_inds(r,0);
        int k = 0;
        for(j=0; j<blocksize; j++){
            if( all_inds[j] == resp_inds[k] ){
                which_resp_inds[k] = j;
                k++;
            }
        }

        // get the covariance matrix
        // be careful about indexing here!
        // which_resp_inds uses 0-indexing
        // all_inds and resp_inds use 1-indexing
        double locsub[blocksize][dim];
        double ysub[blocksize];

        for(j=0; j<blocksize; j++){
            for(k=0; k<dim; k++){
                locsub[j][k] = locs( all_inds[j]-1, k );
            }
            ysub[j] = y[ all_inds[j]-1 ];
        }

        // initialize variables needed to compute L
        double L[blocksize][blocksize];
        //L = new double* [blocksize];
        //for(j=0; j<blocksize; j++){ L[j] = new double [j+1]; }

        double g[blocksize];
        double sig[blocksize];
//        for(k=0;k<blocksize;k++){ for(j=0;j<blocksize;j++){ Li[k][j] = 0.0; }}

        // get 1,1 entry
        L[0][0] = pow( MaternFunction(0.0, cparms) + nugget, 0.5 );

        // get the rest of the entries
        for(j=2; j<blocksize+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<blocksize+1; k++){ g[k-1] = 0.0; }

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){
                d = 0.0;
                for(el=0; el< dim; el++){
                    d+= pow( locsub[k-1][el] - locsub[j-1][el], 2 );
                }
                d = pow( d, 0.5 );

                sig[k-1] = MaternFunction(d,cparms);
                //g[k-1] = 0.0;
                g[k-1] = sig[k-1];
                //for(el=1; el<k+1; el++){
                if(k>1){
                    for(el=1; el<k; el++){
                        //g[k-1] += Li[k-1][el-1]*sig[el-1];
                        g[k-1] -= L[k-1][el-1]*g[el-1];
                    }
                }
                g[k-1] = g[k-1]/L[k-1][k-1];
                L[j-1][k-1] = g[k-1];
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            L[j-1][j-1] = pow( MaternFunction(0,cparms) + nugget - d, 0.5 );

        }

        // now we have L and we can compute
        // contribution to likelihood

        // get g = L^{-1}y
        g[0] = ysub[0]/L[0][0];
        for(j=1; j<blocksize; j++){
            g[j] = ysub[j];
            for(k=0; k<j; k++){
                g[j] -= L[j][k]*g[k];
            }
            g[j] = g[j]/L[j][j];
        }

        // loop over the responses
        for(j=0; j<r; j++){
            k = which_resp_inds[j];
            ll(0) += -g[k]*g[k]/2 - log( L[k][k] );
        }

        // deallocate memory
        //delete[] L;

    }


    ll(0) += -n*log(2*M_PI)/2;
    return ll;
}
