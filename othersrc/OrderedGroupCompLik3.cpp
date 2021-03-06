#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include "covfuns.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector OrderedGroupCompLik3(NumericVector covparms, NumericVector y,
                                  NumericMatrix locs, List NNlist ) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    int blocksize;
    double cparms[4] = {covparms[0], covparms[1], covparms[2], covparms[3]};

    NumericVector ll(1);
    int n = y.length();
    int dim = locs.ncol();

    // number of neighbors + 1
    double d;

    int nblock = NNlist.length();

    // figure out what the biggest block size is
    IntegerMatrix NNarray = as<IntegerMatrix>(NNlist[0]);
    int m = NNarray.ncol();
    int rmax = NNarray.nrow();
    int r = rmax;
    for(i=0; i<nblock; i++){
        IntegerMatrix NNarray = as<IntegerMatrix>(NNlist[i]);
        r = NNarray.nrow();
        if(r > rmax){ rmax = r; }
    }

    double locsub[rmax*m][dim];
    double ysub[rmax*m];
    double Li[rmax*m*(rmax*m+1)/2];
    double g[rmax*m];
    double sig[rmax*m];

    for(i=0; i<nblock; i++){

        // extract the block that we need
        IntegerMatrix NNarray = as<IntegerMatrix>(NNlist[i]);
        r = NNarray.nrow();

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

        blocksize = all_inds.size();

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
        for(j=0; j<blocksize; j++){
            for(k=0; k<dim; k++){
                locsub[j][k] = locs( all_inds[j]-1, k );
            }
            ysub[j] = y[ all_inds[j]-1 ];
        }

        // initialize variables needed to compute Linverse
        for(k=0;k<blocksize*(blocksize+1)/2;k++){ Li[k] = 0.0; }

        // get 1,1 entry
        Li[0] = pow( covfun(0, cparms), -0.5 );
        // get the rest of the entries
        for(j=2; j<blocksize+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<blocksize+1; k++){ g[k-1] = 0.0; }

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){
                d = pow( pow(locsub[k-1][0] - locsub[j-1][0],2) +
                             pow(locsub[k-1][1] - locsub[j-1][1],2), 0.5 );
                sig[k-1] = covfun(d,cparms);
                g[k-1] = 0.0;
                for(el=1; el<k+1; el++){
                    g[k-1] += Li[(k-1)*k/2+el-1]*sig[el-1];
                }
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            Li[(j-1)*j/2+j-1] = pow( covfun(0,cparms) - d, -0.5 );

            // now get first j-1 entries jth row of Linverse
            for(k=1; k<j; k++){
                for(el=1; el<j+1; el++){
                    Li[(j-1)*j/2+k-1] += g[el-1]*Li[(el-1)*el/2+k-1];
                }
                Li[(j-1)*j/2+k-1] = -Li[(j-1)*j/2+k-1]*Li[(j-1)*j/2+j-1];
            }

        }
/**
        // now we have Linverse, and we can compute
        // contribution to likelihood

        // loop over the responses
        for(j=0; j<r; j++){
            d = 0.0;
            k = which_resp_inds[j];
            for(el=0; el<k+1; el++){
                d+= Li[k*(k+1)/2+el]*ysub[el];
            }
            ll(0) += -d*d/2 + log( Li[k*(k+1)/2+k] );
        }
 **/

    }

    ll(0) += -n*log(2*M_PI)/2;
    return ll;
}
