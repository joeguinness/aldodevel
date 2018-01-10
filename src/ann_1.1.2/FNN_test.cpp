#include <Rcpp.h>

using namespace Rcpp;



// [[Rcpp::export]]
NumericMatrix get_KNNX_kd_test(NumericMatrix data, NumericMatrix query,
                 IntegerVector k, IntegerVector dim,
                 IntegerVector n_pts, IntegerVector m_pts,
                 IntegerMatrix nn_idx, NumericMatrix nn_dist)
{
    const int				d=*dim;		// Actual Dimension
    const int				n=*n_pts;	// Number of Base Data points
    const int				m=*m_pts;	// Number of Query Data points
    const int				K = *k;		// Max. num of NN
    const double		error_bound = 0.0;

    ANNidxArray		index = new ANNidx[K];		// Near neighbor indices
    ANNdistArray	dist = new ANNdist[K];		// Near neighbor squared distance
    //	ANNpointArray	data_pts  = annAllocPts(n, d);	// Base Data points
    //	ANNpointArray	query_pts  = annAllocPts(m, d);	// Query Data points
    ANNpointArray	data_pts  = new ANNpoint[n];
    ANNpointArray	query_pts  = new ANNpoint[m];

    if(data_pts==NULL) error("Cannot allocate memroy for data matrix!\n");
    if(query_pts==NULL) error("Cannot allocate memroy for query data matrix!\n");

    //copy data into  matrix from a vector column by column
    Rvector2ANNarray(data_pts, data, n, d);
    Rvector2ANNarray(query_pts, query, m, d);

    ANNkd_tree	*kd_tree	= new ANNkd_tree(data_pts, n, d);

    register int ptr = 0;
    for(int i = 0; i < m; i++)	// read query points
    {
        kd_tree->annkSearch(	// search
                query_pts[i],		// query point
                         K,				// number of near neighbors
                         index,				// nearest neighbors index
                         dist,				// squared distance
                         error_bound);		// error bound

        for (int j = 0; j < K; j++)      //return result row by row
        {
            nn_dist[ptr] = sqrt(dist[j]);	// unsquare distance
            nn_idx[ptr]  = index[j] + 1;	// put indexes in returned array
            ptr++;
        } // end inner for
    } // end for

    delete[] index;
    delete[] dist;
    //	delete kd_tree;

    //	annDeallocPts(data_pts);
    //	annDeallocPts(query_pts);
    delete[] data_pts;
    delete[] query_pts;
    delete kd_tree;


    annClose();

    return nn_dist;

}//end get_KNNX_kd
