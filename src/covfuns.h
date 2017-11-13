
#ifndef COVFUNS_H
#define COVFUNS_H

inline double Matern_from_dist(double d, double *cparms){

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
            Rf_bessel_k(d/cparms[1],cparms[2],1.0)*cparms[3];
        }
    }
    return d;
}

inline double Matern_from_locs( std::vector<double> loc1, std::vector<double> loc2, double *cparms){

    int dim = loc1.size();
    double d = 0.0;
    int i;
    for(i=0; i<dim; i++){
        d += pow( loc1[i] - loc2[i], 2);
    }
    d = pow(d, 0.5);

    d = Matern_from_dist(d,cparms);
    return d;


}


#endif

