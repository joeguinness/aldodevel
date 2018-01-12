
#ifndef COVFUNS_H
#define COVFUNS_H


inline double MaternFunction(double d, double *cparms){

    // has special cases for 1/2 and 3/2
    if( d == 0.0 ){
        d = cparms[0];
    } else {
        if( cparms[2] == 0.5 ){
            d = cparms[0]*exp(-d/cparms[1]);
        } else if( cparms[2] == 1.5 ){
            d = cparms[0]*(1+d/cparms[1])*exp(-d/cparms[1]);
        } else {
            double normcon = cparms[0]/(pow(2.0,cparms[2]-1)*Rf_gammafn(cparms[2]));
            d = normcon*pow( d/cparms[1], cparms[2] )*
            Rf_bessel_k(d/cparms[1],cparms[2],1.0);
        }
    }
    return d;
}

inline double MaternIsotropic( const std::vector<double>* loc1, const std::vector<double>* loc2, double* cparms){

    // inputs two vectors loc1 and loc2 of length dim and returns
    // isotropic matern covariance using euclidean distance
    // between loc1 and loc2
    int dim = (*loc1).size();
    double d = 0.0;
    int i;
    for(i=0; i<dim; i++){
        d += pow( (*loc1)[i] - (*loc2)[i], 2);
    }
    d = pow(d, 0.5);

    d = MaternFunction(d,cparms);
    return d;

}


inline double MaternSphere( const std::vector<double>* lonlat1, const std::vector<double>* lonlat2, double* cparms){

    // lonlat1, lonlat2: vectors of length 2 holding longitude and latitude
    // lon \in (-180,180), lat \in (-90,90)
    double d;
    double lonrad1 = 2*M_PI*(*lonlat1)[0]/360;
    double lonrad2 = 2*M_PI*(*lonlat2)[0]/360;
    double latrad1 = 2*M_PI*((*lonlat1)[1]+90)/360;
    double latrad2 = 2*M_PI*((*lonlat2)[1]+90)/360;
    double x1 = sin(latrad1)*cos(lonrad1);
    double x2 = sin(latrad2)*cos(lonrad2);
    double y1 = sin(latrad1)*sin(lonrad1);
    double y2 = sin(latrad2)*sin(lonrad2);
    double z1 = cos(latrad1);
    double z2 = cos(latrad2);
    /*
    std::vector<double> loc1(3);
    std::vector<double> loc2(3);
    loc1[0] = x1;
    loc1[1] = y1;
    loc1[2] = z1;
    loc2[0] = x2;
    loc2[1] = y2;
    loc2[2] = z2;
    d = MaternIsotropic( &loc1, &loc2, cparms );
    */
    d = pow( x1-x2, 2) + pow(y1 - y2, 2) + pow( z1-z2, 2);
    d = pow(d, 0.5);
    d = MaternFunction(d, cparms);
    return d;
}



#endif

