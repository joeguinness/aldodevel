
#ifndef COVFUNS_H
#define COVFUNS_H

inline double covfun(double d, double *cparms){

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

#endif

