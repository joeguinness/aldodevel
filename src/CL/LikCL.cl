

    // first, fill in ysub and locsub in reverse order
    // will need to pass smaller versions of NN and y
    for(j=m-1; j>=0; j--){
        locsub[m-1-j][0] = locs( NN[i-1,j]-1, 0 );
        locsub[m-1-j][1] = locs( NN[i-1,j]-1, 1 );
        ysub[m-1-j] = y[ NN(i-1,j)-1 ];
    }

    // initialize Li
    double Li[m][m];
    for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

    Li[1-1][1-1] = pow( covfun(0, cparms), -0.5 );

    for(j=2; j<m+1; j++){  // j = row of Li

        // initialize g

        for(k=1; k<m+1; k++){
            g[k-1] = 0.0;
        }


        // get first j-1 entries of jth row of L (not Linverse!)
        for(k=1; k<j; k++){
            d = pow( pow(locsub[k-1][0] - locsub[j-1][0],2) +
            pow(locsub[k-1][1] - locsub[j-1][1],2), 0.5 );
            sig[k-1] = covfun(d,cparms);
            g[k-1] = 0.0;
            for(el=1; el<k+1; el++){
                g[k-1] += Li[k-1][el-1]*sig[el-1];
            }
        }

        // get diagonal entry
        d = 0.0;
        for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
        Li[j-1][j-1] = pow( covfun(0,cparms) - d, -0.5 );

        // now get first j-1 entries jth row of Linverse
        for(k=1; k<j; k++){
            for(el=1; el<j+1; el++){
                Li[j-1][k-1] += g[el-1]*Li[el-1][k-1];
            }
            Li[j-1][k-1] = -Li[j-1][k-1]*Li[j-1][j-1];
        }
    }

    d = 0.0;
    for(k=1;k<m+1;k++){
        d += Li[m-1][k-1]*ysub[k-1];
    }
    ll(0) += -d*d/2 + log( Li[m-1][m-1] );
