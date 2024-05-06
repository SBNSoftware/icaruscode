#ifndef ICARUSCODE_CRT_CRTUTILS_TRIPLEMATCHINGALG_H
#define ICARUSCODE_CRT_CRTUTILS_TRIPLEMATCHINGALG_H


//#include "icaruscode/IcarusObj/CRTPMTMatching.h"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "canvas/Persistency/Common/Ptr.h" 

namespace icarus::triplematching {

struct Direction {
    double dirx; // Direction of the Track: X 
    double diry; // Direction of the Track: Y 
    double dirz; // Direction of the Track: Z 
    double meanx; // Mean Point of the Track: X 
    double meany; // Mean Point of the Track: Y 
    double meanz; // Mean Point of the Track: Z 
};

Direction PCAfit (std::vector<float> x, std::vector<float> y, std::vector<float> z, float cutmin, float cutmax){
	int min=cutmin*x.size();
    int max=cutmax*x.size();
    int size=max-min;
    double xavg=0, yavg=0, zavg=0;
    for(int k=min; k<max; k++) {
        xavg+=x[k]; yavg+=y[k]; zavg+=z[k];
    }
    xavg=xavg/size;
    yavg=yavg/size;
    zavg=zavg/size;
    float x2=0, y2=0, z2=0, xiy=0, xiz=0, yiz=0;
    for(int k=min; k<max; k++){
        float xadj = x[k] - xavg;
        float yadj = y[k] - yavg;
        float zadj = z[k] - zavg;
        x2+=xadj*xadj;
        y2+=yadj*yadj;
        z2+=zadj*zadj;
        xiy+=xadj*yadj;
        xiz+=xadj*zadj;
        yiz+= yadj*zadj;
    }//end loop calculating covariance matrix elements
    x2=x2/size;
    y2=y2/size;
    z2=z2/size; 
    xiy=xiy/size;
    xiz=xiz/size;
    yiz=yiz/size;
    //Now covariance elements are ready, make covariance matrix and do PCA analysis
    TMatrixD covmat(3,3);
    covmat[0][0] = x2; covmat[1][1]=y2; covmat[2][2]=z2;
    covmat[0][1] = xiy; covmat[1][0] = xiy;
    covmat[0][2] = xiz; covmat[2][0] = xiz;
    covmat[2][1] = yiz; covmat[1][2] = yiz;
    TMatrixDEigen thiscov(covmat);
    TMatrixD this_eval = thiscov.GetEigenValues();
    TMatrixD this_evec = thiscov.GetEigenVectors();
    if(this_eval.GetNrows()!=3 || this_eval.GetNcols()!=3 || this_evec.GetNrows()!=3 || this_evec.GetNcols()!=3){ 
        std::cout << "evals/evects wrong size! continuing....\n";
    }//end if evals/evects aren't 3x3 matrices as expected

    int max_eval = -999999; int maxevalpos = -1;
    for(int k = 0; k < 3; k++){
        if(this_eval[k][k]>max_eval){
            max_eval = this_eval[k][k];
            maxevalpos = k;
        }//end if(this_eval[k][k]>max_eval)
    }//end loop looking for best eval
    		
    TripleMatchingUtils::Direction Direction= {this_evec[0][maxevalpos] , this_evec[1][maxevalpos], this_evec[2][maxevalpos], xavg, yavg, zavg};
    return Direction;
};

}