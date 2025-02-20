#ifndef __SPEChargeLogisticResponse_H
#define __SPEChargeLogisticResponse_H

#include "TMath.h"

namespace pmtcalo{

  class  SPEChargeLogisticResponse : public SPEChargeResponse {
  
    public:

      SPEChargeLogisticResponse(){};
      SPEChargeLogisticResponse( int nstart, int nsum, bool useExpPedestal ) 
        : SPEChargeResponse(nstart,nsum,useExpPedestal) {};

      double operator() (double *x, double *par) const override 
      {
        double mu = par[0];
        double q = par[1];
        double sigma = par[2];
        double amplitude = par[3];

        double val = 0;

        for( int n=m_nstart; n<m_nstart+m_nsum; n++ )
          val+=poissonGauss( x[0], amplitude, n, mu, q, sigma );

        double x0 = par[4];
        double k = par[5];

        val *= logisticFunction( x[0], x0, k );
  
        return val;
      };
  };

}

#endif
