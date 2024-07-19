#ifndef __SPEChargeIdealResponse_H
#define __SPEChargeIdealResponse_H

#include "TMath.h"

namespace pmtcalo{

  class  SPEChargeIdealResponse : public SPEChargeResponse {
  
    public:

      SPEChargeIdealResponse(){};
      SPEChargeIdealResponse( int nstart, int nsum, bool useExpPedestal ) 
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

        if( m_useExpPedestal )
        {  
          double a0 = par[4];
  	  double c0 = par[5];
	  val += pedestal( x[0], a0, c0 );
        }

        return val;
      };
  };

}
#endif
