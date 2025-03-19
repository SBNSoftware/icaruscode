#ifndef __SPEChargeResponse_H
#define __SPEChargeResponse_H

#include "TMath.h"

namespace pmtcalo{

  class  SPEChargeResponse {
 
    public:
    
      SPEChargeResponse(){};
      SPEChargeResponse(int nstart, int nsum, bool useExpPedestal) 
        : m_nstart(nstart),
          m_nsum(nsum), 
          m_useExpPedestal(useExpPedestal)
        {};

      virtual ~SPEChargeResponse() {};
      virtual double operator()(double *x, double *par) const = 0; // overridden

    protected: 

      int m_nstart=1;
      int m_nsum=2;
      bool m_useExpPedestal=false;

      double pedestal( double x, double a0, double c0 ) const
      {
        return a0*TMath::Exp( x*c0 );
      }

      double poissonGauss( double x, double amplitude, double n, double mu, double q, double sigma ) const
      {
        double poisson = TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n);
        double gaus = TMath::Exp(-1.0*(x-q*n)*(x-q*n)/(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*TMath::Pi()*n));
        return amplitude*poisson*gaus;  
      }

      double logisticFunction( double x, double x0, double k) const
      {
        return 1./(1.+TMath::Exp(-k*(x-x0)));
      }
  };

}

#endif
