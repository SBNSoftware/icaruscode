#ifndef  __IdeaPmtResponse_H
#define __IdeaPmtResponse_H

#include "TMath.h"

namespace pmtcalo{
	class IdealPmtResponse;
}

class  IdealPmtResponse {
 
	public:

		IdealPmtResponse(){};

		IdealPmtResponse( int nstart, int nsum, bool useExpPedestal ) 
		: m_nstart(nstart),
		  m_nsum(nsum), 
		  m_useExpPedestal(useExpPedestal)
		{};


		double operator() (double *x, double *par) {

			double mu = par[0];
    		double q = par[1];
    		double sigma = par[2];
    		double amplitude = par[3];

			double val = 0;

			for( int n=m_nstart; n<m_nstart+m_nsum; n++ )
				val+=poissonGauss( x[0], amplitude, n, mu, q, sigma );

			if( m_useExpPedestal ){

				double a0 = par[4];
  	  			double c0 = par[5];
				val += pedestal( x[0], a0, c0 );
			}

			return val;
      
		};

	private: 

		int m_nstart=1;
		int m_nsum=4;
		bool m_useExpPedestal=false;

		double pedestal( double x, double a0, double c0 ){

			return a0*TMath::Exp( x*c0 );
             
		}

		double poissonGauss( double x, double amplitude, double n, double mu, double q, double sigma ){

			return amplitude*(TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)
                *TMath::Exp(-1.0*(x-q*n)*(x-q*n)/(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*TMath::Pi()*n))) ;

		};

};

#endif
