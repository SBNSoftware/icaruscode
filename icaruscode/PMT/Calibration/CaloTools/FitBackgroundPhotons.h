#ifndef __FitBackgroundPhotons_H
#define __FitBackgroundPhotons_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TF1.h"
#include "SPEChargeResponse.h"

namespace pmtcalo{

  class FitBackgroundPhotons {

    public: 

      FitBackgroundPhotons() {};
      FitBackgroundPhotons( unsigned int nparameters, pmtcalo::SPEChargeResponse *fitf );

      // preparing and executing the fit
      void setFitRange( float low, float high ){ m_fitrange[0]=low; m_fitrange[1]=high; };
      void setFitParameters(std::vector<double> params, float qmin=0.1, float qmax=1.0);
      void fixFitParameter(int num, float param){ m_fitf->FixParameter(num, param); };
      void fitHistogram( TH1D *hist, std::string fitoption="RQ0+" );

      // extracting the fit results
      float getParameter(int num){ return m_vresults[num]; };
      float getParError(int num){ return m_verrors[num]; };
      float getChi2(){ return m_fitf->GetChisquare(); };
      int getNDF(){ return m_fitf->GetNDF(); };
      int getFitStatus(){ return m_fitstatus; };
      void getFitRange( float &low, float &high ){ low=m_fitrange[0]; high=m_fitrange[1]; };

      // add fit to pdf
      void SaveToPdf( TH1D *hist, TH1D *histcut, std::string name, int cut=1);

      private: 

        unsigned int m_nparameters=6;
        int m_maxrefit=4;
        int m_fitstatus=-99;

        float m_fitrange[2] = { 0.0, 2 };

        TF1 *m_fitf;

        std::vector<float> m_vresults; //(m_nparameters, 1.);
        std::vector<float> m_verrors;// (m_nparameters, 1.);

  };
}
#endif
