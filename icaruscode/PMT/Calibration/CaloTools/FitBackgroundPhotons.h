#ifndef  __FitBackgroundPhotons_H
#define __FitBackgroundPhotons_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "IdealPmtResponse.h"


//namespace pmtcalo{
//	class FitBackgroundPhotons;
//}

class FitBackgroundPhotons {

	public: 

		FitBackgroundPhotons();
		FitBackgroundPhotons( unsigned int nparameters, IdealPmtResponse fitf );

		void fitHistogram( TH1D *hist, std::string fitoption="RQ0+" );

		float getParameter(int num){ return m_vresults[num]; };
		float getParError(int num){ return m_verrors[num]; };
		float getChi2(){ return m_fitf->GetChisquare(); };
		int getNDF(){ return m_fitf->GetNDF(); };
		int getFitStatus(){ return m_fitstatus; };

		void setFitRange( float low, float high ){ m_fitrange[0]=low; m_fitrange[1]=high; };
		void getFitRange( float &low, float &high ){ low=m_fitrange[0]; high=m_fitrange[1]; };

		void saveToCanvas( int pmt, bool m_exponential, TCanvas & canvas );

	private: 

		unsigned int m_nparameters=6;
		int m_maxrefit=4;
		int m_fitstatus=-99;

		float m_fitrange[2] = { 0.1, 2 };

		TH1D * m_hist;
		TF1 *m_fitf;

		//float m_start[m_nparameters]
		std::vector<float> m_vresults; //(m_nparameters, 1.);
		std::vector<float> m_verrors;// (m_nparameters, 1.);

};

#endif