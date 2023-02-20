#include "FitBackgroundPhotons.h"

FitBackgroundPhotons::FitBackgroundPhotons(){}

FitBackgroundPhotons::FitBackgroundPhotons( unsigned int nparameters, IdealPmtResponse fitf )
{
	m_nparameters=nparameters;
	m_fitf = new TF1( "fitFunction", fitf, m_fitrange[0], m_fitrange[1], m_nparameters);

}

void FitBackgroundPhotons::fitHistogram( TH1D *hist, std::string fitoption )
{

	m_fitf->SetParameters( 0.1, 0.8, 0.3, hist->Integral()*0.2, hist->Integral()*0.8, -3. );
	m_fitf->SetParLimits( 0, 0.0, 2.0 );
	m_fitf->SetParLimits( 1, 0.0, 2.0 );
	m_fitf->SetParLimits( 2, 0.0, 2.0 );
	m_fitf->SetParLimits( 3, 0.0, hist->Integral()*10 );
	m_fitf->SetParLimits( 4, 0.0, hist->Integral()*10 );
	m_fitf->SetParLimits( 5, -10, 0.0 );

	int fitstatus = hist->Fit("fitFunction", fitoption.c_str(), "", m_fitrange[0], m_fitrange[1]);

	int nrefit = 0;

	m_vresults.resize(m_nparameters);
	m_verrors.resize(m_nparameters);

	while( fitstatus !=0 && nrefit < m_maxrefit ){

		for(unsigned int i=0; i<m_nparameters; i++){
			m_vresults[i] = m_fitf->GetParameter(i);
			m_verrors[i] = m_fitf->GetParError(i);
		}

		fitstatus = hist->Fit("fitFunction", fitoption.c_str(), "", m_fitrange[0], m_fitrange[1]);

		nrefit++;

	}

	// Final assigment of the parameters and parameters errors
	for(unsigned int i=0; i<m_nparameters; i++){
		m_vresults[i] = m_fitf->GetParameter(i);
		m_verrors[i] = m_fitf->GetParError(i);
	}

	m_fitstatus=fitstatus;

}