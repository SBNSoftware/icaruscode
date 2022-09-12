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
	m_hist=hist;

}

void FitBackgroundPhotons::saveToCanvas( int pmt, bool m_exponential, TCanvas & canvas )
{


	m_hist->GetXaxis()->SetRangeUser( m_fitrange[0]-0.5, m_fitrange[1]+1 );
	m_hist->GetYaxis()->SetTitle( "Number of pulses" );
	m_hist->Draw("hist E1");

	m_fitf->SetLineStyle(2);
	m_fitf->SetLineColor(kRed);
	m_fitf->Draw("same");

	char gname[100];
  	TPaveText* ptGainPar = new TPaveText(0.6, 0.55, 0.9, 0.9,"brNDC");
  	sprintf(gname,"PMT: %d", pmt);
  	ptGainPar->AddText(gname);
  	sprintf( gname,"#chi^2 / ndf = %.0f/%d", m_fitf->GetChisquare(), m_fitf->GetNDF() );
  	ptGainPar->AddText(gname);
 	sprintf(gname,"#mu = %.2e #pm %.2e", m_vresults[0], m_verrors[0] );
  	ptGainPar->AddText(gname);
  	sprintf(gname,"q = %.2e #pm %.2e", m_vresults[1], m_verrors[1] );
 	ptGainPar->AddText(gname);
  	sprintf(gname,"#sigma = %.2e #pm %.2e", m_vresults[2], m_verrors[2] );
  	ptGainPar->AddText(gname);
  	sprintf(gname,"A = %.2e #pm %.2e", m_vresults[3], m_verrors[3] );
  	ptGainPar->AddText(gname);
  		
  	if( m_exponential )
  	{
  		sprintf(gname,"A_{0} = %.2e #pm %.2e", m_vresults[4], m_verrors[4] );
  		ptGainPar->AddText(gname);
  		sprintf(gname,"c_{0} = %.2e #pm %.2e", m_vresults[5], m_verrors[5] );
  		ptGainPar->AddText(gname);
  	}

  	ptGainPar->SetBorderSize(0);
  	ptGainPar->SetFillColor(0);
  	ptGainPar->SetTextFont(42);
  	ptGainPar->SetTextAlign(12);
  	ptGainPar->Draw();


  	return;

}