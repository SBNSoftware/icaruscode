#include "FitBackgroundPhotons.h"
#include "TCanvas.h"

pmtcalo::FitBackgroundPhotons::FitBackgroundPhotons( unsigned int nparameters, SPEChargeResponse *fitf )
{
  m_nparameters=nparameters;
  m_fitf = new TF1( "fitFunction", fitf, m_fitrange[0], m_fitrange[1], m_nparameters);
}

void pmtcalo::FitBackgroundPhotons::setFitParameters( std::vector<double> params )
{

  // set parameters
  for(unsigned int i=0; i<params.size(); i++ )
    m_fitf->SetParameter( i, params[i] );
	
  // set reasonable limits for the poissonGauss
  m_fitf->SetParLimits( 0, 0.0, 2.0 );
  m_fitf->SetParLimits( 1, 0.0, 2.0 );
  m_fitf->SetParLimits( 2, 0.0, 2.0 );
  m_fitf->SetParLimits( 3, 0.0, params[3]*12.5);

}

void pmtcalo::FitBackgroundPhotons::fitHistogram( TH1D *hist, std::string fitoption )
{

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

  // final assigment of the parameters and parameters errors
  for(unsigned int i=0; i<m_nparameters; i++){
    m_vresults[i] = m_fitf->GetParameter(i);
    m_verrors[i] = m_fitf->GetParError(i);
  }

  m_fitstatus=fitstatus;

}

void pmtcalo::FitBackgroundPhotons::SaveToPdf( TH1D *hist, TH1D* histcut, std::string name, int cut )
{

  char cname[200];
  sprintf( cname, "c_%s", hist->GetName() );
  TCanvas *c = new TCanvas(cname,"fit",800.,600.); 	
        
  hist->Draw();
  if (cut == 1)
  {
    histcut->SetLineColor(kMagenta);
    histcut->Draw("SAME");
  }
  
  m_fitf->Draw("SAME");
  
  // add canvas to pdf
  c->Print(name.c_str(),"pdf");

}
