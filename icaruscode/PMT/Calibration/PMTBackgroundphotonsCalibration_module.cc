////////////////////////////////////////////////////////////////////////
// Class:       PMTBackgroundphotonsCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTBackgroundphotonsCalibration_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
//
//  mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////


#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Utilities/Exception.h"

#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

#include "icaruscode/PMT/Calibration/CaloTools/Utils.h"

#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"

namespace pmtcalo {
  class PMTBackgroundphotonsCalibration;
  class Gaussian;
}


  // Fit function for the pulse
class pmtcalo::Gaussian {
    public:
      // use constructor to customize the function object
      double operator() (double* x, double * par){
        
        double A = par[0];
        double mu = par[1];
        double sigma = par[2];

        double arg = (( x[0]-mu )*( x[0]-mu )) / (sigma*sigma) ;
        double val = A/( 2*sqrt(2*TMath::Pi()) )*exp( -arg/2.0 );

        return val;
      }
};


class pmtcalo::PMTBackgroundphotonsCalibration : public art::EDAnalyzer {

public:

  explicit PMTBackgroundphotonsCalibration(fhicl::ParameterSet const& pset);

  PMTBackgroundphotonsCalibration(PMTBackgroundphotonsCalibration const&) = delete;
  PMTBackgroundphotonsCalibration(PMTBackgroundphotonsCalibration&&) = delete;
  PMTBackgroundphotonsCalibration& operator=(PMTBackgroundphotonsCalibration const&) = delete;
  PMTBackgroundphotonsCalibration& operator=(PMTBackgroundphotonsCalibration&&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  void analyze(art::Event const& event) override;

  //void doFitGauss(TH1D* hist, int pmt, std::vector<double>& m_params);

  void clean();

private:

  art::InputTag m_ophit_label;


  TTree *m_ophit_ttree;

  std::map<unsigned int, TH1D*> hintegral;
  std::map<unsigned int, TH1D*> hamplitude;
  std::map<unsigned int, TH1D*> hpulses;

  int m_run;
  int m_subrun;
  int m_event;
  uint32_t m_timestamp;

  double adc_to_mV;
  double adc_to_pC;
  double m_threshold;
  std::vector< unsigned int > m_channel_mask;

  TH1D *hequalization;

  std::vector<unsigned int> *m_channel_id = NULL;
  std::vector<float> *m_tstart = NULL;
  std::vector<float> *m_tmax = NULL;
  std::vector<float> *m_amplitude = NULL;
  std::vector<float> *m_integral = NULL;

  art::ServiceHandle<art::TFileService> tfs;

  std::map<unsigned int, unsigned int> m_pulses_count;


  double echarge = 1.602176634; // In units of 10^-7 pC

};


//------------------------------------------------------------------------------


pmtcalo::PMTBackgroundphotonsCalibration::PMTBackgroundphotonsCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{


   m_ophit_label = pset.get<art::InputTag>("OpHitModule", "ophit");

   m_threshold = pset.get<double>("AmplitudeThreshold");

   adc_to_mV = pset.get<double>("ADCmV");  // Conversion from ADC to mV
   adc_to_pC = pset.get<double>("ADCpC");  // Conversion from ADC to pC 

   m_channel_mask = pset.get< std::vector< unsigned int > >
                            ("ChannelMasks", std::vector< unsigned int >());

  // Add histogram bins and range in form of a triplet (hist low, hist high, binsize)
  // ADD fit Ranges 
      
}


//------------------------------------------------------------------------------


void pmtcalo::PMTBackgroundphotonsCalibration::beginJob()
{

  m_ophit_ttree = tfs->make<TTree>("pulsetree","OpHit TTree");
  
  m_ophit_ttree->Branch("run", &m_run, "run/I" );
  m_ophit_ttree->Branch("subrun", &m_subrun, "subrun/I" );
  m_ophit_ttree->Branch("event", &m_event, "event/I" );
  m_ophit_ttree->Branch("timestamp", &m_timestamp, "timestamp/I" );
  m_ophit_ttree->Branch("channel_id", &m_channel_id);
  m_ophit_ttree->Branch("tstart", &m_tstart );
  m_ophit_ttree->Branch("tmax", &m_tmax );
  m_ophit_ttree->Branch("amplitude", &m_amplitude );
  m_ophit_ttree->Branch("integral", &m_integral );


  auto const geop = lar::providerFrom<geo::Geometry>();
  const unsigned int nPMTs = geop->NOpChannels();

  char histname[100]; 
  char histtitle[100];
  
  for(unsigned int opch=0; opch<nPMTs; ++opch)
  {

    // Create a histogram for any valid PMT 
    if( hasChannel(opch, m_channel_mask) )
      continue;

    sprintf(histname, "hintegral%u", opch);
    sprintf(histtitle, "PMT: %u;Pulse charge [10^7 electrons]", opch);
    hintegral[opch] = tfs->make<TH1D>(histname, histtitle, 200, 0, 8);

    sprintf(histname, "hamplitude%u", opch);
    sprintf(histtitle, "PMT: %u;Pulse amplitude [mV]", opch);
    hamplitude[opch] = tfs->make<TH1D>( histname, histtitle, 200, 0, 120 );

    sprintf(histname, "hpulses%u", opch);
    sprintf(histtitle, "PMT: %u;", opch);
    hpulses[opch] = tfs->make<TH1D>( histname, histtitle, 20, 0, 20 );

    m_pulses_count[opch] = 0;

  }


  sprintf(histname, "hequalization");
  sprintf(histtitle, "Average Pulse charge [10^7 electrons]; Num of PMTs");
  hequalization = tfs->make<TH1D>(histname, histtitle, 80, 0, 2.0);

  
}


//-----------------------------------------------------------------------------


void pmtcalo::PMTBackgroundphotonsCalibration::analyze(art::Event const& event)
{

   // Timestamp of the first event 

   m_run = event.id().run();
   m_subrun = event.id().subRun();
   m_event = event.id().event();
   m_timestamp = event.time().timeHigh(); // We just need precision at the s level

   // First thing we sort the ophit in their respective channels
   art::Handle< std::vector< recob::OpHit > > ophitHandle;
   event.getByLabel(m_ophit_label, ophitHandle);

   
   for( auto const& ophit : (*ophitHandle) )
   {
      
      unsigned int opch = ophit.OpChannel();

      if( hasChannel(opch, m_channel_mask) )
        continue;

      m_channel_id->push_back( opch );
      m_tstart->push_back( ophit.PeakTimeAbs() );
      m_tmax->push_back( ophit.PeakTime() ); 
      m_integral->push_back( ophit.Area() );
      m_amplitude->push_back( ophit.Amplitude() );

      hintegral[opch]->Fill( ophit.Area()*adc_to_pC/echarge );
      hamplitude[opch]->Fill( ophit.Amplitude()*adc_to_mV );

      if( ophit.Amplitude()*adc_to_mV >= m_threshold )
        m_pulses_count[opch]++;

   }

   for( auto counts : m_pulses_count )
   {
    unsigned int opch = counts.first;
    hpulses[opch]->Fill( counts.second );
   }

   m_ophit_ttree->Fill();
   
   clean();
  
} // end analyze


//-----------------------------------------------------------------------------


void pmtcalo::PMTBackgroundphotonsCalibration::clean()
{

  m_channel_id->clear();
  m_tstart->clear();
  m_tmax->clear();
  m_amplitude->clear();
  m_integral->clear();

  m_pulses_count.clear();

}


//-----------------------------------------------------------------------------

/*
void pmtcalo::PMTBackgroundphotonsCalibration::doFitGauss(TH1D* hist, int pmt, std::vector<double>& m_params) 
{

  int npars = 3;

  
  int ncounts = 5;
  int rangeMin = -5;
  int rangeMax = 5;

  double width = hist->GetBinWidth(0);

  int startBin = floor(0.1/width);
  int endBin = floor(3.0/width);

  int binValley=startBin;

  double min=99.; 

  int count=0;
  for( int bin=startBin; bin<endBin; bin++  ) {
    double val = hist->GetBinContent(bin);
    if( val < min ) {
      min = val;
      binValley = bin;
    }

    if(val > min){
      count++;
    }

    if( count >= ncounts ){
      break;
    }

  } 

  int binPeak=binValley;
  double max=0.0;
  for( int bin=binValley; bin<endBin; bin++  ) {
    double val = hist->GetBinContent(bin);
    if( val > max ) {
      max = val;
      binPeak = bin;
    }
  } 

  double pv = float(min)/float(max);

  double startFit = (binPeak+rangeMin)*width;
  double endFit = (binPeak+rangeMax)*width;
  

  Gaussian function_obj;
  TF1* func = new TF1("gaussian", function_obj, 0.0, 8.0, npars);
  func->SetParameters( 100, 0.7, 0.3 );
  func->SetParLimits( 0, 0.0, 1e5 );
  func->SetParLimits( 1, 0.1, 2.0 );
  func->SetParLimits( 2, 0.05, 1.0 );

  int fitstatus = hist->Fit("gaussian", "RNQ+", "", 0.2, 1.0);
  
  m_params.push_back(pmt);
  //m_params.push_back(binValley*width);
  //m_params.push_back(binPeak*width);
  //m_params.push_back(pv);
  m_params.push_back(0);
  m_params.push_back(0);
  m_params.push_back(0);
  for(int num=0; num<npars; num++) {
    m_params.push_back( func->GetParameter(num) );
    m_params.push_back( func->GetParError(num) );
  }

  m_params.push_back(func->GetChisquare());
  m_params.push_back(func->GetNDF());
  m_params.push_back(fitstatus);
  

}
*/


//=============================================================================
//=============================================================================
//=============================================================================


void pmtcalo::PMTBackgroundphotonsCalibration::endJob()
{

  /*
  std::cout << " END OF THE JOB!" << std::endl;
  std::cout << m_timestamp << std::endl;

  std::string line = "pmt,min,max,pv,amp,eamp,mu,emu,sigma,esigma,chi2,ndf,fitstatus\n";
  std::cout << line ;

  for(unsigned int opch=0; opch<hintegral.size(); opch++ )
  {

    if( hasChannel(opch, m_channel_mask) )
      continue;

    
    if( hintegral[opch]->GetEntries() > 0 )
    {

      std::vector<double> m_params;
      doFitGauss(hintegral[opch], opch, m_params);

      hequalization->Fill( m_params[6] );

      for( auto par : m_params )
        std::cout << par << ",";
      std::cout << "\n";

    }

  }

  std::cout << " ALL DONE! " << std::endl;
  */

}


DEFINE_ART_MODULE(pmtcalo::PMTBackgroundphotonsCalibration)
