////////////////////////////////////////////////////////////////////////
// Class:       PMTLaserCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTLaserCalibration_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
// 
//  Data prep for both charge and time calibration using laser pulses
//
//  mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////


#include "larcore/CoreUtils/ServiceUtil.h" 

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

#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "icaruscode/IcarusObj/PMTWaveformTimeCorrection.h"
#include "icaruscode/PMT/Calibration/CaloTools/LaserPulse.h"


#include "TTree.h"

namespace pmtcalo {
  class PMTLaserCalibration;
}


class pmtcalo::PMTLaserCalibration : public art::EDAnalyzer {

public:

  
  explicit PMTLaserCalibration(fhicl::ParameterSet const& pset);

  //PMTLaserCalibration(PMTLaserCalibration const&) = delete;
  //PMTLaserCalibration(PMTLaserCalibration&&) = delete;
  //PMTLaserCalibration& operator=(PMTLaserCalibration const&) = delete;
  //PMTLaserCalibration& operator=(PMTLaserCalibration&&) = delete;

  virtual void beginJob() override;

  void analyze(art::Event const& event) override;

private:

  art::InputTag fOpDetWaveformLabel;

  art::InputTag fTimingCorrectionInputLabel;

  //bool m_filter_noise;

  unsigned int fLaserChannel = 0U;

  fhicl::ParameterSet m_waveform_config;

  bool fDebugMessages;

  //std::size_t fExpectedReadoutSize = 5000;

  icarusDB::IICARUSChannelMap const& fChannelMap;

  detinfo::DetectorClocksData const fClocksData;

  art::ServiceHandle<art::TFileService> tfs;

  std::vector<unsigned int> fIlluminatedChannels;

  LaserPulse *myWaveformAna;

  TTree *m_pulse_ttree;

  int m_run;
  
  int m_event;

  std::vector<int>   *m_channel_id = NULL;
  std::vector<int>   *m_fragment=NULL; 
  //std::vector<double> *m_reftime = NULL;
  std::vector<double> *m_baseline = NULL;
  std::vector<double> *m_rms = NULL;
  std::vector<double> *m_peak_time = NULL;
  std::vector<double> *m_amplitude = NULL;
  std::vector<double> *m_integral = NULL;
  std::vector<double> *m_total_charge = NULL;

  // fitted quantities
  std::vector<double> *m_fit_start_time = NULL;
  std::vector<double> *m_corr_start_time = NULL;
  std::vector<double> *m_error_start_time = NULL;
  std::vector<double> *m_fit_sigma = NULL;
  std::vector<double> *m_error_sigma = NULL;
  std::vector<double> *m_fit_mu = NULL;
  std::vector<double> *m_error_mu = NULL;
  std::vector<double> *m_fit_amplitude = NULL;
  std::vector<double> *m_error_amplitude = NULL;
  std::vector<double> *m_chi2 = NULL;
  std::vector<double> *m_ndf = NULL;
  std::vector<double> *m_fitstatus = NULL; // O:good, >0: bad,  < 0: not working

  bool isIlluminated( unsigned int channel );

  void clean();

};


//------------------------------------------------------------------------------


pmtcalo::PMTLaserCalibration::PMTLaserCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
  , fOpDetWaveformLabel{ pset.get<art::InputTag>("OpDetWaveformLabel", "daqPMTLaser") }
  , fTimingCorrectionInputLabel
    { pset.get<art::InputTag>("PMTWaveformTimingCorrectionLabel", "daqPMT:trgprim") }
  , fLaserChannel{ pset.get<unsigned int>("LaserChannel", 1) }
  , m_waveform_config{ pset.get<fhicl::ParameterSet>("WaveformAnalysis") }
  , fDebugMessages{ pset.get<bool>("DebugMessage", false) }
  , fChannelMap{ *(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}) }
  , fClocksData
    { art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob() }
{ 


  // If the time corrections are already applied at decoder level
  //  fTimingCorrectionInputLabel should remain empty.
  if( fTimingCorrectionInputLabel.empty() ){
    mf::LogError("PMTLaserCalibration")
      << "Timing correction InputLabel not set: timing correction will not be applied"; 
  }


  myWaveformAna = new LaserPulse(m_waveform_config);

}


//------------------------------------------------------------------------------


void pmtcalo::PMTLaserCalibration::beginJob()
{

  //For direct light calibration and timing
  m_pulse_ttree = tfs->make<TTree>("pulsetree","tree with laser pulse characterization");

  m_pulse_ttree->Branch("run", &m_run, "run/I" );
  m_pulse_ttree->Branch("event", &m_event, "event/I" );
  m_pulse_ttree->Branch("fragment", &m_fragment);
  m_pulse_ttree->Branch("channel_id", &m_channel_id );
  m_pulse_ttree->Branch("baseline", &m_baseline );
  m_pulse_ttree->Branch("rms", &m_rms );
  m_pulse_ttree->Branch("peak_time", &m_peak_time );
  m_pulse_ttree->Branch("amplitude", &m_amplitude );
  m_pulse_ttree->Branch("integral", &m_integral );
  m_pulse_ttree->Branch("total_charge", &m_total_charge );

  m_pulse_ttree->Branch("fit_start_time", &m_fit_start_time );
  m_pulse_ttree->Branch("corr_start_time", &m_corr_start_time);
  m_pulse_ttree->Branch("error_start_time", &m_error_start_time );
  m_pulse_ttree->Branch("fit_sigma", &m_fit_sigma);
  m_pulse_ttree->Branch("error_sigma", &m_error_sigma);
  m_pulse_ttree->Branch("fit_mu", &m_fit_mu);
  m_pulse_ttree->Branch("error_mu", &m_error_mu);
  m_pulse_ttree->Branch("fit_amplitude", &m_fit_amplitude);
  m_pulse_ttree->Branch("error_amplitude", &m_error_amplitude);
  m_pulse_ttree->Branch("chi2", &m_chi2);
  m_pulse_ttree->Branch("ndf", &m_ndf);
  m_pulse_ttree->Branch("fitstatus", &m_fitstatus);


  for( std::size_t fragId=0; fragId<24; fragId++ ){
  
    for( icarusDB::PMTChannelInfo_t const & chInfo: fChannelMap.getPMTchannelInfo(fragId) ){
  
      if( chInfo.laserChannelNo == fLaserChannel ){ fIlluminatedChannels.push_back( chInfo.channelID ); };

    }
  
  }

  // Not good if the number of illuminated channels is different than 10 
  assert( fIlluminatedChannels.size() != 10);

}

//-----------------------------------------------------------------------------


bool pmtcalo::PMTLaserCalibration::isIlluminated( unsigned int channel ) {

  auto found = std::find( fIlluminatedChannels.begin(), fIlluminatedChannels.end(), channel );
  return ( found != fIlluminatedChannels.end() );

}



//-----------------------------------------------------------------------------


void pmtcalo::PMTLaserCalibration::analyze(art::Event const& event)
{ 

  m_run = event.id().run();
  
  m_event = event.id().event();

  std::vector<double> startTimeCorrection{360, 0};
  if( !fTimingCorrectionInputLabel.empty() ){

    art::Handle<std::vector<icarus::timing::PMTWaveformTimeCorrection>> pmtTimeCorrectionHandle;
    event.getByLabel(fTimingCorrectionInputLabel , pmtTimeCorrectionHandle);

    if( !pmtTimeCorrectionHandle.isValid() || pmtTimeCorrectionHandle->empty() ) {
      
      for( size_t channelId=0; channelId<(*pmtTimeCorrectionHandle).size(); channelId++ ){

        auto waveformCorr = (*pmtTimeCorrectionHandle)[channelId];
        startTimeCorrection[channelId] = waveformCorr.startTime;

      }

    } else {
      mf::LogError("PMTLaserCalibration")
        << "Invalid input label for product std::vector<icarus::timing::PMTWaveformTimeCorrection>";
      throw; 
    }
     

  } 


  art::Handle< std::vector< raw::OpDetWaveform > > rawWaveformHandle;
  
  if( event.getByLabel(fOpDetWaveformLabel, rawWaveformHandle) ) {


    for( auto const& raw_waveform : (*rawWaveformHandle) ) {

      raw::Channel_t channelId = raw_waveform.ChannelNumber();

      // We are interesed only in the illuminated channels 
      if( !isIlluminated(channelId) ){ continue; }

      m_channel_id->push_back( channelId );
     

      myWaveformAna->loadData( raw_waveform );
      auto pulse = myWaveformAna->getLaserPulse();
      
      // Mostly here we fill up our TTrees
      m_peak_time->push_back( pulse.time_peak );
      
      m_amplitude->push_back( pulse.amplitude );
      
      m_integral->push_back( pulse.integral );
      
      m_total_charge->push_back( myWaveformAna->getTotalCharge() );

      // NB sampling period should be taken from services
      double laser_time = raw_waveform.TimeStamp() + pulse.fit_start_time/1000.; 

      double corr_laser_time = laser_time + startTimeCorrection[channelId];
      
      m_fit_start_time->push_back(laser_time);
      
      // We want the laser position with respect to the trigger time 
      m_corr_start_time->push_back( corr_laser_time );
      
      if ( pulse.amplitude > 100 && fDebugMessages ){

        std::cout << channelId
                  << ", " << pulse.fit_start_time  
                  << ", " << raw_waveform.TimeStamp() 
                  << ", " << fClocksData.TriggerTime() 
                  << ", " << corr_laser_time*1000 
                  << std::endl;

      }

      m_error_start_time->push_back(pulse.error_start_time);

      m_fit_sigma->push_back(pulse.fit_sigma);

      m_error_sigma->push_back(pulse.error_sigma);

      m_fit_mu->push_back(pulse.fit_mu);

      m_error_mu->push_back(pulse.error_mu);

      m_fit_amplitude->push_back(pulse.fit_amplitude);

      m_error_amplitude->push_back(pulse.error_amplitude);

      m_chi2->push_back(pulse.chi2);

      m_ndf->push_back(pulse.ndf);

      m_fitstatus->push_back(pulse.fitstatus);

      // Prepare for the next event
      myWaveformAna->clean();

    } // end loop over pmt channels

    m_pulse_ttree->Fill();
    

    // Cancel the arrays
    clean();

  }

  else {

    mf::LogError("PMTLaserCalibration") 
      << " No OpDetWaveform information with label " << fOpDetWaveformLabel.label() << "\n";
      throw;

  }

} // end analyze


//-----------------------------------------------------------------------------

void pmtcalo::PMTLaserCalibration::clean(){

  m_channel_id->clear();
  m_fragment->clear();
  m_peak_time->clear();
  m_amplitude->clear();
  m_integral->clear();
  m_total_charge->clear();

  m_fit_start_time->clear();
  m_corr_start_time->clear();
  m_error_start_time->clear();
  m_fit_sigma->clear();
  m_error_sigma->clear();
  m_fit_mu->clear();
  m_error_mu->clear();
  m_fit_amplitude->clear();
  m_error_amplitude->clear();
  m_chi2->clear();
  m_ndf->clear();
  m_fitstatus->clear();

}


DEFINE_ART_MODULE(pmtcalo::PMTLaserCalibration)
