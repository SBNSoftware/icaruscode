////////////////////////////////////////////////////////////////////////
// Class:       PMTCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTCalibration_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
//
//  mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////

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

#include "icaruscode/PMT/Calo/CaloTools/Waveform.h"
#include "icaruscode/PMT/Calo/CaloTools/Utils.h"

#include "TTree.h"

namespace pmtcalo {
  class PMTCalibration;
}


class pmtcalo::PMTCalibration : public art::EDAnalyzer {

public:

  explicit PMTCalibration(fhicl::ParameterSet const& pset);

  PMTCalibration(PMTCalibration const&) = delete;
  PMTCalibration(PMTCalibration&&) = delete;
  PMTCalibration& operator=(PMTCalibration const&) = delete;
  PMTCalibration& operator=(PMTCalibration&&) = delete;

  virtual void respondToOpenInputFile(const art::FileBlock& fb) override;
  virtual void beginJob() override;
  virtual void beginSubRun(const art::SubRun &sr) override;


  void analyze(art::Event const& event) override;

private:

  art::InputTag m_data_label;
  std::string m_channel_dbase;
  bool m_filter_noise;
  fhicl::ParameterSet m_waveform_config;


  int m_run;
  int m_subrun;
  int m_event;

  std::map<std::pair<int, int>, int> m_pmtid_map;
  std::map<int, int> m_optch_map;

  TTree *m_charge_ttree;
  //TTree *m_time_ttree;
  TTree *m_geometry_ttree;

  std::vector<float> *m_channel_id = NULL;
  std::vector<float> *m_baseline = NULL;
  std::vector<float> *m_rms = NULL;
  std::vector<float> *m_peak_time = NULL;
  std::vector<float> *m_amplitude = NULL;
  std::vector<float> *m_integral = NULL;
  std::vector<float> *m_total_charge = NULL;

  art::ServiceHandle<art::TFileService> tfs;

  Waveform *myWaveformAna;

};


//------------------------------------------------------------------------------


pmtcalo::PMTCalibration::PMTCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{
   m_data_label = pset.get<art::InputTag>("InputModule", "daq");
   m_filter_noise_label = pset.get<bool>("FilterNoise", false);
   m_waveform_config = pset.get<fhicl::ParameterSet>("WaveformAnalysis");

   myWaveformAna = new Waveform(m_waveform_config);

}


//------------------------------------------------------------------------------


void pmtcalo::PMTCalibration::beginJob()
{

  //For direct light calibration and timing
  m_charge_ttree = tfs->make<TTree>("chargetree","tree for pmt charge");

  m_charge_ttree->Branch("m_run", &m_run, "m_run/I" );
  m_charge_ttree->Branch("m_subrun", &m_subrun, "m_subrun/I" );
  m_charge_ttree->Branch("m_event", &m_event, "m_event/I" );
  m_charge_ttree->Branch("m_channel_id", &m_channel_id );
  m_charge_ttree->Branch("m_baseline", &m_baseline );
  m_charge_ttree->Branch("m_rms", &m_rms );
  m_charge_ttree->Branch("m_peak_time", &m_peak_time );
  m_charge_ttree->Branch("m_amplitude", &m_amplitude );
  m_charge_ttree->Branch("m_integral", &m_integral );
  m_charge_ttree->Branch("m_total_charge", &m_total_charge );


  m_geo_ttree = tfs->make<TTree>("geotree","tree with detector geo info");

}

//------------------------------------------------------------------------------


void pmtcalo::PMTCalibration::respondToOpenInputFile(const art::FileBlock& fb)
{
  m_filename = fb.fileName();
}

//------------------------------------------------------------------------------


 void pmtcalo::PMTCalibration::beginSubRun(const art::SubRun &sr)
 {

   m_run = sr.id().run();
    m_subrun = pmtcalo::fileProgNumber(m_filename);

  } // end beginSubRun


//-----------------------------------------------------------------------------


void pmtcalo::PMTCalibration::analyze(art::Event const& event)
{

   m_event = event.id().event();

   art::Handle< std::vector< raw::OpDetWaveform > > rawHandle;
   event.getByLabel(m_data_label, rawHandle);

   // There is a valid handle per channel
   for( auto const& raw_waveform : (*rawHandle) )
   {

     // Without the correct mapping, we need to set the association with
     // the digitizer id by hand (in future this will be moved to the decoder)
     m_channel_id.push_back( raw_waveform.ChannelNumber() );

     myWaveformAna->loadData( raw_waveform );
     if( m_filter_noise ){ myWaveformAna->filterNoise(); }

     auto pulse = myWaveformAna->getIntegral();

     m_baseline.push_back( _baseline );
     m_rms.push_back( _rms );
     m_peak_time.push_back( pulse.time_peak );
     m_amplitude.push_back( pulse.amplitude );
     m_integral.push_back( pulse.integral );
     m_total_charge.push_back( myWaveformAna->getTotalCharge() );

     // Prepare for the next event
     myWaveformAna->clean();

   } // end loop over pmt channels

   m_charge_ttree->Fill();

   // Cancel the arrays
   clean();


} // end analyze


//-----------------------------------------------------------------------------

void pmtcalo::PMTCalibration::clean(){

  m_channel_id.clear();
  m_baseline.clear();
  m_rms.clear();
  m_peak_time.clear();
  m_amplitude.clear();
  m_integral.clear();
  m_total_charge.clear();

}


DEFINE_ART_MODULE(pmtcalo::PMTCalibration)
