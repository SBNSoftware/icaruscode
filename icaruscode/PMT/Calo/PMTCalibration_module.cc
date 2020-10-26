////////////////////////////////////////////////////////////////////////
// Class:       PMTCalibration
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTCalibration_module.cc
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

  //virtual void respondToOpenInputFile(const art::FileBlock& fb) override;
  virtual void beginJob() override;
  virtual void beginSubRun(const art::SubRun &sr) override;


  void analyze(art::Event const& event) override;

  void clean();

private:

  art::InputTag m_data_label;
  std::string m_channel_dbase;
  bool m_filter_noise;
  fhicl::ParameterSet m_waveform_config;

  int m_run;
  int m_event;

  TTree *m_pulse_ttree;

  std::vector<float> *m_channel_id = NULL;
  std::vector<float> *m_baseline = NULL;
  std::vector<float> *m_rms = NULL;
  std::vector<float> *m_peak_time = NULL;
  std::vector<float> *m_amplitude = NULL;
  std::vector<float> *m_integral = NULL;
  std::vector<float> *m_total_charge = NULL;
  std::vector<float> *m_fragment_timestamp = NULL;

  // fitted quantities
  std::vector<float> *m_fit_start_time = NULL;
  std::vector<float> *m_error_start_time = NULL;
  std::vector<float> *m_fit_sigma = NULL;
  std::vector<float> *m_error_sigma = NULL;
  std::vector<float> *m_fit_mu = NULL;
  std::vector<float> *m_error_mu = NULL;
  std::vector<float> *m_fit_amplitude = NULL;
  std::vector<float> *m_error_amplitude = NULL;
  std::vector<float> *m_chi2 = NULL;
  std::vector<float> *m_ndf = NULL;
  std::vector<float> *m_fitstatus = NULL; // O:good, >0: bad,  < 0: not working

  art::ServiceHandle<art::TFileService> tfs;

  Waveform *myWaveformAna;

};


//------------------------------------------------------------------------------


pmtcalo::PMTCalibration::PMTCalibration(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{
   m_data_label = pset.get<art::InputTag>("InputModule", "daq");
   m_filter_noise = pset.get<bool>("FilterNoise", false);
   m_waveform_config = pset.get<fhicl::ParameterSet>("WaveformAnalysis");

   myWaveformAna = new Waveform(m_waveform_config);

}


//------------------------------------------------------------------------------


void pmtcalo::PMTCalibration::beginJob()
{

  //For direct light calibration and timing
  m_pulse_ttree = tfs->make<TTree>("pulsetree","tree with laser pulse characterization");

  m_pulse_ttree->Branch("run", &m_run, "run/I" );
  //m_pulse_ttree->Branch("subrun", &m_subrun, "subrun/I" );
  m_pulse_ttree->Branch("event", &m_event, "event/I" );
  m_pulse_ttree->Branch("channel_id", &m_channel_id );
  m_pulse_ttree->Branch("baseline", &m_baseline );
  m_pulse_ttree->Branch("rms", &m_rms );
  m_pulse_ttree->Branch("peak_time", &m_peak_time );
  m_pulse_ttree->Branch("amplitude", &m_amplitude );
  m_pulse_ttree->Branch("integral", &m_integral );
  m_pulse_ttree->Branch("total_charge", &m_total_charge );
  m_pulse_ttree->Branch("m_fragment_timestamp", &m_fragment_timestamp );

  m_pulse_ttree->Branch("fit_start_time", &m_fit_start_time );
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

  /*
  m_geo_ttree = tfs->make<TTree>("geotree","tree with detector geo info");

  std::vector<double> pmtX, pmtY, pmtZ;
  std::vector<double> minX, minY, minZ;
  std::vector<double> maxX, maxY, maxZ;
  auto const geop = lar::providerFrom<geo::Geometry>();
  double PMTxyz[3];
  for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {
    geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);
    pmtX.push_back(PMTxyz[0]);
    pmtY.push_back(PMTxyz[1]);
    pmtZ.push_back(PMTxyz[2]);
  }
  for(auto iter=geop->begin_TPC(); iter!=geop->end_TPC(); ++iter) {
    auto const& tpc = (*iter);
    minX.push_back(tpc.BoundingBox().MinX());
    minY.push_back(tpc.BoundingBox().MinY());
    minZ.push_back(tpc.BoundingBox().MinZ());
    maxX.push_back(tpc.BoundingBox().MaxX());
    maxY.push_back(tpc.BoundingBox().MaxY());
    maxZ.push_back(tpc.BoundingBox().MaxZ());
  }
  m_geo_ttree->Branch("pmtX",&pmtX);
  m_geo_ttree->Branch("pmtY",&pmtY);
  m_geo_ttree->Branch("pmtZ",&pmtZ);
  m_geo_ttree->Branch("minX",&minX);
  m_geo_ttree->Branch("minY",&minY);
  m_geo_ttree->Branch("minZ",&minZ);
  m_geo_ttree->Branch("maxX",&maxX);
  m_geo_ttree->Branch("maxY",&maxY);
  m_geo_ttree->Branch("maxZ",&maxZ);
  m_geo_ttree->Fill();
  */
}

//------------------------------------------------------------------------------

/*
void pmtcalo::PMTCalibration::respondToOpenInputFile(const art::FileBlock& fb)
{
  m_filename = fb.fileName();
}
*/

//------------------------------------------------------------------------------


 void pmtcalo::PMTCalibration::beginSubRun(const art::SubRun &sr)
 {

   m_run = sr.id().run();
   //m_subrun = 0;//pmtcalo::fileProgNumber(m_filename);

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
     m_channel_id->push_back( raw_waveform.ChannelNumber() );
     m_fragment_timestamp->push_back( raw_waveform.TimeStamp() );

     myWaveformAna->loadData( raw_waveform );
     if( m_filter_noise ){ myWaveformAna->filterNoise(); }

     auto pulse = myWaveformAna->getIntegral();

     m_baseline->push_back( myWaveformAna->getBaselineMean() );
     m_rms->push_back( myWaveformAna->getBaselineWidth() );
     m_peak_time->push_back( pulse.time_peak );
     m_amplitude->push_back( pulse.amplitude );
     m_integral->push_back( pulse.integral );
     m_total_charge->push_back( myWaveformAna->getTotalCharge() );

     m_fit_start_time->push_back(pulse.fit_start_time);
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


} // end analyze


//-----------------------------------------------------------------------------

void pmtcalo::PMTCalibration::clean(){

  m_channel_id->clear();
  m_baseline->clear();
  m_rms->clear();
  m_peak_time->clear();
  m_amplitude->clear();
  m_integral->clear();
  m_total_charge->clear();
  m_fragment_timestamp->clear();

  m_fit_start_time->clear();
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


DEFINE_ART_MODULE(pmtcalo::PMTCalibration)
