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

  virtual void respondToOpenInputFile(const art::FileBlock& fb) override;
  virtual void beginJob() override;
  virtual void beginSubRun(const art::SubRun &sr) override;


  void analyze(art::Event const& event) override;

  void clean();

private:

  art::InputTag m_data_label;
  std::string m_channel_dbase;
  bool m_filter_noise;
  fhicl::ParameterSet m_waveform_config;
  std::string m_filename;


  int m_run;
  int m_subrun;
  int m_event;
  int m_channel_id;
  float m_baseline;
  float m_rms;
  float m_total_charge;

  TTree *m_wfrm_ttree;
  TTree *m_geo_ttree;

  std::vector<double> *m_waveforom = NULL;
  std::vector<double> *m_ticks = NULL;

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
  m_wfrm_ttree = tfs->make<TTree>("wfrm","tree for pmt charge");

  m_wfrm_ttree->Branch("run", &m_run, "run/I" );
  m_wfrm_ttree->Branch("subrun", &m_subrun, "subrun/I" );
  m_wfrm_ttree->Branch("event", &m_event, "event/I" );
  m_wfrm_ttree->Branch("channel_id", &m_channel_id, "channel_id/I");
  m_wfrm_ttree->Branch("baseline", &m_baseline, "baseline/F");
  m_wfrm_ttree->Branch("rms", &m_rms, "rms/F" );
  m_wfrm_ttree->Branch("total_charge", &m_total_charge, "total_charge/F" );
  m_wfrm_ttree->Branch("waveforom", &m_waveforom );
  m_wfrm_ttree->Branch("ticks", &m_ticks );

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

     m_channel_id = raw_waveform.ChannelNumber();

     myWaveformAna->loadData( raw_waveform );
     if( m_filter_noise ){ myWaveformAna->filterNoise(); }

     m_baseline = myWaveformAna->getBaselineMean();
     m_rms = myWaveformAna->getBaselineWidth();
     m_total_charge = myWaveformAna->getTotalCharge();

     for( auto _adc : myWaveformAna->getWaveform() ){
       m_waveforom->push_back( _adc );
     }

     m_ticks->resize(m_waveforom->size());
     for(size_t i=0; i<m_waveforom->size(); i++){
       m_ticks->at(i) = i*2.0; // in ns
     }

     m_wfrm_ttree->Fill();


     // Prepare for the next event
     myWaveformAna->clean();
     clean();

   } // end loop over pmt channels

} // end analyze


//-----------------------------------------------------------------------------

void pmtcalo::PMTCalibration::clean(){

  m_waveforom->clear();
  m_ticks->clear();

}


DEFINE_ART_MODULE(pmtcalo::PMTCalibration)
