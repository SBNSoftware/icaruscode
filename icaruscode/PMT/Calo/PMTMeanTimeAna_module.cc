////////////////////////////////////////////////////////////////////////
// Class:       PMTMeanTimeAna
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTMeanTimeAna_module.cc
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

#include "TTree.h"

namespace pmtcalo {
  class PMTMeanTimeAna;
}


class pmtcalo::PMTMeanTimeAna : public art::EDAnalyzer {

public:

  explicit PMTMeanTimeAna(fhicl::ParameterSet const& pset);

  PMTMeanTimeAna(PMTMeanTimeAna const&) = delete;
  PMTMeanTimeAna(PMTMeanTimeAna&&) = delete;
  PMTMeanTimeAna& operator=(PMTMeanTimeAna const&) = delete;
  PMTMeanTimeAna& operator=(PMTMeanTimeAna&&) = delete;

  virtual void beginJob() override;

  void analyze(art::Event const& event) override;

  double _findMedian(std::vector<short int> a, size_t n);
  float _getTotalCharge( std::vector<short int> wfrm );
  void clean();

private:

  art::InputTag m_decode_label;
  art::InputTag m_ophit_label;

  fhicl::ParameterSet m_waveform_config;

  TTree *m_ophit_ttree;
  TTree *m_geo_ttree;

  int m_run;
  int m_subrun;
  int m_event;
  
  std::vector<int>   *m_channel_id = NULL;
  std::vector<float> *m_tstart = NULL;
  std::vector<float> *m_tmax = NULL;
  std::vector<float> *m_amplitude = NULL;
  std::vector<float> *m_integral = NULL;
  std::vector<float> *m_total_charge = NULL;

  art::ServiceHandle<art::TFileService> tfs;

};


//------------------------------------------------------------------------------


pmtcalo::PMTMeanTimeAna::PMTMeanTimeAna(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{
   m_decode_label = pset.get<art::InputTag>("DecoderModule", "daqPMT");
   m_ophit_label = pset.get<art::InputTag>("OpHitModule", "ophit");
}


//------------------------------------------------------------------------------


void pmtcalo::PMTMeanTimeAna::beginJob()
{

  m_ophit_ttree = tfs->make<TTree>("ophittree","OpHit TTree");
  
  m_ophit_ttree->Branch("run", &m_run, "run/I" );
  m_ophit_ttree->Branch("subrun", &m_subrun, "subrun/I" );
  m_ophit_ttree->Branch("event", &m_event, "event/I" );
  m_ophit_ttree->Branch("channel_id", &m_channel_id);
  m_ophit_ttree->Branch("tstart", &m_tstart );
  m_ophit_ttree->Branch("tmax", &m_tmax );
  m_ophit_ttree->Branch("amplitude", &m_amplitude );
  m_ophit_ttree->Branch("integral", &m_integral );
  m_ophit_ttree->Branch("total_charge", &m_total_charge );

  //-----------------------------------------------------------------------------

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

  m_geo_ttree->Branch("pmtX",&pmtX);
  m_geo_ttree->Branch("pmtY",&pmtY);
  m_geo_ttree->Branch("pmtZ",&pmtZ);

  m_geo_ttree->Fill();
  
}


//-----------------------------------------------------------------------------


void pmtcalo::PMTMeanTimeAna::analyze(art::Event const& event)
{


  
} // end analyze


//-----------------------------------------------------------------------------

void pmtcalo::PMTMeanTimeAna::clean(){


}


DEFINE_ART_MODULE(pmtcalo::PMTMeanTimeAna)
