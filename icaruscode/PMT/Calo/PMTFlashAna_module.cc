////////////////////////////////////////////////////////////////////////
// Class:       PMTFlashAna
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTFlashAna_module.cc
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "TTree.h"

namespace pmtcalo {
  class PMTFlashAna;
}


class pmtcalo::PMTFlashAna : public art::EDAnalyzer {

public:

  explicit PMTFlashAna(fhicl::ParameterSet const& pset);

  PMTFlashAna(PMTFlashAna const&) = delete;
  PMTFlashAna(PMTFlashAna&&) = delete;
  PMTFlashAna& operator=(PMTFlashAna const&) = delete;
  PMTFlashAna& operator=(PMTFlashAna&&) = delete;

  virtual void beginJob() override;

  void analyze(art::Event const& event) override;
  void clean();

private:

  std::vector<std::string> m_decode_labels;
  std::vector<std::string> m_ophit_labels;
  std::vector<std::string> m_flash_labels;

  std::vector<TTree*> m_flashtrees;
  TTree *m_geo_ttree;

  int m_run;
  int m_subrun;
  int m_event;

  double m_time;
  double m_pe_sum;
  std::vector<double> m_pe_v;
  double m_y;
  double m_z;
  //double m_nphotons;
  
  //std::vector<float> *m_tstart = NULL;
  //std::vector<float> *m_tmax = NULL;
  //std::vector<float> *m_amplitude = NULL;
  //std::vector<float> *m_integral = NULL;
  //std::vector<float> *m_total_charge = NULL;

  art::ServiceHandle<art::TFileService> tfs;

};


//------------------------------------------------------------------------------


pmtcalo::PMTFlashAna::PMTFlashAna(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{
   m_decode_labels = pset.get<std::vector<std::string>>("DecoderModulesList");
   m_ophit_labels = pset.get<std::vector<std::string>>("OpHitModulesList");
   m_flash_labels = pset.get<std::vector<std::string>>("OpFlashModulesList");
}


//------------------------------------------------------------------------------


void pmtcalo::PMTFlashAna::beginJob()
{

for(auto const& label : m_flash_labels) {
 
    std::string name = label + "ttree"; 

    auto flashtree = tfs->make<TTree>(name.c_str(),name.c_str()); 

    flashtree->Branch("run",&m_run,"run/I");
    flashtree->Branch("event",&m_event,"event/I");
    flashtree->Branch("time",&m_time,"time/D");
    flashtree->Branch("pe_v",&m_pe_v);
    flashtree->Branch("pe_sum",&m_pe_sum,"pe_sum/D");
    flashtree->Branch("y",&m_y,"y/D");
    flashtree->Branch("z",&m_z,"z/D");
    
    m_flashtrees.push_back(flashtree);

  }

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


void pmtcalo::PMTFlashAna::analyze(art::Event const& event)
{

   m_run = event.id().run();
   m_subrun = event.id().subRun();
   m_event = event.id().event();

   // Loop over the flash labels
  for( size_t flash_label_id=0; flash_label_id<m_flash_labels.size(); flash_label_id++ ) {

    auto const & flash_label = m_flash_labels[flash_label_id];
    auto const & flashtree = m_flashtrees[flash_label_id];

    auto const flashes = event.getValidHandle<std::vector<recob::OpFlash>>(flash_label);
    const art::FindManyP<recob::OpHit> findOpHits( flashes, event, flash_label );

    for(size_t flash_id=0; flash_id < flashes->size(); flash_id++) { 

      auto const & flash = flashes->at(flash_id);

      m_time = flash.Time();
      //m_time_width = flash.TimeWidth();
      m_pe_v = flash.PEs();
      m_pe_sum = 0.0; //std::accumulate(m_pe_v.begin(),m_pe_v.end());
      m_y = flash.YCenter();
      m_z = flash.ZCenter();


    std::vector<art::Ptr<recob::OpHit>> ophits = findOpHits.at(flash_id);

    std::cout << "OpHit: " << ophits.size() << std::endl;

    flashtree->Fill();
     
    }

  } // end loop over labels
   
} // end analyze


//-----------------------------------------------------------------------------

void pmtcalo::PMTFlashAna::clean(){


}


DEFINE_ART_MODULE(pmtcalo::PMTFlashAna)
