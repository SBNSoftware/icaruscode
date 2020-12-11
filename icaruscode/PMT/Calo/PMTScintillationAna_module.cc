////////////////////////////////////////////////////////////////////////
// Class:       PMTScintillationAna
// Plugin Type: analyzer (art v3_05_00)
// File:        PMTScintillationAna_module.cc
//
// Generated at Mon Sep 21 15:21:37 2020 by Andrea Scarpelli
//
// This analyzer 
//
// mailto:ascarpell@bnl.gov
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
#include "canvas/Persistency/Common/FindManyP.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "icaruscode/PMT/Calo/CaloTools/Waveform.h"
#include "icaruscode/PMT/Calo/CaloTools/Utils.h"

#include "TTree.h"

namespace pmtcalo {
  class PMTScintillationAna;
}


class pmtcalo::PMTScintillationAna : public art::EDAnalyzer {

public:

  explicit PMTScintillationAna(fhicl::ParameterSet const& pset);

  PMTScintillationAna(PMTScintillationAna const&) = delete;
  PMTScintillationAna(PMTScintillationAna&&) = delete;
  PMTScintillationAna& operator=(PMTScintillationAna const&) = delete;
  PMTScintillationAna& operator=(PMTScintillationAna&&) = delete;

  virtual void beginJob() override;

  void cropWaveform( std::vector<double> &tmp_waveform, std::vector<double> & tmp_time, int ch, int binstart );
  void analyze(art::Event const& event) override;

  void clean();


private:

  std::string m_decode_label;
  std::vector<std::string> m_ophit_labels;
  std::vector<std::string> m_flash_labels;
  std::vector<int> m_window;
  std::vector<int> m_pad;
  bool m_dump_waveform;
  fhicl::ParameterSet m_waveform_config;

  TTree* m_scinttree;
  TTree *m_geo_ttree;

  int m_run;
  int m_subrun;
  int m_event;
  int m_wall_id;
  int m_channel_id;
  float m_tstart;
  float m_amplitude;
  float m_integral;
  float m_early_light;
  float m_late_light;
  float m_total_charge;

  std::vector<double>* m_wfrm_crop=NULL;
  std::vector<double>* m_time_crop=NULL;

  art::ServiceHandle<art::TFileService> tfs;

  Waveform *myWaveformAna;

  std::map<int, Waveform::Waveform_t > m_wfrm_map;

};


//------------------------------------------------------------------------------


pmtcalo::PMTScintillationAna::PMTScintillationAna(fhicl::ParameterSet const& pset)
  : art::EDAnalyzer(pset)  // ,
{
   m_decode_label = pset.get<std::string>("DecoderLabel");
   m_ophit_labels = pset.get<std::vector<std::string>>("OpHitModulesList");
   m_flash_labels = pset.get<std::vector<std::string>>("OpFlashModulesList");
   m_window = pset.get<std::vector<int>>("ScintillationWindow");
   m_pad = pset.get<std::vector<int>>("PadWindow");
   m_dump_waveform = pset.get<bool>("DumpWaveform", false);
   m_waveform_config = pset.get<fhicl::ParameterSet>("WaveformAnalysis");


   myWaveformAna = new Waveform(m_waveform_config);


}


//------------------------------------------------------------------------------


void pmtcalo::PMTScintillationAna::beginJob()
{

  m_scinttree = tfs->make<TTree>("scinttree", "scintillation waveform ttree"); 

  m_scinttree->Branch("run", &m_run);
  m_scinttree->Branch("subrun", &m_subrun);
  m_scinttree->Branch("event", &m_event);
  m_scinttree->Branch("wall_id", &m_wall_id);
  m_scinttree->Branch("channel_id", &m_channel_id);
  m_scinttree->Branch("tstart", &m_tstart);
  m_scinttree->Branch("amplitude", &m_amplitude);
  m_scinttree->Branch("integral", &m_integral);
  m_scinttree->Branch("early_light", &m_early_light);
  m_scinttree->Branch("late_light", &m_late_light);
  m_scinttree->Branch("total_charge", &m_total_charge);

  if( m_dump_waveform ) {
    m_scinttree->Branch( "waveform", &m_wfrm_crop );
    m_scinttree->Branch( "time", &m_time_crop );        
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

//----------------------------------------------------------------------------

void pmtcalo::PMTScintillationAna::cropWaveform( std::vector<double> &tmp_waveform, std::vector<double> & tmp_time, int ch, int binstart ) {

 
  double time = m_window[0]*2.0;

  for( int bin=binstart+m_window[0]; bin<binstart+m_window[1]; bin++ ) {
    tmp_waveform.push_back( m_wfrm_map[ch][bin] );
    tmp_time.push_back( time );   
    time += 2.0;
  }

  for( size_t i=0; i<tmp_time.size(); i++ ) {
    if( tmp_time[i] > 90 ){
      m_late_light += tmp_waveform[i];
    }
    else {
      m_early_light += tmp_waveform[i];
    }
    m_total_charge += tmp_waveform[i];
  }

}


//-----------------------------------------------------------------------------


void pmtcalo::PMTScintillationAna::analyze(art::Event const& event)
{

  m_run = event.id().run();
  m_subrun = event.id().subRun();
  m_event = event.id().event();

  art::Handle< std::vector< raw::OpDetWaveform > > rawHandle;
  event.getByLabel(m_decode_label, rawHandle);

  // There is a valid handle per channel
  for( auto const& raw_waveform : (*rawHandle) ) { 

    int channel_id = raw_waveform.ChannelNumber();
    myWaveformAna->loadData( raw_waveform );
    m_wfrm_map[channel_id] = myWaveformAna->getWaveform();
    myWaveformAna->clean();

  }

  // Loop over the flash labels
  for( size_t flash_label_id=0; flash_label_id<m_flash_labels.size(); flash_label_id++ ) {

    m_wall_id = flash_label_id; // for now the label ID is also my wall Id
    auto const & flash_label = m_flash_labels[flash_label_id];
    auto const flashes = event.getValidHandle<std::vector<recob::OpFlash>>(flash_label);
    const art::FindManyP<recob::OpHit> findOpHits( flashes, event, flash_label );

    for(size_t flash_id=0; flash_id < flashes->size(); flash_id++) { 

      // List of the opHit found in the flash
      std::vector<art::Ptr<recob::OpHit>> ophits = findOpHits.at(flash_id);

      for( auto const & ophit : ophits ) {

        m_channel_id =  ophit->OpChannel();
        m_tstart = ophit->OpChannel();
        m_tstart = ophit->PeakTimeAbs()*1000; // in ns
        m_integral = ophit->Area();
        m_amplitude = ophit->Amplitude();

        int bin_start = int( m_tstart/2.0 ); // TODO: use detector services

        if( bin_start < m_pad[0] || bin_start > m_pad[1] ) {
          //std::cout << " Start bin: " << bin_start << 
          //      " outside allowed region [ " << m_pad[0] << ", " << m_pad[1] << " ]" << std::endl;
          continue;
        }

        // Get the cropped waveform around the ROI
        std::vector<double> tmp_waveform;
        std::vector<double> tmp_time;

        m_total_charge = 0; m_early_light = 0; m_late_light = 0;

        cropWaveform( tmp_waveform, tmp_time, m_channel_id, bin_start );

        if ( m_dump_waveform ) {
          m_wfrm_crop = &tmp_waveform;
          m_time_crop = &tmp_time; 
        }

        m_scinttree->Fill();

      }
    }
  }

  m_wfrm_map.clear();
  
} // end analyze


//-----------------------------------------------------------------------------


void pmtcalo::PMTScintillationAna::clean() {



   if( m_wfrm_crop )
      m_wfrm_crop->clear();
   if( m_time_crop )
      m_time_crop->clear();

}


DEFINE_ART_MODULE(pmtcalo::PMTScintillationAna)
