////////////////////////////////////////////////////////////////////////
// Class:       ICARUSFlashAssAna
// Plugin Type: analyzer (art v3_06_03)
// File:        ICARUSFlashAssAna_module.cc
//
// Generated at Tue Jun 29 13:43:54 2021 by Andrea Scarpelli using cetskelgen
// from cetlib version v3_11_01.
//
// Module that dumps the assciation between Flashes and OpHit
//
// mailto:ascarpel@bnl.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "art_root_io/TFileService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

#include "TTree.h"

#include <vector>
#include <map>



namespace opana {
  class ICARUSFlashAssAna;
}


class opana::ICARUSFlashAssAna : public art::EDAnalyzer {

  public:

    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Sequence<art::InputTag> OpHitLabels {
        Name("OpHitLabels"),
        Comment("Tags for the recob::OpHit data products")
      };

      fhicl::Sequence<art::InputTag> FlashLabels {
        Name("FlashLabels"),
        Comment("Tags for the recob::Flashe data products")
      };

      fhicl::Atom<float> PEOpHitThreshold {
        Name("PEOpHitThreshold"),
        Comment("Threshold in PE for an OpHit to be considered in the information calculated for a flash")
      };

    };

    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit ICARUSFlashAssAna(Parameters const& config);

    ICARUSFlashAssAna(ICARUSFlashAssAna const&) = delete;
    ICARUSFlashAssAna(ICARUSFlashAssAna&&) = delete;
    ICARUSFlashAssAna& operator=(ICARUSFlashAssAna const&) = delete;
    ICARUSFlashAssAna& operator=(ICARUSFlashAssAna&&) = delete;

    void analyze(art::Event const& e) override;

    void beginJob() override;
    void endJob() override;

    int getSideByChannel( const int channel );

    void processOpHits( std::vector<art::Ptr<recob::OpHit>> const &ophits, 
                        int &multiplicity_left, int &multiplicity_right, 
                        float &sum_pe_left, float &sum_pe_right, 
                        float *xyz  ); 

  private:

    std::vector<art::InputTag> fOpHitLabel;
    std::vector<art::InputTag> fFlashLabels;
    float fPEOpHitThreshold;


    TTree *fEventTree;
    std::vector<TTree*> fOpFlashTrees;

    int m_run;
    int m_event;
    int m_timestamp;
    int m_nflashes;
    int m_nophit;
    int m_flash_id;
    int m_multiplicity;
    int m_multiplicity_left;
    int m_multiplicity_right;
    float m_sum_pe;
    float m_sum_pe_left;
    float m_sum_pe_right;
    float m_flash_time;
    float m_flash_x;
    float m_flash_y;
    float m_flash_z;

    std::vector<float> m_pmt_x;
    std::vector<float> m_pmt_y;
    std::vector<float> m_pmt_z;

};


opana::ICARUSFlashAssAna::ICARUSFlashAssAna(Parameters const& config)
  : EDAnalyzer(config)
  , fOpHitLabel( config().OpHitLabels() )
  , fFlashLabels( config().FlashLabels() )
  , fPEOpHitThreshold( config().PEOpHitThreshold() )
{ }


void opana::ICARUSFlashAssAna::beginJob() {

  art::ServiceHandle<art::TFileService const> tfs;

  TTree* fGeoTree = tfs->make<TTree>("geotree", "geometry information" );
  fGeoTree->Branch("pmt_x",&m_pmt_x);
  fGeoTree->Branch("pmt_y",&m_pmt_y);
  fGeoTree->Branch("pmt_z",&m_pmt_z);
  
  auto const geop = lar::providerFrom<geo::Geometry>();

  double PMTxyz[3];
  for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {

    geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);

    m_pmt_x.push_back(PMTxyz[0]);
    m_pmt_y.push_back(PMTxyz[1]);
    m_pmt_z.push_back(PMTxyz[2]);

  }

  fGeoTree->Fill();

  fEventTree = tfs->make<TTree>("eventstree", "higher level information on the event" );
  fEventTree->Branch("run", &m_run, "run/I");
  fEventTree->Branch("event", &m_event, "event/I");
  fEventTree->Branch("timestamp", &m_timestamp, "timestamp/I");
  fEventTree->Branch("nflashes", &m_nflashes, "nflashes/I");
  fEventTree->Branch("nophits", &m_nophit, "nophits/I");

  /*
  if ( !fOpHitLabels.empty() ) {

    for( auto const & label : fOpHitLabels ) {

        std::string name = label+"_ophittree"
        std::string info = "Three for the recob::OpHit with label "+label;

        auto ttree = tfs->make<TTree>(name.c_str(), info.c_str() );
        ttree->Branch("run", &m_run, "run/I");
        ttree->Branch("event" &m_event, "event/I");
        //ttree->Branch("flash_id", &m_flash_id, "flash_id/I");

        fOpHitTrees.push_back( ttree );
    }
  }
  */

  if ( !fFlashLabels.empty() ) {

    for( auto const & label : fFlashLabels ) {

        std::string name = label.label()+"_flashtree";
        std::string info = "Three for the recob::Flashes with label "+label.label();

        TTree* ttree = tfs->make<TTree>(name.c_str(), info.c_str() );
        ttree->Branch("run", &m_run, "run/I");
        ttree->Branch("event", &m_event, "event/I");
        ttree->Branch("timestamp", &m_timestamp, "timestamp/I");
        ttree->Branch("flash_id", &m_flash_id, "flash_id/I");
        ttree->Branch("multiplicity", &m_multiplicity, "multiplicity/I");
        ttree->Branch("multiplicity_right", &m_multiplicity_right, "multiplicity_right/I" );
        ttree->Branch("multiplicity_left", &m_multiplicity_left, "multiplicity_left/I" );
        ttree->Branch("sum_pe", &m_sum_pe, "sum_pe/F");
        //ttree->Branch("sum_pe_right", &m_sum_pe_right, "sum_pe_right/F");
        //ttree->Branch("sum_pe_left", &m_sum_pe_left, "sum_pe_left/F");
        ttree->Branch("flash_time", &m_flash_time, "flash_time/F");
        //ttree->Branch("flash_x", &m_flash_x, "flash_x/F");
        //ttree->Branch("flash_y", &m_flash_y, "flash_y/F");
        //ttree->Branch("flash_z", &m_flash_z, "flash_z/F");

        fOpFlashTrees.push_back( ttree );
    }
  }
}


int opana::ICARUSFlashAssAna::getSideByChannel( const int channel ) {

  /* 
  Channels are numbered from east to west, from North (cryo side) to South (beam side)
  We look in the opposide direction wrt to the beam direction South->North: 

  - Left is the east wall of each cryostat;

  - Right is the west side of each cryostat;

  - [ 0:89 ] and [180:269] are on the left, 
    the return value of the function is 0;

  - [ 90-179 ] and [ 270:359 ] are on the right,
    the return value of the function is 1;
  */


  int side = channel / 90; // always round down

  return side % 2;
}


void opana::ICARUSFlashAssAna::processOpHits( std::vector<art::Ptr<recob::OpHit>> const &ophits, 
                                              int &multiplicity_left, int &multiplicity_right, 
                                              float &sum_pe_left, float &sum_pe_right, 
                                              float *xyz  ) {


  std::unordered_map<int, float > sumpe_map;

  // We caluclate the total charge clustered in the flash per channel taking part to the flash
  for( auto const ophit : ophits ) {

    if ( ophit->PE() < fPEOpHitThreshold ) { continue; }

    const int channel_id = ophit->OpChannel(); 

    sumpe_map[ channel_id ]+=ophit->PE() ;

    xyz[0] = m_pmt_x[channel_id]*ophit->PE();
    xyz[1] = m_pmt_y[channel_id]*ophit->PE();
    xyz[2] = m_pmt_z[channel_id]*ophit->PE();

  }

  //m_multiplicity_left = sumpe_map[0].size();

  //m_multiplicity_right = sumpe_map[1].size();

  //m_sum_pe_left = std::accumulate( sumpe_map[0].begin(), sumpe_map[0].end(), 0 );
  
  //m_sum_pe_right = std::accumulate( sumpe_map[1].begin(), sumpe_map[1].end(), 0 );  

  for( int i=0; i<3; i++ ){ xyz[i] /= m_sum_pe_left+ m_sum_pe_right; }
  
}


void opana::ICARUSFlashAssAna::endJob() {

}


void opana::ICARUSFlashAssAna::analyze(art::Event const& e) {


  m_run = e.id().run();
  m_event = e.id().event();
  m_timestamp = e.time().timeHigh(); // precision to the second 


  if ( !fFlashLabels.empty() ) {

    for ( size_t i=0; i<fFlashLabels.size(); i++  ) {

      auto const label = fFlashLabels[i];

      art::Handle<std::vector<recob::OpFlash>> flash_handle;
      e.getByLabel( label, flash_handle );

      if( flash_handle.isValid() ) {

        art::FindManyP<recob::OpHit> ophitsPtr( flash_handle, e, label );

        for ( size_t idx=0; idx<flash_handle->size(); idx++ ) {

          m_flash_id = idx;
          auto const & flash = (*flash_handle)[idx];

          m_flash_time = flash.Time();
          m_sum_pe = flash.TotalPE();
            
          auto const & ophits = ophitsPtr.at(idx);

          // Get the multiplicity, the position and the number of PE per Side
          float xyz[3] = {0.0, 0.0, 0.0};
          processOpHits( ophits, 
                         m_multiplicity_left, m_multiplicity_right, 
                         m_sum_pe_left, m_sum_pe_right, xyz );

          //std::cout << "\tflash id: " << idx << ", time: " << m_flash_time;
          //std::cout << ", multiplicity left: " << m_multiplicity_left << ", multiplicity right: " << m_multiplicity_right;
          //std::cout << ", sum pe left: " << m_sum_pe_left << ", sum pe right: " << m_sum_pe_right;
          //std::cout << " coor: [" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << "]";
          //std::cout  << "\n";

          m_multiplicity = m_multiplicity_left+m_multiplicity_right;

          m_flash_x = xyz[0];
          m_flash_y = xyz[1];
          m_flash_z = xyz[2];

          fOpFlashTrees[i]->Fill();

        }

      }

      else {

        mf::LogError("ICARUSFlashAssAna")
           << "Invalid recob::OpFlash with label"+label.label()+"\n";
      }

    } // end for on flash input tags

  }

  else {
    mf::LogError("ICARUSFlashAssAna")
          << "No recob::OpFlash labels selected\n";
  }


  fEventTree->Fill();

}


DEFINE_ART_MODULE(opana::ICARUSFlashAssAna)
