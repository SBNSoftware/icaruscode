/**
 * @file   CRTCalibrationAnalysis_module.cc
 * @brief  Access CRT data and reco products and compare to MCTruth info 
 * @author Chris Hilgenberg (Chris.Hilgenberg@colostate.edu)
 * 
 * The last revision of this code was done in October 2018 with LArSoft v07_06_01.
 */

// LArSoft includes
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TROOT.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <iostream>
#include <utility>
#include <array>

// CRT data products
#include "sbnobj/ICARUS/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;
using std::to_string;

namespace icarus {
namespace crt {

  class CRTCalibrationAnalysis : public art::EDAnalyzer
  {
  public:

    struct Config {
      
      // Save some typing:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> CRTDAQLabel {
        Name("CRTDAQLabel"),
        Comment("tag of the input data product with calibrated CRT data")
        };
   

    }; // Config
    
    using Parameters = art::EDAnalyzer::Table<Config>;
    
  
    /// Constructor: configures the module (see the Config structure above)
    explicit CRTCalibrationAnalysis(Parameters const& config);

    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;


  private:

    art::ServiceHandle<art::TFileService> tfs;

    map<uint8_t,vector<TH1F*>*> macToHistos;

    // The parameters we'll read from the .fcl file.
    art::InputTag fCRTDAQProducerLabel;

    // Other variables that will be shared between different methods.
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event
    CRTCommonUtils* fCrtutils;  
  }; // class CRTCalibrationAnalysis


  //-----------------------------------------------------------------------
   
  CRTCalibrationAnalysis::CRTCalibrationAnalysis(Parameters const& config)
    : EDAnalyzer(config)
    , fCRTDAQProducerLabel(config().CRTDAQLabel())
    , fCrtutils(new CRTCommonUtils())
  {
    // Get a pointer to the geometry service provider.
    fGeometryService = lar::providerFrom<geo::Geometry>();
    // The same for detector TDC clock services.
    // Access to detector properties.
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fTriggerOffset = trigger_offset(clockData);

    for(int i=1; i<94; i++){
      
      macToHistos[i] = new vector<TH1F*>();
      
      for(int ch=0; ch<32; ch++){
	
	string hname = "hadc_"+to_string(i)+"_"+to_string(ch);
	string htitle = "raw charge: mac5 "+to_string(i)+", ch. "+to_string(ch);
	macToHistos[i]->push_back(tfs->make<TH1F>(hname.c_str(),htitle.c_str(),4100,0,4100));
      }
    }
  }
  
  //-----------------------------------------------------------------------
  void CRTCalibrationAnalysis::beginJob()
  {
  }
   
  void CRTCalibrationAnalysis::beginRun(const art::Run& /*run*/)
  {
  }

  //-----------------------------------------------------------------------
  void CRTCalibrationAnalysis::analyze(const art::Event& event) 
  {
    MF_LOG_DEBUG("CRTCalibrationAnalysis") << "beginning analyis" << '\n';


    art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle;
    bool isCRTDAQ = event.getByLabel(fCRTDAQProducerLabel, crtDAQHandle);

    if (isCRTDAQ)  {
      MF_LOG_DEBUG("CRTCalibrationAnalysis") << "about to loop over CRTDAQ entries" << '\n';

     for ( auto const& febdat : (*crtDAQHandle) ) {

        for(int ch=0; ch<32; ch++) {
	  macToHistos[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
        } 

     } //for CRT FEB events
     
    }//if crtdetsim products present

    else 
      mf::LogError("CRTCalibrationAnalysis") << "CRTDAQ products not found!" << std::endl;
    
  } // CRTCalibrationAnalysis::analyze()
  
  DEFINE_ART_MODULE(CRTCalibrationAnalysis)
} // namespace crt
} // namespace icarus

