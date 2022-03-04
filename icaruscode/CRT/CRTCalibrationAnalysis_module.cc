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

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition
  /**
   * @brief Example analyzer
   * 
   * This class extracts information from the generated and reconstructed
   * particles.
   *
   * It produces histograms for the simulated particles in the input file:
   * - PDG ID (flavor) of all particles
   * - momentum of the primary particles selected to have a specific PDG ID
   * - length of the selected particle trajectory
   * 
   * It also produces two ROOT trees.
   *
   * The first ROOT tree contains information on the simulated
   * particles, including "dEdx", a binned histogram of collected
   * charge as function of track range.
   * 
   * Configuration parameters
   * =========================
   * 
   * - *SimulationLabel* (string, default: "largeant"): tag of the input data
   *   product with the detector simulation information (typically an instance
   *   of the LArG4 module)
   *
   */
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
    
    // -------------------------------------------------------------------
    // -------------------------------------------------------------------
    // Standard constructor for an ART module with configuration validation;
    // we don't need a special destructor here.

    /// Constructor: configures the module (see the Config structure above)
    explicit CRTCalibrationAnalysis(Parameters const& config);

    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;
//    virtual void endJob() override;

  private:

    art::ServiceHandle<art::TFileService> tfs;

    TTree* fDAQNtuple;

    map<uint8_t,vector<TH1F*>*> macToHistos;

    void FillFebMap();
 //   TTree *calTree;

    // The parameters we'll read from the .fcl file.
    art::InputTag fCRTDAQProducerLabel;

    static map<int, vector<pair<int,int>>> fFebMap;

    // The n-tuples we'll create.

    // The comment lines with the @ symbols define groups in doxygen. 
    /// @name The variables that will go into both n-tuples.
    /// @{
    int fEvent;        ///< number of the event being processed
    int fRun;          ///< number of the run being processed
    int fSubRun;       ///< number of the sub-run being processed
    /// @}

    /// @name The variables that will go into the CosmicDisplay n-tuple.
    /// @{
    static const int LAR_PROP_DELAY = 1.0/(30.0/1.38); //[ns/cm]

    //CRT data product vars
    //static const int kDetMax = 64;
    int      fDetEvent;
    int      fNChan; ///< number of channels above threshold for this front-end board readout
    int      fEntry; ///< front-end board entry number (reset for each event)
    int      fFEBReg; ///< CRT region for this front-end board
    int      fMac5; ///< Mac5 address for this front-end board
    int      fDetSubSys;
    double   fT0;///< signal time w.r.t. global event time
    double   fT1;///< signal time w.r.t. PPS
    int      fADC[64];///< signal amplitude
    vector<vector<int>> fTrackID;///< track ID(s) of particle that produced the signal
    vector<vector<int>> fDetPDG; /// signal inducing particle(s)' PDG code

    //CRT hit product vars
//    int      fHitEvent;
//    float    fXHit; ///< reconstructed X position of CRT hit (cm)
//    float    fYHit; ///< reconstructed Y position of CRT hit (cm)
//    float    fZHit; ///< reconstructed Z position of CRT hit (cm)
//    float    fXErrHit; ///< stat error of CRT hit reco X (cm)
//    float    fYErrHit; ///< stat error of CRT hit reco Y (cm)
//    float    fZErrHit; ///< stat error of CRT hit reco Z (cm)
//    int32_t    fT0Hit; ///< hit time w.r.t. global event time
//    int32_t    fT1Hit; ///< hit time w.r.t. PPS
//    int       fHitReg; ///< region code of CRT hit
//    int       fHitSubSys;
//    int       fNHit; ///< number of CRT hits for this event
//    int       fHitStrip;
//    int       fHitMod;
//    int       fNHitFeb;
//    float     fHitTotPe;

    
    // Other variables that will be shared between different methods.
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event
    CRTCommonUtils* fCrtutils;  
  }; // class CRTCalibrationAnalysis


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // Constructor
  // 
  // Note that config is a Table<Config>, and to access the Config
  // value we need to use an operator: "config()". In the same way,
  // each element in Config is an Atom<Type>, so to access the type we
  // again use the call operator, e.g. "SimulationLabel()".

  map<int,vector<pair<int,int>>> CRTCalibrationAnalysis::fFebMap;
 
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
//			macToHistos[i]->at(ch)->Write();
		}
	}
  }

  void CRTCalibrationAnalysis::FillFebMap() { 
    if(!this->fFebMap.empty())
        return;
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file("feb_map.txt",fullFileName);
    std::ifstream fin;
    fin.open(fullFileName,std::ios::in);
    if(fin.good()) std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
    else //std::cout << "could not open file 'feb_map.txt' for reading!" << std::endl;
        throw cet::exception("CRTDataAnalysis::FillFebMap") << "Unable to find/open file 'feb_map.txt'" << std::endl;
    std::vector<std::string> row;
    std::string line, word;
    while(getline(fin,line)) {
        row.clear();
        std::stringstream s(line);
        int mod;
        while (std::getline(s, word, ',')) {
            row.push_back(word);
        }
        mod = std::stoi(row[0]);
        (this->fFebMap)[mod].push_back(std::make_pair(std::stoi(row[1]),std::stoi(row[2])));
        if(row.size()>3)
            (this->fFebMap)[mod].push_back(std::make_pair(std::stoi(row[3]),std::stoi(row[4])));
    }
    std::cout << "filled febMap with " << (this->fFebMap).size() << " entries" << std::endl;
    fin.close();
  }

  
  //-----------------------------------------------------------------------
  void CRTCalibrationAnalysis::beginJob()
  {
    std::cout << " starting analysis job" << std::endl;

}
   
  void CRTCalibrationAnalysis::beginRun(const art::Run& /*run*/)
  {
  }

  //-----------------------------------------------------------------------
  void CRTCalibrationAnalysis::analyze(const art::Event& event) 
  {
    MF_LOG_DEBUG("CRT") << "beginning analyis" << '\n';

    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

	std::cout << "event: " << fEvent << ", run " << fRun << ", subrun " << fSubRun << " found." << std::endl;

    FillFebMap();//febMap);

    art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle;
    bool isCRTDAQ = event.getByLabel(fCRTDAQProducerLabel, crtDAQHandle);

    if (isCRTDAQ)  {
     std::cout << "about to loop over CRTDAQ entries" << std::endl;
     for ( auto const& febdat : (*crtDAQHandle) ) {
        fDetEvent       = fEvent;
        fMac5           = febdat.fMac5;
	if(fMac5>93) continue;
	std::cout << "found mac5: " << fMac5 << std::endl;
        fEntry          = febdat.fEntry;
	fFEBReg         = fCrtutils->AuxDetRegionNameToNum(fCrtutils->MacToRegion(fMac5));
        fNChan = 0;
	fDetSubSys = fCrtutils->MacToTypeCode(fMac5);
        fT0 = febdat.fTs0;
        fT1 = febdat.fTs1;

        int maxchan =0;
        if(fDetSubSys!=2) maxchan=32;
        else maxchan = 64;
	std::cout << "maxchan = " << maxchan << std::endl;
        for(int ch=0; ch<maxchan; ch++) {
		std::cout << "channel " << ch << " found!" << std::endl;
                fADC[ch] = febdat.fAdc[ch]; 
		macToHistos[fMac5]->at(ch)->Fill( fADC[ch] );
        } 

//        fDAQNtuple->Fill();

     } //for CRT FEB events

    }//if crtdetsim products present

    else 
	throw cet::exception("CRTSimAnalysis") << "CRTDAQ products not found!" << std::endl;

  } // CRTCalibrationAnalysis::analyze()
 
  DEFINE_ART_MODULE(CRTCalibrationAnalysis)
} // namespace crt
} // namespace icarus

