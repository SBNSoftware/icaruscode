/**
 * @file   CRTAnalysis_module.cc
 * @brief  Access CRT data and reco products and compare to MCTruth info 
 * @author Chris Hilgenberg (Chris.Hilgenberg@rams.colostate.edu)
 * 
 * The last revision of this code was done in August 2018 with LArSoft v07_01_00.
 */

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "TMath.h"

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
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTHit.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;

namespace {

  uint32_t ModToTypeCode(geo::AuxDetGeo const& adgeo); 
  char ModToAuxDetType(geo::AuxDetGeo const& adgeo);
  int GetAuxDetRegion(geo::AuxDetGeo const& adgeo);
  int ProcessToICode(string const& p);
  uint32_t MacToADReg(uint32_t mac);
  //char MacToType(uint32_t mac);
  uint32_t MacToTypeCode(uint32_t mac);
  uint32_t ADToMac(char type, uint32_t adid);
  //uint32_t ADSToFEBChannel(char type, uint32_t adid, uint32_t adsid);

} // local namespace


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
  class CRTAnalysis : public art::EDAnalyzer
  {
  public:

    struct Config {
      
      // Save some typing:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      // One Atom for each parameter
      fhicl::Atom<art::InputTag> SimulationLabel {
        Name("SimulationLabel"),
        Comment("tag of the input data product with the detector simulation information")
        };

      fhicl::Atom<art::InputTag> CRTSimHitLabel {
        Name("CRTSimHitLabel"),
        Comment("tag of the input data product with reconstructed CRT hits")
        };

      fhicl::Atom<art::InputTag> CRTDetSimLabel {
        Name("CRTDetSimLabel"),
        Comment("tag of the input data product with simulated CRT hits")
        };
   
      fhicl::Sequence<int> PDGs {
        Name("PDGs"),
        Comment("consider only these PDGs")
        };

      fhicl::Sequence<float> MinMomenta {
        Name("MinMomenta"),
        Comment("minimum momentum for each PDG selected")
        };

      fhicl::Sequence<float> MaxMomenta {
        Name("MaxMomenta"),
        Comment("maximum momentum for each PDG selected")
        };

    }; // Config
    
    using Parameters = art::EDAnalyzer::Table<Config>;
    
    // -------------------------------------------------------------------
    // -------------------------------------------------------------------
    // Standard constructor for an ART module with configuration validation;
    // we don't need a special destructor here.

    /// Constructor: configures the module (see the Config structure above)
    explicit CRTAnalysis(Parameters const& config);

    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;

  private:

    // The parameters we'll read from the .fcl file.
    art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    art::InputTag fCRTSimHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fCRTDetSimProducerLabel;
    vector<int> fPDGs;                       ///< PDG code of particle we'll focus on
    vector<float> fMinMomenta;
    vector<float> fMaxMomenta;

    // Pointers to the histograms we'll create. 
    //TH1D* fPDGCodeHist;     ///< PDG code of all particles
    //TH1D* fMomentumHist;    ///< momentum [GeV] of all selected particles
    //TH1D* fTrackLengthHist; ///< true length [cm] of all selected particles

    // The n-tuples we'll create.
    TTree* fSimulationNtuple;     ///< tuple with simulated data
    TTree* fRegionsNtuple;
    TTree* fDetSimNtuple;
    TTree* fSimHitNtuple;
    TTree* fCRTTruthNtuple; 
    TTree* fCRTTruthMatchNtuple;

    // The comment lines with the @ symbols define groups in doxygen. 
    /// @name The variables that will go into both n-tuples.
    /// @{
    int fEvent;        ///< number of the event being processed
    int fRun;          ///< number of the run being processed
    int fSubRun;       ///< number of the sub-run being processed
    /// @}

    static const int kMaxAD = 100; //for memory allocation 

    /// @name The variables that will go into the simulation n-tuple.
    /// @{
    uint32_t fSimHits; ///< number of trajectory points for each MCParticle
    float    fTrackLength; ///< total track length for each MCParticle
    int      fSimPDG;       ///< PDG ID of the particle being processed
    int      fSimProcess;   ///< process that created the particle (e.g. brehmstralung)
    int      fSimEndProcess; ///< process the killed the particle (e.g. annihilation)
    int      fSimTrackID;   ///< GEANT ID of the particle being processed
    uint32_t fNAuxDet;   ///< Number of scintillator strips hit
    uint32_t fAuxDetID[kMaxAD];  ///< Global CRT module ID
    uint32_t fAuxDetSensitiveID[kMaxAD]; ///< Strip ID in module
    float    fADEDep[kMaxAD]; ///< Energy deposited in CRT strip (GeV)
    float    fADdEdx[kMaxAD]; ///< average dEdx for particle traversing CRT strip
    float    fADTrackLength[kMaxAD]; ///< Track length in CRT strip (cm)
    uint32_t fAuxDetReg[kMaxAD]; ///< CRT region code
    uint32_t fADMac[kMaxAD]; ///< Mac5 address of the CRT module
    uint32_t  fADType[kMaxAD];

    float   fADEnterXYZT[kMaxAD][4]; ///< 4-position of entry into CRT strip
    float   fADExitXYZT[kMaxAD][4]; ///< 4-position of exit from CRT strip

    int fParentPDG;
    float fParentE;

    // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
    float fStartXYZT[4];///< (x,y,z,t) of the true start of the particle
    float fEndXYZT[4];  ///< (x,y,z,t) of the true end of the particle
    float fStartPE[4];  ///< (Px,Py,Pz,E) at the true start of the particle
    float fEndPE[4];    ///< (Px,Py,Pz,E) at the true end of the particle    

    int fProgenitor; ///< G4 track ID of the primary particle that ultimately led to this one
    int fMother; ///< G4 track ID of mother that directly produced this MCParticle
    int fNDaught; ///< number of daughters belonging to this MCParticle

    //Regions tree vars
    static const size_t kMaxReg = 50;
    uint32_t fRegEvent;
    uint32_t fNReg;
    uint32_t fRegFid;
    uint32_t fRegActive;
    uint32_t fRegInactive;
    uint32_t fRegCryo;
    uint32_t fRegCRTs;
    int      fRegRegions[kMaxReg];
    int      fRegPDG[kMaxReg];
    int      fRegTrkID[kMaxReg];
    float    fRegEDep[kMaxReg];
    float    fRegdL[kMaxReg];
    float    fRegEntryPE[kMaxReg][4];
    float    fRegExitPE[kMaxReg][4];
    float    fRegEntryXYZT[kMaxReg][4];
    float    fRegExitXYZT[kMaxReg][4];
    float    fRegEntrySlope[kMaxReg][3];
    float    fRegExitSlope[kMaxReg][3];

    //CRT data product vars
    static const int kDetMax = 64;
    uint32_t fDetEvent;
    uint32_t fChan[kDetMax]; ///< front-end board channel (0-31 or 0-63)
    int      fT0[kDetMax]; ///< signal time w.r.t. global event time
    int      fT1[kDetMax]; ///< signal time w.r.t. PPS
    uint32_t fADC[kDetMax]; ///< signal amplitude
    int      fTrackID[kDetMax]; ///< track ID of particle that produced the signal
    int      fDetPDG[kDetMax];
    uint32_t fNChan; ///< number of channels above threshold for this front-end board readout
    uint32_t fEntry; ///< front-end board entry number (reset for each event)
    int      fTTrig;      ///< signal time w.r.t. global event time of primary channel providing trigger
    uint32_t fChanTrig; ///< channel which provided the trigger for this front-end board
    uint32_t fFEBReg; ///< CRT region for this front-end board
    uint32_t fMac5; ///< Mac5 address for this front-end board
    uint32_t fTriggerPair[2]; ///< two channels which provided the coincidence (useful for C or D modules)
    uint32_t fMacPair[2]; ///< two front-end boards with pairwise coincidence ( useful for M modules)
    uint32_t  fDetSubSys;
      
    //CRT hit product vars
    float    fHitEvent;
    float    fXHit; ///< reconstructed X position of CRT hit (cm)
    float    fYHit; ///< reconstructed Y position of CRT hit (cm)
    float    fZHit; ///< reconstructed Z position of CRT hit (cm)
    float    fXErrHit; ///< stat error of CRT hit reco X (cm)
    float    fYErrHit; ///< stat error of CRT hit reco Y (cm)
    float    fZErrHit; ///< stat error of CRT hit reco Z (cm)
    float    fT0Hit; ///< hit time w.r.t. global event time
    float    fT1Hit; ///< hit time w.r.t. PPS
    uint32_t fHitReg; ///< region code of CRT hit
    uint32_t fNHit; ///< number of CRT hits for this event
    int      fHitTrk;

    //truth matching stats for CRTDetSim
    uint32_t fNmuTruth;    //true N primary muons entering CRT sensive volumes
    uint32_t fNmuTruthCRTTrig;
    uint32_t fNmuTruthCRTTrigVec;
    uint32_t fNmuTruthCRTTrigC;
    uint32_t fNmuTruthCRTTrigM;
    uint32_t fNmuTruthCRTTrigD;
    uint32_t fNmuTruthMissCRT;
    uint32_t fNmuTruthCRTTag;
    uint32_t fNmuTruthCRTTagC;
    uint32_t fNmuTruthCRTTagD;
    uint32_t fNmuTruthCRTTagM;
    uint32_t fNmuTruthCRTTagVec;
    TH1F* fStripMultHistC; ///< true N strips hit / C-module / muon track 
    TH1F* fStripMultHistM; ///< true N strips hit / M-module / muon track
    TH1F* fStripMultHistD; ///< true N strips hit / D-module / muon track
    TH1F* fModMultHistC;   ///< true N C-modules hit / muon track
    TH1F* fModMultHistM;   ///< true N M-modules hit / muon track
    TH1F* fModMultHistD;   ///< true N D-modules hit / muon track

    uint32_t fNmuDetCRTTrig; 
    uint32_t fNmuDetCRTTrigVec;
    uint32_t fNmuDetCRTTrigC;
    uint32_t fNmuDetCRTTrigM;
    uint32_t fNmuDetCRTTrigD;

    // note that the following four variables are not used so are being commented out here
    //uint32_t fNmuTagC;     //N muon tracks producing >0 CRT triggers in C-subsystem
    //uint32_t fNmuTagM;     //N muon tracks producing >0 CRT triggers in M-subsystem
    //uint32_t fNmuTagD;     //N muon tracks producing >0 CRT triggers in D-subsystem
    //uint32_t fNmuTagTot;
    TH1F* fChanMultHistC;  //N FEB channels > threshold / muon track
    TH1F* fChanMultHistM;  
    TH1F* fChanMultHistD;
    TH1F* fFEBMultHistC;   //N FEBs w/trigger / muon track
    TH1F* fFEBMultHistM;
    TH1F* fFEBMultHistD;

    float fEffC;        //tagging efficiency of C-subsystem
    float fEffM;        //tagging efficiency of M-subsystem
    float fEffD;        //tagging efficiency of D-subsystem
    float fEffTot;      //tagging efficiency of D-subsystem

    //truth matching histos for CRTHits
    /*TH1F* fXResHistC;
    TH1F* fYResHistC;
    TH1F* fZResHistC;
    TH1F* fTResHistC;
    TH1F* fXResHistM;
    TH1F* fYResHistM;
    TH1F* fZResHistM;
    TH1F* fTResHistM;
    TH1F* fXResHistD;
    TH1F* fYResHistD;
    TH1F* fZResHistD;
    TH1F* fTResHistD;
*/
    /// @}
    
    // Other variables that will be shared between different methods.
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event
    
  }; // class CRTAnalysis


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
  // 
  CRTAnalysis::CRTAnalysis(Parameters const& config)
    : EDAnalyzer(config)
    , fSimulationProducerLabel(config().SimulationLabel())
    , fCRTSimHitProducerLabel(config().CRTSimHitLabel())
    , fCRTDetSimProducerLabel(config().CRTDetSimLabel())
    , fPDGs(config().PDGs())
    , fMinMomenta(config().MinMomenta())
    , fMaxMomenta(config().MaxMomenta())
  {
    // Get a pointer to the geometry service provider.
    fGeometryService = lar::providerFrom<geo::Geometry>();
    // The same for detector TDC clock services.
    //fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
    // Access to detector properties.
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fTriggerOffset = detprop->TriggerOffset();
  }
  
  //-----------------------------------------------------------------------
  void CRTAnalysis::beginJob()
  {
    std::cout << " starting analysis job" << std::endl;

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Define our n-tuples
    fSimulationNtuple    = tfs->make<TTree>("SimTree",          "MyCRTSimulation");
    fRegionsNtuple       = tfs->make<TTree>("RegTree",          "Info about particles crossing boundaries");
    fDetSimNtuple        = tfs->make<TTree>("DetTree",          "MyCRTDetSim");
    fSimHitNtuple        = tfs->make<TTree>("HitTree",          "MyCRTSimHit");
    fCRTTruthNtuple      = tfs->make<TTree>("CRTTruthTree",     "CRT Truth values for muons");
    fCRTTruthMatchNtuple = tfs->make<TTree>("CRTTruthMatchTree","same as CRTTruth but drom DetSim stage");

    // Construct truth matching histograms
    fStripMultHistC   = tfs->make<TH1F>("StripMultC",";no. strips hit / module / #mu;",64,0,64);
    fStripMultHistM   = tfs->make<TH1F>("StripMultM",";no. strips hit / module / #mu;",64,0,64);
    fStripMultHistD   = tfs->make<TH1F>("StripMultD",";no. strips hit / module / #mu;",64,0,64);
    fModMultHistC     = tfs->make<TH1F>("ModMultC",";no. modules hit / #mu;",10,0,10);
    fModMultHistM     = tfs->make<TH1F>("ModMultM",";no. modules hit / #mu;",10,0,10);
    fModMultHistD     = tfs->make<TH1F>("ModMultD",";no. modules hit / #mu;",10,0,10);

    fChanMultHistC    =tfs->make<TH1F>("ChanMultC",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fChanMultHistM    =tfs->make<TH1F>("ChanMultD",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fChanMultHistD    =tfs->make<TH1F>("ChanMultM",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fFEBMultHistC     =tfs->make<TH1F>("FEBMultC",";no. FEB triggers / #mu;",64,0,64);
    fFEBMultHistM     =tfs->make<TH1F>("FEBMultD",";no. FEB triggers / #mu;",64,0,64);
    fFEBMultHistD     =tfs->make<TH1F>("FEBMultM",";no. FEB triggers / #mu;",64,0,64);


    /*fXResHistC        = tfs->make<TH1F>("XResC",";X#_{true}-X#_{reco} (cm);",200-10000,10000); //C
    fYResHistC        = tfs->make<TH1F>("YResC",";Y#_{true}-Y#_{reco} (cm);",200-10000,10000);
    fZResHistC        = tfs->make<TH1F>("ZResC",";Z#_{true}-Z#_{reco} (cm);",200-10000,10000);
    fTResHistC        = tfs->make<TH1F>("TResC",";T#_{true}-T#_{reco} (ns);",200-10000,10000);

    fXResHistM        = tfs->make<TH1F>("XResM",";X#_{true}-X#_{reco} (cm);",200-10000,10000); //M
    fYResHistM        = tfs->make<TH1F>("YResM",";Y#_{true}-Y#_{reco} (cm);",200-10000,10000);
    fZResHistM        = tfs->make<TH1F>("ZResM",";Z#_{true}-Z#_{reco} (cm);",200-10000,10000);
    fTResHistM        = tfs->make<TH1F>("TResM",";T#_{true}-T#_{reco} (ns);",200-10000,10000);

    fXResHistD        = tfs->make<TH1F>("XResD",";X#_{true}-X#_{reco} (cm);",200-10000,10000); //D
    fYResHistD        = tfs->make<TH1F>("YResD",";Y#_{true}-Y#_{reco} (cm);",200-10000,10000);
    fZResHistD        = tfs->make<TH1F>("ZResD",";Z#_{true}-Z#_{reco} (cm);",200-10000,10000);
    fTResHistD        = tfs->make<TH1F>("TResD",";T#_{true}-T#_{reco} (ns);",200-10000,10000);
*/
    // Define the branches of our simulation n-tuple
    fSimulationNtuple->Branch("event",             &fEvent,             "event/I");
    fSimulationNtuple->Branch("nPoints" ,          &fSimHits,           "nPoints/I");
    fSimulationNtuple->Branch("subRun",            &fSubRun,            "subRun/I");
    fSimulationNtuple->Branch("run",               &fRun,               "run/I");
    fSimulationNtuple->Branch("trackID",           &fSimTrackID,        "trackID/I");
    fSimulationNtuple->Branch("pdg",               &fSimPDG,            "pdg/I");
    fSimulationNtuple->Branch("trackLength",       &fTrackLength,       "trackLenth/F");
    fSimulationNtuple->Branch("process",           &fSimProcess,        "process/I");
    fSimulationNtuple->Branch("endProcess",        &fSimEndProcess,     "endProcess/I");
    fSimulationNtuple->Branch("parentPDG",         &fParentPDG,         "parentPDG/I");
    fSimulationNtuple->Branch("parentE",           &fParentE,           "parentE/F");
    fSimulationNtuple->Branch("progenitor",        &fProgenitor,        "progenitor/I");

    // CRT hits
    fSimulationNtuple->Branch("auxDetSensitiveID", fAuxDetSensitiveID, "auxDetSensitiveID[100]/I");
    fSimulationNtuple->Branch("auxDetID",          fAuxDetID,          "auxDetID[100]/I");
    fSimulationNtuple->Branch("auxDetEDep",        fADEDep,            "auxDetEDep[100]/F");
    fSimulationNtuple->Branch("auxDetdEdx",        fADdEdx,            "auxDetdEdx[100]/F");
    fSimulationNtuple->Branch("auxDetTrackLength", fADTrackLength,     "auxDetTrackLength[100]/F");
    fSimulationNtuple->Branch("auxDetEnterXYZT",   fADEnterXYZT,       "auxDetEnterXYZT[100][4]/F");
    fSimulationNtuple->Branch("auxDetExitXYZT",    fADExitXYZT,        "auxDetExitXYZT[100][4]/F");
    fSimulationNtuple->Branch("auxDetRegion",      fAuxDetReg,         "auxDetRegion[100]/I");
    fSimulationNtuple->Branch("mac5",              fADMac,             "adMac[100]/I");
    fSimulationNtuple->Branch("adType",            fADType,            "adType[100]/I");

    fSimulationNtuple->Branch("startXYZT",         fStartXYZT,          "startXYZT[4]/F");
    fSimulationNtuple->Branch("endXYZT",           fEndXYZT,            "endXYZT[4]/F");
    fSimulationNtuple->Branch("startPE",           fStartPE,            "startPE[4]/F");
    fSimulationNtuple->Branch("endPE",             fEndPE,              "endPE[4]/F");
    fSimulationNtuple->Branch("nChan",             &fNAuxDet,           "nChan/I");
    fSimulationNtuple->Branch("mother",            &fMother,            "mother/I");
    fSimulationNtuple->Branch("nDaught",           &fNDaught,           "nDaught/I");

    //regions tree
    fRegionsNtuple->Branch("event",                &fRegEvent,           "event/I");
    fRegionsNtuple->Branch("nReg",                 &fNReg,               "nReg/I");
    fRegionsNtuple->Branch("fid",                  &fRegFid,             "fid/I");
    fRegionsNtuple->Branch("active",               &fRegActive,          "active/I");
    fRegionsNtuple->Branch("inactive",             &fRegInactive,        "inactive/I");
    fRegionsNtuple->Branch("cryo",                 &fRegCryo,            "cryo/I");
    fRegionsNtuple->Branch("crts",                 &fRegCRTs,            "crts/I");
    fRegionsNtuple->Branch("regions",              fRegRegions,          "regions[50]/I");
    fRegionsNtuple->Branch("pdgs",                 fRegPDG,              "pdgs[50]/I");
    fRegionsNtuple->Branch("trackIDs",             fRegTrkID,            "trackIDs[50]/I");
    fRegionsNtuple->Branch("eDep",                 fRegEDep,             "eDep[50]/F");
    fRegionsNtuple->Branch("dL",                   fRegdL,               "dL[50]/F");
    fRegionsNtuple->Branch("entryPE",              fRegEntryPE,          "entryPE[50][4]/F");
    fRegionsNtuple->Branch("entryPE",              fRegExitPE,           "exitPE[50][4]/F");
    fRegionsNtuple->Branch("entryXYZT",            fRegEntryXYZT,        "entryXYZT[50][4]/F");
    fRegionsNtuple->Branch("exitXYZT",             fRegExitXYZT,         "exitXYZT[50][4]/F");
    fRegionsNtuple->Branch("entrySlope",           fRegEntrySlope,       "entrySlope[50][3]/F");
    fRegionsNtuple->Branch("exitSlope",            fRegExitSlope,        "exitSlope[50][3]/F");

    fCRTTruthNtuple->Branch("Event",             &fEvent,             "Event/I");
    fCRTTruthNtuple->Branch("NMu",               &fNmuTruth,          "NMu/I");
    fCRTTruthNtuple->Branch("NMuCRTTag",         &fNmuTruthCRTTag,    "NMuCRT/I");
    fCRTTruthNtuple->Branch("NMuCRTTagVec",      &fNmuTruthCRTTagVec, "NMuCRTVec/I");
    fCRTTruthNtuple->Branch("NMuCRTMiss",        &fNmuTruthMissCRT,   "NMuCRTMiss/I");
    fCRTTruthNtuple->Branch("NMuCRTTagC",        &fNmuTruthCRTTagC,   "NMuCRTTagC/I");
    fCRTTruthNtuple->Branch("NMuCRTTagM",        &fNmuTruthCRTTagM,   "NMuCRTTagM/I");
    fCRTTruthNtuple->Branch("NMuCRTTagD",        &fNmuTruthCRTTagD,   "NMuCRTTagD/I");
    fCRTTruthNtuple->Branch("NMuCRTTrig",        &fNmuTruthCRTTrig,   "NMuCRTTrig/I");
    fCRTTruthNtuple->Branch("NMuCRTTrigVec",     &fNmuTruthCRTTrigVec,"NMuCRTTrigVec/I");
    fCRTTruthNtuple->Branch("NMuCRTTrigC",       &fNmuTruthCRTTrigC,  "NMuCRTTrigC/I");
    fCRTTruthNtuple->Branch("NMuCRTTrigM",       &fNmuTruthCRTTrigM,  "NMuCRTTrigM/I");
    fCRTTruthNtuple->Branch("NMuCRTTrigD",       &fNmuTruthCRTTrigD,  "NMuCRTTrigD/I");

    fCRTTruthMatchNtuple->Branch("Event",             &fEvent,             "Event/I");
    fCRTTruthMatchNtuple->Branch("NMuCRTTrig",        &fNmuDetCRTTrig,   "NMuCRTTrig/I");
    fCRTTruthMatchNtuple->Branch("NMuCRTTrigVec",     &fNmuDetCRTTrigVec,"NMuCRTTrigVec/I");
    fCRTTruthMatchNtuple->Branch("NMuCRTTrigC",       &fNmuDetCRTTrigC,  "NMuCRTTrigC/I");
    fCRTTruthMatchNtuple->Branch("NMuCRTTrigM",       &fNmuDetCRTTrigM,  "NMuCRTTrigM/I");
    fCRTTruthMatchNtuple->Branch("NMuCRTTrigD",       &fNmuDetCRTTrigD,  "NMuCRTTrigD/I");

    // Define the branches of our DetSim n-tuple 
    fDetSimNtuple->Branch("event",                 &fDetEvent,          "event/I");
    fDetSimNtuple->Branch("nChan",                 &fNChan,             "nChan/I");
    fDetSimNtuple->Branch("channel",               fChan,               "channel[64]/I");
    fDetSimNtuple->Branch("t0",                    fT0,                 "t0[64]/I");
    fDetSimNtuple->Branch("t1",                    fT1,                 "t1[64]/I");
    fDetSimNtuple->Branch("adc",                   fADC,                "adc[64]/I");
    fDetSimNtuple->Branch("trackID",               fTrackID,            "trackID[64]/I");
    fDetSimNtuple->Branch("detPDG",                fDetPDG,             "detPDG[64]/I");
    fDetSimNtuple->Branch("entry",                 &fEntry,             "entry/I");
    fDetSimNtuple->Branch("mac5",                  &fMac5,              "mac5/I");
    fDetSimNtuple->Branch("region",                &fFEBReg,            "region/I");
    fDetSimNtuple->Branch("triggerPair",           fTriggerPair,        "triggerPair[2]/I");
    fDetSimNtuple->Branch("macPair",               fMacPair,            "macPair[2]/I");
    fDetSimNtuple->Branch("chanTrig",              &fChanTrig,          "chanTrig/I");
    fDetSimNtuple->Branch("tTrig",                 &fTTrig,             "tTrig/I");
    fDetSimNtuple->Branch("subSys",                &fDetSubSys,         "subSys/I");

    // Define the branches of our SimHit n-tuple
    fSimHitNtuple->Branch("event",       &fHitEvent,    "event/I");
    fSimHitNtuple->Branch("nHit",        &fNHit,        "nHit/I");
    fSimHitNtuple->Branch("x",           &fXHit,        "x/F");
    fSimHitNtuple->Branch("y",           &fYHit,        "y/F");
    fSimHitNtuple->Branch("z",           &fZHit,        "z/F");
    fSimHitNtuple->Branch("xErr",        &fXErrHit,     "xErr/F");
    fSimHitNtuple->Branch("yErr",        &fYErrHit,     "yErr/F");
    fSimHitNtuple->Branch("zErr",        &fZErrHit,     "zErr/F");
    fSimHitNtuple->Branch("t0",          &fT0Hit,       "t0/F");
    fSimHitNtuple->Branch("t1",          &fT1Hit,       "t1/F");
    fSimHitNtuple->Branch("region",      &fHitReg,      "region/I");  
    fSimHitNtuple->Branch("trackID",     &fHitTrk,      "trackID/I");
}
   
  void CRTAnalysis::beginRun(const art::Run& /*run*/)
  {
    //art::ServiceHandle<sim::LArG4Parameters> larParameters;
    //fElectronsToGeV = 1./larParameters->GeVToElectrons();
    std::cout << "beginning run" << std::endl;
  }

  //-----------------------------------------------------------------------
  void CRTAnalysis::analyze(const art::Event& event) 
  {
    //std::cout << "Start analysis" << std::endl;
    //LOG_DEBUG("CRT") << "beginning analyis" << '\n';
    if (fPDGs.size() != fMinMomenta.size() || fPDGs.size() != fMaxMomenta.size())
        throw cet::exception("CRTAnalysis")
          << " PDG/Momtenta values not set correctly in fhicl - lists have different sizes"
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;

    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    // Define a "handle" to point to a vector of the objects.
    art::Handle< vector<simb::MCParticle> > particleHandle;
    map< int, const simb::MCParticle*> particleMap; //particleMap.clear();

    vector<int> muTrigTrkid;
    vector<int> muTrigTrkidC;
    vector<int> muTrigTrkidM;
    vector<int> muTrigTrkidD;

    fNmuTruthCRTTrig = 0;
    fNmuTruthCRTTrigC = 0;
    fNmuTruthCRTTrigM = 0;
    fNmuTruthCRTTrigD = 0;
    fNmuTruthCRTTrigVec = 0;
    fNmuTruthMissCRT = 0;
    fNmuTruth = 0;
    fNmuTruthCRTTag = 0;
    fNmuTruthCRTTagC = 0;
    fNmuTruthCRTTagM = 0;
    fNmuTruthCRTTagD = 0;
    fNmuTruthCRTTagVec = 0;

    fNmuDetCRTTrig = 0;
    fNmuDetCRTTrigC = 0;
    fNmuDetCRTTrigM = 0;
    fNmuDetCRTTrigD = 0;
    fNmuDetCRTTrigVec = 0;

    // If there aren't any simb::MCParticle object art will 
    // display a "ProductNotFound" exception message and may skip
    // all processing for the rest of this event or stop the execution.
    if (!event.getByLabel(fSimulationProducerLabel, particleHandle)) 
      {
	// If we have no MCParticles at all in an event, then we're in
	// big trouble. Throw an exception.
	throw cet::exception("CRTAnalysis") 
	  << " No simb::MCParticle objects in this event - "
	  << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

    art::Handle<vector<sim::AuxDetSimChannel> > auxDetSimChannelHandle;
    if (!event.getByLabel(fSimulationProducerLabel, auxDetSimChannelHandle)) {
        throw cet::exception("CRTAnalysis")
          << " No sim::AuxDetSimChannel objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    //required to access TPC SimChannels
    //auto simChannelHandle =
    //  event.getValidHandle<vector<sim::SimChannel>>(fSimulationProducerLabel);

    // The MCParticle objects are not necessarily in any particular
    // order. Since we may have to search the list of particles, let's
    // put them into a map. To save both space and time, the map
    // will not contain a copy of the MCParticle, but a pointer to it.
    for ( auto const& particle : (*particleHandle) )
    {
        // Add the address of the MCParticle to the map, with the
        // track ID as the key.
        particleMap.insert(std::make_pair(particle.TrackId(),&particle));
    }

    //get TPC objects
    geo::CryostatGeo const& cryo0 = fGeometryService->Cryostat(0);
    geo::CryostatGeo const& cryo1 = fGeometryService->Cryostat(1);

    geo::TPCGeo const& tpc00 = cryo0.TPC(0);
    geo::TPCGeo const& tpc01 = cryo0.TPC(1);
    geo::TPCGeo const& tpc10 = cryo1.TPC(0);
    geo::TPCGeo const& tpc11 = cryo1.TPC(1);


    //loop over MCParticles
    for ( auto const& particle : (*particleHandle) )
    {
        fSimPDG = particle.PdgCode();
        vector<int>::iterator it = fPDGs.begin(); //iterator to list of interested PDGs from FHiCL
        const TLorentzVector& momentumStart = particle.Momentum(0);//initial momentum
        const double p = (momentumStart.Vect()).Mag();
        size_t index = 0; //index in PDG list

        //if ( (particle.Process() != "primary"  && 
        //find matching PDG code given current MCParticle PDG
        // use for momentum cuts
        while (it!=fPDGs.end()) {
          if (*it==fSimPDG) {
              index = (size_t)(it - fPDGs.begin());
              break;
          }
          it++;
        }
        //if PDG not included in interest list, skip to next
        if (!(fPDGs.size()==1 && fPDGs[0]==0) && it == fPDGs.end()) continue;
        //check momentum is within region of interest
        if ( fMinMomenta[index] != 0 && p < fMinMomenta[index]) continue;
        if ( fMaxMomenta[index] != 0 && p > fMaxMomenta[index]) continue;
        //if ( fMinMomenta.size()==1 && fMinMomenta[0]!=0 && fPDGs.size()==1 && fPDGs[0]==0 && p < fMinMomenta[0] ) continue;
        //if ( fMaxMomenta.size()==1 && fMaxMomenta[0]!=0 && fPDGs.size()==1 && fPDGs[0]==0 && p < fMaxMomenta[0] ) continue;
        //if ( abs(fSimPDG)==11 && particle.Process()!="compt" && particle.Process()!="conv" )
        //  continue;

        //count total number of muons present in event
        if (abs(fSimPDG)==13){
            fNmuTruth++;
        }
	fSimTrackID = particle.TrackId();

        // the following bit attempts to establish the ancestry
        //   of each MCParticle of interest
        // this is useful for matching gammas produced by muons
        //   for evaluation of removal algorithms
        fMother = particle.Mother();
        fNDaught = particle.NumberDaughters();
        fSimProcess = ProcessToICode(particle.Process());
        fSimEndProcess = ProcessToICode(particle.EndProcess());
       
        if(fMother!=0){ //if not primary
            map<int,const simb::MCParticle*>::const_iterator it = particleMap.find(fMother);
            if(it==particleMap.end()){
                //std::cout << "mother MCParticle object not found!" << std::endl;
                fParentPDG=INT_MAX;
                fParentE = FLT_MAX;
                fProgenitor = INT_MAX;
            }//if mother not found in particle list
            else{ //otherwise try to find primary muon 
                fParentPDG = it->second->PdgCode();
                fParentE = it->second->E(0);
                int tmpID=it->second->Mother();
                size_t ctr=0;
                map<int,const simb::MCParticle*>::iterator it2 = particleMap.begin();
		
		if(fParentPDG==13||fParentPDG==-13) fProgenitor = it->second->TrackId();
                else while(it2!=particleMap.end()&&ctr<particleMap.size()){
                    it2=particleMap.find(tmpID);
                    if(it2!=particleMap.end()){
                        if(it2->second->PdgCode()==13||it2->second->PdgCode()==-13){
                            fProgenitor=tmpID;
                            break;
                        }
                        tmpID = it2->second->Mother();
                    }

                    ctr++;
                    if(ctr==particleMap.size()) std::cout<<"manual break!"<<std::endl;
                }
                //std::cout<<"out of progenitor search loop" << std::endl;
            }//else mother is found
        }//if particle has mother (not a primary)
        else{ //if current particle in loop is primary
            fParentPDG=0;
            fProgenitor=-10;
            fParentE=-1.0;
        }
        //end of ancestry matching
        //now get some other useful info about the trajectory

	// A particle has a trajectory, consisting of a set of
	// 4-positions and 4-mommenta.
	fSimHits = particle.NumberTrajectoryPoints();

	// For trajectories, as for vectors and arrays, the first
	// point is #0, not #1.
	const int last = fSimHits - 1;
	const TLorentzVector& positionStart = particle.Position(0);
	const TLorentzVector& positionEnd   = particle.Position(last);
	//const TLorentzVector& momentumStart = particle.Momentum(0);
	const TLorentzVector& momentumEnd   = particle.Momentum(last);

	// Fill arrays with the 4-values.
	positionStart.GetXYZT( fStartXYZT );
	positionEnd.GetXYZT( fEndXYZT );
	momentumStart.GetXYZT( fStartPE );
	momentumEnd.GetXYZT( fEndPE );
        fTrackLength = ( positionEnd - positionStart ).Rho();

        fRegEvent = fEvent;
        fNReg = 0;
        fRegFid = 0;
        fRegActive = 0;
        fRegInactive = 0;
        fRegCryo = 0;
        fRegCRTs = 0;
        //std::cout << "initializing region arrays" << std::endl;

        //initialize arrays
        for (int i=0; i<(int)kMaxReg; i++){
                fRegRegions[i] = -INT_MAX;
                fRegPDG[i] = -INT_MAX;
                fRegTrkID[i] = -INT_MAX;
                fRegEDep[i] = -FLT_MAX;
                for (int j=0; j<4; j++){
                        fRegEntryXYZT[i][j] = -INT_MAX;
                        fRegExitXYZT[i][j] = -INT_MAX;
                        fRegEntryPE[i][j] = -INT_MAX;
                        fRegExitPE[i][j] = -INT_MAX;
                        if (j<3) {
                            fRegEntrySlope[i][j] = -INT_MAX;
                            fRegExitSlope[i][j] = -INT_MAX;
                        }
                }
        }

        int oldreg = -1, oldcryo = -1;
        double tmp4_cryo[4] = {}, tmp4_active[4] = {};
        double tmp4pe_cryo[4] = {}, tmp4pe_active[4] = {};

        //std::cout << "about to loop over trajectory points" << std::endl;
        //loop over trajectory points and extract first and last points in TPC
        for (unsigned int i=0; i<fSimHits; i++){
                const TLorentzVector& pos = particle.Position(i);
                const TLorentzVector& posnext = particle.Position(i+1);
                const TLorentzVector& mom = particle.Momentum(i);
                const double point[3] = {pos(0),pos(1),pos(2)};
                const double pointnext[3] = {posnext(0),posnext(1),posnext(2)};
                //double pointlocal[3] = {0.,0.,0.};

                if(fNReg==kMaxReg) std::cout << "about to seg fault..need more NReg!" << std::endl;

                if(cryo0.ContainsPosition(point)) {
                        if (oldcryo!=10) {
                            pos.GetXYZT(tmp4_cryo);
                            mom.GetXYZT(tmp4pe_cryo);
                            oldcryo = 10;
                        }
                        if (oldreg!=5&&oldreg!=6&&oldreg!=10) oldreg=10;
                        if (!cryo0.ContainsPosition(pointnext)||i==fSimHits-1){
                                for (int j=0; j<4; j++){
                                    fRegEntryXYZT[fNReg][j] = tmp4_cryo[j];
                                    fRegEntryPE[fNReg][j]   = tmp4pe_cryo[j];
                                }
                                pos.GetXYZT(fRegExitXYZT[fNReg]);
                                mom.GetXYZT(fRegExitPE[fNReg]);
                                fRegRegions[fNReg] = 10;
                                fRegdL[fNReg] = TMath::Sqrt(TMath::Power(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
                                fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
                                fNReg++; fRegInactive++;
                        }
                }
                if(cryo1.ContainsPosition(point)) {
                        if (oldcryo!=12){
                            pos.GetXYZT(tmp4_cryo);
                            mom.GetXYZT(tmp4pe_cryo);
                            oldcryo = 12;
                        }
                        if (oldreg!=7&&oldreg!=8&&oldreg!=12) oldreg=12;
                        if (!cryo1.ContainsPosition(pointnext)||i==fSimHits-1){
                                for (int j=0; j<4; j++){
                                    fRegEntryXYZT[fNReg][j] = tmp4_cryo[j];
                                    fRegEntryPE[fNReg][j]   = tmp4pe_cryo[j];
                                }
                                pos.GetXYZT(fRegExitXYZT[fNReg]);
                                mom.GetXYZT(fRegExitPE[fNReg]);
                                fRegRegions[fNReg] = 12;
                                fRegdL[fNReg] = TMath::Sqrt(TMath::Power(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
                                fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
                                fNReg++; fRegInactive++;
                        }
                }

                if(tpc00.ContainsPosition(point)) {
                        if (oldreg!=5) {
                            pos.GetXYZT(tmp4_active);
                            mom.GetXYZT(tmp4pe_active);
                            oldreg = 5;
                        }
                        if(tpc00.InFiducialX(point[0],25,0)&&tpc00.InFiducialY(point[1],25,25)
                              &&tpc00.InFiducialZ(point[2],30,50)) fRegFid++;
                        if (!tpc00.ContainsPosition(pointnext)||i==fSimHits-1){
                                for (int j=0; j<4; j++){
                                    fRegEntryXYZT[fNReg][j] = tmp4_active[j];
                                    fRegEntryPE[fNReg][j]   = tmp4pe_active[j];
                                }
                                pos.GetXYZT(fRegExitXYZT[fNReg]);
                                mom.GetXYZT(fRegExitPE[fNReg]);
                                fRegRegions[fNReg] = 5;
                                fRegdL[fNReg] = TMath::Sqrt(TMath::Power(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
                                fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
                                fNReg++; fRegActive++;
                         }
                }
                if(tpc01.ContainsPosition(point)) {
                        if (oldreg!=6){
                            pos.GetXYZT(tmp4_active);
                            mom.GetXYZT(tmp4pe_active);
                            oldreg = 6;
                        }
                        if(tpc01.InFiducialX(point[0],0,25)&&tpc01.InFiducialY(point[1],25,25)
                              &&tpc01.InFiducialZ(point[2],30,50)) fRegFid++;
                        if (!tpc01.ContainsPosition(pointnext)||i==fSimHits-1){
                                for (int j=0; j<4; j++){
                                    fRegEntryXYZT[fNReg][j] = tmp4_active[j];
                                    fRegEntryPE[fNReg][j]   = tmp4pe_active[j];
                                }
                                pos.GetXYZT(fRegExitXYZT[fNReg]);
                                mom.GetXYZT(fRegExitPE[fNReg]);
                                fRegRegions[fNReg] = 6;
                                fRegdL[fNReg] = TMath::Sqrt(TMath::Power(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
                                fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
                                fNReg++; fRegActive++;
                         }
                }
                if(tpc10.ContainsPosition(point)) {
                        if (oldreg!=7){
                            pos.GetXYZT(tmp4_active);
                            mom.GetXYZT(tmp4pe_active);
                            oldreg = 7;
                        }
                        if(tpc10.InFiducialX(point[0],25,0)&&tpc10.InFiducialY(point[1],25,25)
                              &&tpc10.InFiducialZ(point[2],30,50)) fRegFid++;
                        if (!tpc10.ContainsPosition(pointnext)||i==fSimHits-1){
                                for (int j=0; j<4; j++){
                                    fRegEntryXYZT[fNReg][j] = tmp4_active[j];
                                    fRegEntryPE[fNReg][j]   = tmp4pe_active[j];
                                }
                                pos.GetXYZT(fRegExitXYZT[fNReg]);
                                mom.GetXYZT(fRegExitPE[fNReg]);
                                fRegRegions[fNReg] = 7;
                                fRegdL[fNReg] = TMath::Sqrt(TMath::Power(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
                                fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
                                fNReg++; fRegActive++;
                         }
                }
                if(tpc11.ContainsPosition(point)) {
                        if (oldreg!=8){
                            pos.GetXYZT(tmp4_active);
                            mom.GetXYZT(tmp4pe_active);
                            oldreg = 8;
                        }
                        if(tpc11.InFiducialX(point[0],0,25)&&tpc11.InFiducialY(point[1],25,25)
                              &&tpc11.InFiducialZ(point[2],30,50)) fRegFid++;
                        if (!tpc11.ContainsPosition(pointnext)||i==fSimHits-1){
                                for (int j=0; j<4; j++){
                                    fRegEntryXYZT[fNReg][j] = tmp4_active[j];
                                    fRegEntryPE[fNReg][j]   = tmp4pe_active[j];
                                }
                                pos.GetXYZT(fRegExitXYZT[fNReg]);
                                mom.GetXYZT(fRegExitPE[fNReg]);
                                fRegRegions[fNReg] = 8;
                                fRegdL[fNReg] = TMath::Sqrt(TMath::Power(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                               +TMath::Power(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
                                fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
                                fNReg++; fRegActive++;
                         }
                }
        }//for trajectory points

        fRegionsNtuple->Fill();

        //std::cout << "out of point loop" << std::endl;

	/*LOG_DEBUG("CRTAnalysis")
	  << "track ID=" << fSimTrackID 
	  << " (PDG ID: " << fSimPDG << ") "
	  << trackLength << " cm long, momentum " 
	  << momentumStart.P() << " GeV/c, has " 
	  << fSimHits << " trajectory points";
	   */ 

        //map module IDs to strip IDs hit by muons
        map< uint16_t,set<uint8_t>* > muHitMapC; //hits in C modules only
        map< uint16_t,set<uint8_t>* > muHitMapM; //hits in M modules only
        map< uint16_t,set<uint8_t>* > muHitMapD; //hits in D modules only
        map< uint16_t,set<uint8_t>* > muHitMap; //all hits

        //reinitialize ADChannel vars
	fNAuxDet = 0;
        for (size_t i=0; i<kMaxAD; i++) {
            fADTrackLength[i] = FLT_MAX;
            fADEDep[i] = FLT_MAX;
            fADdEdx[i] = FLT_MAX;
            fAuxDetID[i] = UINT32_MAX;
            fAuxDetSensitiveID[i] = UINT32_MAX;
            fAuxDetReg[i] = UINT32_MAX;
            fADMac[i] = UINT32_MAX;
            fADType[i] = UINT32_MAX;
            for (size_t j=0; j<4; j++) {
                fADEnterXYZT[i][j] = FLT_MAX;
                fADExitXYZT[i][j] = FLT_MAX;
            }
        }

	// To look at the energy deposited by this particle's track,
	// we loop over the AuxDetSimChannel objects in the event. 
	// Note all volumes are included, not just ones with energy deps
	for ( auto const& channel : (*auxDetSimChannelHandle) )
	{  
            auto const& adGeo = fGeometryService->AuxDet(channel.AuxDetID());

	    // Get vector of hits in this AuxDet channel
	    auto const& auxDetIDEs = channel.AuxDetIDEs();

	    // For every hit in this channel:
	    for ( auto const& ide : auxDetIDEs )
	    {
		    // Check if the track that deposited the
		    // energy matches the track of the MCParticle.
		    if ( ide.trackID != fSimTrackID ) continue;
		    if ( ide.energyDeposited * 1.0e6 < 50 ) continue; 
                    // Ignore strips w/ID=0 to get around bug (will be fixed soon)
                    if ( channel.AuxDetSensitiveID() == 0 ) continue;


                    //calculate track length in strip
		    double dx = ide.entryX-ide.exitX;
		    double dy = ide.entryY-ide.exitY;
		    double dz = ide.entryZ-ide.exitZ;
                    double adlength = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
                    if ( adlength < 0.0001)  continue;

		    fADTrackLength[fNAuxDet] = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
	            fADEDep[fNAuxDet] = ide.energyDeposited;
	            fADdEdx[fNAuxDet] = ide.energyDeposited/fADTrackLength[fNAuxDet];
		    fAuxDetID[fNAuxDet] = channel.AuxDetID();
		    fAuxDetSensitiveID[fNAuxDet] = channel.AuxDetSensitiveID();
		    fADEnterXYZT[fNAuxDet][0] = ide.entryX;
	            fADEnterXYZT[fNAuxDet][1] = ide.entryY;
	            fADEnterXYZT[fNAuxDet][2] = ide.entryZ;
	            fADEnterXYZT[fNAuxDet][3] = ide.entryT;
                    fADExitXYZT[fNAuxDet][0] = ide.exitX;
                    fADExitXYZT[fNAuxDet][1] = ide.exitY;
                    fADExitXYZT[fNAuxDet][2] = ide.exitZ;
                    fADExitXYZT[fNAuxDet][3] = ide.exitT;
                    fAuxDetReg[fNAuxDet] = GetAuxDetRegion(fGeometryService->AuxDet(channel.AuxDetID()));
                    fADMac[fNAuxDet] = ADToMac(ModToAuxDetType(adGeo),channel.AuxDetID());
                    //fADType[fNAuxDet] = ModToAuxDetType(adGeo);
                    fADType[fNAuxDet] = ModToTypeCode(adGeo);
		    fNAuxDet++;

                    if (abs(fSimPDG)==13) {
                        switch (ModToAuxDetType(adGeo)) {

                          case 'c' : 
                            if( muHitMapC.find(channel.AuxDetID())==muHitMapC.end() )
                              muHitMapC[channel.AuxDetID()] = new set<uint8_t>();
                            muHitMapC[channel.AuxDetID()]->insert(channel.AuxDetSensitiveID());
                            break;

                          case 'm' : 
                            if( muHitMapM.find(channel.AuxDetID())==muHitMapM.end() )
                              muHitMapM[channel.AuxDetID()] = new set<uint8_t>();
                            muHitMapM[channel.AuxDetID()]->insert(channel.AuxDetSensitiveID());
                            break;

                          case 'd' :
                            if( muHitMapD.find(channel.AuxDetID())==muHitMapD.end() )
                              muHitMapD[channel.AuxDetID()] = new set<uint8_t>();
                            muHitMapD[channel.AuxDetID()]->insert(channel.AuxDetSensitiveID());
                            break;
                        }//switch

                        if( muHitMap.find(channel.AuxDetID())==muHitMap.end() )
                          muHitMap[channel.AuxDetID()] = new set<uint8_t>();
                        muHitMap[channel.AuxDetID()]->insert(channel.AuxDetSensitiveID());
                    }//if muon
	    } // For each IDE (strip hit by muon)
              
	} // For each SimChannel (module)

        // write values to tree for this event and particle
        fSimulationNtuple->Fill();

        set<uint16_t> modsTrigPoss;
        set<uint16_t> modsTrigPossC;
        set<uint16_t> modsTrigPossM;
        set<uint16_t> modsTrigPossD;

        if( abs(fSimPDG)==13 ) {

          if (muHitMapC.size()>0) {
              fNmuTruthCRTTagC++;
              fModMultHistC->Fill(muHitMapC.size());
              for( auto const& strips : muHitMapC ) {
                 fStripMultHistC->Fill(strips.second->size());
                 if (strips.second->size()>1) {
                   modsTrigPossC.insert(strips.first);
                   modsTrigPoss.insert(strips.first);
                 }//if more than 1 strip hit
              }//for modules in muHitMap
              if (modsTrigPossC.size()>0) {
                fNmuTruthCRTTrigC++;
                muTrigTrkidC.push_back(fSimTrackID);
              }
          } //if c hit

          if (muHitMapM.size()>0) {
              fNmuTruthCRTTagM++;
              fModMultHistM->Fill(muHitMapM.size());
              for( auto const& strips : muHitMapM ) {
                 fStripMultHistM->Fill(strips.second->size());
                 if (muHitMapM.size()>0) {
                   modsTrigPossM.insert(strips.first);
                   modsTrigPoss.insert(strips.first);
                 }//if more than 1 module hit
              }//for modules in muHitMap
              if (modsTrigPossM.size()>1) {
                fNmuTruthCRTTrigM++;
                muTrigTrkidM.push_back(fSimTrackID);
              }
          } //if m hit

          if (muHitMapD.size()>0) {
              fNmuTruthCRTTagD++;
              fModMultHistD->Fill(muHitMapD.size());
              for( auto const& strips : muHitMapC ) {
                 fStripMultHistD->Fill(strips.second->size());
                 if (strips.second->size()>1) {
                   modsTrigPossD.insert(strips.first);
                   modsTrigPoss.insert(strips.first);
                 }//if more than 1 strip hit
              }//for modules in muHitMap
              if (modsTrigPossD.size()>0) {
                fNmuTruthCRTTrigD++;
                muTrigTrkidD.push_back(fSimTrackID);
              }
          } //if d hit

          if (muHitMap.size()>0) {
            fNmuTruthCRTTag++;
          }
          else
            fNmuTruthMissCRT++;

          if (muHitMap.size()>1)
            fNmuTruthCRTTagVec++;

          //fix me! doesn't check coincidence condition!
          if (modsTrigPoss.size()>0){
            fNmuTruthCRTTrig++;
            muTrigTrkid.push_back(fSimTrackID);
          }
          if (modsTrigPoss.size()>1)
            fNmuTruthCRTTrigVec++;

        }//if muon

    } // loop over all particles in the event. 

    fCRTTruthNtuple->Fill();

    art::Handle<vector<icarus::crt::CRTData>> crtDetSimHandle;
    bool isCRTDetSim = event.getByLabel(fCRTDetSimProducerLabel, crtDetSimHandle);
    set<uint8_t> febidsC;
    set<uint8_t> febidsM;
    set<uint8_t> febidsD;
    set<uint8_t> taggedTrksC;
    set<uint8_t> taggedTrksM;
    set<uint8_t> taggedTrksD;
    set<uint8_t> taggedTrks;

    if (isCRTDetSim)  {
     std::cout << "about to loop over detsim entries" << std::endl;
     for ( auto const& febdat : (*crtDetSimHandle) ) {
        fDetEvent       = febdat.Event();
        fMac5           = febdat.Mac5();
        fChanTrig       = febdat.ChanTrig();
        fEntry          = febdat.Entry();
        fTTrig          = febdat.TTrig();
        pair<uint32_t,uint32_t> tmpPair = febdat.TrigPair();
        fTriggerPair[0] = tmpPair.first;
        fTriggerPair[1] = tmpPair.second;
        tmpPair         = febdat.MacPair();
        fMacPair[0]     = tmpPair.first;
        fMacPair[1]     = tmpPair.second;
        fFEBReg         = MacToADReg(fMac5);
        fNChan = 0;
        //fDetSubSys = MacToType(fMac5);
        fDetSubSys = MacToTypeCode(fMac5);

        for (int i=0; i<kDetMax; i++) {
          fChan[i] = INT_MAX;
          fT0[i]   = INT_MAX;
          fT1[i]   = INT_MAX;
          fADC[i]  = INT_MAX;
          fTrackID[i] = INT_MAX;
          fDetPDG[i] = INT_MAX;
        }

        set<uint8_t> chanidsC;
        set<uint8_t> chanidsM;
        set<uint8_t> chanidsD;
 
        //std::cout << "loop over chandata" << std::endl;
        for ( auto const chandat : febdat.ChanData()) {
          //DetSim tree contains all entries (not neccessarily from muons)
          fChan[fNChan]    = chandat.Channel();
          fT0[fNChan]      = chandat.T0();
          fT1[fNChan]      = chandat.T1();
          fADC[fNChan]     = chandat.ADC();
          fTrackID[fNChan] = chandat.TrackID()[0]; 
          //fNChan++;

          //std::cout << "check if muon" << std::endl;
          //if channel hit came from muon loop over all associated hits
          if( particleMap.find(chandat.TrackID()[0])!=particleMap.end()) {
              fDetPDG[fNChan]  = particleMap[chandat.TrackID()[0]]->PdgCode();

              if (abs(particleMap[chandat.TrackID()[0]]->PdgCode())==13 ) {

              std::cout << "found muon!" << std::endl;
              taggedTrks.insert(chandat.TrackID()[0]);
              for (auto const chandat2 : febdat.ChanData()) {
                  if( chandat.TrackID()[0] != chandat2.TrackID()[0] ) continue;
                  switch( fDetSubSys ){
                    case 'c' :
                      febidsC.insert(fMac5);
                      chanidsC.insert(chandat2.Channel());
                      taggedTrksC.insert(chandat.TrackID()[0]);
                      break;
                    case 'm' :
                      febidsM.insert(fMac5);
                      chanidsM.insert(chandat2.Channel());
                      taggedTrksM.insert(chandat.TrackID()[0]);
                      break;
                    case 'd' :
                      febidsD.insert(fMac5);
                      chanidsD.insert(chandat2.Channel());
                      taggedTrksD.insert(chandat.TrackID()[0]);
                      break;
                  }//switch
              }//inner chandat loop
              switch (fDetSubSys) {
                case 'c' :
                  fChanMultHistC->Fill(chanidsC.size());
                  chanidsC.clear();
                  break;
                case 'm' :
                  fChanMultHistM->Fill(chanidsM.size());
                  chanidsM.clear();
                  break;
                case 'd' :
                  fChanMultHistD->Fill(chanidsD.size());
                  chanidsD.clear();
                  break;
              }
            }//if muon 
         }
         else {
             std::cout << "trackID from DetSim not found in ParticleMap!" << std::endl; 
             fDetPDG[fNChan] = INT_MAX;
         }
         fNChan++;
          //std::cout << "end of loop over chandat" << std::endl;
        }//outer chandat loop

        //std::cout << "about to fill detsimtree" << std::endl;
        fDetSimNtuple->Fill();

     } //for CRT FEB events

     std::cout << '\n'
        << "muTrigTrkid size: " << muTrigTrkid.size() << '\n'
        << "muTrigTrkidC size: " << muTrigTrkidC.size() << '\n'
        << "muTrigTrkidM size: " << muTrigTrkidM.size() << '\n'
        << "muTrigTrkidD size: " << muTrigTrkidD.size() << '\n'
        << '\n'
        << "taggedTrks size: " << taggedTrks.size() << '\n'
        << "taggedTrksC size: " << taggedTrksC.size() << '\n'
        << "taggedTrksM size: " << taggedTrksM.size() << '\n'
        << "taggedTrksD size: " << taggedTrksD.size() << '\n'
     << std::endl;

     for (auto id : muTrigTrkid) {
        for (auto detrks : taggedTrks) {
           if (id==detrks) {
              fNmuDetCRTTrig++;
              break;
           }
        }
     }   

     for (auto id : muTrigTrkidC) {
        for (auto detrks : taggedTrksC) {
           if (id==detrks) {
              fNmuDetCRTTrigC++;
              break;
           }
        }
     }

     for (auto id : muTrigTrkidM) {
        for (auto detrks : taggedTrksM) {
           if (id==detrks) {
              fNmuDetCRTTrigM++;
              break;
           }
        }
     }

     for (auto id : muTrigTrkidD) {
        for (auto detrks : taggedTrksD) {
           if (id==detrks) {
              fNmuDetCRTTrigD++;
              break;
           }
        }
     }

     fCRTTruthMatchNtuple->Fill();

     fFEBMultHistC->Fill(febidsC.size());
     fFEBMultHistM->Fill(febidsM.size());
     fFEBMultHistD->Fill(febidsD.size());

     fEffTot = 1.0*fNmuDetCRTTrig/fNmuTruthCRTTrig;
     fEffC = 1.0*fNmuDetCRTTrigC/fNmuTruthCRTTrigC;
     fEffM = 1.0*fNmuDetCRTTrigM/fNmuTruthCRTTrigM;
     fEffD = 1.0*fNmuDetCRTTrigD/fNmuTruthCRTTrigD;

     mf::LogInfo("CRT") << '\n'
       << " Total muon tracks entering ADS: " << fNmuTruthCRTTag << '\n'
       << " Total muon tracks in C ADS: " << fNmuTruthCRTTagC << '\n'
       << " Total muon tracks in M ADS: " << fNmuTruthCRTTagM << '\n'
       << " Total muon tracks in D ADS: " << fNmuTruthCRTTagD << '\n'
       << " Total muon tracks w/coinc.: " << fNmuTruthCRTTrig << " (" 
          << 1.0*fNmuTruthCRTTrig/fNmuTruthCRTTag << ")" << '\n'
       << " Total muon tracks in C w/coinc.: " << fNmuTruthCRTTrigC << " ("
          << 1.0*fNmuTruthCRTTrigC/fNmuTruthCRTTagC << ")" << '\n'
       << " Total muon tracks in M w/coinc.: " << fNmuTruthCRTTrigM << " ("
          << 1.0*fNmuTruthCRTTrigM/fNmuTruthCRTTagM << ")" << '\n'
       << " Total muon tracks in D w/coinc.: " << fNmuTruthCRTTrigD << " ("
          << 1.0*fNmuTruthCRTTrigD/fNmuTruthCRTTagD << ")" << '\n'
       << " EffC: " << fEffC << '\n'
       << " EffM: " << fEffM << '\n'
       << " EffD: " << fEffD << '\n'
       << " EffTot: " << fEffTot << '\n';

    }//if crtdetsim products present

    else std::cout << "CRTData products not found! (expected if gen/G4 step)" << std::endl;

    art::Handle<vector<icarus::crt::CRTHit>> crtSimHitHandle;
    bool isCRTSimHit = event.getByLabel(fCRTSimHitProducerLabel, crtSimHitHandle);
    if (isCRTSimHit) std::cout << "Found " << crtSimHitHandle->size() << " CRT Hits" << std::endl;

    fNHit = 0;

    if (isCRTSimHit) {
        for ( auto const& hit : (*crtSimHitHandle) )
        {
            fNHit++;
            fXHit    = hit.X();
            fYHit    = hit.Y();
            fZHit    = hit.Z();
            fXErrHit = hit.XErr();
            fYErrHit = hit.YErr();
            fZErrHit = hit.ZErr();
            fT0Hit   = hit.T0();
            fT1Hit   = hit.T1();
            fHitReg  = hit.Region();
            fHitTrk  = hit.TrackID();

            fSimHitNtuple->Fill();
        }//for CRT Hits
    }//if CRT Hits

    else std::cout << "CRTHit products not found! (expected if gen/G4/detsim step)" << std::endl;


  } // CRTAnalysis::analyze()
  
  
  DEFINE_ART_MODULE(CRTAnalysis)

} // namespace crt
} // namespace icarus


// Back to our local namespace.
namespace {

  char ModToAuxDetType(geo::AuxDetGeo const& adgeo) {
    size_t nstrips = adgeo.NSensitiveVolume();
    if (nstrips==16) return 'c'; 
    if (nstrips==20) return 'm';
    if (nstrips==64) return 'd';
    return 'e';
  }

  uint32_t ModToTypeCode(geo::AuxDetGeo const& adgeo) {
    size_t nstrips = adgeo.NSensitiveVolume();
    if (nstrips==16) return 0; //'c'
    if (nstrips==20) return 1; //'m'
    if (nstrips==64) return 2; //'d'
    return UINT32_MAX;
  }

  /*size_t ModToAuxDetRegion(size_t mod)
  {
    if ( mod >= 0  && mod < 54  ) return 50; //left
    if ( mod > 53  && mod < 108 ) return 54; //right
    if ( mod > 107 && mod < 128 ) return 44; //front
    if ( mod > 127 && mod < 148 ) return 42; //back
    if ( mod > 147 && mod < 162 ) return 58; //bottom
    if ( mod > 161 && mod < 246 ) return 38; //top
    if ( mod > 277 && mod < 284 ) return 46; //slope back
    return UINT_MAX;
  }*/

  int GetAuxDetRegion(geo::AuxDetGeo const& adgeo)
  {
    char type = ModToAuxDetType(adgeo);
    string base = "volAuxDet_";//module_xxx_";
    string base2 = "_module_xxx_";
    string volName(adgeo.TotalVolume()->GetName());

    if (type == 'm') base += "MINOS" + base2;
    if (type == 'd') return 58;
    if (type == 'c') base += "CERN" + base2;

    string reg  = volName.substr(base.size());

    if(reg == "Top")        return 38;
    if(reg == "SlopeLeft")  return 52;
    if(reg == "SlopeRight") return 56;
    if(reg == "SlopeFront") return 48;
    if(reg == "SlopeBack")  return 46;
    if(reg == "Left")       return 50;
    if(reg == "Right")      return 54;
    if(reg == "Front")      return 44;
    if(reg == "Back")       return 42;
    if(reg == "Bottom")     return 58;
    return -1;
  }


  int ProcessToICode(string const& p)
  {
    int icode = -1;

    if(p=="primary")               icode=0;

    //gamma (PDG 22)
    if(p=="eBrem")                 icode=1; 
    if(p=="muBrems")               icode=2; 
    if(p=="annihil")               icode=3;
    if(p=="muMinusCaptureAtRest")  icode=4;
    if(p=="nCapture")              icode=5;
    if(p=="hBertiniCaptureAtRest") icode=6;
    if(p=="photonNuclear")         icode=7;
    if(p=="muonNuclear")           icode=8;
    if(p=="neutronInelastic")      icode=9; 
    if(p=="protonInelastic")       icode=10;
    if(p=="pi-Inelastic")          icode=11;
    if(p=="Decay")                 icode=12;

    //e (PDG +/- 11)
    if (p=="muIoni")               icode=13;
    if (p=="eIoni")                icode=14;
    if (p=="hIoni")                icode=15;
    if (p=="compt")                icode=16;
    if (p=="conv")                 icode=17;
    if (p=="muPairProd")           icode=18;
    if (p=="phot")                 icode=19;

    if (p=="hadElastic")           icode=20;

    if (p=="LArVoxelReadoutScoringProcess") icode=30;
    if (p=="CoupledTransportation")         icode=31;

    return icode;

  }

  uint32_t MacToADReg(uint32_t mac) {

      if(mac>=162 && mac<=245) return 38; //top
      if(mac>=246 && mac<=258) return 52; //slope left
      if(mac>=259 && mac<=271) return 56; //slope right
      if(mac>=272 && mac<=277) return 48; //slope front
      if(mac>=278 && mac<=283) return 46; //slope back
      if(            mac<=17 ) return 50; //left  NB removed the test mac >= 0 since always true
      if(mac>=50  && mac<=67 ) return 50; //left
      if(mac>=18  && mac<=35 ) return 54; //right
      if(mac>=68  && mac<=85 ) return 54; //right
      if(mac>=36  && mac<=42 ) return 44; //front
      if(mac>=86  && mac<=92 ) return 44; //front
      if(mac>=43  && mac<=49 ) return 42; //back
      if(mac>=93  && mac<=99 ) return 42; //back
      if(mac>=148 && mac<=161) return 58; //bottom

      return 0;
  }

  /*char MacToType(uint32_t mac) {

      uint32_t reg = MacToADReg(mac);

      if( reg==38 || reg==52 || reg==56 || reg==48 || reg==46 )
        return 'c';
      if( reg==50 || reg==54 || reg==44 || reg==42 )
        return 'm';
      if( reg==58)
        return 'd';

      return 'e';
  }*/

  uint32_t MacToTypeCode(uint32_t mac) {

      uint32_t reg = MacToADReg(mac);

      if( reg==38 || reg==52 || reg==56 || reg==48 || reg==46 )
        return 0; //'c';
      if( reg==50 || reg==54 || reg==44 || reg==42 )
        return 1; //'m';
      if( reg==58)
        return 2; //'d';

      return UINT32_MAX;//'e';
  }

  //for C- and D-modules, mac address is same as AD ID
  //three M-modules / FEB, each modules readout at both ends
  //  numbering convention is module from FEB i 
  //  is readout on the opposite end by FEB i+50
  //  return FEB i
  uint32_t ADToMac(char type, uint32_t adid) {

      switch (type){
          case 'c' :
              return adid;
              //break;
          case 'd' :
              return adid;
              //break;
          case 'm' :
              return adid/3;
              //break;
      }
      return UINT32_MAX;
  }

  // for C-modules, 2 channels per strip, return lowest val chan
  // for M-modules, 2 strips per channel, readout at both ends
  // fir D-modules, 1 channel per strip, adsid same as channel
  /*uint32_t ADSToFEBChannel (char type, uint32_t adid, uint32_t adsid) {

      switch (type){
          case 'c' :
              return 2*adsid;
          case 'd' :
              return adsid;
          case 'm' :
              return adsid/2 + 10*(adid % 3);
      }
      return UINT32_MAX;
  }*/

}//local namespace
