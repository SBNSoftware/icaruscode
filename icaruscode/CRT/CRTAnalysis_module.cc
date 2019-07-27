/**
 * @file   CRTAnalysis_module.cc
 * @brief  Access CRT data and reco products and compare to MCTruth info 
 * @author Chris Hilgenberg (Chris.Hilgenberg@colostate.edu)
 * 
 * The last revision of this code was done in October 2018 with LArSoft v07_06_01.
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
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTTruthMatchUtils.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;

//using cmath::pow;
//using cmath::sqrt;

namespace {
  void FillFebMap(map<int,vector<pair<int,int>>>& m);
  uint32_t ModToTypeCode(geo::AuxDetGeo const& adgeo); 
  char ModToAuxDetType(geo::AuxDetGeo const& adgeo);
  int GetAuxDetRegion(geo::AuxDetGeo const& adgeo);
  int ProcessToICode(string const& p);
  uint32_t MacToADReg(uint32_t mac);
  char MacToType(uint32_t mac);
  uint32_t MacToTypeCode(uint32_t mac);
  std::pair<uint32_t,uint32_t> ADToMac(const map<int,vector<pair<int,int>>>& febMap, uint32_t adid);
  //uint32_t ADSToFEBChannel(char type, uint32_t adid, uint32_t adsid);
  int MacToAuxDetID(const map<int,vector<pair<int,int>>>& febMap, int mac, int chan);
  int ChannelToAuxDetSensitiveID(int mac, int chan);

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

    // The n-tuples we'll create.
    TTree* fCosmicDisplayNtuple;  ///< for ROOT based event display
    TTree* fGenNtuple;
    TTree* fSimulationNtuple;     ///< tuple with simulated data
    TTree* fRegionsNtuple;
    TTree* fDetSimNtuple;
    TTree* fSimHitNtuple;
    TTree* fTrueCRTHitNtuple;

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
    static const int kMaxSeg = 4000;
    int        fCDEvent;
    int        fCDTrackID;
    int        fNCD;
    int        fCDRegions[kMaxSeg];
    int        fCDpdg;
    double     fCDSlopes[kMaxSeg][3]; //direction cosines
    double     fCDpe[kMaxSeg][4];   //4-momentum
    double     fCDxyzt[kMaxSeg][4];   //4-position
    /// @}

    /// @name The variables that will go into the Gen n-tuple.
    /// @{
    static const int kMaxGen = 500;
    int      fNGen;
    int      fGenTrack[kMaxGen];
    int      fGenPDG[kMaxGen];
    double   fGenStartXYZT[kMaxGen][4];
    double   fGenEndXYZT[kMaxGen][4];
    double   fGenStartPE[kMaxGen][4];
    double   fGenEndPE[kMaxGen][4];
    /// @}

    /// @name The variables that will go into the Simulation n-tuple.
    /// @{
    static const int kMaxAD = 100; //for memory allocation
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
    static const size_t kMaxReg = 500;
    int      fRegEvent;
    int      fNReg;
    int      fRegFid;
    int      fRegActive;
    int      fRegInactive;
    int      fRegCRTs;
    int      fRegTrkID;
    int      fRegPDG;
    int      fRegRegions[kMaxReg];
    double   fRegEDep[kMaxReg];
    double   fRegdL[kMaxReg];
    double   fRegEntryPE[kMaxReg][4];
    double   fRegExitPE[kMaxReg][4];
    double   fRegEntryXYZT[kMaxReg][4];
    double   fRegExitXYZT[kMaxReg][4];
    double   fRegEntrySlope[kMaxReg][3];
    double   fRegExitSlope[kMaxReg][3];
    double   fRegDistToOpDet[kMaxReg];
    double   fRegOpDetXYZT[kMaxReg][4];
    int      fRegOpDetID[kMaxReg];

    //CRT data product vars
    static const int kDetMax = 64;
    int      fDetEvent;
    int      fChan[kDetMax]; ///< front-end board channel (0-31 or 0-63)
    double      fT0[kDetMax]; ///< signal time w.r.t. global event time
    double      fT1[kDetMax]; ///< signal time w.r.t. PPS
    int      fADC[kDetMax]; ///< signal amplitude
    int      fTrackID[kDetMax]; ///< track ID of particle that produced the signal
    int      fDetPDG[kDetMax];
    int      fNChan; ///< number of channels above threshold for this front-end board readout
    int      fEntry; ///< front-end board entry number (reset for each event)
    double      fTTrig;      ///< signal time w.r.t. global event time of primary channel providing trigger
    int      fChanTrig; ///< channel which provided the trigger for this front-end board
    int      fFEBReg; ///< CRT region for this front-end board
    int      fMac5; ///< Mac5 address for this front-end board
    int      fTriggerPair[2]; ///< two channels which provided the coincidence (useful for C or D modules)
    int      fMacPair[2]; ///< two front-end boards with pairwise coincidence ( useful for M modules)
    int      fDetSubSys;
      
    //CRT hit product vars
    int       fHitEvent;
    float    fXHit; ///< reconstructed X position of CRT hit (cm)
    float    fYHit; ///< reconstructed Y position of CRT hit (cm)
    float    fZHit; ///< reconstructed Z position of CRT hit (cm)
    float    fXErrHit; ///< stat error of CRT hit reco X (cm)
    float    fYErrHit; ///< stat error of CRT hit reco Y (cm)
    float    fZErrHit; ///< stat error of CRT hit reco Z (cm)
    int32_t    fT0Hit; ///< hit time w.r.t. global event time
    int32_t    fT1Hit; ///< hit time w.r.t. PPS
    //double    fT0CorrHit;
    //double    fT1CorrHit;
    int       fHitReg; ///< region code of CRT hit
    int       fHitSubSys;
    int       fNHit; ///< number of CRT hits for this event
    int       fHitTrk[64];
    int       fHitPDG[64];
    int       fHitStrip;
    int       fHitMod;

    // truth CRT hit vars
    int       fTrueHitEvent;
    double    fTrueXHit; ///< reconstructed X position of CRT hit (cm)
    double    fTrueYHit; ///< reconstructed Y position of CRT hit (cm)
    double    fTrueZHit; ///< reconstructed Z position of CRT hit (cm)
    double    fTrueXHitErr; ///< stat error of CRT hit reco X (cm)
    double    fTrueYHitErr; ///< stat error of CRT hit reco Y (cm)
    double    fTrueZHitErr; ///< stat error of CRT hit reco Z (cm)
    double    fTrueTHit; ///< hit time w.r.t. global event time
    double    fTrueTHitErr; ///< hit time w.r.t. PPS
    int       fTrueHitReg; ///< region code of CRT hit
    int       fTrueHitSubSys;
    //uint32_t fTrueNHit; ///< number of CRT hits for this event
    int       fTrueHitTrk;
    int       fTrueHitPDG;
    int       fTrueHitModID[2]; // one for c or d type, 2 entries for m type


    TH1F* fModMultHistC;   ///< true N C-modules hit / muon track
    TH1F* fModMultHistM;   ///< true N M-modules hit / muon track
    TH1F* fModMultHistD;   ///< true N D-modules hit / muon track

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
    fCosmicDisplayNtuple = tfs->make<TTree>("DisplayTree",      "track information for ROOT event display");
    fGenNtuple           = tfs->make<TTree>("GenTree",          "truth information from the generator");
    fSimulationNtuple    = tfs->make<TTree>("SimTree",          "MyCRTSimulation");
    fRegionsNtuple       = tfs->make<TTree>("RegTree",          "Info about particles crossing boundaries");
    fDetSimNtuple        = tfs->make<TTree>("DetTree",          "MyCRTDetSim");
    fSimHitNtuple        = tfs->make<TTree>("HitTree",          "MyCRTSimHit");
    fTrueCRTHitNtuple    = tfs->make<TTree>("TrueCRTHitTree",   "CRT hits from truth info");

    // Construct truth matching histograms
    //fStripMultHistC   = tfs->make<TH1F>("StripMultC",";no. strips hit / module / #mu;",64,0,64);
    //fStripMultHistM   = tfs->make<TH1F>("StripMultM",";no. strips hit / module / #mu;",64,0,64);
    //fStripMultHistD   = tfs->make<TH1F>("StripMultD",";no. strips hit / module / #mu;",64,0,64);
    fModMultHistC     = tfs->make<TH1F>("ModMultC",";no. modules hit / #mu;",10,0,10);
    fModMultHistM     = tfs->make<TH1F>("ModMultM",";no. modules hit / #mu;",10,0,10);
    fModMultHistD     = tfs->make<TH1F>("ModMultD",";no. modules hit / #mu;",10,0,10);

    fChanMultHistC    =tfs->make<TH1F>("ChanMultC",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fChanMultHistM    =tfs->make<TH1F>("ChanMultD",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fChanMultHistD    =tfs->make<TH1F>("ChanMultM",";no. FEB channels > threshold / FEB / #mu;",64,0,64);
    fFEBMultHistC     =tfs->make<TH1F>("FEBMultC",";no. FEB triggers / #mu;",64,0,64);
    fFEBMultHistM     =tfs->make<TH1F>("FEBMultD",";no. FEB triggers / #mu;",64,0,64);
    fFEBMultHistD     =tfs->make<TH1F>("FEBMultM",";no. FEB triggers / #mu;",64,0,64);

    // Define the branches of our event display n-tuple
    fCosmicDisplayNtuple->Branch("event",             &fCDEvent,             "event/I");
    fCosmicDisplayNtuple->Branch("trackID",           &fCDTrackID,           "trackID/I");
    fCosmicDisplayNtuple->Branch("nSeg",              &fNCD,                 "nSeg/I");
    fCosmicDisplayNtuple->Branch("regions",           fCDRegions,            "regions[4000]/I");
    fCosmicDisplayNtuple->Branch("pdg",               &fCDpdg,               "pdg/I");
    fCosmicDisplayNtuple->Branch("slopes",            fCDSlopes,             "slopes[4000][3]/D");
    fCosmicDisplayNtuple->Branch("pe",                fCDpe,                 "pe[4000][4]/D");
    fCosmicDisplayNtuple->Branch("xyzt",              fCDxyzt,               "xyzt[4000][4]/D");

    // Define the branches of our Gen n-tuple
    fGenNtuple->Branch("event",        &fEvent,         "event/I");
    fGenNtuple->Branch("nGen",         &fNGen,          "nGen/I");
    fGenNtuple->Branch("trackID",      fGenTrack,       "trackID[500]/I");
    fGenNtuple->Branch("pdg",          fGenPDG,         "pdg[500]/I");
    fGenNtuple->Branch("startXYZT",    fGenStartXYZT,   "startXYZT[500][4]/D");
    fGenNtuple->Branch("endXYZT",      fGenEndXYZT,     "endXYZT[500][4]/D");
    fGenNtuple->Branch("startPE",      fGenStartPE,     "startPE[500][4]/D");
    fGenNtuple->Branch("endPE",        fGenEndPE,       "endPE[500][4]/D");

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
    fRegionsNtuple->Branch("fiducial",             &fRegFid,             "fiducial/I");
    fRegionsNtuple->Branch("active",               &fRegActive,          "active/I");
    fRegionsNtuple->Branch("inactive",             &fRegInactive,        "inactive/I");
    fRegionsNtuple->Branch("crts",                 &fRegCRTs,            "crts/I");
    fRegionsNtuple->Branch("regions",              fRegRegions,          "regions[500]/I");
    fRegionsNtuple->Branch("pdg",                  &fRegPDG,             "pdg/I");
    fRegionsNtuple->Branch("trackID",              &fRegTrkID,           "trackID/I");
    fRegionsNtuple->Branch("eDep",                 fRegEDep,             "eDep[500]/D");
    fRegionsNtuple->Branch("dL",                   fRegdL,               "dL[500]/D");
    fRegionsNtuple->Branch("opDetID",              fRegOpDetID,          "opDetID[500]/I");
    fRegionsNtuple->Branch("distToOpDet",          fRegDistToOpDet,      "distToOpDet[500]/D");
    fRegionsNtuple->Branch("opDetXYZT",            fRegOpDetXYZT,        "opDetXYZT[500][4]/D");
    fRegionsNtuple->Branch("entryPE",              fRegEntryPE,          "entryPE[500][4]/D");
    fRegionsNtuple->Branch("exitPE",               fRegExitPE,           "exitPE[500][4]/D");
    fRegionsNtuple->Branch("entryXYZT",            fRegEntryXYZT,        "entryXYZT[500][4]/D");
    fRegionsNtuple->Branch("exitXYZT",             fRegExitXYZT,         "exitXYZT[500][4]/D");
    fRegionsNtuple->Branch("entrySlope",           fRegEntrySlope,       "entrySlope[500][3]/D");
    fRegionsNtuple->Branch("exitSlope",            fRegExitSlope,        "exitSlope[500][3]/D");

    // Define the branches of our DetSim n-tuple 
    fDetSimNtuple->Branch("event",                 &fDetEvent,          "event/I");
    fDetSimNtuple->Branch("nChan",                 &fNChan,             "nChan/I");
    fDetSimNtuple->Branch("channel",               fChan,               "channel[64]/I");
    fDetSimNtuple->Branch("t0",                    fT0,                 "t0[64]/D");
    fDetSimNtuple->Branch("t1",                    fT1,                 "t1[64]/D");
    fDetSimNtuple->Branch("adc",                   fADC,                "adc[64]/I");
    fDetSimNtuple->Branch("trackID",               fTrackID,            "trackID[64]/I");
    fDetSimNtuple->Branch("detPDG",                fDetPDG,             "detPDG[64]/I");
    fDetSimNtuple->Branch("entry",                 &fEntry,             "entry/I");
    fDetSimNtuple->Branch("mac5",                  &fMac5,              "mac5/I");
    fDetSimNtuple->Branch("region",                &fFEBReg,            "region/I");
    fDetSimNtuple->Branch("triggerPair",           fTriggerPair,        "triggerPair[2]/I");
    fDetSimNtuple->Branch("macPair",               fMacPair,            "macPair[2]/I");
    fDetSimNtuple->Branch("chanTrig",              &fChanTrig,          "chanTrig/I");
    fDetSimNtuple->Branch("tTrig",                 &fTTrig,             "tTrig/D");
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
    fSimHitNtuple->Branch("t0",          &fT0Hit,       "t0/I");
    fSimHitNtuple->Branch("t1",          &fT1Hit,       "t1/I");
    //fSimHitNtuple->Branch("t0Corr",      &fT0CorrHit,   "t0/D");
    //fSimHitNtuple->Branch("t1Corr",      &fT1CorrHit,   "t1/D");
    fSimHitNtuple->Branch("region",      &fHitReg,      "region/I");  
    fSimHitNtuple->Branch("subSys",      &fHitSubSys,   "subSys/I");
    fSimHitNtuple->Branch("trackID",     &fHitTrk,      "trackID[64]/I");
    fSimHitNtuple->Branch("pdg",         &fHitPDG,      "pdg[64]/I");
    fSimHitNtuple->Branch("modID",       &fHitMod,      "modID/I");
    fSimHitNtuple->Branch("stripID",     &fHitStrip,    "stripID/I");

    // Define the branches of our SimTrueHit n-tuple
    fTrueCRTHitNtuple->Branch("event",       &fTrueHitEvent,    "event/I");
    //fSimHitNtuple->Branch("nHit",        &fTrueNHit,        "nHit/I");
    fTrueCRTHitNtuple->Branch("x",           &fTrueXHit,        "x/D");
    fTrueCRTHitNtuple->Branch("y",           &fTrueYHit,        "y/D");
    fTrueCRTHitNtuple->Branch("z",           &fTrueZHit,        "z/D");
    fTrueCRTHitNtuple->Branch("xErr",        &fTrueXHitErr,     "xErr/D");
    fTrueCRTHitNtuple->Branch("yErr",        &fTrueYHitErr,     "yErr/D");
    fTrueCRTHitNtuple->Branch("zErr",        &fTrueZHitErr,     "zErr/D");
    fTrueCRTHitNtuple->Branch("t",           &fTrueTHit,        "t/D");
    fTrueCRTHitNtuple->Branch("tErr",        &fTrueTHitErr,     "tErr/D");
    fTrueCRTHitNtuple->Branch("region",      &fTrueHitReg,      "region/I");
    fTrueCRTHitNtuple->Branch("subSys",      &fTrueHitSubSys,   "subSys/I");
    fTrueCRTHitNtuple->Branch("trackID",     &fTrueHitTrk,      "trackID/I");
    fTrueCRTHitNtuple->Branch("pdg",         &fTrueHitPDG,      "pdg/I");
    fTrueCRTHitNtuple->Branch("modID",       fTrueHitModID,     "modID[2]/I");

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
    MF_LOG_DEBUG("CRT") << "beginning analyis" << '\n';

    // Check that lists for momenta limits is same size as last of PDGs from FHiCL
    if (fPDGs.size() != fMinMomenta.size() || fPDGs.size() != fMaxMomenta.size())
        throw cet::exception("CRTAnalysis")
          << " PDG/Momtenta values not set correctly in fhicl - lists have different sizes"
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;


    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    // Define "handle" to Generator level MCTruth objects
    art::Handle< vector<simb::MCTruth>> genHandle;

    // Define a "handle" to point to a vector of MCParticle objects.
    art::Handle< vector<simb::MCParticle> > particleHandle;
    map< int, const simb::MCParticle*> particleMap;

    if (!event.getByLabel("generator", genHandle)) {
        std::cout << "could not get handle to gen objects!!!" << std::endl;
    }

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

    // Handle to AuxDetSimChannel (CRT module) objects generated by LArG4
    art::Handle<vector<sim::AuxDetSimChannel> > auxDetSimChannelHandle;
    if (!event.getByLabel(fSimulationProducerLabel, auxDetSimChannelHandle)) {
        throw cet::exception("CRTAnalysis")
          << " No sim::AuxDetSimChannel objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    map<int,vector<pair<int,int>>> febMap;
    FillFebMap(febMap);
    std::cout << "in analyze: filled febMap with " << febMap.size() << " modules!" << std::endl;

    for ( auto const& truth : (*genHandle) )
    {
        fNGen = 0;
        for ( int i=0; i< kMaxGen; i++ ) {
            fGenTrack[i] = -0.5*INT_MAX;
            fGenPDG[i] = -0.5*INT_MAX;
            for (int j=0; j<4; j++) {
                fGenStartXYZT[i][j] = -0.5*DBL_MAX;
                fGenEndXYZT[i][j] = -0.5*DBL_MAX;
                fGenStartPE[i][j] = -0.5*DBL_MAX;
                fGenEndPE[i][j] = -0.5*DBL_MAX;
            }
        }
        fNGen = truth.NParticles();
        //std::cout << "loop over truth for event : " << fEvent << " with " << fNGen << " particles" << std::endl;
        for ( int i=0; i<fNGen; i++ )
        {
            if (fNGen==kMaxGen) {break; std::cout << "max out gen!" << std::endl;}

            auto const& part = truth.GetParticle(i); //simb::MCParticle
            fGenTrack[i] = part.TrackId();
            fGenPDG[i] = part.PdgCode();

            //std::cout <<" i: " << i << ", trackID: " << part.TrackId() << ", pdg: " << part.PdgCode() << std::endl;

            const TLorentzVector startPos = part.Position(0);
            const TLorentzVector endPos = part.EndPosition();
            startPos.GetXYZT(fGenStartXYZT[i]);
            endPos.GetXYZT(fGenEndXYZT[i]);

            const TLorentzVector startMom = part.Momentum(0);
            const TLorentzVector endMom = part.EndMomentum();
            startMom.GetXYZT(fGenStartPE[i]);
            endMom.GetXYZT(fGenEndPE[i]);
        }

        fGenNtuple->Fill();
    }

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

    std::cout << "event " << fEvent << " with " << particleMap.size() << " MCParticles" << std::endl;

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
        //if (abs(fSimPDG)==13){
        //    fNmuTruth++;
       // }
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

        fCDEvent = fEvent;
        fNCD = 0;
        fCDpdg = fSimPDG;
        fCDTrackID = fSimTrackID;

        fRegEvent = fEvent;
        fNReg = 0;
        fRegFid = 0;
        fRegActive = 0;
        fRegInactive = 0;
        fRegCRTs = 0;
        fRegPDG = fSimPDG;;
        fRegTrkID = fSimTrackID;
        //std::cout << "initializing region arrays" << std::endl;

        for (size_t i=0; i<kMaxSeg; i++){
            fCDRegions[i] = -0.5*INT_MAX;
            for (int j=0; j<4; j++) {
                if (j<3) fCDSlopes[i][j] = -0.5*DBL_MAX;
                fCDpe[i][j] = -0.5*DBL_MAX;
                fCDxyzt[i][j] = -0.5*DBL_MAX;
            }
        }

        //initialize arrays
        for (size_t i=0; i<kMaxReg; i++){
                fRegRegions[i] = -0.5*INT_MAX;
                fRegEDep[i] = -0.5*DBL_MAX;
                fRegDistToOpDet[i] = -0.5*DBL_MAX;
                fRegOpDetID[i] = -0.5*INT_MAX;
                for (int j=0; j<4; j++){
                        fRegEntryXYZT[i][j] = -0.5*DBL_MAX;
                        fRegExitXYZT[i][j] = -0.5*DBL_MAX;
                        fRegEntryPE[i][j] = -0.5*DBL_MAX;
                        fRegExitPE[i][j] = -0.5*DBL_MAX;
                        fRegOpDetXYZT[i][j] = -0.5*DBL_MAX;
                        if (j<3) {
                            fRegEntrySlope[i][j] = -0.5*DBL_MAX;
                            fRegExitSlope[i][j] = -0.5*DBL_MAX;
                        }
                }
        }

        int oldreg = -1;

        MF_LOG_DEBUG("CRT") << "about to loop over trajectory points" << '\n';

        //loop over trajectory points
        for (unsigned int i=0; i<fSimHits; i++){
                const TLorentzVector& pos = particle.Position(i); // 4-position in World coordinates
                const TLorentzVector& posnext = particle.Position(i+1); // problem for last point???
                const TLorentzVector& mom = particle.Momentum(i); // 4-momentum
                const double point[3] = {pos.X(),pos.Y(),pos.Z()};
                const double pointnext[3] = {posnext.X(),posnext.Y(),posnext.Z()};
                double opDetPos[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
                double entryPos[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
                double entryT = -FLT_MAX;
                bool active0 = false, active1 = false, activenext0 = false, activenext1 = false;

                // CosmicDisplay info
                pos.GetXYZT(fCDxyzt[fNCD]); // 4-position
                mom.GetXYZT(fCDpe[fNCD]); // 4-momentum
                fCDSlopes[fNCD][0] = mom.Px()/mom.P(); // direction cosines
                fCDSlopes[fNCD][1] = mom.Py()/mom.P(); 
                fCDSlopes[fNCD][2] = mom.Pz()/mom.P();
                fNCD++;

                if(fNReg==kMaxReg) std::cout << "about to seg fault..need more NReg!" << std::endl;

                // Regions info
                // Check if trajectory points are in cryostats (active + inactve LAr ) 
                if(cryo0.ContainsPosition(point)) {

                        active0 = tpc00.ContainsPosition(point);
                        active1 = tpc01.ContainsPosition(point);
                        activenext0 = tpc00.ContainsPosition(pointnext);
                        activenext1 = tpc01.ContainsPosition(pointnext);

                        // if last point was not in this cryostat or is now entering AV
                        if ( (oldreg!=10&&!active0&&!active1) || (active0&&oldreg!=5) || (active1&&oldreg!=6)) {
                            pos.GetXYZT(fRegEntryXYZT[fNReg]);
                            mom.GetXYZT(fRegEntryPE[fNReg]);
                            oldreg = 10;
                            if (active0) oldreg = 5;
                            if (active1) oldreg = 6;
                        }

                        // if next point is outside of this volume or is last traj. point
                        if (!cryo0.ContainsPosition(pointnext) || (oldreg==10&&(activenext0||activenext1))
                            || i==fSimHits-1
                            || (active0 && !activenext0) || (active1&&!activenext1) ){

                                pos.GetXYZT(fRegExitXYZT[fNReg]);
                                mom.GetXYZT(fRegExitPE[fNReg]);
                                if (active0) {
                                    fRegRegions[fNReg] = 5;
                                    fRegActive++;
                                    if(tpc00.InFiducialX(point[0],25,0) && tpc00.InFiducialY(point[1],25,25)
                                      && tpc00.InFiducialZ(point[2],30,50)) 
                                        fRegFid++;

                                }
                                else if (active1) {
                                    fRegRegions[fNReg] = 6;
                                    fRegActive++;
                                    if(tpc01.InFiducialX(point[0],25,0) && tpc01.InFiducialY(point[1],25,25)
                                      && tpc01.InFiducialZ(point[2],30,50))  
                                        fRegFid++;
                                }
                                else {
                                    fRegRegions[fNReg] = 10;
                                    fRegInactive++;
                                }
                                fRegdL[fNReg] = sqrt(pow(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                    +pow(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                    +pow(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
                                fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
                                for (int index=0; index<3; index++) entryPos[index] = fRegEntryXYZT[fNReg][index];
                                entryT = fRegEntryXYZT[fNReg][3];
                                fRegOpDetID[fNReg] = cryo0.GetClosestOpDet(entryPos);
                                geo::OpDetGeo const& opDet0 = cryo0.OpDet(fRegOpDetID[fNReg]);
                                opDet0.GetCenter(opDetPos);
                                fRegDistToOpDet[fNReg] = sqrt(pow(opDetPos[0]-entryPos[0],2)
                                                            + pow(opDetPos[1]-entryPos[1],2)
                                                            + pow(opDetPos[2]-entryPos[2],2));
                                for (int index=0; index<3; index++) fRegOpDetXYZT[fNReg][index] = opDetPos[index];
                                fRegOpDetXYZT[fNReg][3] = entryT + fRegDistToOpDet[fNReg]*LAR_PROP_DELAY;
                                fNReg++;
                        }
                } //if cryo0

                // if this point in the other cryostat
                if(cryo1.ContainsPosition(point)) {
 
                        //check if this or next points are in active volumes
                        active0 = tpc10.ContainsPosition(point);
                        active1 = tpc11.ContainsPosition(point);
                        activenext0 = tpc10.ContainsPosition(pointnext);
                        activenext1 = tpc11.ContainsPosition(pointnext);

                        // if last point was not in this cryostat or is now entering AV
                        if ( (oldreg!=12&&!active0&&!active1) || (active0&&oldreg!=7) || (active1&&oldreg!=8)) {
                            pos.GetXYZT(fRegEntryXYZT[fNReg]);
                            mom.GetXYZT(fRegEntryPE[fNReg]);
                            oldreg = 12;
                            if (active0) oldreg = 7;
                            if (active1) oldreg = 8;
                        }

                        if (!cryo1.ContainsPosition(pointnext) || (oldreg==12&&(activenext0||activenext1))
                            || i==fSimHits-1
                            || (active0 && !activenext0) || (active1&&!activenext1) ){

                                pos.GetXYZT(fRegExitXYZT[fNReg]);
                                mom.GetXYZT(fRegExitPE[fNReg]);
                                if (active0) {
                                    fRegRegions[fNReg] = 5;
                                    fRegActive++;
                                    if(tpc00.InFiducialX(point[0],25,0) && tpc00.InFiducialY(point[1],25,25)
                                      && tpc00.InFiducialZ(point[2],30,50))
                                        fRegFid++;

                                }
                                else if (active1) {
                                    fRegRegions[fNReg] = 6;
                                    fRegActive++;
                                    if(tpc01.InFiducialX(point[0],25,0) && tpc01.InFiducialY(point[1],25,25)
                                      && tpc01.InFiducialZ(point[2],30,50))
                                        fRegFid++;
                                }
                                else {
                                    fRegRegions[fNReg] = 10;
                                    fRegInactive++;
                                }
                                fRegdL[fNReg] = sqrt(pow(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                    +pow(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                    +pow(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
                                fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
                                for (int index=0; index<3; index++) entryPos[index] = fRegEntryXYZT[fNReg][index];
                                entryT = fRegEntryXYZT[fNReg][3];
                                fRegOpDetID[fNReg] = cryo0.GetClosestOpDet(entryPos);
                                geo::OpDetGeo const& opDet0 = cryo0.OpDet(fRegOpDetID[fNReg]);
                                opDet0.GetCenter(opDetPos);
                                fRegDistToOpDet[fNReg] = sqrt(pow(opDetPos[0]-entryPos[0],2)
                                                            + pow(opDetPos[1]-entryPos[1],2)
                                                            + pow(opDetPos[2]-entryPos[2],2));
                                for (int index=0; index<3; index++) fRegOpDetXYZT[fNReg][index] = opDetPos[index];
                                fRegOpDetXYZT[fNReg][3] = entryT + fRegDistToOpDet[fNReg]*LAR_PROP_DELAY;
                                fNReg++;
                        } // if exiting from volume
                } //if cryo1

        }//for trajectory points

        fCosmicDisplayNtuple->Fill();
        //fRegionsNtuple->Fill();

        //map module IDs to strip IDs hit by muons
        //map< uint16_t,set<uint8_t>* > muHitMapC; //hits in C modules only
        //map< uint16_t,set<uint8_t>* > muHitMapM; //hits in M modules only
        //map< uint16_t,set<uint8_t>* > muHitMapD; //hits in D modules only
        //map< uint16_t,set<uint8_t>* > muHitMap; //all hits

        map< int, vector<double> > regCRTEnter, regCRTExit;

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

        struct tagger {
            char type;
            int region;
            std::set<int> layerID;
            std::map<int,int> stripLayer;
            std::pair<int,int> modPair;
            std::map<int,std::vector<double>> xyzt;
        };
           
        std::map<int,tagger> taggers;    

	// To look at the energy deposited by this particle's track,
	// we loop over the AuxDetSimChannel objects in the event. 
	// Note all volumes are included, not just ones with energy deps
	for ( auto const& channel : (*auxDetSimChannelHandle) )
	{  
	    // Get vector of hits in this AuxDet channel
	    auto const& auxDetIDEs = channel.AuxDetIDEs();

	    // For every hit in this channel:
	    for ( auto const& ide : auxDetIDEs )
	    {
		    // Check if the track that deposited the
		    // energy matches the track of the MCParticle.
		    if ( ide.trackID != fSimTrackID ) continue;
		    if ( ide.energyDeposited * 1.0e6 < 50 ) continue; 

                    tagger& tag = taggers[channel.AuxDetID()];
                    auto const& adGeo = fGeometryService->AuxDet(channel.AuxDetID());
                    tag.type = ModToAuxDetType(adGeo);
                    tag.region = GetAuxDetRegion(adGeo);
                    auto const& adsGeo = adGeo.SensitiveVolume(channel.AuxDetSensitiveID());

                    std::set<std::string> volNames = { adsGeo.TotalVolume()->GetName() };
                    std::vector<std::vector<TGeoNode const*> > paths = fGeometryService->FindAllVolumePaths(volNames);

                    std::string path = "";
                    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
                        path += paths.at(0).at(inode)->GetName();
                        if (inode < paths.at(0).size() - 1) {
                            path += "/";
                        }
                    }
                    TGeoManager* manager = fGeometryService->ROOTGeoManager();
                    manager->cd(path.c_str());
                    TGeoNode* nodeStrip = manager->GetCurrentNode();
                    TGeoNode* nodeInner = manager->GetMother(1);
                    TGeoNode* nodeModule = manager->GetMother(2);
                    double origin[3] = {0, 0, 0};
                    int layid = 0.5*INT_MAX; //set to 0 or 1 if layerid determined
                    //int mac5=LONG_MAX;
                    // Module position in parent (tagger) frame
                    double modulePosMother[3]; //position in CRT region volume
                    nodeModule->LocalToMaster(origin, modulePosMother);

                    // strip position in module frame
                    double stripPosMother[3];
                    double stripPosModule[3];
                    nodeStrip->LocalToMaster(origin, stripPosMother);
                    nodeInner->LocalToMaster(stripPosMother,stripPosModule);

                    fADType[fNAuxDet] = ModToTypeCode(adGeo);
                    fAuxDetReg[fNAuxDet] = GetAuxDetRegion(adGeo);

                    if ( fADType[fNAuxDet] == 0 || fADType[fNAuxDet] == 2 )
                        layid = (stripPosModule[1] > 0);

                    // if 'm' type
                    if ( fADType[fNAuxDet] == 1 ) {
                        // if east or west stacks (6 in total)
                        if ( fAuxDetReg[fNAuxDet] >=40 && fAuxDetReg[fNAuxDet] <=45 ) {
                            layid = ( abs(modulePosMother[0]>0) );
                        }
                        // if front or back
                        if ( fAuxDetReg[fNAuxDet] == 46 || fAuxDetReg[fNAuxDet] == 47) {
                            layid = ( modulePosMother[2]> 0 );
                        }
                    }

                    // What is the distance from the hit (centroid of the entry
                    // and exit points) to the readout end?
                    double trueX = (ide.entryX + ide.exitX) / 2.0;
                    double trueY = (ide.entryY + ide.exitY) / 2.0;
                    double trueZ = (ide.entryZ + ide.exitZ) / 2.0;
                    double trueT = (ide.entryT + ide.exitT) / 2.0;
                    double world[3] = {trueX, trueY, trueZ};
                    double svHitPosLocal[3];
                    double modHitPosLocal[3];
                    adsGeo.WorldToLocal(world, svHitPosLocal); //position in strip frame  (origin at center)
                    adGeo.WorldToLocal(world, modHitPosLocal); //position in module frame (origin at center)

                    if (trueX==0) {
                              std::cout << "hit in " << tag.type << "-module " << channel.AuxDetID() << ", strip "
                              << channel.AuxDetSensitiveID() << ", in layer " << layid << '\n'
                              << "  x: " << trueX << ", y: " << trueY << ", z: " << trueZ 
                              << ", t: " << trueT << std::endl;
                    }

                    tag.layerID.insert(layid);
                    tag.stripLayer[channel.AuxDetSensitiveID()] = layid;
                    std::vector<double>& truePos = tag.xyzt[channel.AuxDetSensitiveID()];
                    truePos.push_back(trueX);
                    truePos.push_back(trueY);
                    truePos.push_back(trueZ);
                    truePos.push_back(trueT);
                    //tag.modPair = std::make_pair(channel.AuxDetID(),channel.AuxDetID());

                    //std::cout << " tagger xyzt length: " << truePos.size() << std::endl;

                    //calculate track length in strip
		    double dx = ide.entryX-ide.exitX;
		    double dy = ide.entryY-ide.exitY;
		    double dz = ide.entryZ-ide.exitZ;
                    double adlength = sqrt(dx*dx+dy*dy+dz*dz);
                    if ( adlength < 0.0001)  continue;

		    fADTrackLength[fNAuxDet] = sqrt(dx*dx+dy*dy+dz*dz);
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
                    fADMac[fNAuxDet] = (ADToMac(febMap,channel.AuxDetID())).first;
                    fADType[fNAuxDet] = ModToTypeCode(adGeo);
		    //fNAuxDet++;

                    vector<double> vtmp = {ide.entryX,ide.entryY,ide.entryZ,ide.entryT};
                    if (regCRTEnter.find(fAuxDetReg[fNAuxDet])!=regCRTEnter.end()) {
                        if (regCRTEnter[fAuxDetReg[fNAuxDet]][3] > ide.entryT) {
                            //regCRTEnter[fAuxDetReg[fNAuxDet]].clear();
                            regCRTEnter[fAuxDetReg[fNAuxDet]] = vtmp;
                        }
                    }
                    else regCRTEnter[fAuxDetReg[fNAuxDet]] = vtmp;

                    vtmp = {ide.exitX,ide.exitY,ide.exitZ,ide.exitT};
                    if (regCRTExit.find(fAuxDetReg[fNAuxDet])!=regCRTExit.end()) {
                        if (regCRTExit[fAuxDetReg[fNAuxDet]][3] < ide.exitT)
                            regCRTExit[fAuxDetReg[fNAuxDet]] = vtmp;
                    }
                    else regCRTExit[fAuxDetReg[fNAuxDet]] = vtmp;

                    fNAuxDet++;

	    } // For each IDE (strip hit by muon)
              
	} // For each SimChannel (module)

        // write values to tree for this event and particle
        fSimulationNtuple->Fill();

        fTrueHitEvent = fEvent;
        fTrueHitTrk = fSimTrackID;
        fTrueHitPDG = fSimPDG;
        map <int,tagger> mTaggers;

        // loop over taggers: modID->strip hit info
        for (auto const& tag : taggers) {

            fTrueXHit = 0;
            fTrueYHit = 0;
            fTrueZHit = 0;
            fTrueTHit = 0;
            fTrueXHitErr = 0;
            fTrueYHitErr = 0;
            fTrueZHitErr = 0;
            fTrueTHitErr = 0;
            fTrueHitModID[0] = tag.first;
            fTrueHitModID[1] = tag.first;
            fTrueHitReg = tag.second.region;
            if (tag.second.type == 'c') fTrueHitSubSys = 0;
            if (tag.second.type == 'd') fTrueHitSubSys = 2;
            std::vector<double> xhits, yhits, zhits, thits;

            // if c or d type module
            if (tag.second.type=='c' || tag.second.type=='d') {
                // if "X-Y" coincidence
                if (tag.second.layerID.size()>1) {
                    // loop over module strips map: stripID->pos 4-vec
                    for (auto const& strip : tag.second.xyzt) {
                        if (strip.second[0] == 0
                         && strip.second[1] == 0
                         && strip.second[2] == 0)
                           continue;

                        xhits.push_back(strip.second[0]);
                        yhits.push_back(strip.second[1]);
                        zhits.push_back(strip.second[2]);
                        thits.push_back(strip.second[3]);
                    }
                    if (xhits.size()==0) continue;

                    for ( size_t i=0; i<xhits.size(); i++) {
                        fTrueXHit+=xhits[i];
                        fTrueYHit+=yhits[i];
                        fTrueZHit+=zhits[i];
                        fTrueTHit+=thits[i];
                    }
                    
                    fTrueXHit *= 1.0/tag.second.xyzt.size();
                    fTrueYHit *= 1.0/tag.second.xyzt.size();
                    fTrueZHit *= 1.0/tag.second.xyzt.size();
                    fTrueTHit *= 1.0/tag.second.xyzt.size();
                      
                    for ( size_t i=0; i<xhits.size(); i++) {
                        fTrueXHitErr += pow(xhits[i]-fTrueXHit,2);
                        fTrueYHitErr += pow(yhits[i]-fTrueYHit,2);
                        fTrueZHitErr += pow(zhits[i]-fTrueZHit,2);
                        fTrueTHitErr += pow(thits[i]-fTrueTHit,2);
                    }

                    fTrueXHitErr = sqrt(fTrueXHitErr/(xhits.size()-1));
                    fTrueYHitErr = sqrt(fTrueYHitErr/(yhits.size()-1));
                    fTrueZHitErr = sqrt(fTrueZHitErr/(zhits.size()-1));
                    fTrueTHitErr = sqrt(fTrueTHitErr/(thits.size()-1));
           
                } //if coincidence
                else continue;

                fTrueCRTHitNtuple->Fill();

            }//if c or d type

            if ( tag.second.type=='m' ) {
                mTaggers[tag.first] = tag.second;
            }

        } //loop over taggers

        set <int> mPairs;
        int nmisspair = 0;
        for (auto const& tag : mTaggers) {

            fTrueXHit = 0;
            fTrueYHit = 0;
            fTrueZHit = 0;
            fTrueTHit = 0;
            fTrueXHitErr = 0;
            fTrueYHitErr = 0;
            fTrueZHitErr = 0;
            fTrueTHitErr = 0;
            fTrueHitModID[0] = tag.first;
            fTrueHitModID[1] = tag.first;
            fTrueHitReg = tag.second.region;
            fTrueHitSubSys = 1;
            std::vector<double> xhits, yhits, zhits, thits;
            bool pairFound = false;

            if (mPairs.find(tag.first) != mPairs.end()) continue;
            for (auto const& tag2 : mTaggers) {
                if (tag.first != tag2.first &&
                  mPairs.find(tag2.first) == mPairs.end() &&
                  tag.second.region == tag2.second.region &&
                  tag.second.layerID != tag2.second.layerID ) {

                    mPairs.insert(tag.first);
                    mPairs.insert(tag2.first);
                    fTrueHitModID[0] = tag.first;
                    fTrueHitModID[1] = tag2.first;

                    // loop over module strips map: stripID->pos 4-vec
                    for (auto const& strip : tag.second.xyzt) {
                        if (strip.second[0] == 0
                         && strip.second[1] == 0
                         && strip.second[2] == 0)
                           continue;

                        xhits.push_back(strip.second[0]);
                        yhits.push_back(strip.second[1]);
                        zhits.push_back(strip.second[2]);
                        thits.push_back(strip.second[3]);
                    }// for xyzt in first tagger (module)

                    for (auto const& strip : tag2.second.xyzt) {
                        if (strip.second[0] == 0
                         && strip.second[1] == 0
                         && strip.second[2] == 0)
                           continue;

                        xhits.push_back(strip.second[0]);
                        yhits.push_back(strip.second[1]);
                        zhits.push_back(strip.second[2]);
                        thits.push_back(strip.second[3]);
                    }// for xyzt in second tagger (module)


                    for ( size_t i=0; i<xhits.size(); i++) {
                        fTrueXHit+=xhits[i];
                        fTrueYHit+=yhits[i];
                        fTrueZHit+=zhits[i];
                        fTrueTHit+=thits[i];
                    } // sum together x,y,z,t values repectively

                    //get the centroid
                    fTrueXHit *= 1.0/xhits.size();
                    fTrueYHit *= 1.0/yhits.size();
                    fTrueZHit *= 1.0/zhits.size();
                    fTrueTHit *= 1.0/thits.size();

                    // calculate rms of hit distrubution
                    for ( size_t i=0; i<xhits.size(); i++) {
                        fTrueXHitErr += pow(xhits[i]-fTrueXHit,2);
                        fTrueYHitErr += pow(yhits[i]-fTrueYHit,2);
                        fTrueZHitErr += pow(zhits[i]-fTrueZHit,2);
                        fTrueTHitErr += pow(thits[i]-fTrueTHit,2);
                    }

                    fTrueXHitErr = sqrt(fTrueXHitErr/(xhits.size()-1));
                    fTrueYHitErr = sqrt(fTrueYHitErr/(yhits.size()-1));
                    fTrueZHitErr = sqrt(fTrueZHitErr/(zhits.size()-1));
                    fTrueTHitErr = sqrt(fTrueTHitErr/(thits.size()-1));
 
                    pairFound = true;             
                    break; // first tagger is matched so break from loop
                }//if match found
            }//inner loop over taggers

            if (pairFound) fTrueCRTHitNtuple->Fill();
            else nmisspair++;

        } // outer loop over taggers

        if (nmisspair>0) std::cout << "missed " << nmisspair << " tagger pairs in trueHit reco for M mods" << std::endl;

        //auto tmpNReg = fNReg;
        for( auto it=regCRTEnter.begin(); it!=regCRTEnter.end(); it++) {
            //std::cout << "found CRT region " << it->first << std::endl;
            fRegRegions[fNReg] = it->first;
            fRegEntryXYZT[fNReg][0] = (it->second)[0];
            fRegEntryXYZT[fNReg][1] = (it->second)[1];
            fRegEntryXYZT[fNReg][2] = (it->second)[2];
            fRegEntryXYZT[fNReg][3] = (it->second)[3];
            fRegExitXYZT[fNReg][0] = regCRTExit[it->first][0];
            fRegExitXYZT[fNReg][1] = regCRTExit[it->first][1];
            fRegExitXYZT[fNReg][2] = regCRTExit[it->first][2];
            fRegExitXYZT[fNReg][3] = regCRTExit[it->first][3];
            fRegdL[fNReg] = sqrt(pow(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                       +pow(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                       +pow(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2));
            fRegEDep[fNReg] = fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3];
            fNReg++; fRegCRTs++;
        }
        //std::cout << "added " << fNReg-tmpNReg << " CRT regions to RegTree" << std::endl;

        //sort region tree entries by entry time
        int flag = 1;    // set flag to 1 to start first pass
        int      tempInt;
        double   tempDoub;
        
        for(int i = 1; (i <= (int)fNReg) && flag; i++)
        {
            flag = 0;
            for (int j=0; j < ((int)fNReg -1); j++)
            {
                if (fRegEntryXYZT[j+1][3] < fRegEntryXYZT[j][3]) 
                { 
                    tempInt = fRegRegions[j];             // swap regions
                    fRegRegions[j] = fRegRegions[j+1];
                    fRegRegions[j+1] = tempInt;
                    tempDoub = fRegEDep[j];             // swap EDep
                    fRegEDep[j] = fRegEDep[j+1];
                    fRegEDep[j+1] = tempDoub;
                    tempDoub = fRegdL[j];             // swap dL
                    fRegdL[j] = fRegdL[j+1];
                    fRegdL[j+1] = tempDoub;
                    tempDoub = fRegDistToOpDet[j];    //swap distToOpDet
                    fRegDistToOpDet[j] = fRegDistToOpDet[j+1];
                    fRegDistToOpDet[j+1] = tempDoub;
                    tempInt = fRegOpDetID[j];          //swap opDetID
                    fRegOpDetID[j] = fRegOpDetID[j+1];
                    fRegOpDetID[j+1] = tempInt;
                    for (int k=0; k<4; k++) {
                        tempDoub = fRegEntryPE[j][k];   //swap entryPE
                        fRegEntryPE[j][k] = fRegEntryPE[j+1][k];
                        fRegEntryPE[j+1][k] = tempDoub;
                        tempDoub  = fRegExitPE[j][k]; //swap exitPE
                        fRegExitPE[j][k] = fRegExitPE[j+1][k];
                        fRegExitPE[j+1][k] = tempDoub;
                        tempDoub = fRegEntryXYZT[j][k]; //swap entryXYZT
                        fRegEntryXYZT[j][k] = fRegEntryXYZT[j+1][k];
                        fRegEntryXYZT[j+1][k] = tempDoub;
                        tempDoub = fRegExitXYZT[j][k]; //swap exitXYZT
                        fRegExitXYZT[j][k] = fRegExitXYZT[j+1][k];
                        fRegExitXYZT[j+1][k] = tempDoub;
                        tempDoub = fRegOpDetXYZT[j][k];  //swap opDetXYZT
                        fRegOpDetXYZT[j][k] = fRegOpDetXYZT[j+1][k];
                        fRegOpDetXYZT[j+1][k] = tempDoub;
                        if(k<3) {
                            tempDoub = fRegEntrySlope[j][k];
                            fRegEntrySlope[j][k] = fRegEntrySlope[j+1][k];
                            fRegEntrySlope[j+1][k] = tempDoub;
                            tempDoub = fRegExitSlope[j][k];
                            fRegExitSlope[j][k] = fRegExitSlope[j+1][k];
                            fRegExitSlope[j+1][k] = tempDoub;
                        }
                    }
                    flag = 1;               // indicates that a swap occurred.
                }
            }
        }

        fRegionsNtuple->Fill();

    } // loop over all particles in the event. 


    art::Handle<vector<icarus::crt::CRTData>> crtDetSimHandle;
    bool isCRTDetSim = event.getByLabel(fCRTDetSimProducerLabel, crtDetSimHandle);

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
          fChan[i] = -0.8*INT_MAX;
          fT0[i]   = -0.8*INT_MAX;
          fT1[i]   = -0.8*INT_MAX;
          fADC[i]  = -0.8*INT_MAX;
          fTrackID[i] = -0.8*INT_MAX;
          fDetPDG[i] = -0.8*INT_MAX;
        }

 
        vector<int> missedIDs;
        //std::cout << "loop over chandata" << std::endl;
        for ( auto const chandat : febdat.ChanData()) {
          //DetSim tree contains all entries (not neccessarily from muons)
          fChan[fNChan]    = chandat.Channel();
          fT0[fNChan]      = chandat.T0();
          fT1[fNChan]      = chandat.T1();
          fADC[fNChan]     = chandat.ADC();
          fTrackID[fNChan] = chandat.TrackID()[0]; 
          if (particleMap.find(fTrackID[fNChan]) != particleMap.end() )
              fDetPDG[fNChan]  = particleMap[fTrackID[fNChan]]->PdgCode();
          else {
              fDetPDG[fNChan] *= -1;
              missedIDs.push_back(fTrackID[fNChan]);
          }

         fNChan++;
        }//outer chandat loop


        if ( missedIDs.size() > 0 ) {
            std::cout 
                << " couldn't match " << missedIDs.size() << " trackIDs from DetSim to MCParticle:" 
            << std::endl;
            for (auto const& id : missedIDs) std::cout << "  - " << id << std::endl;
        }

        fDetSimNtuple->Fill();

     } //for CRT FEB events

    }//if crtdetsim products present

    else std::cout << "CRTData products not found! (expected if gen/G4 step)" << std::endl;

    art::Handle<std::vector<icarus::crt::CRTHit>> crtSimHitHandle;
    //std::vector<art::Ptr<icarus::crt::CRTHit>> crtSimHits;
    
    bool isCRTSimHit = event.getByLabel(fCRTSimHitProducerLabel, crtSimHitHandle);
    std::vector<int> ids;
    if (isCRTSimHit) {

        std::cout << "Found " << crtSimHitHandle->size() << " CRT Hits" << std::endl;
        //art::fill_ptr_vector(crtSimHits,crtSimHitHandle);
        art::FindManyP<icarus::crt::CRTData> findManyData(crtSimHitHandle, event, fCRTSimHitProducerLabel);
        std::vector<art::Ptr<icarus::crt::CRTData>> data = findManyData.at(0);

        for(auto const& dat : data){
          for(auto const& chan : dat->ChanData()){
            for(auto const& trkid : chan.TrackID()){
              ids.push_back(trkid);
            }
          }
        }

        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

        //auto trks = icarus::CRTTruthMatchUtils::AllTrueIds(crtSimHitHandle,event,fCRTSimHitProducerLabel,0);
    }

    fNHit = 0;

    if (isCRTSimHit) {
        for ( auto const& hit : *crtSimHitHandle )
        {
            fNHit++;
            fHitEvent = fEvent;
            fXHit    = hit.x_pos; //X();
            fYHit    = hit.y_pos;//Y();
            fZHit    = hit.z_pos;//Z();
            fXErrHit = hit.x_err;//XErr();
            fYErrHit = hit.y_err;//YErr();
            fZErrHit = hit.z_err;//ZErr();
            fT0Hit   = hit.ts0_ns;//T0();
            fT1Hit   = hit.ts1_ns;//T1();
            //fT0CorrHit = hit.T0Corr();
            //fT1CorrHit = hit.T1Corr();
            
            int mactmp = hit.feb_id[0];
            auto ittmp = hit.pesmap.find(mactmp);
            int chantmp = (*ittmp).second[0].first;
            fHitReg  = MacToADReg(mactmp);
            fHitSubSys = -0.8*INT_MAX;
            if (fHitReg==38||fHitReg==52||fHitReg==56||fHitReg==48||fHitReg==46)
                fHitSubSys = 0;
            if (fHitReg==50||fHitReg==54||fHitReg==42||fHitReg==44)
                fHitSubSys = 1;
            if (fHitReg==58)
                fHitSubSys = 2;
            for (int i=0; i<64; i++) {
                fHitTrk[i]  = -0.8*INT_MAX;
                fHitPDG[i]  = -0.8*INT_MAX;
            }
            //auto trks = CRTTruthMatchUtils::AllTrueIds(crtSimHitHandle,event,fCRTSimHitProducerLabel,fNHit-1);//hit.TrackID();
            size_t index = 0;
            //for ( auto i=trks.begin(); i!=trks.end(); i++ ) {
            for ( auto i=ids.begin(); i!=ids.end(); i++ ) {   
                fHitTrk[index]  = *i;
                if ( particleMap.find(fHitTrk[index]) != particleMap.end())
                    fHitPDG[index] = particleMap[fHitTrk[index]]->PdgCode();
                index++;
            }
            fHitMod  = MacToAuxDetID(febMap, mactmp, chantmp);//hit.Module();
            fHitStrip = ChannelToAuxDetSensitiveID(mactmp, chantmp);//hit.Strip();

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

  void FillFebMap(map<int,vector<pair<int,int>>>& m) {
    std::string dir = "/icarus/app/users/chilgenb/dev_areas/v08_22_00_prof/srcs/icaruscode/icaruscode/Geometry/gdml";
    std::ifstream fin;
    fin.open(dir+"feb_map.txt",std::ios::in);
    if(fin.good()) std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
    else std::cout << "could not open file 'feb_map.txt' for reading!" << std::endl;
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
        m[mod].push_back(std::make_pair(std::stoi(row[1]),std::stoi(row[2])));
        if(row.size()>3)
            m[mod].push_back(std::make_pair(std::stoi(row[3]),std::stoi(row[4])));
    }
    std::cout << "filled febMap with " << m.size() << " entries" << std::endl;
    fin.close();
  }

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

  int GetAuxDetRegion(geo::AuxDetGeo const& adgeo)
  {
    char type = ModToAuxDetType(adgeo);
    std::string base = "volAuxDet_";
    switch ( type ) {
      case 'c' : base+= "CERN"; break;
      case 'd' : base+= "DC"; break;
      case 'm' : base+= "MINOS"; break;
      case 'e' :
          std::cout << "error in GetAuxDetRegion: type error" << std::endl;
          return INT_MAX;
    }
    base+="_module_###_";
    std::string volName(adgeo.TotalVolume()->GetName());

    //module name has 2 possible formats
    //  volAuxDet_<subsystem>_module_###_<region>
    //  volAuxDet_<subsystem>_module_###_cut###_<region>

    std::string reg = volName.substr(base.length(),volName.length());
    if( reg.find("_")!=std::string::npos)
        reg = reg.substr(reg.find("_")+1,reg.length());

    if(reg == "Top")        return 30;
    if(reg == "RimWest")    return 31;
    if(reg == "RimEast")    return 32;
    if(reg == "RimSouth")   return 33;
    if(reg == "RimNorth")   return 34;
    if(reg == "WestSouth")  return 40;
    if(reg == "WestCenter") return 41;
    if(reg == "WestNorth")  return 42;
    if(reg == "EastSouth")  return 43;
    if(reg == "EastCenter") return 44;
    if(reg == "EastNorth")  return 45;
    if(reg == "South")      return 46;
    if(reg == "North")      return 47;
    if(reg == "Bottom")     return 50;

    return INT_MAX;
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

      if(mac>=107 && mac<=190) return 30; //top
      if(mac>=191 && mac<=204) return 31; //rim west
      if(mac>=205 && mac<=218) return 32; //rim east
      if(mac>=219 && mac<=224) return 33; //rim south
      if(mac>=225 && mac<=230) return 34; //rim north
      if(            mac<=12 ) return 40; //west side, south stack
      if(mac>=13  && mac<=24 ) return 41; //west side, center stack
      if(mac>=25  && mac<=36 ) return 42; //west side, north stack
      if(mac>=37  && mac<=48 ) return 43; //east side, south stack
      if(mac>=49  && mac<=60 ) return 44; //east side, center stack
      if(mac>=60  && mac<=61 ) return 45; //east side, north stack
      if(mac>=73  && mac<=84 ) return 46; //south
      if(mac>=85  && mac<=92 ) return 47; //north
      if(mac>=93 && mac<=106) return 50; //bottom

      return 0;
  }

  char MacToType(uint32_t mac) {

      uint32_t reg = MacToADReg(mac);

      if( reg>29 && reg<40 )
        return 'c';
      if( reg>39 && reg<50 )
        return 'm';
      if( reg==50 )
        return 'd';

      return 'e';
  }

  uint32_t MacToTypeCode(uint32_t mac) {

      char reg = MacToType(mac);

      if(reg=='c')
        return 0; //'c';
      if(reg=='m')
        return 1; //'m';
      if(reg=='d')
        return 2; //'d';

      return UINT32_MAX;//'e';
  }

  //for C- and D-modules, mac address is same as AD ID
  //three M-modules / FEB, each modules readout at both ends
  //  numbering convention is module from FEB i 
  //  is readout on the opposite end by FEB i+50
  //  return FEB i
  std::pair<uint32_t,uint32_t> ADToMac(const map<int,vector<pair<int,int>>>& febMap, uint32_t adid) {
      for(auto const& p : febMap) {
          if((uint32_t)p.first!=adid)
              continue;
          if(p.second.size()==2)
              return std::make_pair((uint32_t)p.second[0].first,(uint32_t)p.second[1].first);
          else
              return std::make_pair((uint32_t)p.second[0].first,(uint32_t)p.second[0].first);
      }
      return std::make_pair(UINT32_MAX,UINT32_MAX);
  }

  int MacToAuxDetID(const map<int,vector<pair<int,int>>>& febMap, int mac, int chan){
      char type = MacToType(mac);
      int pos=1;
      if(type=='m')
          pos = chan/10 + 1;

      for(auto const& p : febMap) {
          if(p.second[0].first == mac&&p.second[1].second==pos)
              return (uint32_t)p.first;
          if(p.second.size()==2)
              if(p.second[1].first==mac&&p.second[1].second==pos)
                  return (uint32_t)p.first;
      }

    return INT_MAX;
  }

  int ChannelToAuxDetSensitiveID(int mac, int chan) {
    int type = MacToTypeCode(mac);
    if (type==2) return chan; //d
    if (type==0) return chan/2; //c
    if (type==1) return (chan % 10)*2; //m

    return INT_MAX;
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
