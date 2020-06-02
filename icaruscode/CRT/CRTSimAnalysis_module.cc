/**
 * @file   CRTSimAnalysis_module.cc
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
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "icaruscode/CRT/CRTUtils/CRTBackTracker.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;

namespace {
  int ProcessToICode(string const& p);

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
  class CRTSimAnalysis : public art::EDAnalyzer
  {
  public:
    
    // -------------------------------------------------------------------
    // -------------------------------------------------------------------
    // Standard constructor for an ART module with configuration validation;
    // we don't need a special destructor here.

    /// Constructor: configures the module (see the Config structure above)
    explicit CRTSimAnalysis(fhicl::ParameterSet const& p);

    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;

  private:

    // The parameters we'll read from the .fcl file.
    art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    art::InputTag fCRTSimHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fCRTTrueHitProducerLabel;
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
    int        fCDTrackID;
    int        fNCD;
    int        fCDpdg;
    vector<vector<double>> fCDSlopes; ///< direction cosines
    vector<vector<double>> fCDpe;     ///< 4-momentum
    vector<vector<double>> fCDxyzt;   ///< 4-position
    /// @}

    /// @name The variables that will go into the Gen n-tuple.
    /// @{
    int      fNGen;
    vector<int> fGenTrack;
    vector<int> fGenPDG;
    vector<vector<double>>   fGenStartXYZT;
    vector<vector<double>>   fGenEndXYZT;
    vector<vector<double>>   fGenStartPE;
    vector<vector<double>>   fGenEndPE;
    /// @}

    /// @name The variables that will go into the Simulation n-tuple.
    /// @{
    uint32_t fSimHits; ///< number of trajectory points for each MCParticle
    float    fTrackLength; ///< total track length for each MCParticle
    int      fSimPDG;       ///< PDG ID of the particle being processed
    int      fSimProcess;   ///< process that created the particle (e.g. brehmstralung)
    int      fSimEndProcess; ///< process the killed the particle (e.g. annihilation)
    int      fSimTrackID;   ///< GEANT ID of the particle being processed
    uint32_t fNAuxDet;   ///< Number of scintillator strips hit

    vector<uint32_t> fAuxDetID;  ///< Global CRT module ID
    vector<uint32_t> fAuxDetSensitiveID; ///< Strip ID in module
    vector<double>    fADEDep; ///< Energy deposited in CRT strip (GeV)
    vector<double>    fADdEdx; ///< average dEdx for particle traversing CRT strip
    vector<double>    fADTrackLength; ///< Track length in CRT strip (cm)
    vector<uint32_t> fAuxDetReg; ///< CRT region code
    vector<uint32_t> fADMac; ///< Mac5 address of the CRT module
    vector<uint32_t> fADType;
    vector<vector<double>>   fADEnterXYZT; ///< 4-position of entry into CRT strip
    vector<vector<double>>   fADExitXYZT; ///< 4-position of exit from CRT strip
    vector<vector<double>>   fADEnterPE; ///< 4-position of entry into CRT strip
    vector<vector<double>>   fADExitPE; ///< 4-position of exit from CRT strip

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
    int      fNReg;
    int      fRegFid;
    int      fRegActive;
    int      fRegInactive;
    int      fRegCRTs;
    int      fRegTrkID;
    int      fRegPDG;
    vector<int>            fRegRegions;
    vector<double>         fRegEDep;
    vector<double>         fRegdL;
    vector<vector<double>> fRegEntryPE;
    vector<vector<double>> fRegExitPE;
    vector<vector<double>> fRegEntryXYZT;
    vector<vector<double>> fRegExitXYZT;
    vector<vector<double>> fRegEntrySlope;
    vector<vector<double>> fRegExitSlope;
    vector<int>            fRegOpDetID;
    vector<double>         fRegDistToOpDet;
    vector<vector<double>> fRegOpDetXYZT;

    //CRT data product vars
    //static const int kDetMax = 64;
    //int      fNChan; ///< number of channels above threshold for this front-end board readout
    int      fEntry; ///< front-end board entry number (reset for each event)
    int      fFEBReg; ///< CRT region for this front-end board
    int      fMac5; ///< Mac5 address for this front-end board
    int      fDetSubSys;
    double   fT0;///< signal time w.r.t. global event time
    double   fT1;///< signal time w.r.t. PPS
    int      fADC[64];///< signal amplitude
    vector<int> fTrackID;///< track ID(s) of particle that produced the signal
    vector<int> fDetPDG; /// signal inducing particle(s)' PDG code

    // sim CRT hit product vars
    float    fXHit; ///< reconstructed X position of CRT hit (cm)
    float    fYHit; ///< reconstructed Y position of CRT hit (cm)
    float    fZHit; ///< reconstructed Z position of CRT hit (cm)
    float    fXHitErr; ///< stat error of CRT hit reco X (cm)
    float    fYHitErr; ///< stat error of CRT hit reco Y (cm)
    float    fZHitErr; ///< stat error of CRT hit reco Z (cm)
    int32_t    fT0Hit; ///< hit time w.r.t. global event time
    int32_t    fT1Hit; ///< hit time w.r.t. PPS
    int       fHitReg; ///< region code of CRT hit
    int       fHitSubSys;
    int       fNHit; ///< number of CRT hits for this event
    vector<int> fHitTrk;
    vector<int> fHitPDG;
    int       fHitStrip;
    int       fHitMod;
    int       fNHitFeb;
    vector<float> fHitPe;
    float     fHitTotPe;
    float     fHitPeRms;

    // truth CRT hit vars
    double    fTrueXHit; ///< reconstructed X position of CRT hit (cm)
    double    fTrueYHit; ///< reconstructed Y position of CRT hit (cm)
    double    fTrueZHit; ///< reconstructed Z position of CRT hit (cm)
    double    fTrueXHitErr; ///< stat error of CRT hit reco X (cm)
    double    fTrueYHitErr; ///< stat error of CRT hit reco Y (cm)
    double    fTrueZHitErr; ///< stat error of CRT hit reco Z (cm)
    double    fTrueT0Hit; ///< hit time w.r.t. global event time
    double    fTrueT1Hit; ///< hit time w.r.t. global event time
    int       fTrueHitReg; ///< region code of CRT hit
    int       fTrueHitSubSys;
    uint32_t  fTrueNHit; ///< number of CRT hits for this event
    vector<int> fTrueHitTrk;
    vector<int> fTrueHitPDG;
    int       fTrueHitStrip;
    int       fTrueHitMod;
    int       fTrueNHitFeb;
    vector<float> fTrueHitPe;
    float     fTrueHitTotPe;
    float     fTrueHitPeRms;

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
    CRTBackTracker bt;
    CRTCommonUtils* fCrtutils;    

  }; // class CRTSimAnalysis


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
  CRTSimAnalysis::CRTSimAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}
    , fSimulationProducerLabel(p.get<art::InputTag>("SimulationLabel","largeant"))
    , fCRTSimHitProducerLabel(p.get<art::InputTag>("CRTSimHitLabel","crthit"))
    , fCRTTrueHitProducerLabel(p.get<art::InputTag>("CRTTrueHitLabel","crttruehit"))
    , fCRTDetSimProducerLabel(p.get<art::InputTag>("CRTDetSimLabel","crtdaq"))
    , fPDGs(p.get<vector<int>>("PDGs"))
    , fMinMomenta(p.get<vector<float>>("MinMomenta"))
    , fMaxMomenta(p.get<vector<float>>("MaxMomenta"))
    , bt(p.get<fhicl::ParameterSet>("CRTBackTrack"))
    , fCrtutils(new CRTCommonUtils())
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
  void CRTSimAnalysis::beginJob()
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
    fCosmicDisplayNtuple->Branch("event",             &fEvent,               "event/I");
    fCosmicDisplayNtuple->Branch("trackID",           &fCDTrackID,           "trackID/I");
    fCosmicDisplayNtuple->Branch("nSeg",              &fNCD,                 "nSeg/I");
    fCosmicDisplayNtuple->Branch("pdg",               &fCDpdg,               "pdg/I");
    //fCosmicDisplayNtuple->Branch("regions",           &fCDRegions);
    fCosmicDisplayNtuple->Branch("slopes",            &fCDSlopes);
    fCosmicDisplayNtuple->Branch("pe",                &fCDpe);
    fCosmicDisplayNtuple->Branch("xyzt",              &fCDxyzt);

    // Define the branches of our Gen n-tuple
    fGenNtuple->Branch("event",        &fEvent,         "event/I");
    fGenNtuple->Branch("nGen",         &fNGen,          "nGen/I");
    fGenNtuple->Branch("trackID",      &fGenTrack);
    fGenNtuple->Branch("pdg",          &fGenPDG);
    fGenNtuple->Branch("startXYZT",    &fGenStartXYZT);
    fGenNtuple->Branch("endXYZT",      &fGenEndXYZT);
    fGenNtuple->Branch("startPE",      &fGenStartPE);
    fGenNtuple->Branch("endPE",        &fGenEndPE);

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
    fSimulationNtuple->Branch("auxDetSensitiveID", &fAuxDetSensitiveID);
    fSimulationNtuple->Branch("auxDetID",          &fAuxDetID);
    fSimulationNtuple->Branch("auxDetEDep",        &fADEDep);
    fSimulationNtuple->Branch("auxDetdEdx",        &fADdEdx);
    fSimulationNtuple->Branch("auxDetTrackLength", &fADTrackLength);
    fSimulationNtuple->Branch("auxDetEnterXYZT",   &fADEnterXYZT);
    fSimulationNtuple->Branch("auxDetExitXYZT",    &fADExitXYZT);
    fSimulationNtuple->Branch("auxDetEnterPE",     &fADEnterPE);
    fSimulationNtuple->Branch("auxDetExitPE",      &fADExitPE);
    fSimulationNtuple->Branch("auxDetRegion",      &fAuxDetReg);
    fSimulationNtuple->Branch("mac5",              &fADMac);
    fSimulationNtuple->Branch("adType",            &fADType);

    fSimulationNtuple->Branch("startXYZT",         fStartXYZT,          "startXYZT[4]/F");
    fSimulationNtuple->Branch("endXYZT",           fEndXYZT,            "endXYZT[4]/F");
    fSimulationNtuple->Branch("startPE",           fStartPE,            "startPE[4]/F");
    fSimulationNtuple->Branch("endPE",             fEndPE,              "endPE[4]/F");
    fSimulationNtuple->Branch("nChan",             &fNAuxDet,           "nChan/I");
    fSimulationNtuple->Branch("mother",            &fMother,            "mother/I");
    fSimulationNtuple->Branch("nDaught",           &fNDaught,           "nDaught/I");

    //regions tree
    fRegionsNtuple->Branch("event",                &fEvent,              "event/I");
    fRegionsNtuple->Branch("nReg",                 &fNReg,               "nReg/I");
    fRegionsNtuple->Branch("fiducial",             &fRegFid,             "fiducial/I");
    fRegionsNtuple->Branch("active",               &fRegActive,          "active/I");
    fRegionsNtuple->Branch("inactive",             &fRegInactive,        "inactive/I");
    fRegionsNtuple->Branch("crts",                 &fRegCRTs,            "crts/I");
    fRegionsNtuple->Branch("regions",              &fRegRegions);
    fRegionsNtuple->Branch("pdg",                  &fRegPDG,             "pdg/I");
    fRegionsNtuple->Branch("trackID",              &fRegTrkID,           "trackID/I");
    fRegionsNtuple->Branch("eDep",                 &fRegEDep);
    fRegionsNtuple->Branch("dL",                   &fRegdL);
    fRegionsNtuple->Branch("opDetID",              &fRegOpDetID);
    fRegionsNtuple->Branch("distToOpDet",          &fRegDistToOpDet);
    fRegionsNtuple->Branch("opDetXYZT",            &fRegOpDetXYZT);
    fRegionsNtuple->Branch("entryPE",              &fRegEntryPE);
    fRegionsNtuple->Branch("exitPE",               &fRegExitPE);
    fRegionsNtuple->Branch("entryXYZT",            &fRegEntryXYZT);
    fRegionsNtuple->Branch("exitXYZT",             &fRegExitXYZT);
    fRegionsNtuple->Branch("entrySlope",           &fRegEntrySlope);
    fRegionsNtuple->Branch("exitSlope",            &fRegExitSlope);

    // Define the branches of our DetSim n-tuple 
    fDetSimNtuple->Branch("event",                 &fEvent,             "event/I");
    //fDetSimNtuple->Branch("nChan",                 &fNChan,             "nChan/I");
    fDetSimNtuple->Branch("t0",                    &fT0,                "t0/D");
    fDetSimNtuple->Branch("t1",                    &fT1,                "t1/D");
    fDetSimNtuple->Branch("adc",                   fADC);
    fDetSimNtuple->Branch("trackID",               &fTrackID);
    fDetSimNtuple->Branch("detPDG",                &fDetPDG);
    fDetSimNtuple->Branch("entry",                 &fEntry,             "entry/I");
    fDetSimNtuple->Branch("mac5",                  &fMac5,              "mac5/I");
    fDetSimNtuple->Branch("region",                &fFEBReg,            "region/I");
    fDetSimNtuple->Branch("subSys",                &fDetSubSys,         "subSys/I");

    // Define the branches of our SimHit n-tuple
    fSimHitNtuple->Branch("event",       &fEvent,       "event/I");
    fSimHitNtuple->Branch("nHit",        &fNHit,        "nHit/I");
    fSimHitNtuple->Branch("x",           &fXHit,        "x/F");
    fSimHitNtuple->Branch("y",           &fYHit,        "y/F");
    fSimHitNtuple->Branch("z",           &fZHit,        "z/F");
    fSimHitNtuple->Branch("xErr",        &fXHitErr,     "xErr/F");
    fSimHitNtuple->Branch("yErr",        &fYHitErr,     "yErr/F");
    fSimHitNtuple->Branch("zErr",        &fZHitErr,     "zErr/F");
    fSimHitNtuple->Branch("t0",          &fT0Hit,       "t0/I");
    fSimHitNtuple->Branch("t1",          &fT1Hit,       "t1/I");
    fSimHitNtuple->Branch("region",      &fHitReg,      "region/I");  
    fSimHitNtuple->Branch("subSys",      &fHitSubSys,   "subSys/I");
    fSimHitNtuple->Branch("trackID",     &fHitTrk);
    fSimHitNtuple->Branch("pdg",         &fHitPDG);
    fSimHitNtuple->Branch("modID",       &fHitMod,      "modID/I");
    fSimHitNtuple->Branch("stripID",     &fHitStrip,    "stripID/I");
    fSimHitNtuple->Branch("nfeb",        &fNHitFeb,     "nfeb/I");
    fSimHitNtuple->Branch("hitPe",       &fHitPe);
    fSimHitNtuple->Branch("totPe",       &fHitTotPe,    "totPe/F");
    fSimHitNtuple->Branch("rmsPe",       &fHitPeRms,    "rmsPe/F");

    // Define the branches of our SimTrueHit n-tuple
    fTrueCRTHitNtuple->Branch("event",       &fEvent,           "event/I");
    fTrueCRTHitNtuple->Branch("nHit",        &fTrueNHit,        "nHit/I");
    fTrueCRTHitNtuple->Branch("x",           &fTrueXHit,        "x/D");
    fTrueCRTHitNtuple->Branch("y",           &fTrueYHit,        "y/D");
    fTrueCRTHitNtuple->Branch("z",           &fTrueZHit,        "z/D");
    fTrueCRTHitNtuple->Branch("xErr",        &fTrueXHitErr,     "xErr/D");
    fTrueCRTHitNtuple->Branch("yErr",        &fTrueYHitErr,     "yErr/D");
    fTrueCRTHitNtuple->Branch("zErr",        &fTrueZHitErr,     "zErr/D");
    fTrueCRTHitNtuple->Branch("t0",          &fTrueT0Hit,       "t0/D");
    fTrueCRTHitNtuple->Branch("t1",          &fTrueT1Hit,       "t1/D");
    fTrueCRTHitNtuple->Branch("region",      &fTrueHitReg,      "region/I");
    fTrueCRTHitNtuple->Branch("subSys",      &fTrueHitSubSys,   "subSys/I");
    fTrueCRTHitNtuple->Branch("trackID",     &fTrueHitTrk);
    fTrueCRTHitNtuple->Branch("pdg",         &fTrueHitPDG);
    fTrueCRTHitNtuple->Branch("modID",       &fTrueHitMod,      "modID/I");
    fTrueCRTHitNtuple->Branch("stripID",     &fTrueHitStrip,    "stripID/I");
    fTrueCRTHitNtuple->Branch("nfeb",        &fTrueNHitFeb,     "nfeb/I");
    fTrueCRTHitNtuple->Branch("hitPe",       &fTrueHitPe);
    fTrueCRTHitNtuple->Branch("totPe",       &fTrueHitTotPe,    "totPe/F");
    fTrueCRTHitNtuple->Branch("rmsPe",       &fTrueHitPeRms,    "rmsPe/F");

}
   
  void CRTSimAnalysis::beginRun(const art::Run& /*run*/)
  {
    //art::ServiceHandle<sim::LArG4Parameters> larParameters;
    //fElectronsToGeV = 1./larParameters->GeVToElectrons();
    std::cout << "beginning run" << std::endl;
  }

  //-----------------------------------------------------------------------
  void CRTSimAnalysis::analyze(const art::Event& event) 
  {
    MF_LOG_DEBUG("CRT") << "beginning analyis" << '\n';

    // Check that lists for momenta limits is same size as last of PDGs from FHiCL
    if (fPDGs.size() != fMinMomenta.size() || fPDGs.size() != fMaxMomenta.size())
        throw cet::exception("CRTSimAnalysis")
          << " PDG/Momtenta values not set correctly in fhicl - lists have different sizes"
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;


    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    //CRTBackTracker matches CRTProducts to the true trackIDs
    bt.Initialize(event);

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
	throw cet::exception("CRTSimAnalysis") 
	  << " No simb::MCParticle objects in this event - "
	  << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

    // Handle to AuxDetSimChannel (CRT module) objects generated by LArG4
    art::Handle<vector<sim::AuxDetSimChannel> > auxDetSimChannelHandle;
    if (!event.getByLabel(fSimulationProducerLabel, auxDetSimChannelHandle)) {
        throw cet::exception("CRTSimAnalysis")
          << " No sim::AuxDetSimChannel objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if((*genHandle).size()>1) 
        throw cet::exception("CRTSimAnalysis") << "gen stage MCParticle vector has more than 1 entry!" << '\n';

    auto const& truth = (*genHandle)[0];   
    fNGen = truth.NParticles();
    fGenTrack.clear();
    fGenPDG.clear();
    fGenStartXYZT.clear();
    fGenEndXYZT.clear();
    fGenStartPE.clear();
    fGenEndPE.clear();
   
    for ( int i=0; i<fNGen; i++ )
    {
        auto const& part = truth.GetParticle(i); //simb::MCParticle

        fGenTrack.push_back(part.TrackId());
        fGenPDG.push_back(part.PdgCode());

        const TLorentzVector startPos = part.Position(0);
        const TLorentzVector endPos = part.EndPosition();
        const TLorentzVector startMom = part.Momentum(0);
        const TLorentzVector endMom = part.EndMomentum();

        fGenStartXYZT.push_back({startPos.X(),startPos.Y(),startPos.Z(),startPos.T()});
        fGenEndXYZT.push_back({endPos.X(),endPos.Y(),endPos.Z(),endPos.T()});
        fGenStartPE.push_back({startMom.Px(),startMom.Py(),startMom.Pz(),startMom.E()});
        fGenEndPE.push_back({endMom.Px(),endMom.Py(),endMom.Pz(),endMom.E()});
    }

    fGenNtuple->Fill();

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

        fNCD = 0;
        fCDpdg = fSimPDG;
        fCDTrackID = fSimTrackID;
        //fCDRegions.clear();
        fCDSlopes.clear();
        fCDpe.clear();
        fCDxyzt.clear();

        fNReg = 0;
        fRegFid = 0;
        fRegActive = 0;
        fRegInactive = 0;
        fRegCRTs = 0;
        fRegPDG = fSimPDG;;
        fRegTrkID = fSimTrackID;
        fRegRegions.clear();
        fRegEDep.clear();
        fRegDistToOpDet.clear();
        fRegOpDetID.clear();
        fRegEntryXYZT.clear();
        fRegExitXYZT.clear();
        fRegEntryPE.clear();
        fRegExitPE.clear();
        fRegOpDetXYZT.clear();
        fRegEntrySlope.clear();
        fRegExitSlope.clear();

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
                fCDxyzt.push_back({pos.X(),pos.Y(),pos.Z(),pos.T()});
                fCDpe.push_back({mom.Px(),mom.Py(),mom.Pz(),mom.E()});
                fCDSlopes.push_back({mom.Px()/mom.P(), mom.Py()/mom.P(), mom.Pz()/mom.P()});
                fNCD++;

                // Regions info
                // Check if trajectory points are in cryostats (active + inactve LAr ) 
                if(cryo0.ContainsPosition(point)) {

                        active0 = tpc00.ContainsPosition(point);
                        active1 = tpc01.ContainsPosition(point);
                        activenext0 = tpc00.ContainsPosition(pointnext);
                        activenext1 = tpc01.ContainsPosition(pointnext);

                        // if last point was not in this cryostat or is now entering AV
                        if ( (oldreg!=10&&!active0&&!active1) || (active0&&oldreg!=5) || (active1&&oldreg!=6)) {
                            fRegEntryXYZT.push_back({pos.X(),pos.Y(),pos.Z(),pos.T()});
                            fRegEntryPE.push_back({mom.Px(),mom.Py(),mom.Pz(),mom.E()});
                            fRegEntrySlope.push_back({mom.Px()/mom.P(), mom.Py()/mom.P(), mom.Pz()/mom.P()});
                            oldreg = 10;
                            if (active0) oldreg = 5;
                            if (active1) oldreg = 6;
                        }

                        // if next point is outside of this volume or is last traj. point
                        if (!cryo0.ContainsPosition(pointnext) || (oldreg==10&&(activenext0||activenext1))
                            || i==fSimHits-1
                            || (active0 && !activenext0) || (active1&&!activenext1) ){

                                fRegExitXYZT.push_back({pos.X(),pos.Y(),pos.Z(),pos.T()});
                                fRegExitPE.push_back({mom.Px(),mom.Py(),mom.Pz(),mom.E()});
                                fRegExitSlope.push_back({mom.Px()/mom.P(), mom.Py()/mom.P(), mom.Pz()/mom.P()});
                                if (active0) {
                                    fRegRegions.push_back(5);
                                    fRegActive++;
                                    if(tpc00.InFiducialX(point[0],25,0) && tpc00.InFiducialY(point[1],25,25)
                                      && tpc00.InFiducialZ(point[2],30,50)) 
                                        fRegFid++;

                                }
                                else if (active1) {
                                    fRegRegions.push_back(6);
                                    fRegActive++;
                                    if(tpc01.InFiducialX(point[0],25,0) && tpc01.InFiducialY(point[1],25,25)
                                      && tpc01.InFiducialZ(point[2],30,50))  
                                        fRegFid++;
                                }
                                else {
                                    fRegRegions.push_back(10);
                                    fRegInactive++;
                                }
                                fRegdL.push_back(sqrt(pow(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                    +pow(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                    +pow(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2)));
                                fRegEDep.push_back(fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3]);
                                for (int index=0; index<3; index++) entryPos[index] = fRegEntryXYZT[fNReg][index];
                                entryT = fRegEntryXYZT[fNReg][3];
                                fRegOpDetID.push_back(cryo0.GetClosestOpDet(entryPos));
                                geo::OpDetGeo const& opDet = cryo0.OpDet(fRegOpDetID[fNReg]);
                                opDet.GetCenter(opDetPos);
                                fRegDistToOpDet.push_back(sqrt(pow(opDetPos[0]-entryPos[0],2)
                                                            + pow(opDetPos[1]-entryPos[1],2)
                                                            + pow(opDetPos[2]-entryPos[2],2)));
                                fRegOpDetXYZT.push_back({});
                                for (int index=0; index<3; index++) fRegOpDetXYZT[fNReg].push_back(opDetPos[index]);
                                fRegOpDetXYZT[fNReg].push_back(entryT + fRegDistToOpDet[fNReg]*LAR_PROP_DELAY);
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
                            fRegEntryXYZT.push_back({pos.X(),pos.Y(),pos.Z(),pos.T()});
                            fRegEntryPE.push_back({mom.Px(),mom.Py(),mom.Pz(),mom.E()});
                            fRegEntrySlope.push_back({mom.Px()/mom.P(), mom.Py()/mom.P(), mom.Pz()/mom.P()});
                            oldreg = 12;
                            if (active0) oldreg = 7;
                            if (active1) oldreg = 8;
                        }

                        if (!cryo1.ContainsPosition(pointnext) || (oldreg==12&&(activenext0||activenext1))
                            || i==fSimHits-1
                            || (active0 && !activenext0) || (active1&&!activenext1) ){

                                fRegExitXYZT.push_back({pos.X(),pos.Y(),pos.Z(),pos.T()});
                                fRegExitPE.push_back({mom.Px(),mom.Py(),mom.Pz(),mom.E()});
                                fRegExitSlope.push_back({mom.Px()/mom.P(), mom.Py()/mom.P(), mom.Pz()/mom.P()});
                                if (active0) {
                                    fRegRegions.push_back(7);
                                    fRegActive++;
                                    if(tpc10.InFiducialX(point[0],25,0) && tpc10.InFiducialY(point[1],25,25)
                                      && tpc10.InFiducialZ(point[2],30,50))
                                        fRegFid++;

                                }
                                else if (active1) {
                                    fRegRegions.push_back(8);
                                    fRegActive++;
                                    if(tpc11.InFiducialX(point[0],25,0) && tpc11.InFiducialY(point[1],25,25)
                                      && tpc11.InFiducialZ(point[2],30,50))
                                        fRegFid++;
                                }
                                else {
                                    fRegRegions.push_back(12);
                                    fRegInactive++;
                                }
                                fRegdL.push_back(sqrt(pow(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                                    +pow(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                                    +pow(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2)));
                                fRegEDep.push_back(fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3]);
                                for (int index=0; index<3; index++) entryPos[index] = fRegEntryXYZT[fNReg][index];
                                entryT = fRegEntryXYZT[fNReg][3];
                                fRegOpDetID.push_back(cryo1.GetClosestOpDet(entryPos));
                                geo::OpDetGeo const& opDet = cryo1.OpDet(fRegOpDetID[fNReg]);
                                opDet.GetCenter(opDetPos);
                                fRegDistToOpDet.push_back(sqrt(pow(opDetPos[0]-entryPos[0],2)
                                                            + pow(opDetPos[1]-entryPos[1],2)
                                                            + pow(opDetPos[2]-entryPos[2],2)));
                                fRegOpDetXYZT.push_back({});
                                for (int index=0; index<3; index++) fRegOpDetXYZT[fNReg].push_back(opDetPos[index]);
                                fRegOpDetXYZT[fNReg].push_back(entryT + fRegDistToOpDet[fNReg]*LAR_PROP_DELAY);
                                fNReg++;
                        } // if exiting from volume
                } //if cryo1

        }//for trajectory points

        fCosmicDisplayNtuple->Fill();

        //map module IDs to strip IDs hit by muons
        //map< uint16_t,set<uint8_t>* > muHitMapC; //hits in C modules only
        //map< uint16_t,set<uint8_t>* > muHitMapM; //hits in M modules only
        //map< uint16_t,set<uint8_t>* > muHitMapD; //hits in D modules only
        //map< uint16_t,set<uint8_t>* > muHitMap; //all hits

        map< int, vector<double> > regCRTEnter, regCRTExit, regCRTEnterPE, regCRTExitPE;

        //reinitialize ADChannel vars
        fNAuxDet = 0;
        fADType.clear();
        fAuxDetReg.clear();
        fADTrackLength.clear();
        fADEDep.clear();
        fADdEdx.clear();
        fAuxDetID.clear();
        fAuxDetSensitiveID.clear();
        fADEnterXYZT.clear();
        fADExitXYZT.clear();
        fADEnterPE.clear();
        fADExitPE.clear();
        fADMac.clear();

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

                    size_t adid = channel.AuxDetID();
                    fADType.push_back(fCrtutils->GetAuxDetTypeCode(adid));
                    fAuxDetReg.push_back(
                       fCrtutils->AuxDetRegionNameToNum(fCrtutils->GetAuxDetRegion(adid)));

                    // What is the distance from the hit (centroid of the entry
                    // and exit points) to the readout end?
                    //double trueX = (ide.entryX + ide.exitX) / 2.0;
                    //double trueY = (ide.entryY + ide.exitY) / 2.0;
                    //double trueZ = (ide.entryZ + ide.exitZ) / 2.0;
                    //double trueT = (ide.entryT + ide.exitT) / 2.0;

                    //calculate track length in strip
		    double dx = ide.entryX-ide.exitX;
		    double dy = ide.entryY-ide.exitY;
		    double dz = ide.entryZ-ide.exitZ;
                    double adlength = sqrt(dx*dx+dy*dy+dz*dz);
                    if ( adlength < 0.0001)  continue;

                    fADTrackLength.push_back(adlength);
                    fADEDep.push_back(ide.energyDeposited);
                    fADdEdx.push_back(ide.energyDeposited/fADTrackLength.back());
                    fAuxDetID.push_back(channel.AuxDetID());
                    fAuxDetSensitiveID.push_back(channel.AuxDetSensitiveID());
                    fADEnterXYZT.push_back({ide.entryX,ide.entryY,ide.entryZ,ide.entryT});
                    fADExitXYZT.push_back({ide.exitX,ide.exitY,ide.exitZ,ide.exitT});
                    //we dont have entry mom. or total E info on adsc, so very rough approx. E~P
                    double pmag = sqrt(pow(ide.exitMomentumX,2)+pow(ide.exitMomentumY,2)+pow(ide.exitMomentumZ,2));
                    fADExitPE.push_back({ide.exitMomentumX,ide.exitMomentumY,ide.exitMomentumZ,pmag});
                    fADEnterPE.push_back({fADExitPE[fNAuxDet][0]+pmag/3,fADExitPE[fNAuxDet][1]+pmag/3,
                                          fADExitPE[fNAuxDet][2]+pmag/3,pmag+ide.energyDeposited});
                    fADMac.push_back(fCrtutils->ADToMac(channel.AuxDetID()).first);

                    if (regCRTEnter.find(fAuxDetReg[fNAuxDet])!=regCRTEnter.end()) {
                        if (regCRTEnter[fAuxDetReg[fNAuxDet]][3] > ide.entryT) {
                            regCRTEnter[fAuxDetReg[fNAuxDet]] = fADEnterXYZT[fNAuxDet];
                            regCRTEnterPE[fAuxDetReg[fNAuxDet]] = fADEnterPE[fNAuxDet];
                        }
                    }
                    else {
                        regCRTEnter[fAuxDetReg[fNAuxDet]] = fADEnterXYZT[fNAuxDet];
                        regCRTEnterPE[fAuxDetReg[fNAuxDet]] = fADEnterPE[fNAuxDet];
                    }

                    if (regCRTExit.find(fAuxDetReg[fNAuxDet])!=regCRTExit.end()) {
                        if (regCRTExit[fAuxDetReg[fNAuxDet]][3] < ide.exitT){
                            regCRTExit[fAuxDetReg[fNAuxDet]] = fADExitXYZT[fNAuxDet];
                            regCRTExitPE[fAuxDetReg[fNAuxDet]] = fADExitPE[fNAuxDet];
                        }
                    }
                    else {
                        regCRTExit[fAuxDetReg[fNAuxDet]] = fADExitXYZT[fNAuxDet];
                        regCRTExitPE[fAuxDetReg[fNAuxDet]] = fADExitPE[fNAuxDet];
                    }

                    fNAuxDet++;

	    } // For each IDE (strip hit by muon)
	} // For each SimChannel (module)

        // write values to tree for this event and particle
        fSimulationNtuple->Fill();

        for( auto it=regCRTEnter.begin(); it!=regCRTEnter.end(); it++) {
            //std::cout << "found CRT region " << it->first << std::endl;
            fRegRegions.push_back(it->first);
            fRegEntryXYZT.push_back({(it->second)[0],(it->second)[1],(it->second)[2],(it->second)[3]});
            fRegExitXYZT.push_back({regCRTExit[it->first][0],regCRTExit[it->first][1],regCRTExit[it->first][2],
                                    regCRTExit[it->first][3]});
            fRegdL.push_back(sqrt(pow(fRegExitXYZT[fNReg][0]-fRegEntryXYZT[fNReg][0],2)
                                  +pow(fRegExitXYZT[fNReg][1]-fRegEntryXYZT[fNReg][1],2)
                                  +pow(fRegExitXYZT[fNReg][2]-fRegEntryXYZT[fNReg][2],2)));
            fRegEntryPE.push_back({regCRTEnterPE[it->first][0],regCRTEnterPE[it->first][1],
                                   regCRTEnterPE[it->first][2], regCRTEnterPE[it->first][3]});
            fRegExitPE.push_back({regCRTExitPE[it->first][0],regCRTExitPE[it->first][1],
                                   regCRTExitPE[it->first][2], regCRTExitPE[it->first][3]});
            fRegEDep.push_back(fRegEntryPE[fNReg][3] - fRegExitPE[fNReg][3]);
            fRegDistToOpDet.push_back(-1);
            fRegOpDetID.push_back(-1);
            fRegOpDetXYZT.push_back({-1,-1,-1,-1});
            fRegEntrySlope.push_back({-1,-1,-1});
            fRegExitSlope.push_back({-1,-1,-1});
            
            fNReg++; fRegCRTs++;
        }

        //sort region tree entries by entry time
        int      flag = 1;    // set flag to 1 to start first pass
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

    //CRTData
    art::Handle<vector<CRTData>> crtDetSimHandle;
    if (event.getByLabel(fCRTDetSimProducerLabel, crtDetSimHandle))  {

        vector<art::Ptr<CRTData>> febdata;
        art::fill_ptr_vector(febdata,crtDetSimHandle);
        mf::LogPrint("CRTSimAnalysis") << "found " << febdata.size() << " CRTData" << '\n';

        for ( auto const& data : febdata ) {
           fMac5      = data->fMac5;
           fEntry     = data->fEntry;
           fFEBReg    = fCrtutils->AuxDetRegionNameToNum(fCrtutils->MacToRegion(fMac5));
           fDetSubSys = fCrtutils->MacToTypeCode(fMac5);
           fT0        = data->fTs0;
           fT1        = data->fTs1;
           fTrackID.clear();
           fDetPDG.clear();

           //set max channel by subsystem (32 (c,m) or 64 (d))
           int chmax = 0;
           if(fDetSubSys!=2) chmax=32;
           else chmax=64;

           //loop over all FEB channels
           for(int ch=0; ch<chmax; ch++) {
               fADC[ch] = data->fAdc[ch];
               //loop over all track IDs
           }//loop over channels
           //loop over all track IDs
           for( int trk : bt.AllTrueIds(event,*data)) {
               fTrackID.push_back(trk);
               if (particleMap.find(trk) != particleMap.end() )
                   fDetPDG.push_back(particleMap[trk]->PdgCode());
               else 
                   fDetPDG.push_back(INT_MAX);
           }//for trackIDs

           fDetSimNtuple->Fill();

        } //for CRT FEB events

    }//if crtdetsim products present

    else 
        mf::LogWarning("CRTSimAnalysis") << "CRTData products not found! (expected if gen/G4 step)" << '\n';


    //simulted CRTHits
    art::Handle<vector<CRTHit>> crtSimHitHandle;
    fNHit = 0;
    if (event.getByLabel(fCRTSimHitProducerLabel, crtSimHitHandle)) {

        vector<art::Ptr<CRTHit>> crtSimHits;
        art::fill_ptr_vector(crtSimHits,crtSimHitHandle);

        mf::LogPrint("CRTSimAnalysis") << "found " << crtSimHits.size() << " sim CRTHits" << '\n';
        for ( auto const& hit : crtSimHits )
        {
            fNHit++;
            fXHit     = hit->x_pos;
            fYHit     = hit->y_pos;
            fZHit     = hit->z_pos;
            fXHitErr  = hit->x_err;
            fYHitErr  = hit->y_err;
            fZHitErr  = hit->z_err;
            fT0Hit    = hit->ts0_ns;
            fT1Hit    = hit->ts1_ns;
            fNHitFeb  = hit->feb_id.size();
            fHitTotPe = hit->peshit;
            fHitPeRms = 0.;
            fHitPe.clear();
            fHitTrk.clear();
            fHitPDG.clear();

            uint8_t mactmp = hit->feb_id[0];
            fHitReg  = fCrtutils->AuxDetRegionNameToNum(fCrtutils->MacToRegion(mactmp));
            fHitSubSys = fCrtutils->MacToTypeCode(mactmp);

            if (hit->pesmap.find(mactmp)==hit->pesmap.end()) {
                std::cout << "hitreg: " << fHitReg << std::endl;
                std::cout << "fHitSubSys: "<< fHitSubSys << std::endl;
                std::cout << "mactmp = " << mactmp << std::endl;
                std::cout << "could not find mac in pesmap!" << std::endl;
                continue;
            }

            for(auto const& mactopes : hit->pesmap){
                for(auto const& chanpe : mactopes.second) {
                    fHitPe.push_back(chanpe.second);
                    fHitPeRms+=pow(chanpe.second-fHitTotPe,2);
                }
            }
            fHitPeRms = sqrt(fHitPeRms/(fHitPe.size()-1));

            //loop over all track IDs
            for( int trk : bt.AllTrueIds(event,*hit)) {
                fHitTrk.push_back(trk);
                if ( particleMap.find(trk) != particleMap.end())
                    fHitPDG.push_back(particleMap[trk]->PdgCode());
                else
                    fHitPDG.push_back(INT_MAX);
            }
            fHitMod  = fCrtutils->MacToAuxDetID(mactmp, 0);
            fHitStrip = fCrtutils->ChannelToAuxDetSensitiveID(mactmp,0);

            fSimHitNtuple->Fill();
        }//for CRTHits
    }//if sim CRTHits

    else 
        mf::LogWarning("CRTSimAnalysis") 
            << "CRTHit products not found! (expected if gen/G4/detsim step)" << '\n';

    //true CRTHits
    art::Handle<vector<CRTHit>> crtTrueHitHandle;
    fTrueNHit = 0;
    if (event.getByLabel(fCRTTrueHitProducerLabel, crtTrueHitHandle)) {

        vector<art::Ptr<CRTHit>> crtTrueHits;
        art::fill_ptr_vector(crtTrueHits,crtTrueHitHandle);

        mf::LogPrint("CRTSimAnalysis") << "found " << crtTrueHits.size() << " true CRTHits" << '\n';
        for ( auto const& hit : crtTrueHits )
        {
            fTrueNHit++;
            fTrueXHit    = hit->x_pos;
            fTrueYHit    = hit->y_pos;
            fTrueZHit    = hit->z_pos;
            fTrueXHitErr = hit->x_err;
            fTrueYHitErr = hit->y_err;
            fTrueZHitErr = hit->z_err;
            fTrueT0Hit   = hit->ts0_ns;
            fTrueT1Hit   = hit->ts1_ns;
            fTrueNHitFeb  = hit->feb_id.size();
            fTrueHitTotPe = hit->peshit;
            fTrueHitPeRms = 0.;
            fTrueHitPe.clear();
            fTrueHitTrk.clear();
            fTrueHitPDG.clear();

            uint8_t mactmp = hit->feb_id[0];
            fTrueHitReg  = fCrtutils->AuxDetRegionNameToNum(fCrtutils->MacToRegion(mactmp));
            fTrueHitSubSys = fCrtutils->MacToTypeCode(mactmp);

            if (hit->pesmap.find(mactmp)==hit->pesmap.end()) {
                std::cout << "hitreg: " << fTrueHitReg << std::endl;
                std::cout << "fHitSubSys: "<< fTrueHitSubSys << std::endl;
                std::cout << "mactmp = " << mactmp << std::endl;
                std::cout << "could not find mac in pesmap!" << std::endl;
                continue;
            }

            for(auto const& mactopes : hit->pesmap){
                for(auto const& chanpe : mactopes.second) {
                    fTrueHitPe.push_back(chanpe.second);
                    fTrueHitPeRms+=pow(chanpe.second-fHitTotPe,2);
                }
            }
            fTrueHitPeRms = sqrt(fTrueHitPeRms/(fTrueHitPe.size()-1));

            //loop over all track IDs
            for( int trk : bt.AllTrueIds(event,*hit)) {
                fTrueHitTrk.push_back(trk);
                if ( particleMap.find(trk) != particleMap.end())
                    fTrueHitPDG.push_back(particleMap[trk]->PdgCode());
                else
                    fTrueHitPDG.push_back(INT_MAX);
            }
            //more work to be done with PEs map...for now fix at chan 0
            fTrueHitMod  = fCrtutils->MacToAuxDetID(mactmp, 0);
            fTrueHitStrip = fCrtutils->ChannelToAuxDetSensitiveID(mactmp,0);

            fTrueCRTHitNtuple->Fill();
        }//for CRTHits
    }//if true CRTHits

    else
        mf::LogWarning("CRTSimAnalysis") << "true CRTHit products not found" << '\n';



  } // CRTSimAnalysis::analyze()
  
  DEFINE_ART_MODULE(CRTSimAnalysis)

} // namespace crt
} // namespace icarus


// Back to our local namespace.
namespace {

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

}//local namespace
