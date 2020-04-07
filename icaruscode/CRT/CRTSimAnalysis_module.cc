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
  //void FillFebMap();//map<int,vector<pair<int,int>>>& m);
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
  int RegToTypeCode(int reg);

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
    explicit CRTSimAnalysis(Parameters const& config);

    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;

  private:

    void FillFebMap();

    // The parameters we'll read from the .fcl file.
    art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
    art::InputTag fCRTSimHitProducerLabel;        ///< The name of the producer that created hits
    art::InputTag fCRTDetSimProducerLabel;
    vector<int> fPDGs;                       ///< PDG code of particle we'll focus on
    vector<float> fMinMomenta;
    vector<float> fMaxMomenta;

    static map<int, vector<pair<int,int>>> fFebMap;

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
    //static const int kMaxSeg = 4000;
    int        fCDEvent;
    int        fCDTrackID;
    int        fNCD;
    int        fCDpdg;
    //vector<int> fCDRegions;
    vector<vector<double>> fCDSlopes;
    vector<vector<double>> fCDpe;
    vector<vector<double>> fCDxyzt;
    //int        fCDRegions[kMaxSeg];
    //double     fCDSlopes[kMaxSeg][3]; //direction cosines
    //double     fCDpe[kMaxSeg][4];   //4-momentum
    //double     fCDxyzt[kMaxSeg][4];   //4-position
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
    int      fRegEvent;
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
    int      fDetEvent;
    int      fNChan; ///< number of channels above threshold for this front-end board readout
    int      fEntry; ///< front-end board entry number (reset for each event)
    double      fTTrig;      ///< signal time w.r.t. global event time of primary channel providing trigger
    int      fChanTrig; ///< channel which provided the trigger for this front-end board
    int      fFEBReg; ///< CRT region for this front-end board
    int      fMac5; ///< Mac5 address for this front-end board
    int      fTriggerPair[2]; ///< two channels which provided the coincidence (useful for C or D modules)
    int      fMacPair[2]; ///< two front-end boards with pairwise coincidence ( useful for M modules)
    int      fDetSubSys;
    vector<int> fChan; ///< front-end board channel (0-31 or 0-63)
    vector<double> fT0;///< signal time w.r.t. global event time
    vector<double> fT1;///< signal time w.r.t. PPS
    vector<int> fADC;///< signal amplitude
    vector<vector<int>> fTrackID;///< track ID(s) of particle that produced the signal
    vector<vector<int>> fDetPDG; /// signal inducing particle(s)' PDG code
    /*int      fChan[kDetMax]; ///< front-end board channel (0-31 or 0-63)
    double      fT0[kDetMax]; ///< signal time w.r.t. global event time
    double      fT1[kDetMax]; ///< signal time w.r.t. PPS
    int      fADC[kDetMax]; ///< signal amplitude
    int      fTrackID[kDetMax]; ///< track ID of particle that produced the signal
    int      fDetPDG[kDetMax]; */

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
    vector<int> fHitTrk;
    vector<int> fHitPDG;
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

  map<int,vector<pair<int,int>>> CRTSimAnalysis::fFebMap;
 
  CRTSimAnalysis::CRTSimAnalysis(Parameters const& config)
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

  void CRTSimAnalysis::FillFebMap() { //map<int,vector<pair<int,int>>>& m) {
    if(!this->fFebMap.empty())
        return;
    //std::string dir = "/icarus/app/users/chilgenb/dev_areas/v08_22_00_prof/srcs/icaruscode/icaruscode/Geometry/gdml/";
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file("feb_map.txt",fullFileName);
    std::ifstream fin;
    //fin.open(dir+"feb_map.txt",std::ios::in);
    fin.open(fullFileName,std::ios::in);
    if(fin.good()) std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
    else //std::cout << "could not open file 'feb_map.txt' for reading!" << std::endl;
        throw cet::exception("CRTSimAnalysis::FillFebMap") << "Unable to find/open file 'feb_map.txt'" << std::endl;
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
    fCosmicDisplayNtuple->Branch("event",             &fCDEvent,             "event/I");
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
    fRegionsNtuple->Branch("event",                &fRegEvent,           "event/I");
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
    fDetSimNtuple->Branch("event",                 &fDetEvent,          "event/I");
    fDetSimNtuple->Branch("nChan",                 &fNChan,             "nChan/I");
    fDetSimNtuple->Branch("channel",               &fChan);
    fDetSimNtuple->Branch("t0",                    &fT0);
    fDetSimNtuple->Branch("t1",                    &fT1);
    fDetSimNtuple->Branch("adc",                   &fADC);
    fDetSimNtuple->Branch("trackID",               &fTrackID);
    fDetSimNtuple->Branch("detPDG",                &fDetPDG);
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
    fSimHitNtuple->Branch("trackID",     &fHitTrk);
    fSimHitNtuple->Branch("pdg",         &fHitPDG);
    fSimHitNtuple->Branch("modID",       &fHitMod,      "modID/I");
    fSimHitNtuple->Branch("stripID",     &fHitStrip,    "stripID/I");

    // Define the branches of our SimTrueHit n-tuple
    fTrueCRTHitNtuple->Branch("event",       &fTrueHitEvent,    "event/I");
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

    FillFebMap();//febMap);

    if((*genHandle).size()>1) throw cet::exception("CRTSimAnalysis") << "gen stage MCParticle vector has more than 1 entry!" << std::endl;
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
        //fCDRegions.clear();
        fCDSlopes.clear();
        fCDpe.clear();
        fCDxyzt.clear();

        fRegEvent = fEvent;
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

        struct tagger {
            char type;
            int region;
            std::set<int> layerID;
            std::map<int,int> stripLayer;
            std::pair<int,int> modPair;
            std::map<int,std::vector<double>> xyzt;
        };
           
        std::map<int,tagger> taggers;    

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

                    //fADType[fNAuxDet] = ModToTypeCode(adGeo);
                    //fAuxDetReg[fNAuxDet] = GetAuxDetRegion(adGeo);
                    fADType.push_back(ModToTypeCode(adGeo));
                    fAuxDetReg.push_back(GetAuxDetRegion(adGeo));

                    if ( fADType[fNAuxDet] == 0 || fADType[fNAuxDet] == 2 )
                        layid = (stripPosModule[1] > 0);

                    // if 'm' type
                    if ( fADType[fNAuxDet] == 1 ) {
                        // if east or west stacks (6 in total)
                        if ( fAuxDetReg[fNAuxDet] >=40 && fAuxDetReg[fNAuxDet] <=45 ) {
                            layid = ( modulePosMother[0]>0 );
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
                    fADMac.push_back(ADToMac(this->fFebMap, channel.AuxDetID()).first);

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
        fDetSubSys = MacToTypeCode(fMac5);
        fChan.clear();
        fT0.clear();
        fT1.clear();
        fADC.clear();
        fTrackID.clear();
        fDetPDG.clear();
 
        vector<int> missedIDs;
        for ( auto const chandat : febdat.ChanData()) {
          //DetSim tree contains all entries (not neccessarily from muons)
          fChan.push_back(chandat.Channel());
          fT0.push_back(chandat.T0());
          fT1.push_back(chandat.T1());
          fADC.push_back(chandat.ADC());
          fTrackID.push_back({});
          fDetPDG.push_back({});
          for( int trk : chandat.TrackID()) {
              fTrackID[fNChan].push_back(trk);
              if (particleMap.find(trk) != particleMap.end() )
                  fDetPDG[fNChan].push_back(particleMap[trk]->PdgCode());
              else {
                  fDetPDG[fNChan].push_back(0);
                  missedIDs.push_back(trk);
              }
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
    fNHit = 0;
    if (isCRTSimHit) {

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

        std::cout << "looping over sim hits..." << std::endl;
        for ( auto const& hit : *crtSimHitHandle )
        {
            fNHit++;
            fHitEvent = fEvent;
            fXHit    = hit.x_pos;
            fYHit    = hit.y_pos;
            fZHit    = hit.z_pos;
            fXErrHit = hit.x_err;
            fYErrHit = hit.y_err;
            fZErrHit = hit.z_err;
            fT0Hit   = hit.ts0_ns;
            fT1Hit   = hit.ts1_ns;
            //fT0CorrHit = hit.T0Corr();
            //fT1CorrHit = hit.T1Corr();

            int mactmp = hit.feb_id[0];
            fHitReg  = MacToADReg(mactmp);
            fHitSubSys = RegToTypeCode(fHitReg);

            auto ittmp = hit.pesmap.find(mactmp);
            if (ittmp==hit.pesmap.end()) {
                std::cout << "hitreg: " << fHitReg << std::endl;
                std::cout << "fHitSubSys: "<< fHitSubSys << std::endl;
                std::cout << "mactmp = " << mactmp << std::endl;
                std::cout << "could not find mac in pesmap!" << std::endl;
                continue;
            }
            int chantmp = (*ittmp).second[0].first;

            //auto trks = CRTTruthMatchUtils::AllTrueIds(crtSimHitHandle,event,fCRTSimHitProducerLabel,fNHit-1);//hit.TrackID();
            fHitTrk.clear();
            fHitPDG.clear();
            for ( auto i=ids.begin(); i!=ids.end(); i++ ) {   
                fHitTrk.push_back(*i);
                if ( particleMap.find(fHitTrk.back()) != particleMap.end())
                    fHitPDG.push_back(particleMap[fHitTrk.back()]->PdgCode());
                else
                    fHitPDG.push_back(INT_MAX);
            }
            fHitMod  = MacToAuxDetID(this->fFebMap, mactmp, chantmp);
            fHitStrip = ChannelToAuxDetSensitiveID(mactmp, chantmp);

            fSimHitNtuple->Fill();
        }//for CRT Hits
    }//if CRT Hits

    else std::cout << "CRTHit products not found! (expected if gen/G4/detsim step)" << std::endl;


  } // CRTSimAnalysis::analyze()
  
  
  DEFINE_ART_MODULE(CRTSimAnalysis)

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
      if(mac>=61  && mac<=72 ) return 45; //east side, north stack
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
          if(p.second[0].first == mac && p.second[0].second==pos)
              return (uint32_t)p.first;
          if(p.second.size()==2)
              if(p.second[1].first==mac && p.second[1].second==pos)
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

  int RegToTypeCode(int reg){
      if(reg>=30&&reg<40)
          return 0;
      if(reg>=40&&reg<50)
          return 1;
      if(reg==50)
          return 2;
      std::cout << "ERROR in RegToTypeCode: unknown reg code!" << std::endl;
      return -1;
  }

}//local namespace
