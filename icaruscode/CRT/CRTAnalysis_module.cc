/**
 * @file   CRTAnalysis_module.cc
 * @brief  Access CRT data and reco products and compare to MCTruth info 
 * @author Chris Hilgenberg (Chris.Hilgenberg@rams.colostate.edu)
 * 
 * The last revision of this code was done in August 2018 with LArSoft v06_85_00.
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

#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "icaruscode/CRT/CRTProducts/CRTData.hh"
#include "icaruscode/CRT/CRTProducts/CRTHit.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;

namespace {

  
  char ModToAuxDetType(geo::AuxDetGeo const& adgeo);
  int GetAuxDetRegion(geo::AuxDetGeo const& adgeo);
  int ProcessToICode(string const& p);
  uint32_t MacToADReg(uint32_t mac);
  char MacToType(uint32_t mac);
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
    TH1D* fPDGCodeHist;     ///< PDG code of all particles
    TH1D* fMomentumHist;    ///< momentum [GeV] of all selected particles
    TH1D* fTrackLengthHist; ///< true length [cm] of all selected particles

    // The n-tuples we'll create.
    TTree* fSimulationNtuple;     ///< tuple with simulated data
    TTree* fDetSimNtuple;
    TTree* fSimHitNtuple;
    TTree* fTruthMatchNtuple; 

    // The comment lines with the @ symbols define groups in doxygen. 
    /// @name The variables that will go into both n-tuples.
    /// @{
    int fEvent;        ///< number of the event being processed
    int fRun;          ///< number of the run being processed
    int fSubRun;       ///< number of the sub-run being processed
    /// @}

    static const int kMaxAD = 30;

    /// @name The variables that will go into the simulation n-tuple.
    /// @{
    unsigned int fSimHits;
    int fSimPDG;       ///< PDG ID of the particle being processed
    int fSimProcess;
    int fSimEndProcess;
    int fSimTrackID;   ///< GEANT ID of the particle being processed
    int fNAuxDet;   ///< Number of strips hit
    vector<uint32_t> fAuxDetID;  ///< Global CRT module ID
    vector<uint32_t> fAuxDetSensitiveID; ///< Strip ID in module
    vector<float> fADEDep; ///< Energy deposited in CRT strip
    vector<float> fADdEdx; 
    vector<float> fADTrackLength; ///< Track length in strip
    vector<uint32_t> fAuxDetReg; ///< CRT region
    vector<uint32_t> fADMac;

    double fADEnterXYZT[kMaxAD][4];
    double fADExitXYZT[kMaxAD][4];
    //vector<std::array<double,4>> fADEnterXYZT;
    //vector<std::array<double,4>> fADExitXYZT;

    int fParentPDG;
    double fParentE;

    // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
    double fStartXYZT[4];//[kMaxSeg][4]; ///< (x,y,z,t) of the true start of the particle
    double fEndXYZT[4];   ///< (x,y,z,t) of the true end of the particle
    double fStartPE[4];//[kMaxSeg][4];   ///< (Px,Py,Pz,E) at the true start of the particle
    double fEndPE[4];     ///< (Px,Py,Pz,E) at the true end of the particle    

    int fSimNChan;

    int fProgenitor;
    int fMother;
    int fNDaught;    

    //CRT data product vars
    vector<uint32_t> fChan;
    vector<uint32_t> fT0;
    vector<uint32_t> fT1;
    vector<uint32_t> fADC;
    vector<uint32_t> fTrackID;
    uint32_t fNChan;
    uint32_t fEntry;
    uint32_t fTTrig;
    uint32_t fChanTrig;
    uint32_t fFEBReg;
    uint32_t fMac5;
    uint32_t fTriggerPair[2];
    uint32_t fMacPair[2];
      
    //CRT hit product vars
    vector<float> fXHit;
    vector<float> fYHit;
    vector<float> fZHit;
    vector<float> fXErrHit;
    vector<float> fYErrHit;
    vector<float> fZErrHit;
    vector<float> fT0Hit;
    vector<float> fT1Hit;
    vector<uint32_t> fHitReg;
    uint32_t fNHit;

    //truth matching stats for CRTDetSim
    uint32_t fNmuTruth;    //true N primary muons entering CRT sensive volumes
    TH1F* fStripMultHistC; //true N strips hit / C-module / muon track 
    TH1F* fStripMultHistM; //true N strips hit / M-module / muon track
    TH1F* fStripMultHistD; //true N strips hit / D-module / muon track
    TH1F* fModMultHistC;   //true N C-modules hit / muon track
    TH1F* fModMultHistM;   //true N M-modules hit / muon track
    TH1F* fModMultHistD;   //true N D-modules hit / muon track

    uint32_t fNmuTagC;     //N muon tracks producing >0 CRT triggers in C-subsystem
    uint32_t fNmuTagM;     //N muon tracks producing >0 CRT triggers in M-subsystem
    uint32_t fNmuTagD;     //N muon tracks producing >0 CRT triggers in D-subsystem
    uint32_t fNmuTagTot;
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


    //TH1F* fFEBMultHistC;

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
    LOG_DEBUG("CRT") << " starting analysis job" << '\n';

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Define our n-tuples
    fSimulationNtuple = tfs->make<TTree>("SimTree",        "MyCRTSimulation");
    fDetSimNtuple     = tfs->make<TTree>("DetTree",        "MyCRTDetSim");
    fSimHitNtuple     = tfs->make<TTree>("HitTree",        "MyCRTSimHit");

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
    fSimulationNtuple->Branch("Event",             &fEvent,             "Event/I");
    fSimulationNtuple->Branch("Hits" ,             &fSimHits,           "Hits/I");
    fSimulationNtuple->Branch("SubRun",            &fSubRun,            "SubRun/I");
    fSimulationNtuple->Branch("Run",               &fRun,               "Run/I");
    fSimulationNtuple->Branch("TrackID",           &fSimTrackID,        "TrackID/I");
    fSimulationNtuple->Branch("PDG",               &fSimPDG,            "PDG/I");
    fSimulationNtuple->Branch("NHits",             &fSimHits,           "NHits/I");
    fSimulationNtuple->Branch("Process",           &fSimProcess,        "Process/I");
    fSimulationNtuple->Branch("EndProcess",        &fSimEndProcess,     "EndProcess/I");
    fSimulationNtuple->Branch("ParentPDG",         &fParentPDG,         "ParentPDG/I");
    fSimulationNtuple->Branch("ParentE",           &fParentE,           "ParentE/D");
    fSimulationNtuple->Branch("Progenitor",        &fProgenitor,        "Progenitor/I");
    // CRT hits
    fSimulationNtuple->Branch("AuxDetSensitiveID", &fAuxDetSensitiveID);//,  "AuxDetSensitiveID[30]/I");
    fSimulationNtuple->Branch("AuxDetID",          &fAuxDetID);//,           "AuxDetID[30]/I");
    fSimulationNtuple->Branch("AuxDetEDep",        &fADEDep);//,             "AuxDetEDep[30]/D");
    fSimulationNtuple->Branch("AuxDetdEdx",        &fADdEdx);//,             "ADdEdx[30]/D");
    fSimulationNtuple->Branch("AuxDetTrackLength", &fADTrackLength);//,      "AuxDetTrackLength[30]/D");
    fSimulationNtuple->Branch("AuxDetEnterXYZT",   fADEnterXYZT,        "AuxDetEnterXYZT[30][4]/D");
    fSimulationNtuple->Branch("AuxDetExitXYZT",    fADExitXYZT,         "AuxDetExitXYZT[30][4]/D");
    fSimulationNtuple->Branch("AuxDetRegion",      &fAuxDetReg);//,          "AuxDetRegion[30]/I");
    fSimulationNtuple->Branch("Mac5",      &fADMac);

    fSimulationNtuple->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/D");
    fSimulationNtuple->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/D");
    fSimulationNtuple->Branch("StartPE",     fStartPE,         "StartPE[4]/D");
    fSimulationNtuple->Branch("EndPE",       fEndPE,           "EndPE[4]/D");
    fSimulationNtuple->Branch("NChan",       &fSimNChan, "NChan/I");
    fSimulationNtuple->Branch("Mother",            &fMother,            "Mother/I");
    fSimulationNtuple->Branch("NDaught",           &fNDaught,           "NDaught/I");

    // Define the branches of our DetSim n-tuple 
    fDetSimNtuple->Branch("NChan",       &fNChan,       "NChan/I");
    fDetSimNtuple->Branch("Event",       &fEvent,       "Event/I");
    fDetSimNtuple->Branch("Channel",     &fChan);
    fDetSimNtuple->Branch("T0",          &fT0);
    fDetSimNtuple->Branch("T1",          &fT1);
    fDetSimNtuple->Branch("ADC",         &fADC);
    fDetSimNtuple->Branch("TrackID",     &fTrackID);
    fDetSimNtuple->Branch("Entry",       &fEntry,       "Entry/I");
    fDetSimNtuple->Branch("Mac5",        &fMac5,        "Mac5/I");
    fDetSimNtuple->Branch("Region",      &fFEBReg,      "Region/I");
    fDetSimNtuple->Branch("TriggerPair", fTriggerPair,  "TriggerPair[2]/I");
    fDetSimNtuple->Branch("MacPair", fMacPair,  "MacPair[2]/I");
    fDetSimNtuple->Branch("ChanTrig",    &fChanTrig,    "ChanTrig/I");
    fDetSimNtuple->Branch("TTrig",       &fTTrig,       "TTrig/I");


    // Define the branches of our SimHit n-tuple
    fSimHitNtuple->Branch("NHit",        &fNHit,        "NHit/I");
    fSimHitNtuple->Branch("X",           &fXHit);
    fSimHitNtuple->Branch("Y",           &fYHit);
    fSimHitNtuple->Branch("Z",           &fZHit);
    fSimHitNtuple->Branch("XErr",        &fXErrHit);
    fSimHitNtuple->Branch("YErr",        &fYErrHit);
    fSimHitNtuple->Branch("ZErr",        &fZErrHit);
    fSimHitNtuple->Branch("T0",          &fT0Hit);
    fSimHitNtuple->Branch("T1",          &fT1Hit);
    fSimHitNtuple->Branch("Region",      &fHitReg);  
}
   
  void CRTAnalysis::beginRun(const art::Run& /*run*/)
  {
    //art::ServiceHandle<sim::LArG4Parameters> larParameters;
    //fElectronsToGeV = 1./larParameters->GeVToElectrons();
  }

  //-----------------------------------------------------------------------
  void CRTAnalysis::analyze(const art::Event& event) 
  {
    LOG_DEBUG("CRT") << "beginning analyis" << '\n';
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
    //map< int, vector< pair<uint32_t,uint32_t> > > muHitMap;

    fNmuTruth = 0;

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
    //std::cout << "event " << fEvent << " particleMap filled with " << particleMap.size() << " particles" << std::endl;

    LOG_DEBUG("CRT") << "about to loop over MC particles" << '\n';
    //loop over MCParticles
    for ( auto const& particle : (*particleHandle) )
    {
        fSimPDG = particle.PdgCode();
        vector<int>::iterator it = fPDGs.begin();
        const TLorentzVector& momentumStart = particle.Momentum(0);
        const double p = (momentumStart.Vect()).Mag();
        size_t index = 0;
        //if ( (particle.Process() != "primary"  && 
        while (it!=fPDGs.end()) {
          if (*it==fSimPDG) {
              index = (size_t)(it - fPDGs.begin());
              break;
          }
          it++;
        }
          //if (fSimPDG!=13 && fSimPDG!=-13)// && fSimPDG!=22 && fSimPDG!=11 && fSimPDG!=-11)
        //{  //std::cout << " Process: " <<  particle.Process() << " , PDG: " << fSimPDG << std::endl;
        if (!(fPDGs.size()==1 && fPDGs[0]==0) && it == fPDGs.end()) continue;
        if ( fMinMomenta[index] != 0 && p < fMinMomenta[index]) continue;
        if ( fMaxMomenta[index] != 0 && p > fMaxMomenta[index]) continue;
       // }
        //if ( (fSimPDG==11 || fSimPDG==-11) && particle.Process()!="compt" && particle.Process()!="conv" )
        //  continue;

        if (abs(fSimPDG)==13){
            fNmuTruth++;
        }
	fSimTrackID = particle.TrackId();

        set<uint16_t> modidsC;
        set<uint16_t> modidsM;
        set<uint16_t> modidsD;
        set<uint8_t>  stripidsC;
        set<uint8_t>  stripidsM;
        set<uint8_t>  stripidsD;

        fMother = particle.Mother();
        fNDaught = particle.NumberDaughters();
        fSimProcess = ProcessToICode(particle.Process());
        fSimEndProcess = ProcessToICode(particle.EndProcess());
       
        if(fMother!=0){ //if not primary
            map<int,const simb::MCParticle*>::const_iterator it = particleMap.find(fMother);
            if(it==particleMap.end()){
                //std::cout << "mother MCParticle object not found!" << std::endl;
                fParentPDG=INT_MAX;
                fParentE = DBL_MAX;
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

	// Use a polar-coordinate view of the 4-vectors to
	// get the track length.
	//const double trackLength = ( positionEnd - positionStart ).Rho();

	/*LOG_DEBUG("CRTAnalysis")
	  << "track ID=" << fSimTrackID 
	  << " (PDG ID: " << fSimPDG << ") "
	  << trackLength << " cm long, momentum " 
	  << momentumStart.P() << " GeV/c, has " 
	  << fSimHits << " trajectory points";
	   */ 

	fNAuxDet = 0;
        fADTrackLength.clear();
        fADEDep.clear();
        fADdEdx.clear();
        fAuxDetID.clear();
        fAuxDetSensitiveID.clear();
        //fADEnterXYZT.clear();
        //fADExitXYZT.clear();
        fAuxDetReg.clear();
        for (size_t i=0; i<kMaxAD; i++) {
            for (size_t j=0; j<4; j++) {
                fADEnterXYZT[i][j] = -9e9;
                fADExitXYZT[i][j] = -9e9;
            }
        }

	// To look at the energy deposited by this particle's track,
	// we loop over the AuxDetSimChannel objects in the event. 
	// Note all volumes are included, not just ones with energy deps
	for ( auto const& channel : (*auxDetSimChannelHandle) )
	{
            auto adGeo = fGeometryService->AuxDet(channel.AuxDetID());
            fADMac.push_back(ADToMac(ModToAuxDetType(adGeo),channel.AuxDetID()));

	    // Get vector of hits in this AuxDet channel
	    auto const& auxDetIDEs = channel.AuxDetIDEs();
	   
	    // For every hit in this channel:
	    for ( auto const& ide : auxDetIDEs )
	    {
		    // Check if the track that deposited the
		    // energy matches the track of the particle.
		    if ( ide.trackID != fSimTrackID ) continue;

                    LOG_DEBUG("CRTAnalysis")
            	    << " auxdetsimchannel: " << channel.AuxDetID() << '\n'
		    << " module name: " << (fGeometryService->AuxDet(channel.AuxDetID())).TotalVolume()->GetName() << '\n'
                    << " strip name: " << (fGeometryService->AuxDet(channel.AuxDetID()).SensitiveVolume(channel.AuxDetSensitiveID())).TotalVolume()->GetName() << '\n'
		    << " region " << GetAuxDetRegion(fGeometryService->AuxDet(channel.AuxDetID())) << '\n'
            	    << " sensID: " << channel.AuxDetSensitiveID();
 
                    //std::array<double,4> arrenter, arrexit;
		    double dx = ide.entryX-ide.exitX;
		    double dy = ide.entryY-ide.exitY;
		    double dz = ide.entryZ-ide.exitZ;

		    fADTrackLength.push_back(TMath::Sqrt(dx*dx+dy*dy+dz*dz));
	            fADEDep.push_back(ide.energyDeposited);
	            fADdEdx.push_back(ide.energyDeposited/fADTrackLength[fNAuxDet]);
		    fAuxDetID.push_back(channel.AuxDetID());
		    fAuxDetSensitiveID.push_back(channel.AuxDetSensitiveID());
                    /*arrenter[0] = ide.entryX;
                    arrenter[1] = ide.entryY;
                    arrenter[2] = ide.entryZ;
                    arrenter[3] = ide.entryT;
                    arrexit[0] = ide.exitX;
                    arrexit[1] = ide.exitY;
                    arrexit[2] = ide.exitZ;
         /           arrexit[3] = ide.exitT;
                    //fADEnterXYZT.push_back(arrenter);
                    //fADExitXYZT.push_back(arrexit);*/
		    fADEnterXYZT[fNAuxDet][0] = ide.entryX;
	            fADEnterXYZT[fNAuxDet][1] = ide.entryY;
	            fADEnterXYZT[fNAuxDet][2] = ide.entryZ;
	            fADEnterXYZT[fNAuxDet][3] = ide.entryT;
                    fADExitXYZT[fNAuxDet][0] = ide.exitX;
                    fADExitXYZT[fNAuxDet][1] = ide.exitY;
                    fADExitXYZT[fNAuxDet][2] = ide.exitZ;
                    fADExitXYZT[fNAuxDet][3] = ide.exitT;
                    fAuxDetReg.push_back(GetAuxDetRegion(fGeometryService->AuxDet(channel.AuxDetID())));
		    fNAuxDet++;

                    if (abs(fSimPDG)==13) {
                        switch (ModToAuxDetType(adGeo)) {
                          case 'c' : 
                            modidsC.insert(channel.AuxDetID());
                            stripidsC.insert(channel.AuxDetSensitiveID());
                            break;
                          case 'm' : 
                            modidsM.insert(channel.AuxDetID());
                            stripidsM.insert(channel.AuxDetSensitiveID());
                            break;
                          case 'd' :
                            modidsD.insert(channel.AuxDetID());
                            stripidsD.insert(channel.AuxDetSensitiveID());
                            break;
                        }//switch
                    }//if muon
	      } // For each IDE (strip hit by muon)
              
              if( abs(fSimPDG)==13 ) {
                  fStripMultHistC->Fill(stripidsC.size());
                  fStripMultHistM->Fill(stripidsM.size());
                  fStripMultHistD->Fill(stripidsD.size());
                  stripidsC.clear();
                  stripidsM.clear();
                  stripidsD.clear();
              }//if muon

	  } // For each SimChannel (module)

        // write values to tree for this event and particle
        fSimulationNtuple->Fill();

        fModMultHistC->Fill(modidsC.size());
        fModMultHistM->Fill(modidsM.size());
        fModMultHistD->Fill(modidsD.size());

      } // loop over all particles in the event. 

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
     LOG_DEBUG("CRT") << "about to loop over detsim entries" << '\n';
     for ( auto const& febdat : (*crtDetSimHandle) ) {
        fMac5 = febdat.Mac5();
        fChanTrig = febdat.ChanTrig();
        fEntry = febdat.Entry();
        fTTrig = febdat.TTrig();
        std::pair<uint32_t,uint32_t> tmpPair = febdat.TrigPair();
        fTriggerPair[0] = tmpPair.first;
        fTriggerPair[1] = tmpPair.second;
        tmpPair = febdat.MacPair();
        fMacPair[0] = tmpPair.first;
        fMacPair[1] = tmpPair.second;
        fFEBReg = MacToADReg(fMac5);
        fNChan = 0;

        char type = MacToType(fMac5);
        set<uint8_t> chanidsC;
        set<uint8_t> chanidsM;
        set<uint8_t> chanidsD;
 
        for ( auto const chandat : febdat.ChanData()) {
          //DetSim tree contains all entries (not neccessarily from muons)
          fNChan++;
          fChan.push_back(chandat.Channel());
          fT0.push_back(chandat.T0());
          fT1.push_back(chandat.T1());
          fADC.push_back(chandat.ADC());
          fTrackID.push_back(chandat.TrackID()); 

          //if channel hit came from muon loop over all associated hits
          if( abs(particleMap[chandat.TrackID()]->PdgCode())==13 ) {

              taggedTrks.insert(chandat.TrackID());
              for (auto const chandat2 : febdat.ChanData()) {
                  if( chandat.TrackID() != chandat2.TrackID() ) continue;
                  switch( type ){
                    case 'c' :
                      febidsC.insert(fMac5);
                      chanidsC.insert(chandat2.Channel());
                      taggedTrksC.insert(chandat.TrackID());
                      break;
                    case 'm' :
                      febidsM.insert(fMac5);
                      chanidsM.insert(chandat2.Channel());
                      taggedTrksM.insert(chandat.TrackID());
                      break;
                    case 'd' :
                      febidsD.insert(fMac5);
                      chanidsD.insert(chandat2.Channel());
                      taggedTrksD.insert(chandat.TrackID());
                      break;
                  }//switch
              }//inner chandat loop
              switch (type) {
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
        }//outer chandat loop

        fDetSimNtuple->Fill();

     } //for CRT FEB events

     fFEBMultHistC->Fill(febidsC.size());
     fFEBMultHistM->Fill(febidsM.size());
     fFEBMultHistD->Fill(febidsD.size());
     uint32_t nmissC=0, nmissM=0, nmissD=0, nmissTot=0;

     for (auto pmit = particleMap.begin(); pmit!=particleMap.end(); pmit++) {
       if ( abs((*pmit).second->PdgCode())==13 ) {
         if (taggedTrksC.find((*pmit).first)==taggedTrksC.end()) nmissC++;
         if (taggedTrksM.find((*pmit).first)==taggedTrksM.end()) nmissM++;
         if (taggedTrksD.find((*pmit).first)==taggedTrksD.end()) nmissD++;
         if (taggedTrks.find((*pmit).first) ==taggedTrks.end())  nmissTot++; 
       }//if MC particle is muon, get trackID
     }//for MC particles

     fNmuTagC   = fNmuTruth - nmissC;
     fNmuTagC   = fNmuTruth - nmissM;
     fNmuTagC   = fNmuTruth - nmissD;
     fNmuTagTot = fNmuTruth - nmissTot;
     fEffC = 1.0*fNmuTagC/fNmuTruth;
     fEffM = 1.0*fNmuTagM/fNmuTruth;
     fEffD = 1.0*fNmuTagD/fNmuTruth;
     fEffTot = 1.0*fNmuTagTot/fNmuTruth;

     mf::LogInfo("CRT") << '\n'
       << " Total muon tracks: " << fNmuTruth << '\n'
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
            fXHit.push_back(hit.X());
            fYHit.push_back(hit.Y());
            fZHit.push_back(hit.Z());
            fXErrHit.push_back(hit.XErr());
            fYErrHit.push_back(hit.YErr());
            fZErrHit.push_back(hit.ZErr());
            fT0Hit.push_back(hit.T0());
            fT1Hit.push_back(hit.T1());
            fHitReg.push_back(hit.Region());

        }//for CRT Hits

        fSimHitNtuple->Fill();

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
      if(mac>=0   && mac<=17 ) return 50; //left
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

  char MacToType(uint32_t mac) {

      uint32_t reg = MacToADReg(mac);

      if( reg==38 || reg==52 || reg==56 || reg==48 || reg==46 )
        return 'c';
      if( reg==50 || reg==54 || reg==44 || reg==42 )
        return 'm';
      if( reg==58)
        return 'd';

      return 'e';
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
