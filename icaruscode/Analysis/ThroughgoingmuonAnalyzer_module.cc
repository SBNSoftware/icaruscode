/////////////////////////////////////////////////////////////////////////////////
// Class:       ThroughgoingmuonAnalyzer
// Plugin Type: analyzer (art v3_02_06)
// File:        ThroughgoingmuonAnalyzer_module.cc
//
// This is based on the Tracy's TrackHitEfficiency module.
//
// Generated at Thu Nov 14 15:21:00 2019 by Biswaranjan Behera using cetskelgen
// from cetlib version v3_07_02.
/////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"


#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/Simulation/LArG4Parameters.h"

// Eigen
#include <Eigen/Dense>

//Root
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"

#include <cmath>
#include <algorithm>

namespace thrugoingmuon {

class ThroughgoingmuonAnalyzer : public art::EDAnalyzer 
{
public:
  explicit ThroughgoingmuonAnalyzer(fhicl::ParameterSet const& pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ThroughgoingmuonAnalyzer(ThroughgoingmuonAnalyzer const&) = delete;
  ThroughgoingmuonAnalyzer(ThroughgoingmuonAnalyzer&&) = delete;
  ThroughgoingmuonAnalyzer& operator=(ThroughgoingmuonAnalyzer const&) = delete;
  ThroughgoingmuonAnalyzer& operator=(ThroughgoingmuonAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;
  
private:
  
  
  // Declare member data here.
  // Fcl parameters.
  std::vector<art::InputTag>  fRawDigitProducerLabelVec;
  std::vector<art::InputTag>  fWireProducerLabelVec;
  std::vector<art::InputTag>  fHitProducerLabelVec;
  std::vector<art::InputTag>  fTrackProducerLabelVec;
  art::InputTag               fMCParticleProducerLabel;
  art::InputTag               fSimChannelProducerLabel;
  art::InputTag               fSimEnergyProducerLabel;
  art::InputTag               fBadChannelProducerLabel;
  bool                        fUseBadChannelDB;
  bool                        fUseICARUSGeometry;
  std::vector<int>            fOffsetVec;              ///< Allow offsets for each plane
  std::vector<float>          fSigmaVec;               ///< Window size for matching to SimChannels
  int                         fMinAllowedChanStatus;   ///< Don't consider channels with lower status
  float                       ftotalelctroncut;        ///< Threshold for electron (10k)

  //  double                   fElectronsToGeV;    ///< conversion factor
  // TTree variables
  TTree*             fTree;
  int fEvent;        ///< number of the event being processed
  int fRun;          ///< number of the run being processed
  int fSubRun;       ///< number of the sub-run being processed
  /// @}

  /// @name The variables that will go into the simulation n-tuple.
  /// @{
  int fSimPDG;       ///< PDG ID of the particle being processed
  int fSimTrackID;   ///< GEANT ID of the particle being processed
  int fNtracks_reco;
  int fNtracks_primary;

  //mutable std::vector<float> fHitSummedADCVec;
  mutable std::vector<int>   fTPCVec;
  mutable std::vector<int>   fCryoVec;
  mutable std::vector<int>   fPlaneVec;
  mutable std::vector<int>   fWireVec;

  mutable std::vector<float>   fmcpartvx;
  mutable std::vector<float>   fmcpartvy;
  mutable std::vector<float>   fmcpartvz;

  mutable std::vector<float>   fmcpartpx;
  mutable std::vector<float>   fmcpartpy;
  mutable std::vector<float>   fmcpartpz;
  mutable std::vector<float>   fmcparte;

  mutable std::vector<float> fTotalElectronsVec;
  mutable std::vector<float> fTotalElecEnergyVec;
  mutable std::vector<float> fMaxElectronsVec;
  mutable std::vector<int>   fStartTickVec;
  mutable std::vector<int>   fStopTickVec;
  mutable std::vector<int>   fMaxETickVec;
  mutable std::vector<float> fPartDirX;
  mutable std::vector<float> fPartDirY;
  mutable std::vector<float> fPartDirZ;

  mutable std::vector<int>   fNMatchedWires;
  mutable std::vector<int>   fNMatchedHits;

  mutable std::vector<float> fHitPeakTimeVec;
  mutable std::vector<float> fHitPeakAmpVec;
  mutable std::vector<float> fHitPeakRMSVec;
  mutable std::vector<float> fHitBaselinevec;
  mutable std::vector<float> fHitSummedADCVec;
  mutable std::vector<float> fHitIntegralVec;
  mutable std::vector<int>   fHitStartTickVec;
  mutable std::vector<int>   fHitStopTickVec;
  mutable std::vector<int>   fHitMultiplicityVec;
  mutable std::vector<int>   fHitLocalIndexVec;
  mutable std::vector<float> fHitGoodnessVec;
  mutable std::vector<int>   fNumDegreesVec;

  mutable std::vector<float>   fmcstartx;
  mutable std::vector<float>   fmcstarty;
  mutable std::vector<float>   fmcstartz;

  mutable std::vector<float>   fmcendx;
  mutable std::vector<float>   fmcendy;
  mutable std::vector<float>   fmcendz;

  mutable std::vector<float>   fmcstartdirx;
  mutable std::vector<float>   fmcstartdiry;
  mutable std::vector<float>   fmcstartdirz;

  mutable std::vector<float>   fmcenddirx;
  mutable std::vector<float>   fmcenddiry;
  mutable std::vector<float>   fmcenddirz;

  mutable std::vector<float>   fmctstartx;
  mutable std::vector<float>   fmctstarty;
  mutable std::vector<float>   fmctstartz;

  mutable std::vector<float>   fmctendx;
  mutable std::vector<float>   fmctendy;
  mutable std::vector<float>   fmctendz;

  mutable std::vector<float>   fmctstartdirx;
  mutable std::vector<float>   fmctstartdiry;
  mutable std::vector<float>   fmctstartdirz;

  mutable std::vector<float>   fmctenddirx;
  mutable std::vector<float>   fmctenddiry;
  mutable std::vector<float>   fmctenddirz;

  mutable std::vector<float>   frecostartx;
  mutable std::vector<float>   frecostarty;
  mutable std::vector<float>   frecostartz;

  mutable std::vector<float>   frecoendx;
  mutable std::vector<float>   frecoendy;
  mutable std::vector<float>   frecoendz;

  mutable std::vector<float>   frecostartdirx;
  mutable std::vector<float>   frecostartdiry;
  mutable std::vector<float>   frecostartdirz;

  mutable std::vector<float>   frecoenddirx;
  mutable std::vector<float>   frecoenddiry;
  mutable std::vector<float>   frecoenddirz;

  mutable std::vector<float>   fLength;
  mutable std::vector<float>   fThetaXZ;
  mutable std::vector<float>   fThetaYZ;
  mutable std::vector<float>   fTheta;
  mutable std::vector<float>   fPhi;

  mutable std::vector<float>   fmcThetaXZ;
  mutable std::vector<float>   fmcThetaYZ;
  mutable std::vector<float>   fmcTheta;
  mutable std::vector<float>   fmcPhi;

  mutable std::vector<float>   fmctThetaXZ;
  mutable std::vector<float>   fmctThetaYZ;
  mutable std::vector<float>   fmctTheta;
  mutable std::vector<float>   fmctPhi;

  mutable std::vector<float>   fNSimChannelHitsVec;
  mutable std::vector<float>   fNRecobHitVec;
  mutable std::vector<float>   fHitEfficiencyVec;
  mutable std::vector<float>   fNFakeHitVec;

  mutable std::vector<float>   flocx;
  mutable std::vector<float>   flocy;
  mutable std::vector<float>   flocz;

  std::vector<TProfile*>      fHitEffvselecVec;

  // Implementation of required member function here.
  art::ServiceHandle< art::TFileService > tfs;

  // How to convert from number of electrons to GeV. The ultimate
  // source of this conversion factor is
  // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h.
  // But sim::LArG4Parameters might in principle ask a database for it.
  //art::ServiceHandle<sim::LArG4Parameters const> larParameters;


  // Useful services, keep copies for now (we can update during begin run periods)
  const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
  const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
  const detinfo::DetectorClocks*     fClockService;         ///< Detector clocks service

  // Get geometry.
  //  art::ServiceHandle<geo::Geometry> geom;
};


ThroughgoingmuonAnalyzer::ThroughgoingmuonAnalyzer(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}  // ,
// More initializers here.
{ 
  std::cout << "In analyzer constructor" << std::endl;

  // Call appropriate consumes<>() for any products to be retrieved by this module.
    fRawDigitProducerLabelVec = pset.get< std::vector<art::InputTag>>("RawDigitLabelVec",   std::vector<art::InputTag>() = {"rawdigitfilter"});
    fWireProducerLabelVec     = pset.get< std::vector<art::InputTag>>("WireModuleLabelVec", std::vector<art::InputTag>() = {"decon1droi"});
    fHitProducerLabelVec      = pset.get< std::vector<art::InputTag>>("HitModuleLabelVec",  std::vector<art::InputTag>() = {"gaushit"});
    fTrackProducerLabelVec    = pset.get< std::vector<art::InputTag>>("TrackModuleLabelVec",std::vector<art::InputTag>() = {"pandoraTrackGausCryo0"});
    fMCParticleProducerLabel  = pset.get< art::InputTag             >("MCParticleLabel",    "largeant");
    fSimChannelProducerLabel  = pset.get< art::InputTag             >("SimChannelLabel",    "largeant");
    fSimEnergyProducerLabel   = pset.get< art::InputTag             >("SimEnergyLabel",     "ionization");
    fBadChannelProducerLabel  = pset.get< art::InputTag             >("BadChannelLabel",    "simnfspl1:badchannels");
    fUseBadChannelDB          = pset.get< bool                      >("UseBadChannelDB",    false);
    fUseICARUSGeometry        = pset.get< bool                      >("UseICARUSGeometry",  false);
  //  fLocalDirName             = pset.get<std::string                >("LocalDirName",       std::string("wow"));
    fOffsetVec                = pset.get<std::vector<int>           >("OffsetVec",          std::vector<int>()={0,0,0});
    fSigmaVec                 = pset.get<std::vector<float>         >("SigmaVec",           std::vector<float>()={1.,1.,1.});
    fMinAllowedChanStatus     = pset.get< int                       >("MinAllowedChannelStatus");
    ftotalelctroncut          = pset.get<float                      >("totalelctroncut", 10000.);

    std::cout << "UseICARUSGeometry: " << fUseICARUSGeometry << std::endl;
    std::cout << "Read the fhicl parameters, recovering services" << std::endl;

  //  fHitProducerLabelVec      = pset.get< std::vector<art::InputTag>>("HitModuleLabelVec",  std::vector<art::InputTag>() = {"gauss"});

    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fClockService       = lar::providerFrom<detinfo::DetectorClocksService>();

    art::TFileDirectory dir = tfs->mkdir("histos");
    fHitEffvselecVec.resize(fGeometry->Nplanes());
    for(size_t plane = 0; plane < fGeometry->Nplanes(); plane++)
      {
	fHitEffvselecVec.at(plane)           = dir.make<TProfile>  (("HitEffvselec"    + std::to_string(plane)).c_str(), "Hit Efficiency; Total # e-",        200,  0., 200000., 0., 1.);
      }

    std::cout << "Services recovered, setting up output tree" << std::endl;

    fTree     = tfs->make<TTree>("Analysis",    "Analysis tree");

    fTree->Branch("Event",           &fEvent,           "Event/I");
    fTree->Branch("SubRun",          &fSubRun,          "SubRun/I");
    fTree->Branch("Run",             &fRun,             "Run/I");
    fTree->Branch("TrackID",         &fSimTrackID,      "TrackID/I");
    fTree->Branch("PDG",             &fSimPDG,          "PDG/I");
    fTree->Branch("ntracks_reco",    &fNtracks_reco,    "ntracks_reco/I");
    fTree->Branch("ntracks_primary", &fNtracks_primary, "ntracks_primary/I");

    fTree->Branch("CryostataVec",      "std::vector<int>",   &fCryoVec);
    fTree->Branch("TPCVec",            "std::vector<int>",   &fTPCVec);
    fTree->Branch("PlaneVec",          "std::vector<int>",   &fPlaneVec);
    fTree->Branch("WireVec",           "std::vector<int>",   &fWireVec);

    fTree->Branch("mcpartvx",           "std::vector<float>",   &fmcpartvx);
    fTree->Branch("mcpartvy",           "std::vector<float>",   &fmcpartvy);
    fTree->Branch("mcpartvz",           "std::vector<float>",   &fmcpartvz);

    fTree->Branch("mcpartpx",           "std::vector<float>",   &fmcpartpx);
    fTree->Branch("mcpartpy",           "std::vector<float>",   &fmcpartpy);
    fTree->Branch("mcpartpz",           "std::vector<float>",   &fmcpartpz);
    fTree->Branch("mcparte",            "std::vector<float>",   &fmcparte);

    fTree->Branch("TotalElectronsVec", "std::vector<float>", &fTotalElectronsVec);
    fTree->Branch("TotalElecEnergyVec","std::vector<float>", &fTotalElecEnergyVec);
    fTree->Branch("MaxElectronsVec",   "std::vector<float>", &fMaxElectronsVec);
    fTree->Branch("StartTick",         "std::vector<int>",   &fStartTickVec);
    fTree->Branch("StopTick",          "std::vector<int>",   &fStopTickVec);
    fTree->Branch("MaxETick",          "std::vector<int>",   &fMaxETickVec);
    fTree->Branch("PartDirX",          "std::vector<float>", &fPartDirX);
    fTree->Branch("PartDirY",          "std::vector<float>", &fPartDirY);
    fTree->Branch("PartDirZ",          "std::vector<float>", &fPartDirZ);

    fTree->Branch("NMatchedWires",     "std::vector<int>",   &fNMatchedWires);
    fTree->Branch("NMatchedHits",      "std::vector<int>",   &fNMatchedHits);

    fTree->Branch("HitPeakTimeVec",    "std::vector<float>", &fHitPeakTimeVec);
    fTree->Branch("HitPeakAmpVec",     "std::vector<float>", &fHitPeakAmpVec);
    fTree->Branch("HitPeakRMSVec",     "std::vector<float>", &fHitPeakRMSVec);
    fTree->Branch("HitBaselineVec",    "std::vector<float>", &fHitBaselinevec);
    fTree->Branch("HitSummedADCVec",   "std::vector<float>", &fHitSummedADCVec);
    fTree->Branch("HitIntegralVec",    "std::vector<float>", &fHitIntegralVec);
    fTree->Branch("HitStartTickVec",   "std::vector<int>",   &fHitStartTickVec);
    fTree->Branch("HitStopTickVec",    "std::vector<int>",   &fHitStopTickVec);
    fTree->Branch("HitMultiplicity",   "std::vector<int>",   &fHitMultiplicityVec);
    fTree->Branch("HitLocalIndex",     "std::vector<int>",   &fHitLocalIndexVec);
    fTree->Branch("HitGoodness",       "std::vector<float>", &fHitGoodnessVec);
    fTree->Branch("HitNumDegrees",     "std::vector<int>",   &fNumDegreesVec);

    fTree->Branch("NSimChannelHitsVec", "std::vector<float>",   &fNSimChannelHitsVec);
    fTree->Branch("NRecobHitVec",       "std::vector<float>",   &fNRecobHitVec);
    fTree->Branch("HitEfficiencyVec",   "std::vector<float>",   &fHitEfficiencyVec);
    fTree->Branch("NFakeHitVec",        "std::vector<float>",   &fNFakeHitVec);

    fTree->Branch("mcstartx",           "std::vector<float>", &fmcstartx);
    fTree->Branch("mcstarty",           "std::vector<float>", &fmcstarty);
    fTree->Branch("mcstartz",           "std::vector<float>", &fmcstartz);
    fTree->Branch("mcendx",             "std::vector<float>", &fmcendx);
    fTree->Branch("mcendy",             "std::vector<float>", &fmcendy);
    fTree->Branch("mcendz",             "std::vector<float>", &fmcendz);

    fTree->Branch("mcstartdirx",        "std::vector<float>", &fmcstartdirx);
    fTree->Branch("mcstartdiry",        "std::vector<float>", &fmcstartdiry);
    fTree->Branch("mcstartdirz",        "std::vector<float>", &fmcstartdirz);
    fTree->Branch("mcenddirx",          "std::vector<float>", &fmcenddirx);
    fTree->Branch("mcenddiry",          "std::vector<float>", &fmcenddiry);
    fTree->Branch("mcenddirz",          "std::vector<float>", &fmcenddirz);

    fTree->Branch("mctstartx",          "std::vector<float>", &fmctstartx);
    fTree->Branch("mctstarty",          "std::vector<float>", &fmctstarty);
    fTree->Branch("mctstartz",          "std::vector<float>", &fmctstartz);
    fTree->Branch("mctendx",            "std::vector<float>", &fmctendx);
    fTree->Branch("mctendy",            "std::vector<float>", &fmctendy);
    fTree->Branch("mctendz",            "std::vector<float>", &fmctendz);

    fTree->Branch("mctstartdirx",       "std::vector<float>", &fmctstartdirx);
    fTree->Branch("mctstartdiry",       "std::vector<float>", &fmctstartdiry);
    fTree->Branch("mctstartdirz",       "std::vector<float>", &fmctstartdirz);
    fTree->Branch("mctenddirx",         "std::vector<float>", &fmctenddirx);
    fTree->Branch("mctenddiry",         "std::vector<float>", &fmctenddiry);
    fTree->Branch("mctenddirz",         "std::vector<float>", &fmctenddirz);

    fTree->Branch("recostartx",         "std::vector<float>", &frecostartx);
    fTree->Branch("recostarty",         "std::vector<float>", &frecostarty);
    fTree->Branch("recostartz",         "std::vector<float>", &frecostartz);
    fTree->Branch("recoendx",           "std::vector<float>", &frecoendx);
    fTree->Branch("recoendy",           "std::vector<float>", &frecoendy);
    fTree->Branch("recoendz",           "std::vector<float>", &frecoendz);

    fTree->Branch("recostartdirx",      "std::vector<float>", &frecostartdirx);
    fTree->Branch("recostartdiry",      "std::vector<float>", &frecostartdiry);
    fTree->Branch("recostartdirz",      "std::vector<float>", &frecostartdirz);
    fTree->Branch("recoenddirx",        "std::vector<float>", &frecoenddirx);
    fTree->Branch("recoenddiry",        "std::vector<float>", &frecoenddiry);
    fTree->Branch("recoenddirz",        "std::vector<float>", &frecoenddirz);

    fTree->Branch("length",             "std::vector<float>", &fLength);
    fTree->Branch("thetaxz",            "std::vector<float>", &fThetaXZ);
    fTree->Branch("thetayz",            "std::vector<float>", &fThetaYZ);
    fTree->Branch("theta",              "std::vector<float>", &fTheta);
    fTree->Branch("phi",                "std::vector<float>", &fPhi);
    fTree->Branch("mctheta",            "std::vector<float>", &fmcTheta);
    fTree->Branch("mcphi",              "std::vector<float>", &fmcPhi);
    fTree->Branch("mcthetaxz",          "std::vector<float>", &fmcThetaXZ);
    fTree->Branch("mcthetayz",          "std::vector<float>", &fmcThetaYZ);
    fTree->Branch("mcttheta",           "std::vector<float>", &fmctTheta);
    fTree->Branch("mctphi",             "std::vector<float>", &fmctPhi);
    fTree->Branch("mctthetaxz",         "std::vector<float>", &fmctThetaXZ);
    fTree->Branch("mctthetayz",         "std::vector<float>", &fmctThetaYZ);

    fTree->Branch("locx",         "std::vector<float>", &flocx);
    fTree->Branch("locy",         "std::vector<float>", &flocy);
    fTree->Branch("locz",         "std::vector<float>", &flocz);

    std::cout << "Returning" << std::endl;
    return;
}

void ThroughgoingmuonAnalyzer::analyze(art::Event const& evt)
{
    // Implementation of required member function here.
    // Always clear the tuple
    // clear();

    fNtracks_reco    = 0;
    fNtracks_primary = 0;
    
    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();
    fWireVec.clear();

    fmcpartvx.clear();
    fmcpartvy.clear();
    fmcpartvz.clear();

    fmcpartpx.clear();
    fmcpartpy.clear();
    fmcpartpz.clear();
    fmcparte.clear();

    fTotalElectronsVec. clear();
    fTotalElecEnergyVec.clear();
    fMaxElectronsVec.   clear();
    fStartTickVec.      clear();
    fStopTickVec.       clear();
    fMaxETickVec.       clear();
    fPartDirX.          clear();
    fPartDirY.          clear();
    fPartDirZ.          clear();

    fNMatchedWires.     clear();
    fNMatchedHits.      clear();

    fHitPeakTimeVec.    clear();
    fHitPeakAmpVec.     clear();
    fHitPeakRMSVec.     clear();
    fHitBaselinevec.    clear();
    fHitSummedADCVec.   clear();
    fHitIntegralVec.    clear();
    fHitStartTickVec.   clear();
    fHitStopTickVec.    clear();
    fHitMultiplicityVec.clear();
    fHitLocalIndexVec.  clear();
    fHitGoodnessVec.    clear();
    fNumDegreesVec.     clear();
    fmcstartx.          clear();
    fmcstarty.          clear();
    fmcstartz.          clear();
    fmcendx.            clear();
    fmcendy.            clear();
    fmcendz.            clear();
    fmcstartdirx.       clear();
    fmcstartdiry.       clear();
    fmcstartdirz.       clear();
    fmctstartx.         clear();
    fmctstarty.         clear();
    fmctstartz.         clear();
    fmctendx.           clear();
    fmctendy.           clear();
    fmctendz.           clear();
    fmctstartdirx.      clear();
    fmctstartdiry.      clear();
    fmctstartdirz.      clear();

    flocx.     clear();
    flocy.     clear();
    flocz.     clear();

    fmcenddirx.clear();
    fmcenddiry.clear();
    fmcenddirz.clear();
    frecostartx.clear();
    frecostarty.clear();
    frecostartz.clear();
    frecoendx.clear();
    frecoendy.clear();
    frecoendz.clear();
    frecostartdirx.clear();
    frecostartdiry.clear();
    frecostartdirz.clear();
    frecoenddirx.  clear();
    frecoenddiry.  clear();
    frecoenddirz.  clear();
    fLength. clear();
    fThetaXZ.clear();
    fThetaYZ.clear();
    fTheta.  clear();
    fPhi.    clear();
    fmcThetaXZ.clear();
    fmcThetaYZ.clear();
    fmcTheta.  clear();
    fmcPhi.    clear();
    fmctThetaXZ.clear();
    fmctThetaYZ.clear();
    fmctTheta.  clear();
    fmctPhi.    clear();

    fNSimChannelHitsVec.clear();
    fNRecobHitVec.      clear();
    fHitEfficiencyVec.  clear();
    fNFakeHitVec.       clear();

    fHitEffvselecVec.   clear();

  //*/
  // Start by fetching some basic event information for our n-tuple.
    fEvent  = evt.id().event();
    fRun    = evt.run();
    fSubRun = evt.subRun();

    int ntrk_reco = 0;

    art::Handle< std::vector<sim::SimChannel>> simChannelHandle;
    evt.getByLabel(fSimChannelProducerLabel, simChannelHandle);

    
    art::Handle< std::vector<simb::MCParticle>> mcParticleHandle;
    evt.getByLabel(fMCParticleProducerLabel, mcParticleHandle);

    art::Handle<std::vector<sim::SimEnergyDeposit>> simEnergyHandle;
    evt.getByLabel(fSimEnergyProducerLabel, simEnergyHandle);

    //fElectronsToGeV = 1./larParameters->GeVToElectrons();

    // If there is no sim channel informaton then exit
    if (!simChannelHandle.isValid() || simChannelHandle->empty() || 
        !simEnergyHandle.isValid()  || simEnergyHandle->empty()  ||
        !mcParticleHandle.isValid() ) return;

    /*  std::cout << " Event: " << fEvent
	      << " Run: "   << fRun
	      << " SubRun: "<< fSubRun<< std::endl;
    */
  
    // There are several things going on here... for each channel we have particles (track id's) depositing energy in a range to ticks
    // So... for each channel we want to build a structure that relates particles to tdc ranges and deposited energy (or electrons)
    // Here is a complicated structure:
    using TDCToIDEMap             = std::map<unsigned short, sim::IDE>; // We need this one in order
    using ChanToTDCToIDEMap       = std::map<raw::ChannelID_t, TDCToIDEMap>;
    using PartToChanToTDCToIDEMap = std::unordered_map<int, ChanToTDCToIDEMap>;

    PartToChanToTDCToIDEMap partToChanToTDCToIDEMap;
  
    // Build out the above data structure
    for(const auto& simChannel : *simChannelHandle)
    {
        for(const auto& tdcide : simChannel.TDCIDEMap())
        {
            for(const auto& ide : tdcide.second) partToChanToTDCToIDEMap[ide.trackID][simChannel.Channel()][tdcide.first] = ide;
        }
    }

    // Here we make a map between track ID and associatied SimEnergyDeposit objects
    // We'll need this for sorting out the track direction at each hit
    using SimEnergyDepositVec = std::vector<const sim::SimEnergyDeposit*>;
    using PartToSimEnergyMap  = std::unordered_map<int, SimEnergyDepositVec>;

    PartToSimEnergyMap partToSimEnergyMap;

    for(const auto& simEnergy : *simEnergyHandle)
    {
        partToSimEnergyMap[simEnergy.TrackID()].push_back(&simEnergy);
    }

    // Now we create a data structure to relate hits to their channel ID
    using ChanToHitVecMap = std::unordered_map<raw::ChannelID_t,std::vector<const recob::Hit*>>;

    ChanToHitVecMap channelToHitVec;

    // And now fill it
    for(const auto& hitLabel : fHitProducerLabelVec)
    {
        art::Handle<std::vector<recob::Hit>> hitHandle;
        evt.getByLabel(hitLabel, hitHandle);

        for(const auto& hit : *hitHandle)
            channelToHitVec[hit.Channel()].push_back(&hit);
    }

    // We need a structure to relate hits to MCParticles and their trajectory ID
    using MCParticleTrajIdxPair = std::pair<const simb::MCParticle*,size_t>;
    using MCPartPointPairIdxVec = std::vector<MCParticleTrajIdxPair>;
    using HitToMCPartToIdxMap   = std::unordered_map<const recob::Hit*,MCPartPointPairIdxVec>;

    HitToMCPartToIdxMap hitToMcPartToIdxMap;
    
    // It is useful to create a mapping between trackID and MCParticle
    using TrackIDToMCParticleMap = std::unordered_map<int, const simb::MCParticle*>;
    
    TrackIDToMCParticleMap trackIDToMCParticleMap;
    
    for(const auto& mcParticle : *mcParticleHandle)
        trackIDToMCParticleMap[mcParticle.TrackId()] = &mcParticle;
    
    const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

    // Look up the list of bad channels
    art::Handle< std::vector<int>> badChannelHandle;
    evt.getByLabel(fBadChannelProducerLabel, badChannelHandle);
    
    std::vector<int> nSimChannelHitVec  = {0,0,0};
    std::vector<int> nRecobHitVec       = {0,0,0};
    std::vector<int> nFakeHitVec        = {0,0,0};
    std::vector<int> nSimulatedWiresVec = {0,0,0};

    std::cout << "***************** EVENT " << fEvent << " ******************" << std::endl;
    std::cout << "-- Looping over channels for hit efficiency, # MC Track IDs: " << partToChanToTDCToIDEMap.size() << std::endl;

    // Initiate loop over MC Track IDs <--> SimChannel information by wire and TDC
    for(const auto& partToChanInfo : partToChanToTDCToIDEMap)
    {
        TrackIDToMCParticleMap::const_iterator trackIDToMCPartItr = trackIDToMCParticleMap.find(partToChanInfo.first);
    
        if (trackIDToMCPartItr == trackIDToMCParticleMap.end()) continue;
    
        const simb::MCParticle* mcParticle = trackIDToMCPartItr->second;

        int         trackPDGCode = mcParticle->PdgCode();
        std::string processName  = mcParticle->Process();

          // Looking for primary muons (e.g. CR Tracks)
          //	  fTotalElectronsVec.push_back(totalElectrons);
        if (fabs(trackPDGCode) != 13 || processName != "primary") continue;
          //      if (fabs(trackPDGCode) != 13) continue;
	//if (fEvent != 44) continue;
        std::cout << ">>>> Loop on MCParticle: " << mcParticle << ", pdg: " << trackPDGCode << ", process: " << processName << std::endl;

        // Let's recover the SimEnergyDeposit vector for this track
//        PartToSimEnergyMap::iterator simEneDepItr = partToSimEnergyMap.find(partToChanInfo.first);

//        if (simEneDepItr == partToSimEnergyMap.end())
//        {
//            std::cout << "No match for SimEnergyDeposit" << std::endl;
//            continue;
//        }

        // We should sort this vector by time, this will facilitate searching for matches to SimChannels later
//        SimEnergyDepositVec& simEnergyDepositVec = simEneDepItr->second;

//        std::sort(simEnergyDepositVec.begin(),simEnergyDepositVec.end(),[](const auto& left,const auto& right){return left->T() < right->T();});

//        std::cout << "  -- Processing track id: " << partToChanInfo.first << ", with " << simEnergyDepositVec.size() << " SimEnergyDeposit objects" << ", # channels: " << partToChanInfo.second.size() << std::endl;
      
        fSimPDG = trackPDGCode;
        fSimTrackID = mcParticle->TrackId();

	fmcpartvx.push_back(mcParticle->Vx());
	fmcpartvy.push_back(mcParticle->Vy());
	fmcpartvz.push_back(mcParticle->Vz());

	fmcpartpx.push_back(mcParticle->Px());
	fmcpartpy.push_back(mcParticle->Py());
	fmcpartpz.push_back(mcParticle->Pz());
	fmcparte.push_back(mcParticle->E());
	//std::cout << "Energy: " << mcParticle->E() << std::endl;
	// Recover particle position and angle information
        Eigen::Vector3f partStartPos(mcParticle->Vx(),mcParticle->Vy(),mcParticle->Vz());
        Eigen::Vector3f partStartDir(mcParticle->Px(),mcParticle->Py(),mcParticle->Pz());
    
        partStartDir.normalize();
    
        Eigen::Vector2f partStartDirVecXZ(partStartDir[0],partStartDir[2]);
    
        partStartDirVecXZ.normalize();
    
          // Assuming the SimChannels contain position information (currently not true for WC produced SimChannels)
          // then we want to keep a running position
        std::vector<Eigen::Vector3f> lastPositionVec = {partStartPos,partStartPos,partStartPos};

        std::cout << "  *** Looping over channels, have " << partToChanInfo.second.size() << " channels" << std::endl;

        for(const auto& chanToTDCToIDEMap : partToChanInfo.second)
        {
    	        // skip bad channels
            if (fUseBadChannelDB)
            {
    	          // This is the "correct" way to check and remove bad channels...
                if( chanFilt.Status(chanToTDCToIDEMap.first) < fMinAllowedChanStatus)
                {
                    std::vector<geo::WireID> wids = fGeometry->ChannelToWire(chanToTDCToIDEMap.first);
                    std::cout << "*** skipping bad channel with status: " << chanFilt.Status(chanToTDCToIDEMap.first) 
                          << " for channel: "                         << chanToTDCToIDEMap.first 
                          << ", plane: "                              << wids[0].Plane 
                          << ", wire: "                               << wids[0].Wire    << std::endl;
                          continue;
                }
            }
      
    	      // Was a list made available?
    	      // If so then we try that
            if (badChannelHandle.isValid())
            {
                std::vector<int>::const_iterator badItr = std::find(badChannelHandle->begin(),badChannelHandle->end(),chanToTDCToIDEMap.first);
    
                if (badItr != badChannelHandle->end()) continue;
            }
     
            TDCToIDEMap    tdcToIDEMap = chanToTDCToIDEMap.second;
            float          totalElectrons(0.);
            float          totalEnergy(0.);
            float          maxElectrons(0.);
            unsigned short maxElectronsTDC(0);
            int            nMatchedWires(0);
            int            nMatchedHits(0);
            
        	     // The below try-catch block may no longer be necessary
        	     // Decode the channel and make sure we have a valid one
            std::vector<geo::WireID> wids = fGeometry->ChannelToWire(chanToTDCToIDEMap.first);
        
        	  // Recover plane and wire in the plane
            unsigned int plane = wids[0].Plane;
            //unsigned int wire  = wids[0].Wire;

            Eigen::Vector3f avePosition(0.,0.,0.);
            
	    nSimulatedWiresVec[plane]++;  // Loop insures channels are unique here
            
            for(const auto& ideVal : tdcToIDEMap)
            {
                totalElectrons += ideVal.second.numElectrons;
                totalEnergy    += ideVal.second.energy;
        
                if (maxElectrons < ideVal.second.numElectrons)
                {
                    maxElectrons    = ideVal.second.numElectrons;
                    maxElectronsTDC = ideVal.first;
                }
        
                avePosition += Eigen::Vector3f(ideVal.second.x,ideVal.second.y,ideVal.second.z);
            }
        
    	      // Get local track direction by using the average position of deposited charge as the current position
    	      // and then subtracting the last position
            avePosition /= float(tdcToIDEMap.size());
    
            Eigen::Vector3f partDirVec = avePosition - lastPositionVec[plane];
    
            partDirVec.normalize();
            
            lastPositionVec[plane] = avePosition;

	    // std::cout << "sim wire: "<<  nSimulatedWiresVec[plane]++ << ", sim ch hit : " <<nSimChannelHitVec[plane]++ << std::endl;
	    // Threshold for electrons
	    if (totalElectrons > ftotalelctroncut)          nSimChannelHitVec[plane]++;
        
	    //  std::cout <<" total electron cut : \t" << ftotalelctroncut << std::endl;     
	    //  totalElectrons = std::min(totalElectrons, float(99900.));
            
            //nSimChannelHitVec[plane]++;
            	  //	  std::cout << "after the check --- line 469   "<< nSimChannelHitVec[plane] << std::endl; 	    
            unsigned short startTDC = tdcToIDEMap.begin()->first;
            unsigned short stopTDC  = tdcToIDEMap.rbegin()->first;
            
            	  // Convert to ticks to get in same units as hits
            unsigned short startTick = fClockService->TPCTDC2Tick(startTDC)        + fOffsetVec[plane];
            unsigned short stopTick  = fClockService->TPCTDC2Tick(stopTDC)         + fOffsetVec[plane];
            unsigned short maxETick  = fClockService->TPCTDC2Tick(maxElectronsTDC) + fOffsetVec[plane];
    
    	      //fSimNumTDCVec[plane]->Fill(stopTick - startTick, 1.);
    	      //fSimNumTDCVec[plane]->Fill(stopTick - startTick, 1.);

            // Let's do the match to the SimChannel now just to see that it is working
//          Eigen::Vector3f simChanPos(ideVal.second.x,ideVal.second.y,ideVal.second.z);
//    
//          float                        minDist = std::numeric_limits<float>::max();
//          const sim::SimEnergyDeposit* simEDepMatch(0);
//    
//          for(const auto& simEDep : simEnergyDepositVec)
//          {
//              Eigen::Vector3f simEnePos(simEDep->MidPointX(),simEDep->MidPointY(),simEDep->MidPointZ());
//    
//              float sepDist = (simEnePos - simChanPos).norm();
//    
//              if (sepDist < minDist)
//              {
//                  minDist      = sepDist;
//                  simEDepMatch = simEDep;
//              }
//              else if (sepDist > 1.5 * minDist) break;
//          }
//
//          if (!simEDepMatch || minDist > 0.15)
//          {
//              std::cout << "Bad match, simEDepMatch: " << simEDepMatch << ", minDist: " << minDist << std::endl;
//              continue;
//          }
//
//          // Recalculate the angles of the track at this point
//          Eigen::Vector3f newPartDirVec(simEDepMatch->EndX()-simEDepMatch->StartX(),
//                                        simEDepMatch->EndY()-simEDepMatch->StartY(),
//                                        simEDepMatch->EndZ()-simEDepMatch->StartZ());
//
//          newPartDirVec.normalize();
//
//          // Overwrite the angles
//          partDirVec     = newPartDirVec;
//          projPairDirVec = Eigen::Vector2f(partDirVec[0],partDirVec[2]);
//
//          projPairDirVec.normalize();

    	      // Set up to extract the "best" parameters in the event of more than one hit for this pulse train
            float          nElectronsTotalBest(0.);
            float          hitSummedADCBest(0.);
            float          hitIntegralBest(0.);
            float          hitPeakTimeBest(0.);
            float          hitPeakAmpBest(-100.);
            float          hitRMSBest(0.);
            int            hitMultiplicityBest(0);
            int            hitLocalIndexBest(0);
            float          hitGoodnessBest(0.);
            int            hitNumDegreesBest(0);
            float          hitBaselineBest(0.);
            	  //float          hitSnippetLenBest(0.);
            unsigned short hitStopTickBest(0);
            unsigned short hitStartTickBest(0);
            
            const recob::Hit* rejectedHit = 0;
            const recob::Hit* bestHit     = 0;
             // The next mission is to recover the hits associated to this Wire
    		    // The easiest way to do this is to simply look up all the hits on this channel and then match
            ChanToHitVecMap::iterator hitIter = channelToHitVec.find(chanToTDCToIDEMap.first);

            if (hitIter != channelToHitVec.end())
            {
    		        // Loop through the hits for this channel and look for matches
    		        // In the event of more than one hit associated to the sim channel range, keep only
    		        // the best match (assuming the nearby hits are "extra")
    		        // Note that assumption breaks down for long pulse trains but worry about that later
                for(const auto& hit : hitIter->second)
                {
		  //std::cout <<" ........... sigma vector: \t " << fSigmaVec[plane] << std::endl;
                    unsigned short hitStartTick = hit->PeakTime() - fSigmaVec[plane] * hit->RMS();
                    unsigned short hitStopTick  = hit->PeakTime() + fSigmaVec[plane] * hit->RMS();
    		           	  // If hit is out of range then skip, it is not related to this particle
                    if (hitStartTick > stopTick || hitStopTick < startTick)
                    {
                        nFakeHitVec[plane]++;
                        rejectedHit = hit;
                        continue;
                    }
               
                    float hitHeight = hit->PeakAmplitude();
               
    		           	  // Use the hit with the largest pulse height as the "best"
                    if (hitHeight < hitPeakAmpBest) continue;
               
                    hitPeakAmpBest   = hitHeight;
                    bestHit          = hit;
                    hitStartTickBest = hitStartTick;
                    hitStopTickBest  = hitStopTick;
                } 
    		        // Find a match?
                if (bestHit)
                {
                    nElectronsTotalBest = 0.;
                    hitPeakTimeBest     = bestHit->PeakTime();
                    hitIntegralBest     = bestHit->Integral();
                    hitSummedADCBest    = bestHit->SummedADC();
                    hitRMSBest          = bestHit->RMS();
                    hitMultiplicityBest = bestHit->Multiplicity();
                    hitLocalIndexBest   = bestHit->LocalIndex();
                    hitGoodnessBest     = bestHit->GoodnessOfFit();
                    hitNumDegreesBest   = bestHit->DegreesOfFreedom();
		    //hitSnippetLenBest   = bestHit->EndTick() - bestHit->StartTick();
		    hitBaselineBest     = 0.;  // To do...
             
		    nMatchedHits++;
		    
		    // Get the number of electrons
		    for(unsigned short tick = hitStartTickBest; tick <= hitStopTickBest; tick++)
		      {
                        unsigned short hitTDC = fClockService->TPCTick2TDC(tick - fOffsetVec[plane]);
			
                        TDCToIDEMap::iterator ideIterator = tdcToIDEMap.find(hitTDC);
			
                        if (ideIterator != tdcToIDEMap.end()) nElectronsTotalBest += ideIterator->second.numElectrons;
		      }
                    // Ok, now we need to figure out which trajectory point this hit is associated to 
                    // Use the associated IDE to get the x,y,z position for this hit
                    Eigen::Vector3f hitIDEPos(avePosition[0],avePosition[1],avePosition[2]);
                    size_t          matchIdx(0);
                    float           bestMatch(std::numeric_limits<float>::max());

                    for (size_t idx = 0; idx < mcParticle->NumberTrajectoryPoints(); idx++) 
                    {
                        const TLorentzVector& pos = mcParticle->Position(idx); // 4-position in World coordinates
                        Eigen::Vector3f trackPos(pos.X(),pos.Y(),pos.Z());
                        float           missDist = (trackPos - hitIDEPos).norm();
                        if (missDist < bestMatch)
                        {
                            bestMatch = missDist;
                            matchIdx  = idx;
                        }
                    }

//                    std::cout << "  --> Hit: " << bestHit << ", traj idx: " << matchIdx << " with dist: " << bestMatch << std::endl;
//                    if (matchIdx > 0)
//                    {
//                        const TLorentzVector& pos  = mcParticle->Position(matchIdx); // 4-position in World coordinates
//                        const TLorentzVector& last = mcParticle->Position(matchIdx-1); // 4-position in World coordinates
//                   float dLastPos = (Eigen::Vector3f(pos.X(),pos.Y(),pos.Z()) - Eigen::Vector3f(last.X(),last.Y(),last.Z())).norm();
//                   std::cout << "      --> distance to last trajectory point: " << dLastPos << std::endl;
//                    }
                    if (bestMatch < std::numeric_limits<float>::max())
                    {
                        hitToMcPartToIdxMap[bestHit].push_back(MCParticleTrajIdxPair(mcParticle,matchIdx));
                        //std::cout << "-hit " << bestHit  << " associated to MC " << mcParticle << " at point " << matchIdx << std::endl;
    
                        partDirVec = Eigen::Vector3f(mcParticle->Px(matchIdx),mcParticle->Py(matchIdx),mcParticle->Pz(matchIdx));
    
                        partDirVec.normalize();
                    }
       	        } // end of besthit
             
                if (nMatchedHits > 0)
                {
    	         	    //float chgRatioADC = hitSummedADCBest > 0. ? totalElectrons / hitSummedADCBest : 0.;
    	         	    //float chgRatioInt = hitIntegralBest  > 0. ? totalElectrons / hitIntegralBest  : 0.;
		  if (totalElectrons > ftotalelctroncut) nRecobHitVec[plane]++;    
		  // nRecobHitVec[plane]++;
                }
                else if (rejectedHit)
                {
    	         	    //unsigned short hitStartTick = rejectedHit->PeakTime() - fSigmaVec[plane] * rejectedHit->RMS();
    	         	    //unsigned short hitStopTick  = rejectedHit->PeakTime() + fSigmaVec[plane] * rejectedHit->RMS();
                    std::cout << "    -Rejected: " << rejectedHit->WireID().Plane << ", ph: " << rejectedHit->PeakAmplitude() << std::endl;
                }
                else
                {
                    mf::LogDebug("ThroughgoingmuonAnalyzer") << "==> No match, TPC/Plane/Wire: " 
                    << "/"               << wids[0].TPC 
                    << "/"               << wids[0].Plane 
                    << "/"               << wids[0].Wire 
                    << ", # electrons: " << totalElectrons 
                    << ", startTick: "   << startTick 
                    << ", stopTick: "    << stopTick << std::endl;
                } // end of matched hit
            } // end of channel to hit map iterator
                
	    float matchHit   = std::min(nMatchedHits,1);
	    fHitEffvselecVec[plane]->Fill(totalElectrons,   matchHit,   1.);

            fTPCVec.push_back(wids[0].TPC);
            fCryoVec.push_back(wids[0].Cryostat);
            fPlaneVec.push_back(wids[0].Plane);
            fWireVec.push_back(wids[0].Wire);
            
            fTotalElectronsVec.push_back(totalElectrons);
            fTotalElecEnergyVec.push_back(totalEnergy);
            fMaxElectronsVec.push_back(maxElectrons);
            fStartTickVec.push_back(startTick);
            fStopTickVec.push_back(stopTick);
            fMaxETickVec.push_back(maxETick);
            fPartDirX.push_back(partDirVec[0]);
            fPartDirY.push_back(partDirVec[1]);
            fPartDirZ.push_back(partDirVec[2]);
            
            fNMatchedWires.push_back(nMatchedWires);
            fNMatchedHits.push_back(nMatchedHits);
            
            fHitPeakTimeVec.push_back(hitPeakTimeBest);
            fHitPeakAmpVec.push_back(hitPeakAmpBest);
            fHitPeakRMSVec.push_back(hitRMSBest);
            fHitBaselinevec.push_back(hitBaselineBest);
            fHitSummedADCVec.push_back(hitSummedADCBest);
            fHitIntegralVec.push_back(hitIntegralBest);
            fHitStartTickVec.push_back(hitStartTickBest);
            fHitStopTickVec.push_back(hitStopTickBest);
            fHitMultiplicityVec.push_back(hitMultiplicityBest);
            fHitLocalIndexVec.push_back(hitLocalIndexBest);
            fHitGoodnessVec.push_back(hitGoodnessBest);
            fNumDegreesVec.push_back(hitNumDegreesBest);
        } //end of partcle to id map info
    } // Done looping over track <--> MCParticle associations

    std::cout << "-- Now starting loop over tracks" << std::endl;
    std::cout << "-- Hit/MCParticle map size: " << hitToMcPartToIdxMap.size() << std::endl;

    // Let's keep track of track to MCParticle and hits. 
    using MCPartToHitVecMap = std::unordered_map<const simb::MCParticle*,std::vector<const recob::Hit*>>;
    using TrackToMCMap      = std::unordered_map<const recob::Track*, MCPartToHitVecMap>;

    TrackToMCMap trackToMCMap;

    std::cout << "** Starting outer loop over all sets of tracks" << std::endl;

    // The game plan for this module is to look at recob::Tracks and objects associated to tracks
    // The first thing we need to do is match tracks to MCParticles which we can do using structures 
    // created in the above loops over hits. 
    for(const auto& trackLabel : fTrackProducerLabelVec)
    {
        art::Handle<std::vector<recob::Track> > trackHandle;
         //std::vector<art::Ptr<recob::Track> > tracklist;
        evt.getByLabel(trackLabel, trackHandle);
              if (!trackHandle.isValid()) return;
     
         //track information
        ntrk_reco     += trackHandle->size();  
        fNtracks_reco  = ntrk_reco;

         // Recover the collection of associations between tracks and hits
        art::FindMany<recob::Hit> trackHitAssns(trackHandle, evt, trackLabel);
     
         // First mission will be to find the tracks associated to our primary MCParticle...
        for(size_t trackIdx = 0; trackIdx < trackHandle->size(); trackIdx++)
        {
            art::Ptr<recob::Track> ptrack(trackHandle,trackIdx);

            // Associated hits
            const std::vector<const recob::Hit*>& trackHitVec = trackHitAssns.at(ptrack.key());

            // Keep track of the number of hits per MCParticle for this track
            MCPartToHitVecMap mcPartCountMap;

            std::cout << "- Track idx " << trackIdx << ", ptr: " << ptrack.get() << " has " << trackHitVec.size() << " hits" << std::endl;

            // Loop through hits and use associations tables to relate to MCParticle
            for(const auto& hit : trackHitVec)
            {
                // Get list of MCParticles/traj point ids
                HitToMCPartToIdxMap::iterator hitToMcPartIdxItr = hitToMcPartToIdxMap.find(hit);

                if (hitToMcPartIdxItr != hitToMcPartToIdxMap.end())
                {
                    for(const auto& mcPartIdxPair : hitToMcPartIdxItr->second)
                    {
                        mcPartCountMap[mcPartIdxPair.first].push_back(hit);
                    }
                }
                else 
                    mcPartCountMap[nullptr].push_back(hit);
            }

            std::cout << "  - associated to " << mcPartCountMap.size() << " MCParticles" << std::endl;

	    for(const auto& mcPair : mcPartCountMap) 
	      std::cout << "    - MCParticle: " << mcPair.first <<  ", # hits: " << mcPair.second.size() << std::endl;

            // Find the best match (simple majority logic)
            const simb::MCParticle* bestMCMatch(nullptr);
            size_t                  bestCount(0);

            for(const auto& pair : mcPartCountMap)
            {
	      //std::cout << "    - MCParticle: " << pair.first <<", pdg: " << pair.first->PdgCode() << ", # hits: " << pair.second.size() << std::endl;
                if (pair.first && pair.second.size() > bestCount)
                {
                    bestMCMatch = pair.first;
                    bestCount   = pair.second.size();
                }
            }

	    //            std::cout << "    ==> bestMCMatch: " << bestMCMatch << ", pdg:" << bestMCMatch->PdgCode()<< ", bestCount: " << bestCount << std::endl;
	    //std::cout << "    ==> bestMCMatch: " << bestMCMatch << ", bestCount: " << bestCount << std::endl;

            if (bestCount > 0)
            {
                MCPartToHitVecMap::iterator bestItr = mcPartCountMap.find(bestMCMatch);

                trackToMCMap[ptrack.get()][bestItr->first] = bestItr->second;

	    }
        }
    }

    std::cout << "====> trackToMCMap size: " << trackToMCMap.size() << std::endl;
    //    std::cout << "====> mcPartCountMap size: " <<  MCPartToHitVecMap.size() << std::endl;

    fNtracks_primary = trackToMCMap.size();
    std::cout << "==================================================================================> " << std::endl;
    if (fNtracks_primary != 1)
        std::cout << "**> tracks associated to primary = " << fNtracks_primary << ", event: " << evt.id().event() << std::endl;


    // Go through the tracks matched to the primary
    for(const auto& trackMCPair : trackToMCMap)
    {
      const recob::Track& track = *trackMCPair.first;

      //	      std::cout << track.ParticleId()<< std::endl;
      unsigned int ntraj = track.NumberTrajectoryPoints();
      //std::cout << "ntraj: \t"<< ntraj << std::endl; 
      
      if(ntraj <= 0) continue;
      

      //  const geo::TPCGeo* tpcGeo = fGeometry->PositionToTPCptr(track.Start());

      //if (!tpcGeo) continue;

      //std::cout << tpcGeo->ActiveBoundingBox().ContainsPosition(track.Start()) << std::endl;
      //std::cout << "track length : "<< track.Length() << std::endl;


      //      const recob::tracking::Point_t& pos
      //	      TVector3 pos = track.Vertex<TVector3>();
      //TVector3 pos = track.Start<TVector3>();
      //TVector3 dir = track.StartDirection<TVector3>();
      //TVector3 end = track.End<TVector3>();
      // TVector3 enddir = track.EndDirection<TVector3>();
      //    double   pstart = track.VertexMomentum();
      //double   pend   = track.EndMomentum();		
      
      //	      double length   = track.Length();
      // double theta_xz = std::atan2(dir.X(), dir.Z());
      //double theta_yz = std::atan2(dir.Y(), dir.Z());
      //	      std::cout << end.X() << "\t" << end.Y() << "\t" << end.Z()<<std::endl;
      //	      TVector3 pos(0,2.43493,851.593);
      //	      if(!cryo0.ContainsPosition(pos) && !cryo1.ContainsPosition(pos)) continue;
      //	      std::cout << trackLabel << "  no.of track:\t" <<trackHandle->size() << "\t"<< ntracks_reco<< std::endl;
      // This seems a redundant test for fit tracks? 
      //  geo::Point_t trackPoint(pos.X(),pos.Y(),pos.Z());
      //      std::cout << track.X() << "\t" << track.Y() << "\t" << track.Z()<<std::endl;      

      //std::cout << "Reco Track--- Before: \t" << pos.X() <<"\t"<< pos.Y() <<"\t"<< pos.Z() << std::endl;      
      // geo::Point_t const trackPoint(track.Start());
      //geo::Point_t const trackEndPoint(track.End());

      //const geo::TPCGeo* tpcGeo = fGeometry->PositionToTPCptr(trackPoint);
      
      //if (!tpcGeo) continue;
      //if (!tpcGeo->ActiveBoundingBox().ContainsPosition(trackPoint)) continue; // out of active volume
      //if (!tpcGeo->ActiveBoundingBox().ContainsPosition(trackEndPoint)) continue; // out of active volume
    
      // std::cout << "Reco Track start point : \t" << trackPoint.X() <<"\t"<< trackPoint.Y() <<"\t"<< trackPoint.Z() << std::endl;
      
     

      //      if (tpcGeo->ActiveBoundingBox().ContainsPosition(trackPoint)) trackEndPoint = trackPoint;
      //else break; // this means that once you are out of the active volume the first time, you are done; think if that's what you want, 
                  //since tracks scatter a lot in LAr
      //      std::cout << "Reco Track start point : \t" << trackPoint.X() <<"\t"<< trackPoint.Y() <<"\t"<< trackPoint.Z() << std::endl;
      //std::cout << "Reco Track end point : \t" << trackEndPoint.X() <<"\t"<< trackEndPoint.Y() <<"\t"<< trackEndPoint.Z() << std::endl;
      //      std::cout <<: \t" << end.X() <<"\t"<< end.Y() <<"\t"<< end.Z() << std::endl;            

      //geo::Point_t  trackstartPoint ;
      //      geo::Point_t  trackendPoint ;
      std::vector<double> posx, posy, posz;
      std::vector<double> dirx, diry, dirz;
      //      const geo::CryostatID* cryo;
      //      std::cout <<      geo::CryostatID(0)<< std::endl;

	//      std::cout << "cryostat: " << cryostat << std::endl;
      const geo::CryostatID C0id { 0 };
      // ...
      //std::cout << "cryostat: " << C0id << std::endl;
      for (unsigned int i= 0; i < ntraj; i++)
	{
	  //geo::Point_t  trackstartPoint = track.LocationAtPoint(i);
	  //geo::Point_t  trackendPoint = track.LocationAtPoint(i);
	  TVector3 mom = track.MomentumVectorAtPoint<TVector3>(i);
	  geo::Point_t trackstartPoint = track.LocationAtPoint(i);
	  //trackendPoint = track.LocationAtPoint(i);
	  const geo::TPCGeo* tpcGeom = fGeometry->PositionToTPCptr(trackstartPoint);

	  if (!tpcGeom) continue;
	  if (tpcGeom->ID() != C0id) continue; // point not in cryostat 0
	  if (!tpcGeom->ActiveBoundingBox().ContainsPosition(trackstartPoint)) continue;
	  //if (tpcGeom->ActiveBoundingBox().ContainsPosition(trackstartPoint)) trackendPoint = trackstartPoint;
	  //else break; 
	  posx.push_back( trackstartPoint.X() );
	  posy.push_back( trackstartPoint.Y() );
	  posz.push_back( trackstartPoint.Z() );

	  dirx.push_back( mom.X() );
          diry.push_back( mom.Y() );
          dirz.push_back( mom.Z() );

	  //  if(i == 10) std::cout<< "track length inside loop: "<< track.Length() << std::endl;
	}

      //std::cout << "Reco Track start point : \t" << posx.front() <<"\t"<< posy.front() <<"\t"<< posz.front() << std::endl;
      //std::cout << "Reco Track end point : \t" << posx.back() <<"\t"<< posy.back() <<"\t"<< posz.back() << std::endl;
      

      // frecostartx.push_back( pos.X() );
      // frecostarty.push_back( pos.Y() );
      // frecostartz.push_back( pos.Z() );
      
      // frecoendx.push_back( end.X() );
      // frecoendy.push_back( end.Y() );
      // frecoendz.push_back( end.Z() );
      if (!posx.empty())
	{
	  frecostartx.push_back( posx.front() );
	  frecostarty.push_back( posy.front() );
	  frecostartz.push_back( posz.front() );
	  
	  frecoendx.push_back( posx.back() );
	  frecoendy.push_back( posy.back() );
	  frecoendz.push_back( posz.back() );

	  frecostartdirx.push_back( dirx.front() );
          frecostartdiry.push_back( diry.front() );
          frecostartdirz.push_back( dirz.front() );

          frecoenddirx.push_back( dirx.back() );
          frecoenddiry.push_back( diry.back() );
          frecoenddirz.push_back( dirz.back() );

          double theta_xz = std::atan2(dirx.front(), dirz.front());
          double theta_yz = std::atan2(diry.front(), dirz.front());

          fThetaXZ.push_back( theta_xz );
          fThetaYZ.push_back( theta_yz );
	}
      //      frecostartdirx.push_back( dir.X() );
      // frecostartdiry.push_back( dir.Y() );
      // frecostartdirz.push_back( dir.Z() );
      
      // frecoenddirx.push_back( enddir.X() );
      //frecoenddiry.push_back( enddir.Y() );
      //frecoenddirz.push_back( enddir.Z() );
      
      //      fLength. push_back( track.Length() );
      // fThetaXZ.push_back( theta_xz );
      //fThetaYZ.push_back( theta_yz );

      fLength. push_back( track.Length() );
      fTheta.  push_back( track.Theta() );
      fPhi.    push_back( track.Phi() );

      posx.clear(); posy.clear(); posz.clear();
      dirx.clear(); diry.clear(); dirz.clear();
      //		}
      //std::cout << pos.X() << "\t" << pos.Y() << "\t" << pos.Z()<<std::endl;
      //	std::cout << end.X() << "\t" << end.Y() << "\t" << end.Z()<<std::endl;
      
      //std::cout << "-- Now starting loop over trajectory points for track " << track.ID() << std::endl;
      
      // Now start going through the associated MCParticles, focus is on the primary for now
      for(const auto& mcPartHitVecPair : trackMCPair.second)
        {
	  const simb::MCParticle* mcParticle = mcPartHitVecPair.first;
	  
	  // Skip the unassociated hits for now
	  if (!mcParticle) continue;
	  

	  geo::Point_t vtxPoint(mcParticle->Vx(), mcParticle->Vy(), mcParticle->Vz());

	  const geo::TPCGeo* tpcGeo = fGeometry->PositionToTPCptr(vtxPoint);

	  if (!tpcGeo) continue;
	  if (tpcGeo->ID() != C0id) continue; // point not in cryostat 0
	  if (!tpcGeo->ActiveBoundingBox().ContainsPosition(vtxPoint)) continue; // out of active volume

	  // Find the extreme trajectory points associated to hits for this track
	  size_t minTrajIdx = std::numeric_limits<size_t>::max();
	  size_t maxTrajIdx = 0;
	  
	  for(const auto& hit : mcPartHitVecPair.second)
            {
	      //  std::cout << " # Hit: "<< mcPartHitVecPair.second.size()<< "; Hit Index: " << &hit - &mcPartHitVecPair.second[0] << std::endl;

	      geo::Point_t loc = track.LocationAtPoint(&hit - &mcPartHitVecPair.second[0]);

	      //      std::cout << "location at x: " << loc.X() << " ; location at y: " << loc.Y()<< " ; location at z: " << loc.Z() << std::endl;	      

	      flocx.push_back(loc.X());
	      flocy.push_back(loc.Y());
	      flocz.push_back(loc.Z());
	      
	      // **THIS NEEDS TO BE CHECKED **
	      // There should be only one entry for the trajectory index the way the code is set up to run now 
	      // (Only looking at the primary)
	      size_t trajIdx = hitToMcPartToIdxMap[hit].front().second;
	      
	      minTrajIdx = std::min(minTrajIdx,trajIdx);
	      maxTrajIdx = std::max(maxTrajIdx,trajIdx);
            }
          
	  //  std::cout << "   # MCParticle trajectory points: " << mcParticle->NumberTrajectoryPoints() << ", max: " << maxTrajIdx << ", min: " << minTrajIdx << std::endl;

	  if (mcParticle->Px(maxTrajIdx) == 0 && mcParticle->Py(maxTrajIdx) == 0 && mcParticle->Pz(maxTrajIdx) == 0) maxTrajIdx -= 1;	  
	  //if (mcParticle->Px(minTrajIdx) == 0 && mcParticle->Py(minTrajIdx) == 0 && mcParticle->Pz(minTrajIdx) == 0){minTrajIdx = minTrajIdx - 1;}	  
	  
	  // Fill MC info for trajectory point limits
	 
	  TVector3 minPointDir(mcParticle->Px(minTrajIdx),mcParticle->Py(minTrajIdx),mcParticle->Pz(minTrajIdx));
	  TVector3 maxPointDir(mcParticle->Px(maxTrajIdx),mcParticle->Py(maxTrajIdx),mcParticle->Pz(maxTrajIdx));
	  
	  // Px,Py,Pz are the momentum of the track, "SetMag(1.)" is meant to normalize so you get a direction vector.
	  // if one of the two trajectory points somehow has zero momentum, avoid those cases.
	  if (!(maxPointDir.Mag() > 0.))  {maxTrajIdx -= 1; maxPointDir=TVector3(mcParticle->Px(maxTrajIdx),mcParticle->Py(maxTrajIdx),mcParticle->Pz(maxTrajIdx));}
	  if (!(minPointDir.Mag() > 0.))  {minTrajIdx -= 1; minPointDir=TVector3(mcParticle->Px(minTrajIdx),mcParticle->Py(minTrajIdx),mcParticle->Pz(minTrajIdx));}
	  //	  if (!(maxPointDir.Mag() > 0.)) maxTrajIdx -= 1;

	  //std:: cout << " minx: " << minPointDir.X() <<", miny: " << minPointDir.Y() <<", minz: " << minPointDir.Z() <<", minMag: " << minPointDir.Mag() 
	  //	     << ", maxx: " << maxPointDir.X() <<", maxy: " << maxPointDir.Y() <<", maxz: " << maxPointDir.Z() <<", maxMag: " << maxPointDir.Mag() << std::endl;	  


	  minPointDir.SetMag(1.);
	  maxPointDir.SetMag(1.);

	  fmctstartx.push_back(mcParticle->Vx(minTrajIdx));
	  fmctstarty.push_back(mcParticle->Vy(minTrajIdx));
	  fmctstartz.push_back(mcParticle->Vz(minTrajIdx));
          
	  fmctendx.push_back(mcParticle->Vx(maxTrajIdx));
	  fmctendy.push_back(mcParticle->Vy(maxTrajIdx));
	  fmctendz.push_back(mcParticle->Vz(maxTrajIdx));
          
	  fmctstartdirx.push_back(minPointDir[0]);
	  fmctstartdiry.push_back(minPointDir[1]);
	  fmctstartdirz.push_back(minPointDir[2]);
          
	  fmctenddirx.push_back(maxPointDir[0]);
	  fmctenddiry.push_back(maxPointDir[1]);
	  fmctenddirz.push_back(maxPointDir[2]);
	  
	  fmctThetaXZ.push_back( std::atan2(minPointDir.X(), minPointDir.Z()) );
	  fmctThetaYZ.push_back( std::atan2(minPointDir.Y(), minPointDir.Z()) );
          
	  fmctTheta.  push_back( minPointDir.Theta() );
	  fmctPhi.    push_back( minPointDir.Phi()   );

	  //fmcstartx.push_back(mcParticle->Vx());
	  //double earliest_time = 1e15;
	  // TVector point;
	  //std::cout << "Before: \t" <<partStartPos[0] <<"\t"<< partStartPos[1] <<"\t"<<  partStartPos[2] << std::endl;
	  std::vector<double> x, y, z, t;
	  std::vector<double> px, py, pz, e;
	  TLorentzVector dir_start;
	  //      std::vector<float> endx, endy, endz, endt;
	  //int last = (int)trackIDToMCPartItr->second->NumberTrajectoryPoints() -1;
	  //std::vector<double> vtime;
	  for (unsigned int i=0; i<mcParticle->NumberTrajectoryPoints(); i++) 
            {
	      const TLorentzVector& pos = mcParticle->Position(i); // 4-position in World coordinates
	      const TLorentzVector& mom = mcParticle->Momentum(i); // 4-position in World coordinates
	      //point(pos.X(),pos.Y(),pos.Z()
	      //const double point[3] = {pos.X(),pos.Y(),pos.Z()};

	      //std::cout << "Mc Part--- Before: \t" << pos.X() <<"\t"<< pos.Y() <<"\t"<< pos.Z() << std::endl;      
	      // if (mcParticle->Vx() < -894.951 && mcParticle->Vx() > 894.951) std::cout <<"---------------->>>>>>>>>>>>>>  "  << mcParticle->Vx() << std::endl;

	      geo::Point_t trackPoint(pos.X(),pos.Y(),pos.Z());
	      
	      const geo::TPCGeo* tpcGeo = fGeometry->PositionToTPCptr(trackPoint);
	      
	      if (!tpcGeo) continue;
	      if (tpcGeo->ID() != C0id) continue; // point not in cryostat 0
	      if (!tpcGeo->ActiveBoundingBox().ContainsPosition(trackPoint)) continue; // out of active volume

	      //  std::cout << "McPart--- After: " << tpcGeo->ID() << " - " << pos.X() <<"\t"<< pos.Y() <<"\t"<< pos.Z() << std::endl;    
	      //std::cout << "Mc Part--- After: \t" << pos.X() <<"\t"<< pos.Y() <<"\t"<< pos.Z() << std::endl;      

	      // x.push_back(mcParticle->Vx());
	      x.push_back(pos.X());
	      y.push_back(pos.Y());
	      z.push_back(pos.Z());
	      t.push_back(pos.T());
              
	      px.push_back(mom.X());
	      py.push_back(mom.Y());
	      pz.push_back(mom.Z());
	      e.push_back(mom.T());
	      
            } // no. of trajectory loop
              //      for (unsigned int i =0; i< x.size(); i++){
              //std::cout << x[i] << std::endl;
              //}
	  
	  // std::cout << "   ==> filling overall MC arrays, size of x: " << x.size() << std::endl;
	  //std::cout << "---- front -------- x: " << x.front() << ", y: " << y.front() << ",z: " << z.front() << std::endl;
	  //std::cout << "---- back ---------- x: " << x.back() << ", y: " << y.back() << ",z: " << z.back() << std::endl;

            // Make sure something got added
	  if (!x.empty())
            {
	      fmcstartx.push_back(x.front());
	      fmcstarty.push_back(y.front());
	      fmcstartz.push_back(z.front());
              
	      fmcendx.push_back(x.back());
	      fmcendy.push_back(y.back());
	      fmcendz.push_back(z.back());
	      
	      fmcstartdirx.push_back(px.front());
	      fmcstartdiry.push_back(py.front());
	      fmcstartdirz.push_back(pz.front());
              
	      fmcenddirx.push_back(px.back());
	      fmcenddiry.push_back(py.back());
	      fmcenddirz.push_back(pz.back());
              
	      dir_start.SetPxPyPzE(px.front(),py.front(),pz.front(),e.front());      
	      
	      fmcThetaXZ.push_back( std::atan2(px.front(), pz.front()) );
	      fmcThetaYZ.push_back( std::atan2(py.front(), pz.front()) );
	      
	      fmcTheta.  push_back( dir_start.Theta() );
	      fmcPhi.    push_back( dir_start.Phi()   );
            }

	  x.clear(); y.clear(); z.clear(); t.clear();
          px.clear(); py.clear(); pz.clear(); e.clear();
        } // Done with loop over MCParticles associted to this track
    }
    
    //    float diffex = std::min(std::abs((*recostartx)[0] - (*mcendx)[0]), std::abs((*recoendx)[0]-(*mcendx)[0]));
    //float diffey = std::min(std::abs(frecostarty[0] - fmcendy[0]), std::abs(frecoendy[0]-fmcendy[0]));
    // float diffez = std::min(std::abs((*recostartz)[0] - (*mcendz)[0]), std::abs((*recoendz)[0]-(*mcendz)[0]));
    //if (diffey < 34 && diffey > 32)
    //std::cout << diffey << std::endl;
    /*
    if (fNtracks_primary == 1){ 
    if (!fmcstartx.empty()){
      if (!frecostartx.empty()){
	//	float diffx = std::min(std::abs(frecostartx[0] - fmcstartx[0]), std::abs(frecoendx[0]-fmcstartx[0]));
    

	float diffxer = std::abs(frecostartx[0] - fmcstartx[0]);
	float diffyer = std::abs(frecostarty[0] - fmcstarty[0]);
	float diffzer = std::abs(frecostartz[0] - fmcstartz[0]);

	float diffxsr = std::abs(frecoendx[0] - fmcendx[0]);
	float diffysr = std::abs(frecoendy[0] - fmcendy[0]);
	float diffzsr = std::abs(frecoendz[0] - fmcendz[0]);

	float diffdr = std::min(  std::sqrt(std::pow(diffxsr,2)+std::pow(diffysr,2)+std::pow(diffzsr,2)) , 
				  std::sqrt(std::pow(diffxer,2)+std::pow(diffyer,2)+std::pow(diffzer,2))  );
	if  (diffdr <= 1.04){
	  if ( frecostartx[0] < -268.49 ||  frecostartx[0] > -171.94){
	    std::cout <<"------> outside fiducial volume ------->mcstartx: " << fmcstartx[0] <<" ; recoendx: "
		      << frecoendx[0] << "; mcendx: "<< fmcendx[0] <<", recostartx: "<< frecostartx[0]
		      << ", (reco - true) startx: "<<  frecostartx[0] - fmcstartx[0] <<" ; mcstarty: " << fmcstarty[0] <<" ; recoendy: "
                      << frecoendy[0] << "; mcendy: "<< fmcendy[0] <<", recostarty: "<< frecostarty[0]
		      << ", (reco - true) starty: "<<  frecostarty[0] - fmcstarty[0] << " ; mcstartz: " << fmcstartz[0] <<" ; recoendz: "
                      << frecoendz[0] << "; mcendz: "<< fmcendz[0] <<", recostartz: "<< frecostartz[0]
                      << ", (reco - true) startz: "<<  frecostartz[0] - fmcstartz[0] << " ; checking minm: "<< diffdr << std::endl;
	  // std::cout <<"-----------------------------> outside fiducial volume ........ mcstartx: " << (*mcstarty)[0] <<" ; recoendx: "
	  //        << (*recoendy)[0] << "; mcendx: "<< (*mcendy)[0] <<", recostartx: "<< (*recostarty)[0]
	  //        << ", (reco - true) startx: "<<  (*recostarty)[0] - (*mcstarty)[0] << "; checking minm: "<< diffdr << std::endl;
	  //    std::cout <<"-----------------------------> outside fiducial volume ........ mcstartx: " << (*mcstartx)[0] <<" ; recoendx: "
	  //  << (*recoendx)[0] << "; mcendx: "<< (*mcendx)[0] <<", recostartx: "<< (*recostartx)[0]
	  //              << ", (reco - true) startx: "<<  (*recostartx)[0] - (*mcstartx)[0] << "; checking minm: "<< diffdr << std::endl;
	  }
	}

	//    std::cout << "***************************************************************> everything work upto here before track loop" << std::endl;
	//    if (!fmcstartx.empty()) //&& fmcendx[0] > 0)    
	//      std::cout <<"........................ mcstartx: " << fmcstartx[0] <<" , mcgenerationvx: " << fmcendx[0] <<", recostartx: "<< frecoendx[0] << std::endl;
	
	//std::cout <<"mcstartx: " << fmcstartx[0] <<" ; recoendx: " << frecoendx[0] << "; mcendx: "<< fmcendx[0] <<", recostartx: "<< frecostartx[0] 
	//	  << ", (reco - true) startx: "<<  frecostartx[0] - fmcstartx[0] << "; checking minm: "<< diffx << std::endl;
	
	//if ( frecostartx[0] < -268.49 || frecostartx[0] > -171.94){
	//std::cout <<"-----------------------------> outside fiducial volume ........ mcstartx: " << fmcstartx[0] <<" ; recoendx: " 
	//	    << frecoendx[0] << "; mcendx: "<< fmcendx[0] <<", recostartx: "<< frecostartx[0]
	//	    << ", (reco - true) startx: "<<  frecostartx[0] - fmcstartx[0] << "; checking minm: "<< diffx << std::endl;
	}
	
	}
	}
	}
    */
    //std::cout << x[0] << "\t" << y[0] << "\t" << z[0] << "\t" << 
    //x[last] << "\t" << y[last] << "\t" << z[last] << std::endl;
    
    std::cout << "-- final quick loop to fill summary hists" << std::endl;
    
    for(size_t idx = 0; idx < fGeometry->Nplanes();idx++)
      {
        if (nSimChannelHitVec[idx] > 10)
	  {
            float hitEfficiency = float(nRecobHitVec[idx]) / float(nSimChannelHitVec[idx]);
	    
	    // fNSimChannelHitsVec.push_back(std::min(nSimChannelHitVec[idx],1999));
            //fNRecobHitVec.push_back(std::min(nRecobHitVec[idx],1999));
	    fNSimChannelHitsVec.push_back(nSimChannelHitVec[idx]);
            fNRecobHitVec.push_back(nRecobHitVec[idx]);
            fNFakeHitVec.push_back(nFakeHitVec[idx]/(float)nSimulatedWiresVec[idx]);
            fHitEfficiencyVec.push_back(hitEfficiency);
	  }  
      }
    
    
    //  std::cout << "# wire:  \t "<< wireno << std::endl;
    //MF_LOG_DEBUG("ThroughgoingmuonAnalyzer")
    //<< "Event" << fEvent
    //<< "Run: " << fRun 
    //<< "SubRun"<< fSubRun;
    
    std::cout << "-- fill the tree" << std::endl;
    
    fTree->Fill();  
    
    std::cout << "-- done " << std::endl;
    
    return;
}
DEFINE_ART_MODULE(ThroughgoingmuonAnalyzer)
} // end of namespace
