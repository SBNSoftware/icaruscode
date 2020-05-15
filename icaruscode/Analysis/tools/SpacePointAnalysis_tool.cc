
#include "icaruscode/Analysis/tools/IHitEfficiencyHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larsim/Simulation/LArVoxelID.h"

// Eigen
#include <Eigen/Dense>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"

#include <cmath>
#include <algorithm>

namespace SpacePointAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       SpacePointAnalysis
    // Module Type: producer
    // File:        SpacePointAnalysis.cc
    //
    //              The intent of this module is to provide methods for
    //              "analyzing" hits on waveforms
    //
    // Configuration parameters:
    //
    // TruncMeanFraction     - the fraction of waveform bins to discard when
    //
    // Created by Tracy Usher (usher@slac.stanford.edu) on February 19, 2016
    //
    ////////////////////////////////////////////////////////////////////////
    
// The following typedefs will, obviously, be useful
using HitPtrVec       = std::vector<art::Ptr<recob::Hit>>;
using ViewHitMap      = std::map<size_t,HitPtrVec>;
using TrackViewHitMap = std::map<int,ViewHitMap>;

class SpacePointAnalysis : virtual public IHitEfficiencyHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit SpacePointAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~SpacePointAnalysis();
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset) override;

    /**
     *  @brief Interface for initializing the histograms to be filled
     *
     *  @param TFileService   handle to the TFile service
     *  @param string         subdirectory to store the hists in
     */
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&) override;

    /**
     *  @brief Interface for initializing the tuple variables
     *
     *  @param TTree          pointer to a TTree object to which to add variables
     */
    void initializeTuple(TTree*) override;

    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    void endJob(int numEvents) override;
    
    /**
     *  @brief Interface for filling histograms
     */
    void fillHistograms(const art::Event&)  const override;
    
private:
    
    // Clear mutable variables
    void clear() const;
    
    // Fcl parameters.
    std::vector<art::InputTag>  fRecobHitLabelVec;
    std::vector<art::InputTag>  fSpacePointLabelVec;
    art::InputTag               fMCParticleProducerLabel;
    art::InputTag               fSimChannelProducerLabel;
    art::InputTag               fSimEnergyProducerLabel;
    art::InputTag               fBadChannelProducerLabel;
    std::vector<int>            fOffsetVec;              ///< Allow offsets for each plane
    float                       fSimChannelMinEnergy;
    float                       fSimEnergyMinEnergy;
    
    // TTree variables
    mutable TTree*             fTree;
    
    mutable std::vector<int>   fTPCVec;
    mutable std::vector<int>   fCryoVec;
    mutable std::vector<int>   fPlaneVec;

    // Output tuples for all SpacePoints
    mutable std::vector<int>   fNumIDEsHit0Vec;
    mutable std::vector<int>   fNumIDEsHit1Vec;
    mutable std::vector<int>   fNumIDEsHit2Vec;
    mutable std::vector<int>   fNumIDEsSpacePointVec;
    
    mutable std::vector<float> fSPQualityVec;
    mutable std::vector<float> fSPTotalChargeVec;
    mutable std::vector<float> fSPAsymmetryVec;
    
    mutable std::vector<float> fSmallestPHVec;
    mutable std::vector<float> fAveragePHVec;
    mutable std::vector<float> fLargestDelTVec;

    // Output looking at matching of hits to simchannel information
    mutable std::vector<int>   fTicksSimChannel0Vec;
    mutable std::vector<int>   fTicksSimChanMost0Vec;
    mutable std::vector<int>   fTicksTotHit0Vec;
    mutable std::vector<int>   fTicksMaxSimRel0Vec;
    mutable std::vector<int>   fTicksDiffSimHit0Vec;
    mutable std::vector<float> fEneTotDepHit0Vec;
    mutable std::vector<float> fEneBestDepHit0Vec;
    mutable std::vector<float> fEneMaxDepHit0Vec;
    mutable std::vector<int>   fNDFHit0Vec;
    mutable std::vector<int>   fMultiplicityHit0Vec;
    mutable std::vector<int>   fLocalIndexHit0Vec;
    mutable std::vector<float> fChiSquareHit0Vec;
    mutable std::vector<float> fChargeHit0Vec;
    
    mutable std::vector<int>   fTicksSimChannel1Vec;
    mutable std::vector<int>   fTicksSimChanMost1Vec;
    mutable std::vector<int>   fTicksTotHit1Vec;
    mutable std::vector<int>   fTicksMaxSimRel1Vec;
    mutable std::vector<int>   fTicksDiffSimHit1Vec;
    mutable std::vector<float> fEneTotDepHit1Vec;
    mutable std::vector<float> fEneBestDepHit1Vec;
    mutable std::vector<float> fEneMaxDepHit1Vec;
    mutable std::vector<int>   fNDFHit1Vec;
    mutable std::vector<int>   fMultiplicityHit1Vec;
    mutable std::vector<int>   fLocalIndexHit1Vec;
    mutable std::vector<float> fChiSquareHit1Vec;
    mutable std::vector<float> fChargeHit1Vec;
    
    mutable std::vector<int>   fTicksSimChannel2Vec;
    mutable std::vector<int>   fTicksSimChanMost2Vec;
    mutable std::vector<int>   fTicksTotHit2Vec;
    mutable std::vector<int>   fTicksMaxSimRel2Vec;
    mutable std::vector<int>   fTicksDiffSimHit2Vec;
    mutable std::vector<float> fEneTotDepHit2Vec;
    mutable std::vector<float> fEneBestDepHit2Vec;
    mutable std::vector<float> fEneMaxDepHit2Vec;
    mutable std::vector<int>   fNDFHit2Vec;
    mutable std::vector<int>   fMultiplicityHit2Vec;
    mutable std::vector<int>   fLocalIndexHit2Vec;
    mutable std::vector<float> fChiSquareHit2Vec;
    mutable std::vector<float> fChargeHit2Vec;
    
    // Output tuples for SpacePoints with one or more hits not matching MC
    mutable std::vector<int>   fNumIDEsHit0NoMVec;
    mutable std::vector<int>   fNumIDEsHit1NoMVec;
    mutable std::vector<int>   fNumIDEsHit2NoMVec;
    
    mutable std::vector<float> fSPQualityNoMVec;
    mutable std::vector<float> fSPTotalChargeNoMVec;
    mutable std::vector<float> fSPAsymmetryNoMVec;
    
    mutable std::vector<float> fSmallestPHNoMVec;
    mutable std::vector<float> fAveragePHNoMVec;
    mutable std::vector<float> fLargestDelTNoMVec;

    // Output tuples for Ghost SpacePoints
    mutable std::vector<int>   fNumIDEsHit0GhostVec;
    mutable std::vector<int>   fNumIDEsHit1GhostVec;
    mutable std::vector<int>   fNumIDEsHit2GhostVec;
    
    mutable std::vector<float> fSPQualityGhostVec;
    mutable std::vector<float> fSPTotalChargeGhostVec;
    mutable std::vector<float> fSPAsymmetryGhostVec;
    
    mutable std::vector<float> fSmallestPHGhostVec;
    mutable std::vector<float> fAveragePHGhostVec;
    mutable std::vector<float> fLargestDelTGhostVec;

    // Output tuples for matched space points
    mutable std::vector<int>   fNumIDEsHit0MatchVec;
    mutable std::vector<int>   fNumIDEsHit1MatchVec;
    mutable std::vector<int>   fNumIDEsHit2MatchVec;
    mutable std::vector<int>   fNumIDEsSpacePointMatchVec;

    mutable std::vector<float> fSPQualityMatchVec;
    mutable std::vector<float> fSPTotalChargeMatchVec;
    mutable std::vector<float> fSPAsymmetryMatchVec;
    
    mutable std::vector<float> fSmallestPHMatchVec;
    mutable std::vector<float> fAveragePHMatchVec;
    mutable std::vector<float> fLargestDelTMatchVec;

    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
    const detinfo::DetectorClocks*     fClockService;         ///< Detector clocks service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
SpacePointAnalysis::SpacePointAnalysis(fhicl::ParameterSet const & pset) : fTree(nullptr)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fClockService       = lar::providerFrom<detinfo::DetectorClocksService>();
    
    configure(pset);
    
    // Report.
    mf::LogInfo("SpacePointAnalysis") << "SpacePointAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
SpacePointAnalysis::~SpacePointAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void SpacePointAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fRecobHitLabelVec         = pset.get< std::vector<art::InputTag>>("HitLabelVec",         std::vector<art::InputTag>() = {"cluster3d"});
    fSpacePointLabelVec       = pset.get< std::vector<art::InputTag>>("SpacePointLabelVec",  std::vector<art::InputTag>() = {"cluster3d"});
    fMCParticleProducerLabel  = pset.get< art::InputTag             >("MCParticleLabel",     "largeant");
    fSimChannelProducerLabel  = pset.get< art::InputTag             >("SimChannelLabel",     "largeant");
    fSimEnergyProducerLabel   = pset.get< art::InputTag             >("SimEnergyLabel",      "largeant");
    fOffsetVec                = pset.get<std::vector<int>           >("OffsetVec",           std::vector<int>()={0,0,0});
    fSimChannelMinEnergy      = pset.get<float                      >("SimChannelMinEnergy", std::numeric_limits<float>::epsilon());
    fSimEnergyMinEnergy       = pset.get<float                      >("SimEnergyMinEnergy",  std::numeric_limits<float>::epsilon());

    return;
}

//----------------------------------------------------------------------------
/// Begin job method.
void SpacePointAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
//    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());

    return;
}
    
void SpacePointAnalysis::initializeTuple(TTree* tree)
{
    fTree = tree;

    fTree->Branch("CryostataVec",       "std::vector<int>",   &fCryoVec);
    fTree->Branch("TPCVec",             "std::vector<int>",   &fTPCVec);
    fTree->Branch("PlaneVec",           "std::vector<int>",   &fPlaneVec);
 
    fTree->Branch("NumIDEsHit0",        "std::vector<int>",   &fNumIDEsHit0Vec);
    fTree->Branch("NumIDEsHit1",        "std::vector<int>",   &fNumIDEsHit1Vec);
    fTree->Branch("NumIDEsHit2",        "std::vector<int>",   &fNumIDEsHit2Vec);
    fTree->Branch("NumIDEsSpacePoint",  "std::vector<int>",   &fNumIDEsSpacePointVec);
    
    fTree->Branch("SPQuality",          "std::vector<float>", &fSPQualityVec);
    fTree->Branch("SPTotalCharge",      "std::vector<float>", &fSPTotalChargeVec);
    fTree->Branch("SPAsymmetry",        "std::vector<float>", &fSPAsymmetryVec);
    fTree->Branch("SmallestPH",         "std::vector<float>", &fSmallestPHVec);
    fTree->Branch("AveragePH",          "std::vector<float>", &fAveragePHVec);
    fTree->Branch("LargestDelT",        "std::vector<float>", &fLargestDelTVec);

    fTree->Branch("TicksSimChannel0",   "std::vector<int>",   &fTicksSimChannel0Vec);
    fTree->Branch("TicksSimChanMost0",  "std::vector<int>",   &fTicksSimChanMost0Vec);
    fTree->Branch("TicksTotHit0",       "std::vector<int>",   &fTicksTotHit0Vec);
    fTree->Branch("TicksMaxSimRel0",    "std::vector<int>",   &fTicksMaxSimRel0Vec);
    fTree->Branch("TicksDiffSimHit0",   "std::vector<int>",   &fTicksDiffSimHit0Vec);
    fTree->Branch("EneTotDepHit0",      "std::vector<float>", &fEneTotDepHit0Vec);
    fTree->Branch("EneBestDepHit0",     "std::vector<float>", &fEneBestDepHit0Vec);
    fTree->Branch("EneMaxDepHit0",      "std::vector<float>", &fEneMaxDepHit0Vec);
    fTree->Branch("NDFHit0",            "std::vector<int>",   &fNDFHit0Vec);
    fTree->Branch("MultiplicityHit0",   "std::vector<int>",   &fMultiplicityHit0Vec);
    fTree->Branch("LocalIndexHit0",     "std::vector<int>",   &fLocalIndexHit0Vec);
    fTree->Branch("ChiSquareHit0",      "std::vector<float>", &fChiSquareHit0Vec);
    fTree->Branch("ChargeHit0",         "std::vector<float>", &fChargeHit0Vec);

    fTree->Branch("TicksSimChannel1",   "std::vector<int>",   &fTicksSimChannel1Vec);
    fTree->Branch("TicksSimChanMost1",  "std::vector<int>",   &fTicksSimChanMost1Vec);
    fTree->Branch("TicksTotHit1",       "std::vector<int>",   &fTicksTotHit1Vec);
    fTree->Branch("TicksMaxSimRel1",    "std::vector<int>",   &fTicksMaxSimRel1Vec);
    fTree->Branch("TicksDiffSimHit1",   "std::vector<int>",   &fTicksDiffSimHit1Vec);
    fTree->Branch("EneTotDepHit1",      "std::vector<float>", &fEneTotDepHit1Vec);
    fTree->Branch("EneBestDepHit1",     "std::vector<float>", &fEneBestDepHit1Vec);
    fTree->Branch("EneMaxDepHit1",      "std::vector<float>", &fEneMaxDepHit1Vec);
    fTree->Branch("NDFHit1",            "std::vector<int>",   &fNDFHit1Vec);
    fTree->Branch("MultiplicityHit1",   "std::vector<int>",   &fMultiplicityHit1Vec);
    fTree->Branch("LocalIndexHit1",     "std::vector<int>",   &fLocalIndexHit1Vec);
    fTree->Branch("ChiSquareHit1",      "std::vector<float>", &fChiSquareHit1Vec);
    fTree->Branch("ChargeHit1",         "std::vector<float>", &fChargeHit1Vec);

    fTree->Branch("TicksSimChannel2",   "std::vector<int>",   &fTicksSimChannel2Vec);
    fTree->Branch("TicksSimChanMost2",  "std::vector<int>",   &fTicksSimChanMost2Vec);
    fTree->Branch("TicksTotHit2",       "std::vector<int>",   &fTicksTotHit2Vec);
    fTree->Branch("TicksMaxSimRel2",    "std::vector<int>",   &fTicksMaxSimRel2Vec);
    fTree->Branch("TicksDiffSimHit2",   "std::vector<int>",   &fTicksDiffSimHit2Vec);
    fTree->Branch("EneTotDepHit2",      "std::vector<float>", &fEneTotDepHit2Vec);
    fTree->Branch("EneBestDepHit2",     "std::vector<float>", &fEneBestDepHit2Vec);
    fTree->Branch("EneMaxDepHit2",      "std::vector<float>", &fEneMaxDepHit2Vec);    
    fTree->Branch("NDFHit2",            "std::vector<int>",   &fNDFHit2Vec);
    fTree->Branch("MultiplicityHit2",   "std::vector<int>",   &fMultiplicityHit2Vec);
    fTree->Branch("LocalIndexHit2",     "std::vector<int>",   &fLocalIndexHit2Vec);
    fTree->Branch("ChiSquareHit2",      "std::vector<float>", &fChiSquareHit2Vec);
    fTree->Branch("ChargeHit2",         "std::vector<float>", &fChargeHit2Vec);

    fTree->Branch("NumIDEsHit0NoM",     "std::vector<int>",   &fNumIDEsHit0NoMVec);
    fTree->Branch("NumIDEsHit1NoM",     "std::vector<int>",   &fNumIDEsHit1NoMVec);
    fTree->Branch("NumIDEsHit2NoM",     "std::vector<int>",   &fNumIDEsHit2NoMVec);

    fTree->Branch("SPQualityNoM",       "std::vector<float>", &fSPQualityNoMVec);
    fTree->Branch("SPTotalChargeNoM",   "std::vector<float>", &fSPTotalChargeNoMVec);
    fTree->Branch("SPAsymmetryNoM",     "std::vector<float>", &fSPAsymmetryNoMVec);
    fTree->Branch("SmallestPHNoM",      "std::vector<float>", &fSmallestPHNoMVec);
    fTree->Branch("AveragePHNoM",       "std::vector<float>", &fAveragePHNoMVec);
    fTree->Branch("LargestDelTNoM",     "std::vector<float>", &fLargestDelTNoMVec);

    fTree->Branch("NumIDEsHit0Match",   "std::vector<int>",   &fNumIDEsHit0MatchVec);
    fTree->Branch("NumIDEsHit1Match",   "std::vector<int>",   &fNumIDEsHit1MatchVec);
    fTree->Branch("NumIDEsHit2Match",   "std::vector<int>",   &fNumIDEsHit2MatchVec);
    fTree->Branch("NumIDEsSPMatch",     "std::vector<int>",   &fNumIDEsSpacePointMatchVec);
    fTree->Branch("SmallestPHMatch",    "std::vector<float>", &fSmallestPHMatchVec);
    fTree->Branch("AveragePHMatch",     "std::vector<float>", &fAveragePHMatchVec);
    fTree->Branch("LargestDelTMatch",   "std::vector<float>", &fLargestDelTMatchVec);

    fTree->Branch("SPQualityMatch",     "std::vector<float>", &fSPQualityMatchVec);
    fTree->Branch("SPTotalChargeMatch", "std::vector<float>", &fSPTotalChargeMatchVec);
    fTree->Branch("SPAsymmetryMatch",   "std::vector<float>", &fSPAsymmetryMatchVec);
    
    fTree->Branch("NumIDEsHit0Ghost",   "std::vector<int>",   &fNumIDEsHit0GhostVec);
    fTree->Branch("NumIDEsHit1Ghost",   "std::vector<int>",   &fNumIDEsHit1GhostVec);
    fTree->Branch("NumIDEsHit2Ghost",   "std::vector<int>",   &fNumIDEsHit2GhostVec);

    fTree->Branch("SPQualityGhost",     "std::vector<float>", &fSPQualityGhostVec);
    fTree->Branch("SPTotalChargeGhost", "std::vector<float>", &fSPTotalChargeGhostVec);
    fTree->Branch("SPAsymmetryGhost",   "std::vector<float>", &fSPAsymmetryGhostVec);
    fTree->Branch("SmallestPHGhost",    "std::vector<float>", &fSmallestPHGhostVec);
    fTree->Branch("AveragePHGhost",     "std::vector<float>", &fAveragePHGhostVec);
    fTree->Branch("LargestDelTGhost",   "std::vector<float>", &fLargestDelTGhostVec);

    clear();

    return;
}
    
void SpacePointAnalysis::clear() const
{
    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();

    fTicksSimChannel0Vec.clear();
    fTicksSimChanMost0Vec.clear();
    fTicksTotHit0Vec.clear();
    fTicksMaxSimRel0Vec.clear();
    fTicksDiffSimHit0Vec.clear();
    fEneTotDepHit0Vec.clear();
    fEneBestDepHit0Vec.clear();
    fEneMaxDepHit0Vec.clear();
    fNDFHit0Vec.clear();
    fMultiplicityHit0Vec.clear();
    fLocalIndexHit0Vec.clear();
    fChiSquareHit0Vec.clear();
    fChargeHit0Vec.clear();

    fTicksSimChannel1Vec.clear();
    fTicksSimChanMost1Vec.clear();
    fTicksTotHit1Vec.clear();
    fTicksMaxSimRel1Vec.clear();
    fTicksDiffSimHit1Vec.clear();
    fEneTotDepHit1Vec.clear();
    fEneBestDepHit1Vec.clear();
    fEneMaxDepHit1Vec.clear();
    fNDFHit1Vec.clear();
    fMultiplicityHit1Vec.clear();
    fLocalIndexHit1Vec.clear();
    fChiSquareHit1Vec.clear();
    fChargeHit1Vec.clear();

    fTicksSimChannel2Vec.clear();
    fTicksSimChanMost2Vec.clear();
    fTicksTotHit2Vec.clear();
    fTicksMaxSimRel2Vec.clear();
    fTicksDiffSimHit2Vec.clear();
    fEneTotDepHit2Vec.clear();
    fEneBestDepHit2Vec.clear();
    fEneMaxDepHit2Vec.clear();
    fNDFHit2Vec.clear();
    fMultiplicityHit2Vec.clear();
    fLocalIndexHit2Vec.clear();
    fChiSquareHit2Vec.clear();
    fChargeHit2Vec.clear();

    fNumIDEsHit0Vec.clear();
    fNumIDEsHit1Vec.clear();
    fNumIDEsHit2Vec.clear();
    fNumIDEsSpacePointVec.clear();
    
    fSPQualityVec.clear();
    fSPTotalChargeVec.clear();
    fSPAsymmetryVec.clear();
    fSmallestPHVec.clear();
    fAveragePHVec.clear();
    fLargestDelTVec.clear();

    fNumIDEsHit0NoMVec.clear();
    fNumIDEsHit1NoMVec.clear();
    fNumIDEsHit2NoMVec.clear();

    fSPQualityNoMVec.clear();
    fSPTotalChargeNoMVec.clear();
    fSPAsymmetryNoMVec.clear();
    fSmallestPHNoMVec.clear();
    fAveragePHNoMVec.clear();
    fLargestDelTNoMVec.clear();

    fNumIDEsHit0GhostVec.clear();
    fNumIDEsHit1GhostVec.clear();
    fNumIDEsHit2GhostVec.clear();

    fSPQualityGhostVec.clear();
    fSPTotalChargeGhostVec.clear();
    fSPAsymmetryGhostVec.clear();
    fSmallestPHGhostVec.clear();
    fAveragePHGhostVec.clear();
    fLargestDelTGhostVec.clear();

    fNumIDEsHit0MatchVec.clear();
    fNumIDEsHit1MatchVec.clear();
    fNumIDEsHit2MatchVec.clear();
    fNumIDEsSpacePointMatchVec.clear();

    fSPQualityMatchVec.clear();
    fSPTotalChargeMatchVec.clear();
    fSPAsymmetryMatchVec.clear();
    fSmallestPHMatchVec.clear();
    fAveragePHMatchVec.clear();
    fLargestDelTMatchVec.clear();

    return;
}
   
// Create a struct allowing us to sort IDEs in a set by largest to smallest energy
struct ideCompare
{
    bool operator() (const sim::IDE* left, const sim::IDE* right) const {return left->energy > right->energy;}
};
     
void SpacePointAnalysis::fillHistograms(const art::Event& event) const
{
    // Ok... this is starting to grow too much and get out of control... we will need to break it up directly...
    
    // Always clear the tuple
    clear();
    
    art::Handle<std::vector<sim::SimChannel>> simChannelHandle;
    event.getByLabel(fSimChannelProducerLabel, simChannelHandle);
    
    if (!simChannelHandle.isValid() || simChannelHandle->empty() ) return;
    
    art::Handle<std::vector<sim::SimEnergyDeposit>> simEnergyHandle;
    event.getByLabel(fSimEnergyProducerLabel, simEnergyHandle);
    
    if (!simEnergyHandle.isValid() || simEnergyHandle->empty()) return;
    
    art::Handle<std::vector<simb::MCParticle>> mcParticleHandle;
    event.getByLabel(fMCParticleProducerLabel, mcParticleHandle);

    // If there is no sim channel informaton then exit
    if (!mcParticleHandle.isValid()) return;
    
    mf::LogDebug("SpacePointAnalysis") << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    // First task is to build a map between ides and voxel ids (that we calcualate based on position)
    // and also get the reverse since it will be useful in the end.
    // At the same time should also build a mapping of ides per channel so we can do quick hit lookup
    using SimIDESet          = std::set<const sim::IDE*,ideCompare>;
    using IDEToVoxelIDMap    = std::unordered_map<const sim::IDE*, sim::LArVoxelID>;
    using VoxelIDToIDESetMap = std::map<sim::LArVoxelID, SimIDESet>;
    using TDCToIDEMap        = std::map<unsigned short, SimIDESet>; // We need this one in order
    using ChanToTDCToIDEMap  = std::map<raw::ChannelID_t, TDCToIDEMap>;
    using VoxelIDSet         = std::set<sim::LArVoxelID>;

    IDEToVoxelIDMap    ideToVoxelIDMap;
    VoxelIDToIDESetMap voxelIDToIDEMap;
    ChanToTDCToIDEMap  chanToTDCToIDEMap;
    VoxelIDSet         simChannelVoxelIDSet;

    // Fill the above maps/structures
    for(const auto& simChannel : *simChannelHandle)
    {
        for(const auto& tdcide : simChannel.TDCIDEMap())
        {
            for(const auto& ide : tdcide.second) //chanToTDCToIDEMap[simChannel.Channel()][tdcide.first] = ide;
            {
                if (ide.energy < fSimChannelMinEnergy) continue;
                
                sim::LArVoxelID voxelID(ide.x,ide.y,ide.z,0.);
                
                ideToVoxelIDMap[&ide]    = voxelID;
                voxelIDToIDEMap[voxelID].insert(&ide);
                chanToTDCToIDEMap[simChannel.Channel()][tdcide.first].insert(&ide);
                simChannelVoxelIDSet.insert(voxelID);
                
                if (ide.energy < std::numeric_limits<float>::epsilon()) mf::LogDebug("SpacePointAnalysis") << ">> epsilon simchan deposited energy: " << ide.energy << std::endl;
            }
        }
    }

    // Now we go throught the SimEnergyDeposit objects and try to make similar mappings
    // It is worth noting that in this case there can be multiple SimEnergyDeposit objects per voxel
    // We assume that calculating the voxel ID as above from the mean position of the SimEnergyDeposit objects will
    // result in the correct Voxel ID for relating to SimChannels (to be demonstrated!)
    using SimEnergyToVoxelIDMap    = std::unordered_map<const sim::SimEnergyDeposit*, sim::LArVoxelID>;
    using VoxelIDToSimEnergySetMap = std::map<sim::LArVoxelID, std::set<const sim::SimEnergyDeposit*>>;
    
    SimEnergyToVoxelIDMap    simEnergyToVoxelIDMap;
    VoxelIDToSimEnergySetMap voxelIDToSimEnergySetMap;
    VoxelIDSet               simEnergyVoxelIDSet;

    for(const auto& simEnergy : *simEnergyHandle)
    {
        if (simEnergy.Energy() < fSimEnergyMinEnergy) continue;
        
        geo::Point_t    midPoint = simEnergy.MidPoint();
        sim::LArVoxelID voxelID(midPoint.X(),midPoint.Y(),midPoint.Z(),0.);

        simEnergyToVoxelIDMap[&simEnergy] = voxelID;
        voxelIDToSimEnergySetMap[voxelID].insert(&simEnergy);
        simEnergyVoxelIDSet.insert(voxelID);
        
        if (simEnergy.Energy() < std::numeric_limits<float>::epsilon()) mf::LogDebug("SpacePointAnalysis") << ">> epsilon simenergy deposited energy: " << simEnergy.Energy() << std::endl;
    }

    // Ok, for my next trick I want to build a mapping between hits and voxel IDs. Note that any given hit can be associated to more than one voxel...
    // We do this on the entire hit collection, ultimately we will want to consider SpacePoint efficiency (this could be done in the loop over SpacePoints
    // using the associated hits and would save time/memory)
    using VoxelIDSet           = std::set<sim::LArVoxelID>;
//    using VoxelIDSetVec        = std::vector<VoxelIDSet>;
//    using RecobHitToVoxelIDMap = std::unordered_map<const recob::Hit*, VoxelIDSetVec>;
    using RecobHitToVoxelIDMap = std::unordered_map<const recob::Hit*, VoxelIDSet>;
    
    RecobHitToVoxelIDMap recobHitToVoxelIDMap;
    
    // And now fill it
    for(const auto& hitLabel : fRecobHitLabelVec)
    {
        art::Handle< std::vector<recob::Hit> > hitHandle;
        event.getByLabel(hitLabel, hitHandle);
        
        // Outer loop over hits in this hit collection
        for(const auto& hit : *hitHandle)
        {
            // ********** TEMPORARY CUT ***********
            // Try ignoring hit which are not "good" fits
            if (hit.GoodnessOfFit() < 0.) continue;

            // Recover channel information based on this hit
            ChanToTDCToIDEMap::const_iterator chanToTDCToIDEItr = chanToTDCToIDEMap.find(hit.Channel());

            // This at least weeds out the case where a channel displayed a hit but had no signal otherwise
            if (chanToTDCToIDEItr != chanToTDCToIDEMap.end())
            {
                // Recover hit time range (in ticks), cast a wide net here
                int peakTick  = std::round(hit.PeakTime());
                int startTick = std::floor(hit.PeakTime() - 3. * hit.RMS());
                int endTick   = std::ceil(hit.PeakTime() + 3. * hit.RMS());

                // Initial data structures
                using TickToIDEMap        = std::map<unsigned short, const sim::IDE*>;
                using TrackIDToTickIDEMap = std::map<int,TickToIDEMap>;

                TrackIDToTickIDEMap trackIDToTickIDEMap;
                
                const TDCToIDEMap& tdcToIDEMap = chanToTDCToIDEItr->second;

                // The idea here is to build up the list of all ide's that are potentially assiociated to this hit
                // and to keep track of them by deposited energy
                for(unsigned short tick =startTick; tick <= endTick; tick++)
                {
                    unsigned short hitTDC = fClockService->TPCTick2TDC(tick - fOffsetVec[hit.WireID().Plane]);
                    
                    TDCToIDEMap::const_iterator ideIterator = tdcToIDEMap.find(hitTDC);
                    
                    if (ideIterator != tdcToIDEMap.end())
                        for (const auto& ide : ideIterator->second) trackIDToTickIDEMap[ide->trackID][tick] = ide;
                }

                // Is the map filled? The real test for noise hits
                if (trackIDToTickIDEMap.empty()) continue;

                // Try to determine "the" hit/MC match
                // This is going to select out the track ID that deposits the most energy in the range of the hit
                int   bestTrackID(0);
                int   bestTicks(0);
                float bestDepEne(-1.);
                float totDepEne(0.);

                for(const auto& trackInfo : trackIDToTickIDEMap)
                {
                    float depEne(0.);

                    for(const auto& tickInfo : trackInfo.second)
                    {
                        depEne += tickInfo.second->energy;
                    }

                    totDepEne += depEne;

                    if (depEne > bestDepEne)
                    {
                        bestTrackID = trackInfo.first;
                        bestTicks   = trackInfo.second.size();
                        bestDepEne  = depEne;
                    }
                }

                // Get the selected TickToIDEMap
                const TickToIDEMap& tickToIDEMap = trackIDToTickIDEMap[bestTrackID];

                // Find the max dep ene tick
                unsigned short maxDepTick(0);
                float          maxDepEneTick(-1.);

                for(const auto& tickInfo : tickToIDEMap)
                {
                    if (tickInfo.second->energy > maxDepEneTick)
                    {
                        maxDepTick    = tickInfo.first;
                        maxDepEneTick = tickInfo.second->energy;
                    }
                }

                // One final time through to find sim ticks that "matter"
                // We define this as the collection of IDE's that make up to 90% of the total deposited energy
                std::vector<std::pair<unsigned short,const sim::IDE*>> tickIDEVec;

                for(const auto& tickInfo : tickToIDEMap) tickIDEVec.emplace_back(tickInfo);

                std::sort(tickIDEVec.begin(),tickIDEVec.end(),[](const auto& left,const auto& right){return left.second->energy > right.second->energy;});

                // Grab the voxelID set for this tick
                VoxelIDSet voxelIDSet;

                int   bestTicksGood(0);
                float sumEne(0.);

                for(const auto& tickInfo : tickIDEVec)
                {
                    // At the same time we can keep track of the voxels associated to the best track
                    const sim::LArVoxelID& voxelID = ideToVoxelIDMap[tickInfo.second];

                    sumEne += tickInfo.second->energy;
                    bestTicksGood++;

                    voxelIDSet.insert(voxelID);

                    if (sumEne > 0.9 * bestDepEne) break;
                }
            
                // Finally, grab the voxels from the track leaving the most energy
//                recobHitToVoxelIDMap[&hit].emplace_back(voxelIDSet);
                recobHitToVoxelIDMap[&hit] = voxelIDSet;

                // Fill depending on the plane
                if (hit.WireID().Plane == 0)
                {
                    fTicksSimChannel0Vec.emplace_back(bestTicks);
                    fTicksSimChanMost0Vec.emplace_back(bestTicksGood);
                    fTicksTotHit0Vec.emplace_back(endTick-startTick+1);
                    fTicksMaxSimRel0Vec.emplace_back(maxDepTick-startTick);
                    fTicksDiffSimHit0Vec.emplace_back(peakTick-startTick-bestTicks);
                    fEneTotDepHit0Vec.emplace_back(totDepEne);
                    fEneBestDepHit0Vec.emplace_back(bestDepEne);
                    fEneMaxDepHit0Vec.emplace_back(maxDepEneTick);
                    fNDFHit0Vec.emplace_back(hit.DegreesOfFreedom());
                    fMultiplicityHit0Vec.emplace_back(hit.Multiplicity());
                    fLocalIndexHit0Vec.emplace_back(hit.LocalIndex());
                    fChiSquareHit0Vec.emplace_back(hit.GoodnessOfFit());
                    fChargeHit0Vec.emplace_back(hit.SummedADC());
                }
                else if (hit.WireID().Plane == 1)
                {
                    fTicksSimChannel1Vec.emplace_back(bestTicks);
                    fTicksSimChanMost1Vec.emplace_back(bestTicksGood);
                    fTicksTotHit1Vec.emplace_back(endTick-startTick+1);
                    fTicksMaxSimRel1Vec.emplace_back(maxDepTick-startTick);
                    fTicksDiffSimHit1Vec.emplace_back(peakTick-startTick-bestTicks);
                    fEneTotDepHit1Vec.emplace_back(totDepEne);
                    fEneBestDepHit1Vec.emplace_back(bestDepEne);
                    fEneMaxDepHit1Vec.emplace_back(maxDepEneTick);
                    fNDFHit1Vec.emplace_back(hit.DegreesOfFreedom());
                    fMultiplicityHit1Vec.emplace_back(hit.Multiplicity());
                    fLocalIndexHit1Vec.emplace_back(hit.LocalIndex());
                    fChiSquareHit1Vec.emplace_back(hit.GoodnessOfFit());
                    fChargeHit1Vec.emplace_back(hit.SummedADC());
                }
                else
                {
                    fTicksSimChannel2Vec.emplace_back(bestTicks);
                    fTicksSimChanMost2Vec.emplace_back(bestTicksGood);
                    fTicksTotHit2Vec.emplace_back(endTick-startTick+1);
                    fTicksMaxSimRel2Vec.emplace_back(maxDepTick-startTick);
                    fTicksDiffSimHit2Vec.emplace_back(peakTick-startTick-bestTicks);
                    fEneTotDepHit2Vec.emplace_back(totDepEne);
                    fEneBestDepHit2Vec.emplace_back(bestDepEne);
                    fEneMaxDepHit2Vec.emplace_back(maxDepEneTick);
                    fNDFHit2Vec.emplace_back(hit.DegreesOfFreedom());
                    fMultiplicityHit2Vec.emplace_back(hit.Multiplicity());
                    fLocalIndexHit2Vec.emplace_back(hit.LocalIndex());
                    fChiSquareHit2Vec.emplace_back(hit.GoodnessOfFit());
                    fChargeHit2Vec.emplace_back(hit.SummedADC());
                }
            }
        }
    }

    // Armed with these maps we can now process the SpacePoints...
    if (!recobHitToVoxelIDMap.empty())
    {
        // So now we loop through the various SpacePoint sources
        for(const auto& spacePointLabel : fSpacePointLabelVec)
        {
            art::Handle< std::vector<recob::SpacePoint> > spacePointHandle;
            event.getByLabel(spacePointLabel, spacePointHandle);
            
            if (!spacePointHandle.isValid()) continue;

            // Look up assocations to hits
            art::FindManyP<recob::Hit> spHitAssnVec(spacePointHandle, event, spacePointLabel);

            // And now, without further adieu, here we begin the loop that will actually produce some useful output.
            // Loop over all space points and find out their true nature
            for(size_t idx = 0; idx < spacePointHandle->size(); idx++)
            {
                art::Ptr<recob::SpacePoint> spacePointPtr(spacePointHandle,idx);
                
                std::vector<art::Ptr<recob::Hit>> associatedHits(spHitAssnVec.at(spacePointPtr.key()));
                
                if (associatedHits.size() != 3)
                {
                    mf::LogDebug("SpacePointAnalysis") << "I am certain this cannot happen... but here you go, space point with " << associatedHits.size() << " hits" << std::endl;
                    continue;
                }
                
                // Retrieve the magic numbers from the space point
                float spQuality   = spacePointPtr->Chisq();
                float spCharge    = spacePointPtr->ErrXYZ()[1];
                float spAsymmetry = spacePointPtr->ErrXYZ()[3];
                float smallestPH  = std::numeric_limits<float>::max();
                float averagePH   = 0.;
                float averagePT   = 0.;
                float largestDelT = 0.;
                
                std::vector<int> numIDEsHitVec;
                int              numIDEsSpacePoint(0);

                std::vector<RecobHitToVoxelIDMap::const_iterator> recobHitToVoxelIterVec;
                
                // Now we can use our maps to find out if the hits making up the SpacePoint are truly related...
                for(const auto& hitPtr : associatedHits)
                {
                    RecobHitToVoxelIDMap::iterator hitToVoxelItr = recobHitToVoxelIDMap.find(hitPtr.get());
                    
                    float  peakAmplitude = hitPtr->PeakAmplitude();
                    
                    averagePH += peakAmplitude;
                    averagePT += hitPtr->PeakTime();
                    
                    if (peakAmplitude < smallestPH) smallestPH = peakAmplitude;
                    
                    if (hitToVoxelItr == recobHitToVoxelIDMap.end())
                    {
                        numIDEsHitVec.push_back(0);
                        continue;
                    }
                    
                    recobHitToVoxelIterVec.push_back(hitToVoxelItr);
                    numIDEsHitVec.push_back(hitToVoxelItr->second.size());
                }
                
                averagePH /= 3.;
                averagePT /= 3.;
                
                for(const auto& hitPtr : associatedHits)
                {
                    float delT = hitPtr->PeakTime() - averagePT;
                    
                    if (std::abs(delT) > std::abs(largestDelT)) largestDelT = delT;
                }
                
                // If a SpacePoint is made from "true" MC hits then we will have found the relations to the MC info for all three
                // hits. If this condition is not satisfied it means one or more hits making the SpacePoint are noise hits
                if (recobHitToVoxelIterVec.size() == 3)
                {
    
                    bool ghostHit(true);
                    
                    // Find the intersection of the vectors of IDEs for the first two hits
                    std::vector<sim::LArVoxelID> firstIntersectionVec(recobHitToVoxelIterVec[0]->second.size()+recobHitToVoxelIterVec[1]->second.size());
                    
                    std::vector<sim::LArVoxelID>::iterator firstIntersectionItr = std::set_intersection(recobHitToVoxelIterVec[0]->second.begin(),recobHitToVoxelIterVec[0]->second.end(),
                                                                                                        recobHitToVoxelIterVec[1]->second.begin(),recobHitToVoxelIterVec[1]->second.end(),
                                                                                                        firstIntersectionVec.begin());
                    
                    firstIntersectionVec.resize(firstIntersectionItr - firstIntersectionVec.begin());
                    
                    // No intersection means, of course, the hits did not come from the same MC energy deposit
                    if (!firstIntersectionVec.empty())
                    {
                        // Now find the intersection of the resultant intersection above and the third hit
                        std::vector<sim::LArVoxelID> secondIntersectionVec(firstIntersectionVec.size()+recobHitToVoxelIterVec[2]->second.size());
                        
                        std::vector<sim::LArVoxelID>::iterator secondIntersectionItr = std::set_intersection(firstIntersectionVec.begin(),             firstIntersectionVec.end(),
                                                                                                             recobHitToVoxelIterVec[1]->second.begin(),recobHitToVoxelIterVec[1]->second.end(),
                                                                                                             secondIntersectionVec.begin());
                        
                        secondIntersectionVec.resize(secondIntersectionItr - secondIntersectionVec.begin());
                        
                        // Again, no IDEs in the intersection means it is a ghost space point but, of course, we are hoping
                        // there are common IDEs so we can call it a real SpacePoint
                        if (!secondIntersectionVec.empty())
                        {
                            numIDEsSpacePoint = secondIntersectionVec.size();
                            
                            fNumIDEsHit0MatchVec.push_back(numIDEsHitVec[0]);
                            fNumIDEsHit1MatchVec.push_back(numIDEsHitVec[1]);
                            fNumIDEsHit2MatchVec.push_back(numIDEsHitVec[2]);
                            fNumIDEsSpacePointMatchVec.push_back(numIDEsSpacePoint);

                            fSPQualityMatchVec.push_back(spQuality);
                            fSPTotalChargeMatchVec.push_back(spCharge);
                            fSPAsymmetryMatchVec.push_back(spAsymmetry);
                            fSmallestPHMatchVec.push_back(smallestPH);
                            fAveragePHMatchVec.push_back(averagePH);
                            fLargestDelTMatchVec.push_back(largestDelT);

                            ghostHit = false;
                        }
                    }
                    
                    if (ghostHit)
                    {
                        fNumIDEsHit0GhostVec.push_back(numIDEsHitVec[0]);
                        fNumIDEsHit1GhostVec.push_back(numIDEsHitVec[1]);
                        fNumIDEsHit2GhostVec.push_back(numIDEsHitVec[2]);
                        
                        fSPQualityGhostVec.push_back(spQuality);
                        fSPTotalChargeGhostVec.push_back(spCharge);
                        fSPAsymmetryGhostVec.push_back(spAsymmetry);
                        fSmallestPHGhostVec.push_back(smallestPH);
                        fAveragePHGhostVec.push_back(averagePH);
                        fLargestDelTGhostVec.push_back(largestDelT);
                    }
                }
                else
                {
                    fNumIDEsHit0NoMVec.push_back(numIDEsHitVec[0]);
                    fNumIDEsHit1NoMVec.push_back(numIDEsHitVec[1]);
                    fNumIDEsHit2NoMVec.push_back(numIDEsHitVec[2]);
                    
                    fSPQualityNoMVec.push_back(spQuality);
                    fSPTotalChargeNoMVec.push_back(spCharge);
                    fSPAsymmetryNoMVec.push_back(spAsymmetry);
                    fSmallestPHNoMVec.push_back(smallestPH);
                    fAveragePHNoMVec.push_back(averagePH);
                    fLargestDelTNoMVec.push_back(largestDelT);
                }
                
                // Fill for "all" cases
                fSPQualityVec.push_back(spQuality);
                fSPTotalChargeVec.push_back(spCharge);
                fSPAsymmetryVec.push_back(spAsymmetry);
                fSmallestPHVec.push_back(smallestPH);
                fAveragePHVec.push_back(averagePH);
                fLargestDelTVec.push_back(largestDelT);
                
                fNumIDEsHit0Vec.push_back(numIDEsHitVec[0]);
                fNumIDEsHit1Vec.push_back(numIDEsHitVec[1]);
                fNumIDEsHit2Vec.push_back(numIDEsHitVec[2]);
                fNumIDEsSpacePointVec.push_back(numIDEsSpacePoint);
            }
        }
    }

//            // Store tuple variables
//            fTPCVec.push_back(wids[0].TPC);
//            fCryoVec.push_back(wids[0].Cryostat);
//            fPlaneVec.push_back(wids[0].Plane);
//            fWireVec.push_back(wids[0].Wire);

    return;
}
    
// Useful for normalizing histograms
void SpacePointAnalysis::endJob(int numEvents)
{
    return;
}
    
DEFINE_ART_CLASS_TOOL(SpacePointAnalysis)
}
