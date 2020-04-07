
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
#include <tuple>
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
    mutable TTree*             fSpacePointTree;
    mutable TTree*             fMatchedHitTree;
    
    mutable std::vector<int>   fTPCVec;
    mutable std::vector<int>   fCryoVec;
    mutable std::vector<int>   fPlaneVec;

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
    mutable std::vector<float> fPHHit0Vec;
    
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
    mutable std::vector<float> fPHHit1Vec;
    
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
    mutable std::vector<float> fPHHit2Vec;

    // Output tuples for all SpacePoints
    mutable std::vector<int>   fNumIDEsHit0Vec;
    mutable std::vector<int>   fNumIDEsHit1Vec;
    mutable std::vector<int>   fNumIDEsHit2Vec;
    mutable std::vector<int>   fNumIDEsSpacePointVec;
    
    mutable std::vector<float> fSPQualityVec;
    mutable std::vector<float> fSPTotalChargeVec;
    mutable std::vector<float> fSPAsymmetryVec;
    mutable std::vector<float> fSmallestPHVec;
    mutable std::vector<float> fLargestPHVec;
    mutable std::vector<float> fAveragePHVec;
    mutable std::vector<float> fLargestDelTVec;

    mutable std::vector<int>   fNumLongHitsVec;
    mutable std::vector<int>   fNumPlanesSimMatchVec;
    mutable std::vector<int>   fNumIntersectSetVec;

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
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;

    fTree = tree;

    fTree->Branch("CryostataVec",       "std::vector<int>",   &fCryoVec);
    fTree->Branch("TPCVec",             "std::vector<int>",   &fTPCVec);
    fTree->Branch("PlaneVec",           "std::vector<int>",   &fPlaneVec);
    
    // Set up specific branch for space points
    fSpacePointTree = tfs->makeAndRegister<TTree>("SpacePoint_t","SpacePoint Tuple");
 
    fSpacePointTree->Branch("NumIDEsHit0",        "std::vector<int>",   &fNumIDEsHit0Vec);
    fSpacePointTree->Branch("NumIDEsHit1",        "std::vector<int>",   &fNumIDEsHit1Vec);
    fSpacePointTree->Branch("NumIDEsHit2",        "std::vector<int>",   &fNumIDEsHit2Vec);
    fSpacePointTree->Branch("NumIDEsSpacePoint",  "std::vector<int>",   &fNumIDEsSpacePointVec);
    
    fSpacePointTree->Branch("SPQuality",          "std::vector<float>", &fSPQualityVec);
    fSpacePointTree->Branch("SPTotalCharge",      "std::vector<float>", &fSPTotalChargeVec);
    fSpacePointTree->Branch("SPAsymmetry",        "std::vector<float>", &fSPAsymmetryVec);
    fSpacePointTree->Branch("SmallestPH",         "std::vector<float>", &fSmallestPHVec);
    fSpacePointTree->Branch("LargestPH",          "std::vector<float>", &fLargestPHVec);
    fSpacePointTree->Branch("AveragePH",          "std::vector<float>", &fAveragePHVec);
    fSpacePointTree->Branch("LargestDelT",        "std::vector<float>", &fLargestDelTVec);

    fSpacePointTree->Branch("NumLongHitsSP",      "std::vector<int>",   &fNumLongHitsVec);
    fSpacePointTree->Branch("NumPlanesSimMatch",  "std::vector<int>",   &fNumPlanesSimMatchVec);
    fSpacePointTree->Branch("NumIntersectSet",    "std::vector<int>",   &fNumIntersectSetVec);
    
    // Set up specific branch for space points
    fMatchedHitTree = tfs->makeAndRegister<TTree>("MatchedHits_t","Matched Hits Tuple");

    fMatchedHitTree->Branch("TicksSimChannel0",   "std::vector<int>",   &fTicksSimChannel0Vec);
    fMatchedHitTree->Branch("TicksSimChanMost0",  "std::vector<int>",   &fTicksSimChanMost0Vec);
    fMatchedHitTree->Branch("TicksTotHit0",       "std::vector<int>",   &fTicksTotHit0Vec);
    fMatchedHitTree->Branch("TicksMaxSimRel0",    "std::vector<int>",   &fTicksMaxSimRel0Vec);
    fMatchedHitTree->Branch("TicksDiffSimHit0",   "std::vector<int>",   &fTicksDiffSimHit0Vec);
    fMatchedHitTree->Branch("EneTotDepHit0",      "std::vector<float>", &fEneTotDepHit0Vec);
    fMatchedHitTree->Branch("EneBestDepHit0",     "std::vector<float>", &fEneBestDepHit0Vec);
    fMatchedHitTree->Branch("EneMaxDepHit0",      "std::vector<float>", &fEneMaxDepHit0Vec);
    fMatchedHitTree->Branch("NDFHit0",            "std::vector<int>",   &fNDFHit0Vec);
    fMatchedHitTree->Branch("MultiplicityHit0",   "std::vector<int>",   &fMultiplicityHit0Vec);
    fMatchedHitTree->Branch("LocalIndexHit0",     "std::vector<int>",   &fLocalIndexHit0Vec);
    fMatchedHitTree->Branch("ChiSquareHit0",      "std::vector<float>", &fChiSquareHit0Vec);
    fMatchedHitTree->Branch("ChargeHit0",         "std::vector<float>", &fChargeHit0Vec);
    fMatchedHitTree->Branch("PulseHeightHit0",    "std::vector<float>", &fPHHit0Vec);

    fMatchedHitTree->Branch("TicksSimChannel1",   "std::vector<int>",   &fTicksSimChannel1Vec);
    fMatchedHitTree->Branch("TicksSimChanMost1",  "std::vector<int>",   &fTicksSimChanMost1Vec);
    fMatchedHitTree->Branch("TicksTotHit1",       "std::vector<int>",   &fTicksTotHit1Vec);
    fMatchedHitTree->Branch("TicksMaxSimRel1",    "std::vector<int>",   &fTicksMaxSimRel1Vec);
    fMatchedHitTree->Branch("TicksDiffSimHit1",   "std::vector<int>",   &fTicksDiffSimHit1Vec);
    fMatchedHitTree->Branch("EneTotDepHit1",      "std::vector<float>", &fEneTotDepHit1Vec);
    fMatchedHitTree->Branch("EneBestDepHit1",     "std::vector<float>", &fEneBestDepHit1Vec);
    fMatchedHitTree->Branch("EneMaxDepHit1",      "std::vector<float>", &fEneMaxDepHit1Vec);
    fMatchedHitTree->Branch("NDFHit1",            "std::vector<int>",   &fNDFHit1Vec);
    fMatchedHitTree->Branch("MultiplicityHit1",   "std::vector<int>",   &fMultiplicityHit1Vec);
    fMatchedHitTree->Branch("LocalIndexHit1",     "std::vector<int>",   &fLocalIndexHit1Vec);
    fMatchedHitTree->Branch("ChiSquareHit1",      "std::vector<float>", &fChiSquareHit1Vec);
    fMatchedHitTree->Branch("ChargeHit1",         "std::vector<float>", &fChargeHit1Vec);
    fMatchedHitTree->Branch("PulseHeightHit1",    "std::vector<float>", &fPHHit1Vec);

    fMatchedHitTree->Branch("TicksSimChannel2",   "std::vector<int>",   &fTicksSimChannel2Vec);
    fMatchedHitTree->Branch("TicksSimChanMost2",  "std::vector<int>",   &fTicksSimChanMost2Vec);
    fMatchedHitTree->Branch("TicksTotHit2",       "std::vector<int>",   &fTicksTotHit2Vec);
    fMatchedHitTree->Branch("TicksMaxSimRel2",    "std::vector<int>",   &fTicksMaxSimRel2Vec);
    fMatchedHitTree->Branch("TicksDiffSimHit2",   "std::vector<int>",   &fTicksDiffSimHit2Vec);
    fMatchedHitTree->Branch("EneTotDepHit2",      "std::vector<float>", &fEneTotDepHit2Vec);
    fMatchedHitTree->Branch("EneBestDepHit2",     "std::vector<float>", &fEneBestDepHit2Vec);
    fMatchedHitTree->Branch("EneMaxDepHit2",      "std::vector<float>", &fEneMaxDepHit2Vec);    
    fMatchedHitTree->Branch("NDFHit2",            "std::vector<int>",   &fNDFHit2Vec);
    fMatchedHitTree->Branch("MultiplicityHit2",   "std::vector<int>",   &fMultiplicityHit2Vec);
    fMatchedHitTree->Branch("LocalIndexHit2",     "std::vector<int>",   &fLocalIndexHit2Vec);
    fMatchedHitTree->Branch("ChiSquareHit2",      "std::vector<float>", &fChiSquareHit2Vec);
    fMatchedHitTree->Branch("ChargeHit2",         "std::vector<float>", &fChargeHit2Vec);
    fMatchedHitTree->Branch("PulseHeightHit2",    "std::vector<float>", &fPHHit2Vec);

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
    fPHHit0Vec.clear();

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
    fPHHit1Vec.clear();

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
    fPHHit2Vec.clear();

    fNumIDEsHit0Vec.clear();
    fNumIDEsHit1Vec.clear();
    fNumIDEsHit2Vec.clear();
    fNumIDEsSpacePointVec.clear();
    
    fSPQualityVec.clear();
    fSPTotalChargeVec.clear();
    fSPAsymmetryVec.clear();
    fSmallestPHVec.clear();
    fLargestPHVec.clear();
    fAveragePHVec.clear();
    fLargestDelTVec.clear();

    fNumLongHitsVec.clear();
    fNumPlanesSimMatchVec.clear();
    fNumIntersectSetVec.clear();

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

    // The following creates a trackID mapping 
    using TDCIDEPair             = std::pair<unsigned short, const sim::IDE*>;
    using TickTDCIDEVec          = std::vector<TDCIDEPair>;
    using ChanToTDCIDEMap        = std::unordered_map<raw::ChannelID_t,TickTDCIDEVec>;
    using TrackIDChanToTDCIDEMap = std::unordered_map<int,ChanToTDCIDEMap>;

    IDEToVoxelIDMap    ideToVoxelIDMap;
    VoxelIDToIDESetMap voxelIDToIDEMap;
    ChanToTDCToIDEMap  chanToTDCToIDEMap;
    VoxelIDSet         simChannelVoxelIDSet;

    TrackIDChanToTDCIDEMap trackIDChanToTDCIDEMap;

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

                trackIDChanToTDCIDEMap[ide.trackID][simChannel.Channel()].emplace_back(tdcide.first,&ide);
                
                if (ide.energy < std::numeric_limits<float>::epsilon()) mf::LogDebug("SpacePointAnalysis") << ">> epsilon simchan deposited energy: " << ide.energy << std::endl;
            }
        }
    }

    // More data structures, here we want to keep track of the start/peak/end of the charge deposit along a wire for a given track
    using ChargeDeposit        = std::tuple<TDCIDEPair,TDCIDEPair,TDCIDEPair,float>;
    using ChargeDepositVec     = std::vector<ChargeDeposit>;
    using ChanToChargeMap      = std::map<raw::ChannelID_t,ChargeDepositVec>;
    using TrackToChanChargeMap = std::unordered_map<int,ChanToChargeMap>;

    TrackToChanChargeMap trackToChanChargeMap;

    // Go through the list of track to ides and get the total deposited energy per track
    float bestTotDepEne(0.);
    int   bestTrackID(0);

    for(const auto& trackIDEPair : trackIDChanToTDCIDEMap)
    {
        ChanToChargeMap& chanToChargeMap = trackToChanChargeMap[trackIDEPair.first];

        float trackTotDepE(0.);

        for(const auto& chanTDCIDEPair : trackIDEPair.second)
        {
            ChargeDepositVec& chargeDepositVec = chanToChargeMap[chanTDCIDEPair.first];

            // Keep track of first,peak,last/ene
            TDCIDEPair firstPair = chanTDCIDEPair.second.front();
            TDCIDEPair peakPair  = firstPair;
            TDCIDEPair lastPair  = chanTDCIDEPair.second.back();

            // Keep watch for gaps
            TDCIDEPair prevPair  = firstPair;

            // Keep track of deposited energy on a snippet
            float snippetDepEne(0.);

            for(const auto& tdcIDEPair : chanTDCIDEPair.second)
            {
                float depEne = tdcIDEPair.second->energy;

                trackTotDepE += depEne;

                // Watch for a gap...
                if (tdcIDEPair.first - prevPair.first > 1)
                {
                    chargeDepositVec.emplace_back(firstPair,peakPair,prevPair,snippetDepEne);

                    firstPair     = tdcIDEPair;
                    peakPair      = firstPair;
                    snippetDepEne = 0.;
                }

                snippetDepEne += depEne;

                if (depEne > peakPair.second->energy) peakPair = tdcIDEPair;

                prevPair = tdcIDEPair;
            }

            chargeDepositVec.emplace_back(firstPair,peakPair,lastPair,snippetDepEne);
        }

        if (trackTotDepE > bestTotDepEne) 
        {
            bestTrackID   = trackIDEPair.first;
            bestTotDepEne = trackTotDepE;
        }
    }

    // Ok, for my next trick I want to build a mapping between hits and voxel IDs. Note that any given hit can be associated to more than one voxel...
    // We do this on the entire hit collection, ultimately we will want to consider SpacePoint efficiency (this could be done in the loop over SpacePoints
    // using the associated hits and would save time/memory)
    using VoxelIDSet           = std::set<sim::LArVoxelID>;
//    using VoxelIDSetVec        = std::vector<VoxelIDSet>;
//    using RecobHitToVoxelIDMap = std::unordered_map<const recob::Hit*, VoxelIDSetVec>;
    using RecobHitToVoxelIDMap = std::unordered_map<const recob::Hit*, VoxelIDSet>;
    
    RecobHitToVoxelIDMap recobHitToVoxelIDMap;

    // Recover the "best" track info to start
    TrackToChanChargeMap::const_iterator chanToChargeMap = trackToChanChargeMap.find(bestTrackID);
    
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
            //if (hit.GoodnessOfFit() < 0.) continue;

            // Recover channel information based on this hit
            ChanToChargeMap::const_iterator chanToChargeItr = chanToChargeMap->second.find(hit.Channel());

            // This at least weeds out the case where a channel displayed a hit but had no signal otherwise
            if (chanToChargeItr != chanToChargeMap->second.end())
            {
                // Recover hit time range (in ticks), cast a wide net here
                int peakTick  = std::round(hit.PeakTime());
                int startTick = std::max(   0,int(std::floor(hit.PeakTime() - 3. * hit.RMS())));
                int endTick   = std::min(4096,int(std::ceil(hit.PeakTime() + 3. * hit.RMS())));

                int startTDC = fClockService->TPCTick2TDC(startTick - fOffsetVec[hit.WireID().Plane]);
                int peakTDC  = fClockService->TPCTick2TDC(peakTick  - fOffsetVec[hit.WireID().Plane]);
                int endTDC   = fClockService->TPCTick2TDC(endTick   - fOffsetVec[hit.WireID().Plane]);

                // If we have a match then this iterator gets set to the matching values
                ChargeDepositVec::const_iterator chargeMatchItr = chanToChargeItr->second.end();

                int bestPeakDiff = std::numeric_limits<int>::max();

                // Match the hit (if there is one)
                for(ChargeDepositVec::const_iterator chargeInfoItr = chanToChargeItr->second.begin(); chargeInfoItr != chanToChargeItr->second.end(); chargeInfoItr++)
                {
                    // Require some amount of overlap between the hit and the sim info
                    if (endTDC > std::get<0>(*chargeInfoItr).first && startTDC < std::get<2>(*chargeInfoItr).first) 
                    {
                        const TDCIDEPair& peakTDCIDE = std::get<1>(*chargeInfoItr);

                        int peakDiff = peakTDC - int(peakTDCIDE.first);

                        if (std::abs(peakDiff) < std::abs(bestPeakDiff))
                        {
                            bestPeakDiff   = peakDiff;
                            chargeMatchItr = chargeInfoItr;
                        }
                    }
                }

                // If no match then skip
                if (chargeMatchItr == chanToChargeItr->second.end()) continue;

                // Find the max dep ene tick
                const ChargeDeposit& chargeDeposit = *chargeMatchItr;

                int   firstSimTick(std::get<0>(chargeDeposit).first);
                int   lastSimTick(std::get<2>(chargeDeposit).first);
                int   maxDepTick(std::get<1>(chargeDeposit).first);
                float maxDepEneTick(std::get<1>(chargeDeposit).second->energy);
                float bestDepEne(std::get<3>(chargeDeposit));
                float totDepEne(bestDepEne);                                                             //<**** fix this
                int   bestTicks(lastSimTick - firstSimTick + 1);

                // One final time through to find sim ticks that "matter"
                // We define this as the collection of IDE's that make up to 90% of the total deposited energy
                const TickTDCIDEVec& tickToTDCIDEVec = trackIDChanToTDCIDEMap[bestTrackID][hit.Channel()];
                TickTDCIDEVec        tickIDEVec;

                for(const auto& tickInfo : tickToTDCIDEVec)
                {
                    //if (tickInfo.first >= firstSimTick && tickInfo.first <= lastSimTick) tickIDEVec.emplace_back(tickInfo);
                    if (tickInfo.first >= startTDC && tickInfo.first <= endTDC) tickIDEVec.emplace_back(tickInfo);
                }

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
                recobHitToVoxelIDMap[&hit] = voxelIDSet;

                // Fill depending on the plane
                if (hit.WireID().Plane == 0)
                {
                    fTicksSimChannel0Vec.emplace_back(bestTicks);
                    fTicksSimChanMost0Vec.emplace_back(bestTicksGood);
                    fTicksTotHit0Vec.emplace_back(endTick-startTick+1);
                    fTicksMaxSimRel0Vec.emplace_back(maxDepTick-startTDC);
                    fTicksDiffSimHit0Vec.emplace_back(peakTick-startTick-bestTicks);
                    fEneTotDepHit0Vec.emplace_back(totDepEne);
                    fEneBestDepHit0Vec.emplace_back(bestDepEne);
                    fEneMaxDepHit0Vec.emplace_back(maxDepEneTick);
                    fNDFHit0Vec.emplace_back(hit.DegreesOfFreedom());
                    fMultiplicityHit0Vec.emplace_back(hit.Multiplicity());
                    fLocalIndexHit0Vec.emplace_back(hit.LocalIndex());
                    fChiSquareHit0Vec.emplace_back(hit.GoodnessOfFit());
                    fChargeHit0Vec.emplace_back(hit.SummedADC());
                    fPHHit0Vec.emplace_back(hit.PeakAmplitude());
                }
                else if (hit.WireID().Plane == 1)
                {
                    fTicksSimChannel1Vec.emplace_back(bestTicks);
                    fTicksSimChanMost1Vec.emplace_back(bestTicksGood);
                    fTicksTotHit1Vec.emplace_back(endTick-startTick+1);
                    fTicksMaxSimRel1Vec.emplace_back(maxDepTick-startTDC);
                    fTicksDiffSimHit1Vec.emplace_back(peakTick-startTick-bestTicks);
                    fEneTotDepHit1Vec.emplace_back(totDepEne);
                    fEneBestDepHit1Vec.emplace_back(bestDepEne);
                    fEneMaxDepHit1Vec.emplace_back(maxDepEneTick);
                    fNDFHit1Vec.emplace_back(hit.DegreesOfFreedom());
                    fMultiplicityHit1Vec.emplace_back(hit.Multiplicity());
                    fLocalIndexHit1Vec.emplace_back(hit.LocalIndex());
                    fChiSquareHit1Vec.emplace_back(hit.GoodnessOfFit());
                    fChargeHit1Vec.emplace_back(hit.SummedADC());
                    fPHHit1Vec.emplace_back(hit.PeakAmplitude());
                }
                else
                {
                    fTicksSimChannel2Vec.emplace_back(bestTicks);
                    fTicksSimChanMost2Vec.emplace_back(bestTicksGood);
                    fTicksTotHit2Vec.emplace_back(endTick-startTick+1);
                    fTicksMaxSimRel2Vec.emplace_back(maxDepTick-startTDC);
                    fTicksDiffSimHit2Vec.emplace_back(peakTick-startTick-bestTicks);
                    fEneTotDepHit2Vec.emplace_back(totDepEne);
                    fEneBestDepHit2Vec.emplace_back(bestDepEne);
                    fEneMaxDepHit2Vec.emplace_back(maxDepEneTick);
                    fNDFHit2Vec.emplace_back(hit.DegreesOfFreedom());
                    fMultiplicityHit2Vec.emplace_back(hit.Multiplicity());
                    fLocalIndexHit2Vec.emplace_back(hit.LocalIndex());
                    fChiSquareHit2Vec.emplace_back(hit.GoodnessOfFit());
                    fChargeHit2Vec.emplace_back(hit.SummedADC());
                    fPHHit2Vec.emplace_back(hit.PeakAmplitude());
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
                float largestPH   = 0.;
                int   numHits     = 0;
                float averagePH   = 0.;
                float averagePT   = 0.;
                float largestDelT = 0.;
                
                std::vector<int> numIDEsHitVec;
                int              numIDEsSpacePoint(0);
                int              numLongHits(0);
                int              numIntersections(0);

                std::vector<RecobHitToVoxelIDMap::const_iterator> recobHitToVoxelIterVec;
                
                // Now we can use our maps to find out if the hits making up the SpacePoint are truly related...
                for(const auto& hitPtr : associatedHits)
                {
                    RecobHitToVoxelIDMap::iterator hitToVoxelItr = recobHitToVoxelIDMap.find(hitPtr.get());
                    
                    float  peakAmplitude = hitPtr->PeakAmplitude();
                    
                    numHits++;
                    averagePH += peakAmplitude;
                    averagePT += hitPtr->PeakTime();

                    smallestPH = std::min(peakAmplitude,smallestPH);
                    largestPH  = std::max(peakAmplitude,largestPH);

                    if (hitPtr->DegreesOfFreedom() < 2) numLongHits++;
                    
                    if (hitToVoxelItr == recobHitToVoxelIDMap.end())
                    {
                        numIDEsHitVec.push_back(0);
                        continue;
                    }
                    
                    recobHitToVoxelIterVec.push_back(hitToVoxelItr);
                    numIDEsHitVec.push_back(hitToVoxelItr->second.size());
                }
                
                averagePH /= float(numHits);
                averagePT /= float(numHits);
                
                for(const auto& hitPtr : associatedHits)
                {
                    float delT = hitPtr->PeakTime() - averagePT;
                    
                    if (std::abs(delT) > std::abs(largestDelT)) largestDelT = delT;
                }
                
                // If a SpacePoint is made from "true" MC hits then we will have found the relations to the MC info for all three
                // hits. If this condition is not satisfied it means one or more hits making the SpacePoint are noise hits
                if (recobHitToVoxelIterVec.size() == 3)
                {
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

                        numIntersections++;
                        
                        // Again, no IDEs in the intersection means it is a ghost space point but, of course, we are hoping
                        // there are common IDEs so we can call it a real SpacePoint
                        if (!secondIntersectionVec.empty())
                        {
                            numIDEsSpacePoint = secondIntersectionVec.size();

                            numIntersections++;
                        }
                    }
                }
                
                // Fill for "all" cases
                fSPQualityVec.push_back(spQuality);
                fSPTotalChargeVec.push_back(spCharge);
                fSPAsymmetryVec.push_back(spAsymmetry);
                fSmallestPHVec.push_back(smallestPH);
                fLargestPHVec.push_back(largestPH);
                fAveragePHVec.push_back(averagePH);
                fLargestDelTVec.push_back(largestDelT);
                
                fNumIDEsHit0Vec.push_back(numIDEsHitVec[0]);
                fNumIDEsHit1Vec.push_back(numIDEsHitVec[1]);
                fNumIDEsHit2Vec.push_back(numIDEsHitVec[2]);
                fNumIDEsSpacePointVec.push_back(numIDEsSpacePoint);

                fNumLongHitsVec.emplace_back(numLongHits);
                fNumPlanesSimMatchVec.emplace_back(recobHitToVoxelIterVec.size());
                fNumIntersectSetVec.emplace_back(numIntersections);
            }
        }
    }

//            // Store tuple variables
//            fTPCVec.push_back(wids[0].TPC);
//            fCryoVec.push_back(wids[0].Cryostat);
//            fPlaneVec.push_back(wids[0].Plane);
//            fWireVec.push_back(wids[0].Wire);

    fSpacePointTree->Fill();
    fMatchedHitTree->Fill();

    return;
}
    
// Useful for normalizing histograms
void SpacePointAnalysis::endJob(int numEvents)
{
    return;
}
    
DEFINE_ART_CLASS_TOOL(SpacePointAnalysis)
}
