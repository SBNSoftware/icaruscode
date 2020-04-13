
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

// Define object to keep track of hit/spacepoint related items
class HitSpacePointObj
{
public:
    HitSpacePointObj() : fTree(nullptr) {}

    void setBranches(TTree* tree)
    {
        tree->Branch("NumIDEsHit0",        "std::vector<int>",   &fNumIDEsHit0Vec);
        tree->Branch("NumIDEsHit1",        "std::vector<int>",   &fNumIDEsHit1Vec);
        tree->Branch("NumIDEsHit2",        "std::vector<int>",   &fNumIDEsHit2Vec);
        tree->Branch("NumIDEsSpacePoint",  "std::vector<int>",   &fNumIDEsSpacePointVec);

        tree->Branch("SPQuality",          "std::vector<float>", &fSPQualityVec);
        tree->Branch("SPTotalCharge",      "std::vector<float>", &fSPTotalChargeVec);
        tree->Branch("SPAsymmetry",        "std::vector<float>", &fSPAsymmetryVec);
        tree->Branch("SmallestPH",         "std::vector<float>", &fSmallestPHVec);
        tree->Branch("LargestPH",          "std::vector<float>", &fLargestPHVec);
        tree->Branch("AveragePH",          "std::vector<float>", &fAveragePHVec);
        tree->Branch("LargestDelT",        "std::vector<float>", &fLargestDelTVec);

        tree->Branch("NumLongHitsSP",      "std::vector<int>",   &fNumLongHitsVec);
        tree->Branch("NumPlanesSimMatch",  "std::vector<int>",   &fNumPlanesSimMatchVec);
        tree->Branch("NumIntersectSet",    "std::vector<int>",   &fNumIntersectSetVec);

        fTree = tree;
    }

    void fill()
    {
        if (fTree) fTree->Fill();
    }

    void clear()
    {
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
    }

    // Define tuple vars, make public for direct access
    std::vector<int>   fNumIDEsHit0Vec;
    std::vector<int>   fNumIDEsHit1Vec;
    std::vector<int>   fNumIDEsHit2Vec;
    std::vector<int>   fNumIDEsSpacePointVec;

    std::vector<float> fSPQualityVec;
    std::vector<float> fSPTotalChargeVec;
    std::vector<float> fSPAsymmetryVec;
    std::vector<float> fSmallestPHVec;
    std::vector<float> fLargestPHVec;
    std::vector<float> fAveragePHVec;
    std::vector<float> fLargestDelTVec;

    std::vector<int>   fNumLongHitsVec;
    std::vector<int>   fNumPlanesSimMatchVec;
    std::vector<int>   fNumIntersectSetVec;

private:
    TTree* fTree;
};


// Define object to keep track of hit/sim related tuple items 
class HitSimulationTupleObj
{
public:
    HitSimulationTupleObj() : fTree(nullptr) {}

    void setBranches(TTree* tree) 
    {
        tree->Branch("TicksSimChannel",   "std::vector<int>",   &fTicksSimChannelVec);
        tree->Branch("TicksSimChanMost",  "std::vector<int>",   &fTicksSimChanMostVec);
        tree->Branch("TicksTotHit",       "std::vector<int>",   &fTicksTotHitVec);
        tree->Branch("TicksMaxSimRel",    "std::vector<int>",   &fTicksMaxSimRelVec);
        tree->Branch("TicksDiffSimHit",   "std::vector<int>",   &fTicksDiffSimHitVec);
        tree->Branch("EneTotDepHit",      "std::vector<float>", &fEneTotDepHitVec);
        tree->Branch("NElecTotalHit",     "std::vector<float>", &fNElecTotHitVec);           //< Total number elecrons (all sources) for hit
        tree->Branch("EneBestDepHit",     "std::vector<float>", &fEneBestDepHitVec);
        tree->Branch("NElecBestHit",      "std::vector<float>", &fNElecBestHitVec);          //< # electrons from primary track
        tree->Branch("EneMaxDepHit",      "std::vector<float>", &fEneMaxDepHitVec);
        tree->Branch("NDF",               "std::vector<int>",   &fNDFHitVec);                //< Number of degrees of freedom of hit fit
        tree->Branch("Multiplicity",      "std::vector<int>",   &fMultiplicityHitVec);       //< Multiplicity of the snippet the hit is on
        tree->Branch("LocalIndex",        "std::vector<int>",   &fLocalIndexHitVec);         //< The index of the hit within the snippet
        tree->Branch("TimeOrder",         "std::vector<int>",   &fTimeOrderHitVec);          //< Time order of the hit (selection variable)
        tree->Branch("ChiSquare",         "std::vector<float>", &fChiSquareHitVec);          //< Chi square of fit 
        tree->Branch("SummedADC",         "std::vector<float>", &fSummedADCHitVec);          //< Sum of all ADC values start/end of snippet
        tree->Branch("Integral",          "std::vector<float>", &fIntegralHitVec);           //< Integrated charge +/- n sigma about peak center
        tree->Branch("PulseHeight",       "std::vector<float>", &fPHHitVec);                 //< Pulse height of hit
        tree->Branch("RMS",               "std::vector<float>", &fRMSHitVec);                //< RMS of hit (from fit)

        tree->Branch("PulseHeightOrder",  "std::vector<int>",   &fPHOrderHitVec);            //< Local index ordered by pulse height

        fTree = tree;

        return;
    }

    void fill()
    {
        if (fTree) fTree->Fill();
    }

    void clear()
    {
        fTicksSimChannelVec.clear();
        fTicksSimChanMostVec.clear();
        fTicksTotHitVec.clear();
        fTicksMaxSimRelVec.clear();
        fTicksDiffSimHitVec.clear();
        fEneTotDepHitVec.clear();
        fNElecTotHitVec.clear();
        fEneBestDepHitVec.clear();
        fNElecBestHitVec.clear();
        fEneMaxDepHitVec.clear();
        fNDFHitVec.clear();
        fMultiplicityHitVec.clear();
        fLocalIndexHitVec.clear();
        fTimeOrderHitVec.clear();
        fChiSquareHitVec.clear();
        fSummedADCHitVec.clear();
        fIntegralHitVec.clear();
        fPHHitVec.clear();
        fRMSHitVec.clear();

        fPHOrderHitVec.clear();
    }

    void fillSimInfo(int   ticksSimChannel,
                     int   ticksSimChanMost,
                     float totDepEne,
                     float totNumElectrons,
                     float bestDepEne,
                     float bestNumElectrons,
                     float maxDepEneTick
                     )
    {
        fTicksSimChannelVec.emplace_back(ticksSimChannel);
        fTicksSimChanMostVec.emplace_back(ticksSimChanMost);
        fEneTotDepHitVec.emplace_back(totDepEne);
        fNElecTotHitVec.emplace_back(totNumElectrons);
        fEneBestDepHitVec.emplace_back(bestDepEne);
        fNElecBestHitVec.emplace_back(bestNumElectrons);
        fEneMaxDepHitVec.emplace_back(maxDepEneTick);
    }

    void fillMixedInfo(int hitWidth,
                       int ticksToMax,
                       int deltaTicks)
    {
        fTicksTotHitVec.emplace_back(hitWidth);
        fTicksMaxSimRelVec.emplace_back(ticksToMax);
        fTicksDiffSimHitVec.emplace_back(deltaTicks);
    }

    void fillHitInfo(const recob::Hit* hit,int hitOrder)
    {
        fNDFHitVec.emplace_back(hit->DegreesOfFreedom());
        fMultiplicityHitVec.emplace_back(hit->Multiplicity());
        fLocalIndexHitVec.emplace_back(hit->LocalIndex());
        fTimeOrderHitVec.emplace_back(hitOrder);
        fChiSquareHitVec.emplace_back(hit->GoodnessOfFit());
        fSummedADCHitVec.emplace_back(hit->SummedADC());
        fIntegralHitVec.emplace_back(hit->Integral());
        fPHHitVec.emplace_back(hit->PeakAmplitude());
        fRMSHitVec.emplace_back(hit->RMS());
    }

    // Define tuple values, these are public so can be diretly accessed for filling
    std::vector<int>   fTicksSimChannelVec;
    std::vector<int>   fTicksSimChanMostVec;
    std::vector<int>   fTicksTotHitVec;
    std::vector<int>   fTicksMaxSimRelVec;
    std::vector<int>   fTicksDiffSimHitVec;
    std::vector<float> fEneTotDepHitVec;
    std::vector<float> fNElecTotHitVec;
    std::vector<float> fEneBestDepHitVec;
    std::vector<float> fNElecBestHitVec;
    std::vector<float> fEneMaxDepHitVec;
    std::vector<int>   fNDFHitVec;
    std::vector<int>   fMultiplicityHitVec;
    std::vector<int>   fLocalIndexHitVec;
    std::vector<int>   fTimeOrderHitVec;
    std::vector<float> fChiSquareHitVec;
    std::vector<float> fSummedADCHitVec;
    std::vector<float> fIntegralHitVec;
    std::vector<float> fPHHitVec;
    std::vector<float> fRMSHitVec;

    std::vector<int>   fPHOrderHitVec;

private:
    TTree* fTree;
};


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
   
    // Create a struct allowing us to sort IDEs in a set by largest to smallest energy
    struct ideCompare
    {
        bool operator() (const sim::IDE* left, const sim::IDE* right) const {return left->energy > right->energy;}
    };

    // Define structures for relating SimChannel to Voxels
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

    // More data structures, here we want to keep track of the start/peak/end of the charge deposit along a wire for a given track
    using ChargeDeposit        = std::tuple<TDCIDEPair,TDCIDEPair,TDCIDEPair,float,float>;
    using ChargeDepositVec     = std::vector<ChargeDeposit>;
    using ChanToChargeMap      = std::map<raw::ChannelID_t,ChargeDepositVec>;
    using TrackToChanChargeMap = std::unordered_map<int,ChanToChargeMap>;

    // Define a function to map IDE's from SimChannel objects to Track IDs
    void makeTrackToChanChargeMap(const TrackIDChanToTDCIDEMap&, TrackToChanChargeMap&, float&, int&) const;

    // Relate hits to voxels
    using HitPointerVec        = std::vector<const recob::Hit*>;
    using RecobHitToVoxelIDMap = std::unordered_map<const recob::Hit*, VoxelIDSet>;

    void compareHitsToSim(const art::Event&, const ChanToTDCToIDEMap&, const ChanToChargeMap&, const ChanToTDCIDEMap&, const IDEToVoxelIDMap&, RecobHitToVoxelIDMap&) const;

    void matchHitSim(const HitPointerVec&, const ChanToTDCToIDEMap&, const ChargeDepositVec&, const ChanToTDCIDEMap&, const IDEToVoxelIDMap&, RecobHitToVoxelIDMap&) const;

    void compareSpacePointsToSim(const art::Event&, const RecobHitToVoxelIDMap&) const;

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
    mutable TTree*              fTree;
    
    mutable std::vector<int>    fTPCVec;
    mutable std::vector<int>    fCryoVec;
    mutable std::vector<int>    fPlaneVec;

    using HitSimObjVec = std::vector<HitSimulationTupleObj>;

    mutable HitSimObjVec        fHitSimObjVec;
    mutable HitSpacePointObj    fHitSpacePointObj;

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
    TTree* locTree = tfs->makeAndRegister<TTree>("SpacePoint_t","SpacePoint Tuple");

    fHitSpacePointObj.setBranches(locTree);

    fHitSimObjVec.resize(fGeometry->Nplanes());

    for(size_t plane = 0; plane < fGeometry->Nplanes(); plane++)
    {
        // Set up specific branch for space points
        locTree = tfs->makeAndRegister<TTree>("MatchedHits_P"+std::to_string(plane),"Matched Hits Tuple plane "+std::to_string(plane));

        fHitSimObjVec[plane].setBranches(locTree);
    }

    clear();

    return;
}
    
void SpacePointAnalysis::clear() const
{
    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();

    fHitSpacePointObj.clear();

    for(auto& hitObj : fHitSimObjVec) hitObj.clear();

    return;
}
     
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
    TrackToChanChargeMap trackToChanChargeMap;

    // Go through the list of track to ides and get the total deposited energy per track
    float bestTotDepEne(0.);
    int   bestTrackID(0);

    makeTrackToChanChargeMap(trackIDChanToTDCIDEMap, trackToChanChargeMap, bestTotDepEne, bestTrackID);

    // Ok, for my next trick I want to build a mapping between hits and voxel IDs. Note that any given hit can be associated to more than one voxel...
    // We do this on the entire hit collection, ultimately we will want to consider SpacePoint efficiency (this could be done in the loop over SpacePoints
    // using the associated hits and would save time/memory)
    using VoxelIDSet           = std::set<sim::LArVoxelID>;
    using RecobHitToVoxelIDMap = std::unordered_map<const recob::Hit*, VoxelIDSet>;
    
    RecobHitToVoxelIDMap recobHitToVoxelIDMap;

    ChanToTDCIDEMap& chanToTDCIDEMap = trackIDChanToTDCIDEMap[bestTrackID];

    // Recover the "best" track info to start
    TrackToChanChargeMap::const_iterator chanToChargeMapItr = trackToChanChargeMap.find(bestTrackID);

    // Process the hit/simulation 
    compareHitsToSim(event, chanToTDCToIDEMap, chanToChargeMapItr->second, chanToTDCIDEMap, ideToVoxelIDMap, recobHitToVoxelIDMap);

    // Now do the space points
    compareSpacePointsToSim(event, recobHitToVoxelIDMap);
    
    // Make sure the output tuples are filled
    fHitSpacePointObj.fill();

    for(auto& hitObj : fHitSimObjVec) hitObj.fill();

    return;
}

void SpacePointAnalysis::makeTrackToChanChargeMap(const TrackIDChanToTDCIDEMap& trackIDChanToTDCIDEMap, 
                                                  TrackToChanChargeMap&         trackToChanChargeMap,
                                                  float&                        bestTotDepEne,
                                                  int&                          bestTrackID) const
{
    // Pretty straightforward looping here...
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
            float snippetNumElectrons(0.);

            for(const auto& tdcIDEPair : chanTDCIDEPair.second)
            {
                float depEne = tdcIDEPair.second->energy;

                trackTotDepE += depEne;

                // Watch for a gap...
                if (tdcIDEPair.first - prevPair.first > 1)
                {
                    chargeDepositVec.emplace_back(firstPair,peakPair,prevPair,snippetDepEne,snippetNumElectrons);

                    firstPair           = tdcIDEPair;
                    peakPair            = firstPair;
                    snippetDepEne       = 0.;
                    snippetNumElectrons = 0.;
                }

                snippetDepEne       += depEne;
                snippetNumElectrons += tdcIDEPair.second->numElectrons;

                if (depEne > peakPair.second->energy) peakPair = tdcIDEPair;

                prevPair = tdcIDEPair;
            }

            chargeDepositVec.emplace_back(firstPair,peakPair,lastPair,snippetDepEne,snippetNumElectrons);
        }

        if (trackTotDepE > bestTotDepEne) 
        {
            bestTrackID   = trackIDEPair.first;
            bestTotDepEne = trackTotDepE;
        }
    }

    return;
}

void SpacePointAnalysis::compareHitsToSim(const art::Event&        event,                          // For recovering data from event store
                                          const ChanToTDCToIDEMap& chanToTDCToIDEMap,              // This gives us ability to retrieve total charge deposits
                                          const ChanToChargeMap&   chanToChargeMap,                // Charge deposit for specific track 
                                          const ChanToTDCIDEMap&   chanToTDCIDEMap,                // Charge deposit for specific track
                                          const IDEToVoxelIDMap&   ideToVoxelIDMap,                // Mapping of ide info to voxels
                                          RecobHitToVoxelIDMap&    recobHitToVoxelIDMap) const     // The output info
{
    // We start by building a mapping between channels and lists of hits on that channel (ordered by time)
    using HitPointerVec   = std::vector<const recob::Hit*>;
    using ChanToHitVecMap = std::map<raw::ChannelID_t,HitPointerVec>;

    ChanToHitVecMap chanToHitVecMap;

    // And now fill it
    for(const auto& hitLabel : fRecobHitLabelVec)
    {
        art::Handle< std::vector<recob::Hit> > hitHandle;
        event.getByLabel(hitLabel, hitHandle);

        // If no hits then skip
        if ((*hitHandle).empty()) continue;

        for(const auto& hit : *hitHandle) chanToHitVecMap[hit.Channel()].emplace_back(&hit);
    }

    // Now go through and order each vector of hits by time
    for(auto& chanToHitPair : chanToHitVecMap)
    {
        HitPointerVec& hitPtrVec = chanToHitPair.second;

        std::sort(hitPtrVec.begin(),
                  hitPtrVec.end(),
                  [](const auto& left, const auto& right){return left->Channel() == right->Channel() ? left->PeakTime() < right->PeakTime() : left->Channel() < right->Channel();});
    }

    // The idea is to loop over the input sim information so we can look at efficiency as well as resolution issues
    for(const auto& chanToChargePair : chanToChargeMap)
    {
        // Recover the channel 
        raw::ChannelID_t channel = chanToChargePair.first;

        // Look up the hits associated to this channel
        ChanToHitVecMap::const_iterator chanToHitVecItr = chanToHitVecMap.find(channel);

        // For now we simply punt...
        if (chanToHitVecItr == chanToHitVecMap.end()) continue;

        // Recover channel information based on this hit
        const ChargeDepositVec& chargeDepositVec = chanToChargePair.second;

        // Get the hits... 
        const HitPointerVec& hitPtrVec = chanToHitVecItr->second;

        short int lastSnippetStart(-1);

        HitPointerVec hitVec;
        
        // Outer loop over hits in this hit collection
        for(const auto& hitPtr : hitPtrVec)
        {
            // We want to collect together the hits that are on the same snippet. Hits will come grouped and in order 
            // along the snippet, so we simply keep them in a local vector until we hit the end of the snippet...
            //
            // ** It is worth noting this scheme as implemented will miss the last snippet of hits... so think about 
            // that for the future
            short int snippetStart = hitPtr->StartTick();

            if (snippetStart == lastSnippetStart)
            {
                lastSnippetStart = snippetStart;

                hitVec.emplace_back(hitPtr);
                
                continue;
            }

            // Process the current list of hits (which will be on the same snippet)
            matchHitSim(hitVec, chanToTDCToIDEMap, chargeDepositVec, chanToTDCIDEMap, ideToVoxelIDMap, recobHitToVoxelIDMap);

            hitVec.clear();
            hitVec.emplace_back(hitPtr);
            lastSnippetStart = snippetStart;
        }

        // Make sure to catch the last set of hits in the group
        if (!hitVec.empty()) matchHitSim(hitVec, chanToTDCToIDEMap, chargeDepositVec, chanToTDCIDEMap, ideToVoxelIDMap, recobHitToVoxelIDMap);
    }

    return;
}

void SpacePointAnalysis::matchHitSim(const HitPointerVec&     hitPointerVec,                  // Hits to match to simulation
                                     const ChanToTDCToIDEMap& chanToTDCToIDEMap,              // This gives us ability to retrieve total charge deposits
                                     const ChargeDepositVec&  chargeDepositVec,               // Charge deposit for specific track 
                                     const ChanToTDCIDEMap&   chanToTDCIDEMap,                // Charge deposit for specific track
                                     const IDEToVoxelIDMap&   ideToVoxelIDMap,                // Mapping of ide info to voxels
                                     RecobHitToVoxelIDMap&    recobHitToVoxelIDMap) const     // The output info
{
    // Data structure to allow ordering of multiple hits in a snippet
    using HitPeakTimeChargeTuple = std::tuple<int,const recob::Hit*,ChargeDepositVec::const_iterator>;
    using HitPeakTimeChargeVec   = std::vector<HitPeakTimeChargeTuple>;

    HitPeakTimeChargeVec hitPeakTimeChargeVec;

    // If here then we are on to the next hit, so we need to process our current list
    for(const auto& hit : hitPointerVec)
    {
        // Recover hit time range (in ticks), cast a wide net here
        int peakTick  = std::round(hit->PeakTime());
        int startTick = std::max(   0,int(std::floor(hit->PeakTime() - 3. * hit->RMS())));
        int endTick   = std::min(4096,int(std::ceil(hit->PeakTime() + 3. * hit->RMS())));

        int startTDC = fClockService->TPCTick2TDC(startTick - fOffsetVec[hit->WireID().Plane]);
        int peakTDC  = fClockService->TPCTick2TDC(peakTick  - fOffsetVec[hit->WireID().Plane]);
        int endTDC   = fClockService->TPCTick2TDC(endTick   - fOffsetVec[hit->WireID().Plane]);

        // If we have a match then this iterator gets set to the matching values
        ChargeDepositVec::const_iterator chargeMatchItr = chargeDepositVec.end();

        int bestPeakDiff = std::numeric_limits<int>::max();

        // Match the hit (if there is one)
        for(ChargeDepositVec::const_iterator chargeInfoItr = chargeDepositVec.begin(); chargeInfoItr != chargeDepositVec.end(); chargeInfoItr++)
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
        if (chargeMatchItr == chargeDepositVec.end()) continue;

        hitPeakTimeChargeVec.emplace_back(std::make_tuple(bestPeakDiff,hit,chargeMatchItr));
    }

    if (!hitPeakTimeChargeVec.empty())
    {
        // Ok, now we sort this vector by smallest peak time
        std::sort(hitPeakTimeChargeVec.begin(),hitPeakTimeChargeVec.end(),[](const auto& left,const auto& right){return std::abs(std::get<0>(left)) < std::abs(std::get<0>(right));});

        // Keep track of hit ordering on this snippet 
        int hitOrder(0);

        HitSimulationTupleObj& hitObj = fHitSimObjVec[std::get<1>(hitPeakTimeChargeVec.front())->WireID().Plane];

        // Now loop through
        for(const auto& hitPeakCharge : hitPeakTimeChargeVec)
        {
            const recob::Hit*    hit           = std::get<1>(hitPeakCharge);
            const ChargeDeposit& chargeDeposit = *std::get<2>(hitPeakCharge);

            // Recover hit time range (in ticks), cast a wide net here
            int   peakTick  = std::round(hit->PeakTime());
            int   startTick = std::max(   0,int(std::floor(hit->PeakTime() - 3. * hit->RMS())));
            int   endTick   = std::min(4096,int(std::ceil(hit->PeakTime() + 3. * hit->RMS())));

            int   startTDC = fClockService->TPCTick2TDC(startTick - fOffsetVec[hit->WireID().Plane]);
            int   peakTDC  = fClockService->TPCTick2TDC(peakTick  - fOffsetVec[hit->WireID().Plane]);
            int   endTDC   = fClockService->TPCTick2TDC(endTick   - fOffsetVec[hit->WireID().Plane]);

            int   firstSimTick(std::get<0>(chargeDeposit).first);
            int   lastSimTick(std::get<2>(chargeDeposit).first);
            int   maxDepTick(std::get<1>(chargeDeposit).first);
            float maxDepEneTick(std::get<1>(chargeDeposit).second->energy);
            float bestNumElectrons(std::get<4>(chargeDeposit));
            float bestDepEne(std::get<3>(chargeDeposit));
            float totDepEne(0.); 
            float totNumElectrons(0.);
            int   bestTicks(lastSimTick - firstSimTick + 1);

            // We want to get the total energy deposit from all particles in the ticks for this hit
            const TDCToIDEMap& tdcToIDEMap = chanToTDCToIDEMap.find(hit->Channel())->second;
            for(const auto& tdcToIDEPair : tdcToIDEMap)
            {
                for(const auto& ide : tdcToIDEPair.second) 
                {
                    totDepEne       += ide->energy;
                    totNumElectrons += ide->numElectrons;
                }
            }

            // One final time through to find sim ticks that "matter"
            // We define this as the collection of IDE's that make up to 90% of the total deposit
            ChanToTDCIDEMap::const_iterator tickToTDCIDEVecItr = chanToTDCIDEMap.find(hit->Channel());

            if (tickToTDCIDEVecItr == chanToTDCIDEMap.end()) continue;

            const TickTDCIDEVec& tickToTDCIDEVec = tickToTDCIDEVecItr->second;
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
                IDEToVoxelIDMap::const_iterator ideToVoxelIDMapItr = ideToVoxelIDMap.find(tickInfo.second);

                if (ideToVoxelIDMapItr == ideToVoxelIDMap.end()) continue;

                const sim::LArVoxelID& voxelID = ideToVoxelIDMapItr->second;

                sumEne += tickInfo.second->energy;
                bestTicksGood++;

                voxelIDSet.insert(voxelID);

                if (sumEne > 0.9 * bestDepEne) break;
            }

            // Finally, grab the voxels from the track leaving the most energy
            recobHitToVoxelIDMap[hit] = voxelIDSet;

            hitObj.fillSimInfo(bestTicks, bestTicksGood, totDepEne, totNumElectrons, bestDepEne, bestNumElectrons, maxDepEneTick);
            hitObj.fillMixedInfo(endTick-startTick+1, maxDepTick-startTDC, peakTDC-maxDepTick);
            hitObj.fillHitInfo(hit,hitOrder++);
        }

        // Resort in pulse height order (largest to smallest)
        std::sort(hitPeakTimeChargeVec.begin(),hitPeakTimeChargeVec.end(),[](const auto& left,const auto& right){return std::get<1>(left)->PeakAmplitude() > std::get<1>(right)->PeakAmplitude();});

        // Now loop through
        for(const auto& hitPeakCharge : hitPeakTimeChargeVec) hitObj.fPHOrderHitVec.emplace_back(std::get<1>(hitPeakCharge)->LocalIndex());
    }

    return;
}

void SpacePointAnalysis::compareSpacePointsToSim(const art::Event& event, const RecobHitToVoxelIDMap& recobHitToVoxelIDMap) const
{
    // Armed with these maps we can now process the SpacePoints...
    if (!recobHitToVoxelIDMap.empty())
    {
        // So now we loop through the various SpacePoint sources
        for(const auto& spacePointLabel : fSpacePointLabelVec)
        {
            art::Handle< std::vector<recob::SpacePoint>> spacePointHandle;
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
                    RecobHitToVoxelIDMap::const_iterator hitToVoxelItr = recobHitToVoxelIDMap.find(hitPtr.get());
                    
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
                fHitSpacePointObj.fSPQualityVec.push_back(spQuality);
                fHitSpacePointObj.fSPTotalChargeVec.push_back(spCharge);
                fHitSpacePointObj.fSPAsymmetryVec.push_back(spAsymmetry);
                fHitSpacePointObj.fSmallestPHVec.push_back(smallestPH);
                fHitSpacePointObj.fLargestPHVec.push_back(largestPH);
                fHitSpacePointObj.fAveragePHVec.push_back(averagePH);
                fHitSpacePointObj.fLargestDelTVec.push_back(largestDelT);

                fHitSpacePointObj.fNumIDEsHit0Vec.push_back(numIDEsHitVec[0]);
                fHitSpacePointObj.fNumIDEsHit1Vec.push_back(numIDEsHitVec[1]);
                fHitSpacePointObj.fNumIDEsHit2Vec.push_back(numIDEsHitVec[2]);
                fHitSpacePointObj.fNumIDEsSpacePointVec.push_back(numIDEsSpacePoint);

                fHitSpacePointObj.fNumLongHitsVec.emplace_back(numLongHits);
                fHitSpacePointObj.fNumPlanesSimMatchVec.emplace_back(recobHitToVoxelIterVec.size());
                fHitSpacePointObj.fNumIntersectSetVec.emplace_back(numIntersections);
            }
        }
    }

    return;
}

    
// Useful for normalizing histograms
void SpacePointAnalysis::endJob(int numEvents)
{
    return;
}
    
DEFINE_ART_CLASS_TOOL(SpacePointAnalysis)
}
