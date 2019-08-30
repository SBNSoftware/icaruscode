
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

    // Conversion factors
    double                      fPositionToVoxelIDFactor;
    std::vector<double>         fVolumeOffsetsLow;
    std::vector<double>         fVolumeOffsetsHigh;
    std::vector<int>            fNumVoxelsByAxis;
    
    // TTree variables
    mutable TTree*             fTree;
    
    mutable std::vector<int>   fTPCVec;
    mutable std::vector<int>   fCryoVec;
    mutable std::vector<int>   fPlaneVec;
    
    // Look at matching of SimChannel and SimEnergyDeposit
    mutable std::vector<float> fSEDeltaX;
    mutable std::vector<float> fSEDeltaY;
    mutable std::vector<float> fSEDeltaZ;
    mutable std::vector<float> fSEDistance;

    // Output tuples for all SpacePoints
    mutable std::vector<int>   fNumIDEsHit0Vec;
    mutable std::vector<int>   fNumIDEsHit1Vec;
    mutable std::vector<int>   fNumIDEsHit2Vec;
    mutable std::vector<int>   fNumIDEsSpacePointVec;
    
    mutable std::vector<float> fSPQualityVec;
    mutable std::vector<float> fSPTotalChargeVec;
    mutable std::vector<float> fSPAsymmetryVec;
    
    mutable std::vector<float> fSmallestPHVec;
    
    // Output tuples for SpacePoints with one or more hits not matching MC
    mutable std::vector<int>   fNumIDEsHit0NoMVec;
    mutable std::vector<int>   fNumIDEsHit1NoMVec;
    mutable std::vector<int>   fNumIDEsHit2NoMVec;
    
    mutable std::vector<float> fSPQualityNoMVec;
    mutable std::vector<float> fSPTotalChargeNoMVec;
    mutable std::vector<float> fSPAsymmetryNoMVec;
    
    mutable std::vector<float> fSmallestPHNoMVec;

    // Output tuples for Ghost SpacePoints
    mutable std::vector<int>   fNumIDEsHit0GhostVec;
    mutable std::vector<int>   fNumIDEsHit1GhostVec;
    mutable std::vector<int>   fNumIDEsHit2GhostVec;
    
    mutable std::vector<float> fSPQualityGhostVec;
    mutable std::vector<float> fSPTotalChargeGhostVec;
    mutable std::vector<float> fSPAsymmetryGhostVec;
    
    mutable std::vector<float> fSmallestPHGhostVec;

    // Output tuples for matched space points
    mutable std::vector<int>   fNumIDEsHit0MatchVec;
    mutable std::vector<int>   fNumIDEsHit1MatchVec;
    mutable std::vector<int>   fNumIDEsHit2MatchVec;
    mutable std::vector<int>   fNumIDEsSpacePointMatchVec;

    mutable std::vector<float> fSPQualityMatchVec;
    mutable std::vector<float> fSPTotalChargeMatchVec;
    mutable std::vector<float> fSPAsymmetryMatchVec;
    
    mutable std::vector<float> fSmallestPHMatchVec;

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

    // Set up the voxel mapping parameters (note wire pitch is the same for all planes in ICARUS)
    fPositionToVoxelIDFactor = 10. / fGeometry->WirePitch();    // voxel len/wid/height is 1/10 wire pitch
    
    // Recover world volume offsets
    fVolumeOffsetsLow.resize(3,0.);
    fVolumeOffsetsHigh.resize(3,0.);
    
    fGeometry->WorldBox(&fVolumeOffsetsLow[0],&fVolumeOffsetsHigh[0],&fVolumeOffsetsLow[1],&fVolumeOffsetsHigh[1],&fVolumeOffsetsLow[2],&fVolumeOffsetsHigh[2]);
    
    // Get number of voxels in each direction
    fNumVoxelsByAxis.resize(3,0);
    
    fNumVoxelsByAxis[0] = int(fPositionToVoxelIDFactor * (fVolumeOffsetsHigh[0] - fVolumeOffsetsLow[0])) + 1;
    fNumVoxelsByAxis[1] = int(fPositionToVoxelIDFactor * (fVolumeOffsetsHigh[1] - fVolumeOffsetsLow[1])) + 1;
    fNumVoxelsByAxis[2] = int(fPositionToVoxelIDFactor * (fVolumeOffsetsHigh[2] - fVolumeOffsetsLow[2])) + 1;

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

    fTree->Branch("SEDeltaX",          "std::vector<float>", &fSEDeltaX);
    fTree->Branch("SEDeltaY",          "std::vector<float>", &fSEDeltaY);
    fTree->Branch("SEDeltaZ",          "std::vector<float>", &fSEDeltaZ);
    fTree->Branch("SEDistance",        "std::vector<float>", &fSEDistance);

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
    
    fTree->Branch("NumIDEsHit0NoM",     "std::vector<int>",   &fNumIDEsHit0NoMVec);
    fTree->Branch("NumIDEsHit1NoM",     "std::vector<int>",   &fNumIDEsHit1NoMVec);
    fTree->Branch("NumIDEsHit2NoM",     "std::vector<int>",   &fNumIDEsHit2NoMVec);

    fTree->Branch("SPQualityNoM",       "std::vector<float>", &fSPQualityNoMVec);
    fTree->Branch("SPTotalChargeNoM",   "std::vector<float>", &fSPTotalChargeNoMVec);
    fTree->Branch("SPAsymmetryNoM",     "std::vector<float>", &fSPAsymmetryNoMVec);
    fTree->Branch("SmallestPHNoM",      "std::vector<float>", &fSmallestPHNoMVec);
    
    fTree->Branch("NumIDEsHit0Match",   "std::vector<int>",   &fNumIDEsHit0MatchVec);
    fTree->Branch("NumIDEsHit1Match",   "std::vector<int>",   &fNumIDEsHit1MatchVec);
    fTree->Branch("NumIDEsHit2Match",   "std::vector<int>",   &fNumIDEsHit2MatchVec);
    fTree->Branch("NumIDEsSPMatch",     "std::vector<int>",   &fNumIDEsSpacePointMatchVec);
    fTree->Branch("SmallestPHMatch",    "std::vector<float>", &fSmallestPHMatchVec);
    
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
    
    clear();

    return;
}
    
void SpacePointAnalysis::clear() const
{
    fSEDeltaX.clear();
    fSEDeltaY.clear();
    fSEDeltaZ.clear();
    fSEDistance.clear();

    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();
    
    fNumIDEsHit0Vec.clear();
    fNumIDEsHit1Vec.clear();
    fNumIDEsHit2Vec.clear();
    fNumIDEsSpacePointVec.clear();
    
    fSPQualityVec.clear();
    fSPTotalChargeVec.clear();
    fSPAsymmetryVec.clear();
    
    fNumIDEsHit0NoMVec.clear();
    fNumIDEsHit1NoMVec.clear();
    fNumIDEsHit2NoMVec.clear();

    fSPQualityNoMVec.clear();
    fSPTotalChargeNoMVec.clear();
    fSPAsymmetryNoMVec.clear();
    
    fNumIDEsHit0GhostVec.clear();
    fNumIDEsHit1GhostVec.clear();
    fNumIDEsHit2GhostVec.clear();

    fSPQualityGhostVec.clear();
    fSPTotalChargeGhostVec.clear();
    fSPAsymmetryGhostVec.clear();
    
    fNumIDEsHit0MatchVec.clear();
    fNumIDEsHit1MatchVec.clear();
    fNumIDEsHit2MatchVec.clear();
    fNumIDEsSpacePointMatchVec.clear();

    fSPQualityMatchVec.clear();
    fSPTotalChargeMatchVec.clear();
    fSPAsymmetryMatchVec.clear();

    return;
}

void SpacePointAnalysis::fillHistograms(const art::Event& event) const
{
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
    
    // First task is to build a map between ides and voxel ids (that we calcualate based on position)
    // and also get the reverse since it will be useful in the end.
    // At the same time should also build a mapping of ides per channel so we can do quick hit lookup
    using IDEToVoxelIDMap   = std::unordered_map<const sim::IDE*, sim::LArVoxelID>;
    using VoxelIDToIDEMap   = std::map<sim::LArVoxelID, const sim::IDE*>;
    using TDCToIDEMap       = std::map<unsigned short, std::unordered_set<const sim::IDE*>>; // We need this one in order
    using ChanToTDCToIDEMap = std::map<raw::ChannelID_t, TDCToIDEMap>;
    using VoxelIDSet        = std::set<sim::LArVoxelID>;

    IDEToVoxelIDMap   ideToVoxelIDMap;
    VoxelIDToIDEMap   voxelIDToIDEMap;
    ChanToTDCToIDEMap chanToTDCToIDEMap;
    VoxelIDSet        simChannelVoxelIDSet;

    for(const auto& simChannel : *simChannelHandle)
    {
        for(const auto& tdcide : simChannel.TDCIDEMap())
        {
            double timeTDC  = 0.; // tdcide.first;  ==> Note this includes the drift time of the particle
            
            for(const auto& ide : tdcide.second) //chanToTDCToIDEMap[simChannel.Channel()][tdcide.first] = ide;
            {
                sim::LArVoxelID voxelID(ide.x,ide.y,ide.z,timeTDC);
                
                ideToVoxelIDMap[&ide]    = voxelID;
                voxelIDToIDEMap[voxelID] = &ide;
                chanToTDCToIDEMap[simChannel.Channel()][tdcide.first].insert(&ide);
                simChannelVoxelIDSet.insert(voxelID);
            }
        }
    }
    
    // Now we go throught the SimEnergyDeposit objects and try to make similar mappings
    // It is worth noting that in this case there can be multiple SimEnergyDeposit objects per voxel
    // We assume that calculating the voxel ID as above from the mean position of the SimEnergyDeposit objects will
    // result in the correct Voxel ID for relating to SimChannels (to be demonstrated!)
    using SimEnergyToVoxelIDMap    = std::unordered_map<const sim::SimEnergyDeposit*, sim::LArVoxelID>;
    using VoxelIDToSimEnergyVecMap = std::map<sim::LArVoxelID, std::vector<const sim::SimEnergyDeposit*>>;
    
    SimEnergyToVoxelIDMap    simEnergyToVoxelIDMap;
    VoxelIDToSimEnergyVecMap voxelIDToSimEnergyVecMap;
    VoxelIDSet               simEnergyVoxelIDSet;

    for(const auto& simEnergy : *simEnergyHandle)
    {
        geo::Point_t    midPoint = simEnergy.MidPoint();
        //double          timeTDC  = fClockService->TPCClock().Ticks(fClockService->G4ToElecTime(simEnergy.Time()));  // SimEnergy time is G4 time
        double          timeTDC  = 0.; // fClockService->TPCG4Time2Tick(simEnergy.Time()); ==> This does not include the drift time...
        sim::LArVoxelID voxelID(midPoint.X(),midPoint.Y(),midPoint.Z(),timeTDC);

        simEnergyToVoxelIDMap[&simEnergy] = voxelID;
        voxelIDToSimEnergyVecMap[voxelID].push_back(&simEnergy);
        simEnergyVoxelIDSet.insert(voxelID);
    }

    // Look at the difference, this will be the SimEnergyDeposit objects with no corresponding SimChannels...
    std::vector<sim::LArVoxelID>           energyDiffVec(simEnergyVoxelIDSet.size() + simChannelVoxelIDSet.size());
    std::vector<sim::LArVoxelID>::iterator energyDiffVecItr = std::set_difference(simEnergyVoxelIDSet.begin(),simEnergyVoxelIDSet.end(),simChannelVoxelIDSet.begin(),simChannelVoxelIDSet.end(),energyDiffVec.begin());
    
    if (energyDiffVecItr != energyDiffVec.begin())
    {
        energyDiffVec.resize(energyDiffVecItr - energyDiffVec.begin());
        
        std::cout << "==> SimEnergy size: " << simEnergyVoxelIDSet.size() << ", simChannel size: " << simChannelVoxelIDSet.size() << ", # in SimEnergy not in SimChannel: " << energyDiffVec.size() << std::endl;
        
//        // In this direction the diff vec contains voxel IDs in SimEnergyVoxelIDSet but not in SimChannelVoxelIDSet
//        for(const auto& voxelID : energyDiffVec)
//        {
//            std::cout << "  - voxel ID: " << voxelID;
//
//            VoxelIDToSimEnergyVecMap::iterator voxToEnergyItr = voxelIDToSimEnergyVecMap.find(voxelID);
//
//            if (voxToEnergyItr != voxelIDToSimEnergyVecMap.end())
//            {
//                for(const auto& simEnergy : voxToEnergyItr->second)
//                {
//                    const geo::Point_t midPoint = simEnergy->MidPoint();
//
//                    std::cout << ", SimEnergy pos: " << midPoint.X() << "/" << midPoint.Y() << "/" << midPoint.Z();
//                }
//            }
//
//            std::cout << std::endl;
//       }
    }
    else
        std::cout << "==> SimEnergy size: " << simEnergyVoxelIDSet.size() << ", simChannel size: " << simChannelVoxelIDSet.size() << ", No differences" << std::endl;
    
    
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    
    
    // Now go the other direction
    std::vector<sim::LArVoxelID>           simChanDiffVec(simEnergyVoxelIDSet.size() + simChannelVoxelIDSet.size());
    std::vector<sim::LArVoxelID>::iterator simChanDiffVecItr = std::set_difference(simChannelVoxelIDSet.begin(),simChannelVoxelIDSet.end(),simEnergyVoxelIDSet.begin(),simEnergyVoxelIDSet.end(),simChanDiffVec.begin());
    
    if (simChanDiffVecItr != simChanDiffVec.begin())
    {
        simChanDiffVec.resize(simChanDiffVecItr - simChanDiffVec.begin());
        
        std::cout << "==> SimEnergy size: " << simEnergyVoxelIDSet.size() << ", simChannel size: " << simChannelVoxelIDSet.size() << ", # in SimChannel not in SimEnergy: " << simChanDiffVec.size() << std::endl;
        
        // In this direction, diff vec contains voxel IDs in SimChannelVoxelIDSet but not in SimEnergyVoxelIDSet
        for(const auto& voxelID : simChanDiffVec)
        {
            VoxelIDSet::iterator lowerItr = std::lower_bound(simEnergyVoxelIDSet.begin(),simEnergyVoxelIDSet.end(),voxelID);
            VoxelIDSet::iterator upperItr = std::upper_bound(simEnergyVoxelIDSet.begin(),simEnergyVoxelIDSet.end(),voxelID);
            
            if (lowerItr != simEnergyVoxelIDSet.end() && upperItr != simEnergyVoxelIDSet.end())
            {
                TVector3 vecToLower = TVector3(*lowerItr) - TVector3(voxelID);
//                TVector3 vecToUpper = TVector3(*upperItr) - TVector3(voxelID);
                double   manDist    = std::abs(vecToLower.X()) + std::abs(vecToLower.Y()) + std::abs(vecToLower.Z());
                
                fSEDeltaX.push_back(vecToLower.X());
                fSEDeltaY.push_back(vecToLower.Y());
                fSEDeltaZ.push_back(vecToLower.Z());
                fSEDistance.push_back(manDist);
            }
        }
    }
    else
        std::cout << "==> SimEnergy size: " << simEnergyVoxelIDSet.size() << ", simChannel size: " << simChannelVoxelIDSet.size() << ", No differences" << std::endl;

    // Ok, for my next trick I want to build a mapping between hits and voxel IDs. Note that any given hit can be associated to more than one voxel...
    // We do this on the entire hit collection, ultimately we will want to consider SpacePoint efficiency (this could be done in the loop over SpacePoints
    // using the associated hits and would save time/memory)
    using RecobHitToVoxelIDMap = std::unordered_map<const recob::Hit*, std::set<sim::LArVoxelID>>;
    
    RecobHitToVoxelIDMap recobHitToVoxelIDMap;
    
    // And now fill it
    for(const auto& hitLabel : fRecobHitLabelVec)
    {
        art::Handle< std::vector<recob::Hit> > hitHandle;
        event.getByLabel(hitLabel, hitHandle);
        
        for(const auto& hit : *hitHandle)
        {
            // Recover channel information based on this hit
            ChanToTDCToIDEMap::const_iterator chanToTDCToIDEItr = chanToTDCToIDEMap.find(hit.Channel());
            
            if (chanToTDCToIDEItr != chanToTDCToIDEMap.end())
            {
                // Recover hit time range (in ticks)
                int startTick = hit.PeakTime() - 2. * hit.RMS();
                int endTick   = hit.PeakTime() + 2. * hit.RMS() + 1.;
                
                const TDCToIDEMap& tdcToIDEMap = chanToTDCToIDEItr->second;

                // Get the number of electrons
                for(unsigned short tick =startTick; tick <= endTick; tick++)
                {
                    unsigned short hitTDC = fClockService->TPCTick2TDC(tick - fOffsetVec[hit.WireID().Plane]);
                    
                    TDCToIDEMap::const_iterator ideIterator = tdcToIDEMap.find(hitTDC);
                    
                    if (ideIterator != tdcToIDEMap.end())
                    {
                        for (const auto& ide : ideIterator->second)
                        {
                            sim::LArVoxelID voxelID = ideToVoxelIDMap[ide];
                            
                            recobHitToVoxelIDMap[&hit].insert(voxelID);
                        }
                    }
                }
                
                // Why no match I wonder?
                if (recobHitToVoxelIDMap.find(&hit) == recobHitToVoxelIDMap.end())
                {
                    std::cout << "--> Hit C/T/P/W " << hit.WireID().Cryostat << "/" << hit.WireID().TPC << "/" << hit.WireID().Plane << "/" << hit.WireID().Wire << ", idx/Mult: " << hit.LocalIndex() << "/" << hit.Multiplicity() << ", PH: " << hit.PeakAmplitude() << " with " << endTick - startTick << " ticks: " << ", first/last IDE: " << tdcToIDEMap.begin()->first << "/" << tdcToIDEMap.rbegin()->first << ", start/end tick: " << fClockService->TPCTick2TDC(startTick) << "/" << fClockService->TPCTick2TDC(endTick) << " (" << startTick << "/" << endTick << ")" << std::endl;
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
                    std::cout << "I am certain this cannot happen... but here you go, space point with " << associatedHits.size() << " hits" << std::endl;
                    continue;
                }
                
                // Retrieve the magic numbers from the space point
                float spQuality   = spacePointPtr->Chisq();
                float spCharge    = spacePointPtr->ErrXYZ()[1];
                float spAsymmetry = spacePointPtr->ErrXYZ()[3];
                float smallestPH  = std::numeric_limits<float>::max();
                
                std::vector<int> numIDEsHitVec;
                int              numIDEsSpacePoint(0);

                std::vector<RecobHitToVoxelIDMap::const_iterator> recobHitToVoxelIterVec;
                
                // Now we can use our maps to find out if the hits making up the SpacePoint are truly related...
                for(const auto& hitPtr : associatedHits)
                {
                    RecobHitToVoxelIDMap::iterator hitToVoxelItr = recobHitToVoxelIDMap.find(hitPtr.get());
                    
                    if (hitPtr->PeakAmplitude() < smallestPH) smallestPH = hitPtr->PeakAmplitude();
                    
                    if (hitToVoxelItr == recobHitToVoxelIDMap.end())
                    {
                        numIDEsHitVec.push_back(0);
                        continue;
                    }
                    
                    recobHitToVoxelIterVec.push_back(hitToVoxelItr);
                    numIDEsHitVec.push_back(hitToVoxelItr->second.size());
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
                }
                
                // Fill for "all" cases
                fSPQualityVec.push_back(spQuality);
                fSPTotalChargeVec.push_back(spCharge);
                fSPAsymmetryVec.push_back(spAsymmetry);
                fSmallestPHVec.push_back(smallestPH);
                
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
