
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
#include "nusimdata/SimulationBase/MCParticle.h"

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
    
    mutable std::vector<int>   fNumIDEsHit0Vec;
    mutable std::vector<int>   fNumIDEsHit1Vec;
    mutable std::vector<int>   fNumIDEsHit2Vec;
    mutable std::vector<int>   fNumIDEsSpacePointVec;
    
    mutable std::vector<float> fSPQualityAllVec;
    mutable std::vector<float> fSPTotalChargeAllVec;
    mutable std::vector<float> fSPAsymmetryAllVec;
    
    mutable std::vector<float> fSPQualityMatchVec;
    mutable std::vector<float> fSPTotalChargeMatchVec;
    mutable std::vector<float> fSPAsymmetryMatchVec;

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
    
    fTree->Branch("CryostataVec",       "std::vector<int>",   &fCryoVec);
    fTree->Branch("TPCVec",             "std::vector<int>",   &fTPCVec);
    fTree->Branch("PlaneVec",           "std::vector<int>",   &fPlaneVec);
 
    fTree->Branch("NumIDEsHit0",        "std::vector<int>",   &fNumIDEsHit0Vec);
    fTree->Branch("NumIDEsHit1",        "std::vector<int>",   &fNumIDEsHit1Vec);
    fTree->Branch("NumIDEsHit2",        "std::vector<int>",   &fNumIDEsHit2Vec);
    fTree->Branch("NumIDEsSpacePoint",  "std::vector<int>",   &fNumIDEsSpacePointVec);
    
    fTree->Branch("SPQualityAll",       "std::vector<float>", &fSPQualityAllVec);
    fTree->Branch("SPTotalChargeAll",   "std::vector<float>", &fSPTotalChargeAllVec);
    fTree->Branch("SPAsymmetryAll",     "std::vector<float>", &fSPAsymmetryAllVec);
    
    fTree->Branch("SPQualityMatch",     "std::vector<float>", &fSPQualityMatchVec);
    fTree->Branch("SPTotalChargeMatch", "std::vector<float>", &fSPTotalChargeMatchVec);
    fTree->Branch("SPAsymmetryMatch",   "std::vector<float>", &fSPAsymmetryMatchVec);

    clear();

    return;
}
    
void SpacePointAnalysis::clear() const
{
    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();
    
    fNumIDEsHit0Vec.clear();
    fNumIDEsHit1Vec.clear();
    fNumIDEsHit2Vec.clear();
    fNumIDEsSpacePointVec.clear();
    
    fSPQualityAllVec.clear();
    fSPTotalChargeAllVec.clear();
    fSPAsymmetryAllVec.clear();
    
    fSPQualityMatchVec.clear();
    fSPTotalChargeMatchVec.clear();
    fSPAsymmetryMatchVec.clear();

    return;
}

void SpacePointAnalysis::fillHistograms(const art::Event& event) const
{
    // Always clear the tuple
    clear();
    
    art::Handle< std::vector<sim::SimChannel>> simChannelHandle;
    event.getByLabel(fSimChannelProducerLabel, simChannelHandle);
    
    art::Handle< std::vector<simb::MCParticle>> mcParticleHandle;
    event.getByLabel(fMCParticleProducerLabel, mcParticleHandle);

    // If there is no sim channel informaton then exit
    if (!simChannelHandle.isValid() || simChannelHandle->empty() || !mcParticleHandle.isValid()) return;
    
    // First task is to build a map between ides and voxel ids (that we calcualate based on position)
    // and also get the reverse since it will be useful in the end.
    // At the same time should also build a mapping of ides per channel so we can do quick hit lookup
    using IDEToVoxelIDMap   = std::unordered_map<const sim::IDE*, size_t>;
    using VoxelIDToIDEMap   = std::unordered_map<size_t, const sim::IDE*>;
    using TDCToIDEMap       = std::map<unsigned short, const sim::IDE*>; // We need this one in order
    using ChanToTDCToIDEMap = std::map<raw::ChannelID_t, TDCToIDEMap>;

    IDEToVoxelIDMap   ideToVoxelIDMap;
    VoxelIDToIDEMap   voxelIDToIDEMap;
    ChanToTDCToIDEMap chanToTDCToIDEMap;

    for(const auto& simChannel : *simChannelHandle)
    {
        for(const auto& tdcide : simChannel.TDCIDEMap())
        {
            for(const auto& ide : tdcide.second) //chanToTDCToIDEMap[simChannel.Channel()][tdcide.first] = ide;
            {
                size_t voxelID =                        int(fPositionToVoxelIDFactor * (ide.x - fVolumeOffsetsLow[0]) + 1.)
                               + fNumVoxelsByAxis[0] * (int(fPositionToVoxelIDFactor * (ide.y - fVolumeOffsetsLow[1]) + 1.)
                               + fNumVoxelsByAxis[1] * (int(fPositionToVoxelIDFactor * (ide.z - fVolumeOffsetsLow[2]) + 1.)));
                
                ideToVoxelIDMap[&ide]                                 = voxelID;
                voxelIDToIDEMap[voxelID]                              = &ide;
                chanToTDCToIDEMap[simChannel.Channel()][tdcide.first] = &ide;
            }
        }
    }
    
    // Ok, for my next trick I want to build a mapping between hits and voxel IDs. Note that any given hit can be associated to more than one voxel...
    // We do this on the entire hit collection, ultimately we will want to consider SpacePoint efficiency (this could be done in the loop over SpacePoints
    // using the associated hits and would save time/memory)
    using RecobHitToVoxelIDMap = std::unordered_map<const recob::Hit*, std::vector<size_t>>;
    
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
                        const sim::IDE* ide = ideIterator->second;
                        
                        recobHitToVoxelIDMap[&hit].push_back(ideToVoxelIDMap[ide]);
                    }
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

            // And now, without further ado, here we beging the loop that will actually produce some useful output. 
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
                float spQualityAll   = spacePointPtr->Chisq();
                float spChargeAll    = spacePointPtr->ErrXYZ()[1];
                float spAsymmetryAll = spacePointPtr->ErrXYZ()[3];
                
                float spQualityMatch(-1.);
                float spChargeMatch(-1.);
                float spAsymmetryMatch(-2.);
                
                std::vector<int> numIDEsHitVec;
                int              numIDEsSpacePoint(0);

                std::vector<RecobHitToVoxelIDMap::const_iterator> recobHitToVoxelIterVec;
                
                // Now we can use our maps to find out if the hits making up the SpacePoint are truly related...
                for(const auto& hitPtr : associatedHits)
                {
                    RecobHitToVoxelIDMap::iterator hitToVoxelItr = recobHitToVoxelIDMap.find(hitPtr.get());
                    
                    if (hitToVoxelItr == recobHitToVoxelIDMap.end())
                    {
                        numIDEsHitVec.push_back(0);
                        continue;
                    }
                    
                    // Need the list sorted in order to use the set intersection method
                    if (!hitToVoxelItr->second.empty()) std::sort(hitToVoxelItr->second.begin(),hitToVoxelItr->second.end());
                    
                    recobHitToVoxelIterVec.push_back(hitToVoxelItr);
                    numIDEsHitVec.push_back(hitToVoxelItr->second.size());
                }
                
                if (recobHitToVoxelIterVec.size() == 3)
                {
                    std::vector<size_t> firstIntersectionVec(recobHitToVoxelIterVec[0]->second.size(),recobHitToVoxelIterVec[1]->second.size());
                    
                    std::vector<size_t>::iterator firstIntersectionItr = std::set_intersection(recobHitToVoxelIterVec[0]->second.begin(),recobHitToVoxelIterVec[0]->second.end(),
                                                                                               recobHitToVoxelIterVec[1]->second.begin(),recobHitToVoxelIterVec[1]->second.end(),
                                                                                               firstIntersectionVec.begin());
                    
                    firstIntersectionVec.resize(firstIntersectionItr - firstIntersectionVec.begin());
                    
                    if (!firstIntersectionVec.empty())
                    {
                        std::vector<size_t> secondIntersectionVec(firstIntersectionVec.size(),recobHitToVoxelIterVec[2]->second.size());
                        
                        std::vector<size_t>::iterator secondIntersectionItr = std::set_intersection(firstIntersectionVec.begin(),             firstIntersectionVec.end(),
                                                                                                    recobHitToVoxelIterVec[1]->second.begin(),recobHitToVoxelIterVec[1]->second.end(),
                                                                                                    secondIntersectionVec.begin());
                        
                        secondIntersectionVec.resize(secondIntersectionItr - secondIntersectionVec.begin());
                        
                        if (!secondIntersectionVec.empty())
                        {
                            numIDEsSpacePoint = secondIntersectionVec.size();
                            spQualityMatch    = spQualityAll;
                            spChargeMatch     = spChargeAll;
                            spAsymmetryMatch  = spAsymmetryAll;
                        }
                    }
                }
                
                
                fSPQualityAllVec.push_back(spQualityAll);
                fSPTotalChargeAllVec.push_back(spChargeAll);
                fSPAsymmetryAllVec.push_back(spAsymmetryAll);
                
                fNumIDEsHit0Vec.push_back(numIDEsHitVec[0]);
                fNumIDEsHit1Vec.push_back(numIDEsHitVec[1]);
                fNumIDEsHit2Vec.push_back(numIDEsHitVec[2]);
                fNumIDEsSpacePointVec.push_back(numIDEsSpacePoint);
                
                fSPQualityMatchVec.push_back(spQualityMatch);
                fSPTotalChargeMatchVec.push_back(spChargeMatch);
                fSPAsymmetryMatchVec.push_back(spAsymmetryMatch);
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
