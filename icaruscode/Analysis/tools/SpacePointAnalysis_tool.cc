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
#include "canvas/Persistency/Common/FindOneP.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"

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

// Define object to keep track of hit related tuple items
class HitTupleObj
{
public:
    HitTupleObj() : fTree(nullptr) {}

    void setBranches(TTree* tree)
    {
        tree->Branch("Cryostat",          "std::vector<int>",   &fCryostatVec);
        tree->Branch("TPC",               "std::vector<int>",   &fTPCVec);
        tree->Branch("TicksTotHit",       "std::vector<int>",   &fTicksTotHitVec);
        tree->Branch("Tick",              "std::vector<int>",   &fTicksVec);
        tree->Branch("NDF",               "std::vector<int>",   &fNDFHitVec);                //< Number of degrees of freedom of hit fit
        tree->Branch("Multiplicity",      "std::vector<int>",   &fMultiplicityHitVec);       //< Multiplicity of the snippet the hit is on
        tree->Branch("LocalIndex",        "std::vector<int>",   &fLocalIndexHitVec);         //< The index of the hit within the snippet
        tree->Branch("ChiSquare",         "std::vector<float>", &fChiSquareHitVec);          //< Chi square of fit
        tree->Branch("SummedADC",         "std::vector<float>", &fSummedADCHitVec);          //< Sum of all ADC values start/end of snippet
        tree->Branch("Integral",          "std::vector<float>", &fIntegralHitVec);           //< Integrated charge +/- n sigma about peak center
        tree->Branch("PulseHeight",       "std::vector<float>", &fPHHitVec);                 //< Pulse height of hit
        tree->Branch("RMS",               "std::vector<float>", &fRMSHitVec);                //< RMS of hit (from fit)
        tree->Branch("ClusterSize",       "std::vector<int>",   &fClusterSizeVec);           //< Number space points in cluster for this hit

        fTree = tree;

        return;
    }

    void fill()
    {
        if (fTree) fTree->Fill();
    }

    void clear()
    {
        fCryostatVec.clear();
        fTPCVec.clear();
        fTicksTotHitVec.clear();
        fTicksVec.clear();
        fNDFHitVec.clear();
        fMultiplicityHitVec.clear();
        fLocalIndexHitVec.clear();
        fChiSquareHitVec.clear();
        fSummedADCHitVec.clear();
        fIntegralHitVec.clear();
        fPHHitVec.clear();
        fRMSHitVec.clear();
        fPHOrderHitVec.clear();
        fClusterSizeVec.clear();
    }

    void fillHitInfo(const recob::Hit* hit, int hitWidth, int clusterSize)
    {
        fCryostatVec.emplace_back(hit->WireID().Cryostat);
        fTPCVec.emplace_back(hit->WireID().TPC);
        fTicksTotHitVec.emplace_back(hitWidth);
        fTicksVec.emplace_back(hit->PeakTime());
        fNDFHitVec.emplace_back(hit->DegreesOfFreedom());
        fMultiplicityHitVec.emplace_back(hit->Multiplicity());
        fLocalIndexHitVec.emplace_back(hit->LocalIndex());
        fChiSquareHitVec.emplace_back(hit->GoodnessOfFit());
        fSummedADCHitVec.emplace_back(hit->SummedADC());
        fIntegralHitVec.emplace_back(hit->Integral());
        fPHHitVec.emplace_back(hit->PeakAmplitude());
        fRMSHitVec.emplace_back(hit->RMS());
        fClusterSizeVec.emplace_back(clusterSize);
    }

    // Define tuple values, these are public so can be diretly accessed for filling
    std::vector<int>   fCryostatVec;
    std::vector<int>   fTPCVec;
    std::vector<int>   fTicksTotHitVec;
    std::vector<int>   fTicksVec;
    std::vector<int>   fNDFHitVec;
    std::vector<int>   fMultiplicityHitVec;
    std::vector<int>   fLocalIndexHitVec;
    std::vector<float> fChiSquareHitVec;
    std::vector<float> fSummedADCHitVec;
    std::vector<float> fIntegralHitVec;
    std::vector<float> fPHHitVec;
    std::vector<float> fRMSHitVec;

    std::vector<int>   fPHOrderHitVec;
    std::vector<int>   fClusterSizeVec;

private:
    TTree* fTree;
};

// Define object to keep track of hit/spacepoint related items
class HitSpacePointObj
{
public:
    HitSpacePointObj() : fTree(nullptr) {}

    void setBranches(TTree* tree)
    {
        tree->Branch("SPCryostat",         "std::vector<int>",   &fSPCryostatVec);
        tree->Branch("SPTPC",              "std::vector<int>",   &fSPTPCVec);
        tree->Branch("SPQuality",          "std::vector<float>", &fSPQualityVec);
        tree->Branch("SPTotalCharge",      "std::vector<float>", &fSPTotalChargeVec);
        tree->Branch("SPAsymmetry",        "std::vector<float>", &fSPAsymmetryVec);
        tree->Branch("SmallestPH",         "std::vector<float>", &fSmallestPHVec);
        tree->Branch("LargestPH",          "std::vector<float>", &fLargestPHVec);
        tree->Branch("AveragePH",          "std::vector<float>", &fAveragePHVec);
        tree->Branch("LargestDelT",        "std::vector<float>", &fLargestDelTVec);
        tree->Branch("SmallestDelT",       "std::vector<float>", &fSmallestDelTVec);

        tree->Branch("SP_x",               "std::vector<float>", &fSP_x);
        tree->Branch("SP_y",               "std::vector<float>", &fSP_y);
        tree->Branch("SP_z",               "std::vector<float>", &fSP_z);

        tree->Branch("Num2DHits",          "std::vector<int>",   &fNum2DHitsVec);
        tree->Branch("HitPlane",           "std::vector<int>",   &fHitPlaneVec);
        tree->Branch("NumLongHitsSP",      "std::vector<int>",   &fNumLongHitsVec);
        tree->Branch("NumIntersectSet",    "std::vector<int>",   &fNumIntersectSetVec);
        tree->Branch("ClusterNSP",         "std::vector<int>",   &fClusterNSPVec);

        tree->Branch("HitDeltaT10",        "std::vector<float>", &fHitDelta10Vec);
        tree->Branch("HitSigmaT10",        "std::vector<float>", &fHitSigma10Vec);
        tree->Branch("HitDeltaT21",        "std::vector<float>", &fHitDelta21Vec);
        tree->Branch("HitSigmaT21",        "std::vector<float>", &fHitSigma21Vec);
        tree->Branch("HitDeltaT20",        "std::vector<float>", &fHitDelta20Vec);
        tree->Branch("HitSigmaT20",        "std::vector<float>", &fHitSigma20Vec);
        tree->Branch("HitMultProduct",     "std::vector<int>",   &fHitMultProductVec);

        fTree = tree;
    }

    void fill()
    {
        if (fTree) fTree->Fill();
    }

    void clear()
    {
        fSPCryostatVec.clear();
        fSPTPCVec.clear();

        fSPQualityVec.clear();
        fSPTotalChargeVec.clear();
        fSPAsymmetryVec.clear();
        fSmallestPHVec.clear();
        fLargestPHVec.clear();
        fAveragePHVec.clear();
        fLargestDelTVec.clear();
        fSmallestDelTVec.clear();

        fSP_x.clear();
        fSP_y.clear();
        fSP_z.clear();

        fNum2DHitsVec.clear();
        fHitPlaneVec.clear();
        fNumLongHitsVec.clear();
        fNumIntersectSetVec.clear();
        fClusterNSPVec.clear();

        fHitDelta10Vec.clear();
        fHitSigma10Vec.clear();
        fHitDelta21Vec.clear();
        fHitSigma21Vec.clear();
        fHitDelta20Vec.clear();
        fHitSigma20Vec.clear();
        fHitMultProductVec.clear();
    }

    // Define tuple vars, make public for direct access
    std::vector<int>   fSPCryostatVec;
    std::vector<int>   fSPTPCVec;

    std::vector<float> fSPQualityVec;
    std::vector<float> fSPTotalChargeVec;
    std::vector<float> fSPAsymmetryVec;
    std::vector<float> fSmallestPHVec;
    std::vector<float> fLargestPHVec;
    std::vector<float> fAveragePHVec;
    std::vector<float> fLargestDelTVec;
    std::vector<float> fSmallestDelTVec;

    std::vector<float> fSP_x;
    std::vector<float> fSP_y;
    std::vector<float> fSP_z;

    std::vector<int>   fNum2DHitsVec;
    std::vector<int>   fHitPlaneVec;
    std::vector<int>   fNumLongHitsVec;
    std::vector<int>   fNumIntersectSetVec;
    std::vector<int>   fClusterNSPVec;

    std::vector<float> fHitDelta10Vec;
    std::vector<float> fHitDelta21Vec;
    std::vector<float> fHitDelta20Vec;
    std::vector<float> fHitSigma10Vec;
    std::vector<float> fHitSigma21Vec;
    std::vector<float> fHitSigma20Vec;
    std::vector<int>   fHitMultProductVec;

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

    void processSpacePoints(const art::Event&, const detinfo::DetectorClocksData&) const;

    // Relate hits to voxels
    using HitPointerVec        = std::vector<const recob::Hit*>;

    // Fcl parameters.
    using PlaneToT0OffsetMap = std::map<geo::PlaneID,float>;

    std::vector<art::InputTag>  fSpacePointLabelVec;
    art::InputTag               fBadChannelProducerLabel;
    bool                        fUseT0Offsets;

    // TTree variables
    mutable TTree*              fTree;

    mutable std::vector<int>    fTPCVec;
    mutable std::vector<int>    fCryoVec;
    mutable std::vector<int>    fPlaneVec;

    mutable PlaneToT0OffsetMap  fPlaneToT0OffsetMap;

    using HitTuplebjVec = std::vector<HitTupleObj>;

    mutable HitTuplebjVec       fHitTupleObjVec;
    mutable HitSpacePointObj    fHitSpacePointObj;

    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*    fGeometry;             ///< pointer to Geometry service
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
    fGeometry = lar::providerFrom<geo::Geometry>();

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
    fSpacePointLabelVec  = pset.get< std::vector<art::InputTag>>("SpacePointLabelVec",  {"cluster3d"});
    fUseT0Offsets        = pset.get< bool                      >("UseT0Offsets",        false);

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

    fHitTupleObjVec.resize(fGeometry->Nplanes());

    for(size_t plane = 0; plane < fGeometry->Nplanes(); plane++)
    {
        // Set up specific branch for space points
        locTree = tfs->makeAndRegister<TTree>("MatchedHits_P"+std::to_string(plane),"Matched Hits Tuple plane "+std::to_string(plane));

        fHitTupleObjVec[plane].setBranches(locTree);
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

    for(auto& hitObj : fHitTupleObjVec) hitObj.clear();

    return;
}

void SpacePointAnalysis::fillHistograms(const art::Event& event) const
{
    // Ok... this is starting to grow too much and get out of control... we will need to break it up directly...

    // Always clear the tuple
    clear();

    mf::LogDebug("SpacePointAnalysis") << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    // Now do the space points
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    processSpacePoints(event, clockData);

    // Make sure the output tuples are filled
    fHitSpacePointObj.fill();

    for(auto& hitObj : fHitTupleObjVec) hitObj.fill();

    return;
}

void SpacePointAnalysis::processSpacePoints(const art::Event&                  event,
                                            const detinfo::DetectorClocksData& clockData) const
{
    // One time initialization done?

    // Do the one time initialization of the tick offsets. 
    if (fPlaneToT0OffsetMap.empty())
    {
        // Need the detector properties which needs the clocks 
        auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

        for(size_t cryoIdx = 0; cryoIdx < fGeometry->Ncryostats(); cryoIdx++)
        {
            for(size_t tpcIdx = 0; tpcIdx < fGeometry->NTPC(); tpcIdx++)
            {
                for(size_t planeIdx = 0; planeIdx < fGeometry->Nplanes(); planeIdx++)
                {
                    geo::PlaneID planeID(cryoIdx,tpcIdx,planeIdx);

//                    if (fUseT0Offsets) fPlaneToT0OffsetMap[planeID] = int(det_prop.GetXTicksOffset(planeID) - det_prop.GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,0)));
                    if (fUseT0Offsets) fPlaneToT0OffsetMap[planeID] = det_prop.GetXTicksOffset(planeID) - det_prop.GetXTicksOffset(geo::PlaneID(cryoIdx,tpcIdx,0));
                    else               fPlaneToT0OffsetMap[planeID] = 0.;

                    std::cout << "--PlaneID: " << planeID << ", has T0 offset: " << fPlaneToT0OffsetMap.find(planeID)->second << std::endl;
                }
            }
        }   
    }

    // Diagnostics
    using Triplet    = std::tuple<const recob::Hit*, const recob::Hit*, const recob::Hit*>;
    using TripletMap = std::map<Triplet, std::vector<const recob::SpacePoint*>>;

    TripletMap tripletMap;

    // detector properties
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();

    // So now we loop through the various SpacePoint sources
    for(const auto& collectionLabel : fSpacePointLabelVec)
    {
        art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
        event.getByLabel(collectionLabel, pfParticleHandle);

        if (!pfParticleHandle.isValid()) continue;

        art::Handle< std::vector<recob::SpacePoint>> spacePointHandle;
        event.getByLabel(collectionLabel, spacePointHandle);

        if (!spacePointHandle.isValid()) continue;

        // Recover the collection of associations between tracks and hits
        art::FindManyP<recob::SpacePoint> pfPartSpacePointAssns(pfParticleHandle, event, collectionLabel);

        // Look up assocations between pfparticles and space points
        art::FindManyP<recob::Hit> spHitAssnVec(spacePointHandle, event, collectionLabel);

        // Use this to build a map between PFParticles and the number of associated space points 
        using PFParticleToNumSPMap = std::map<const recob::PFParticle*,int>;

        PFParticleToNumSPMap pfParticleToNumSPMap;

        for(size_t idx = 0; idx < pfParticleHandle->size(); idx++)
        {
            art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,idx);

            std::vector<art::Ptr<recob::SpacePoint>> spacePointVec(pfPartSpacePointAssns.at(pfParticle.key()));

            pfParticleToNumSPMap[pfParticle.get()] = spacePointVec.size();

//            if (spacePointVec.size() > 50000) std::cout << "SpacePointAnalysis finds PFParticle with " << spacePointVec.size() << " associated space points, run/event: " << event.id() << std::endl;
        }

        // Ok now we want the reverse look up
        art::FindManyP<recob::PFParticle> spacePointPFPartAssns(spacePointHandle, event, collectionLabel);

        std::unordered_map<const recob::Hit*,int> uniqueHitMap;

        // And now, without further adieu, here we begin the loop that will actually produce some useful output.
        // Loop over all space points and find out their true nature
        for(size_t idx = 0; idx < spacePointHandle->size(); idx++)
        {
            // Recover space point
            art::Ptr<recob::SpacePoint> spacePointPtr(spacePointHandle,idx);

            std::vector<art::Ptr<recob::Hit>> associatedHits(spHitAssnVec.at(spacePointPtr.key()));

            if (associatedHits.size() < 2)
            {
                mf::LogDebug("SpacePointAnalysis") << "I am certain this cannot happen... but here you go, space point with " << associatedHits.size() << " hits" << std::endl;
                continue;
            }

            // Recover the PFParticle associated to this space point and get the number of associated hits
            int nSpacePointsInPFParticle(0);

            std::vector<art::Ptr<recob::PFParticle>> pfParticleVec(spacePointPFPartAssns.at(spacePointPtr.key()));

            if (pfParticleVec.size() == 1) nSpacePointsInPFParticle = pfParticleToNumSPMap.at(pfParticleVec[0].get());

            mf::LogDebug("SpacePointAnalysis") << "==> pfPartVec size: " << pfParticleVec.size() << ", # space points: " << nSpacePointsInPFParticle << std::endl;

            // Retrieve the magic numbers from the space point
            float              spQuality       = spacePointPtr->Chisq();
            float              spCharge        = spacePointPtr->ErrXYZ()[1];
            float              spAsymmetry     = spacePointPtr->ErrXYZ()[3];
            const Double32_t*  spPosition      = spacePointPtr->XYZ();
            float              smallestPH      = std::numeric_limits<float>::max();
            float              largestPH       = 0.;
            int                numHits         = 0;
            int                hitPlane        = 0;
            float              averagePH       = 0.;
            float              averagePT       = 0.;
            float              largestDelT     = 0.;
            float              smallestDelT    = 100000.;
            std::vector<float> hitPeakTimeVec  = {-10000.,-20000.,-30000.};
            std::vector<float> hitPeakRMSVec   = {1000.,1000.,1000.};
            int                hitMultProduct  = 1;
            int                numLongHits(0);
            int                numIntersections(0);
            int                cryostat(-1);
            int                tpc(-1);

            std::vector<const recob::Hit*> recobHitVec(3,nullptr);

            std::vector<float> peakAmpVec;

            // Now we can use our maps to find out if the hits making up the SpacePoint are truly related...
            for(const auto& hitPtr : associatedHits)
            {
                float peakAmplitude = hitPtr->PeakAmplitude();
                float peakTime      = hitPtr->PeakTime();
                float rms           = hitPtr->RMS();
                int   plane         = hitPtr->WireID().Plane;

                peakAmpVec.emplace_back(peakAmplitude);

                // Add to the set
                uniqueHitMap[hitPtr.get()] = nSpacePointsInPFParticle;

                recobHitVec[plane] = hitPtr.get();
                numHits++;
                hitPlane   += 1 << plane;
                averagePH  += peakAmplitude;
                smallestPH  = std::min(peakAmplitude,smallestPH);
                largestPH   = std::max(peakAmplitude,largestPH);

                hitMultProduct *= hitPtr->Multiplicity();

                if (hitPtr->DegreesOfFreedom() < 2) numLongHits++;

                hitPeakTimeVec[plane] = peakTime - fPlaneToT0OffsetMap.find(hitPtr->WireID().planeID())->second;
                hitPeakRMSVec[plane]  = rms;
                averagePT            += hitPeakTimeVec[plane];

                cryostat = hitPtr->WireID().Cryostat;
                tpc      = hitPtr->WireID().TPC;
            }

            Triplet hitTriplet(recobHitVec[0],recobHitVec[1],recobHitVec[2]);

            tripletMap[hitTriplet].emplace_back(spacePointPtr.get());

            averagePH /= float(numHits);
            averagePT /= float(numHits);

            for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
            {
                // Skip if hit missing
                if (hitPeakTimeVec[planeIdx] < 0) continue;

                float delT = hitPeakTimeVec[planeIdx] - averagePT;
                if (std::abs(delT) > std::abs(largestDelT))  largestDelT  = delT;
                if (std::abs(delT) < std::abs(smallestDelT)) smallestDelT = delT;
            }

            // Fill for "all" cases
            fHitSpacePointObj.fSPCryostatVec.emplace_back(cryostat);
            fHitSpacePointObj.fSPTPCVec.emplace_back(tpc);
            fHitSpacePointObj.fSPQualityVec.emplace_back(spQuality);
            fHitSpacePointObj.fSPTotalChargeVec.emplace_back(spCharge);
            fHitSpacePointObj.fSPAsymmetryVec.emplace_back(spAsymmetry);
            fHitSpacePointObj.fSmallestPHVec.emplace_back(smallestPH);
            fHitSpacePointObj.fLargestPHVec.emplace_back(largestPH);
            fHitSpacePointObj.fAveragePHVec.emplace_back(averagePH);
            fHitSpacePointObj.fLargestDelTVec.emplace_back(largestDelT);
            fHitSpacePointObj.fSmallestDelTVec.emplace_back(smallestDelT);
            fHitSpacePointObj.fNum2DHitsVec.emplace_back(numHits);
            fHitSpacePointObj.fHitPlaneVec.emplace_back(hitPlane);
            fHitSpacePointObj.fNumLongHitsVec.emplace_back(numLongHits);
            fHitSpacePointObj.fNumIntersectSetVec.emplace_back(numIntersections);
            fHitSpacePointObj.fClusterNSPVec.emplace_back(nSpacePointsInPFParticle);
            fHitSpacePointObj.fHitDelta10Vec.emplace_back(hitPeakTimeVec[1] - hitPeakTimeVec[0]);
            fHitSpacePointObj.fHitSigma10Vec.emplace_back(std::sqrt(std::pow(hitPeakRMSVec[1],2) + std::pow(hitPeakRMSVec[0],2)));
            fHitSpacePointObj.fHitDelta21Vec.emplace_back(hitPeakTimeVec[2] - hitPeakTimeVec[1]);
            fHitSpacePointObj.fHitSigma21Vec.emplace_back(std::sqrt(std::pow(hitPeakRMSVec[2],2) + std::pow(hitPeakRMSVec[1],2)));
            fHitSpacePointObj.fHitDelta20Vec.emplace_back(hitPeakTimeVec[2] - hitPeakTimeVec[0]);
            fHitSpacePointObj.fHitSigma20Vec.emplace_back(std::sqrt(std::pow(hitPeakRMSVec[2],2) + std::pow(hitPeakRMSVec[0],2)));
            fHitSpacePointObj.fHitMultProductVec.emplace_back(hitMultProduct);

            fHitSpacePointObj.fSP_x.emplace_back(spPosition[0]);
            fHitSpacePointObj.fSP_y.emplace_back(spPosition[1]);
            fHitSpacePointObj.fSP_z.emplace_back(spPosition[2]);
        }

        // Now include hit information for unique hits
        for(const auto& hitItr : uniqueHitMap)
        {
            // Recover hit time range (in ticks), cast a wide net here
            const recob::Hit* hit = hitItr.first;

            float peakTime  = hit->PeakTime();
            float rms       = hit->PeakTime();
            int   startTick = std::max(   0,int(std::floor(peakTime - 3. * rms)));
            int   endTick   = std::min(4096,int(std::ceil( peakTime + 3. * rms)));

            fHitTupleObjVec[hit->WireID().Plane].fillHitInfo(hit,endTick-startTick+1,hitItr.second);

        }
    }
    // Can we check to see if we have duplicates?
    std::vector<int> numSpacePointVec = {0,0,0,0,0};
    for(const auto& tripletPair : tripletMap)
    {
        int numSpacePoints = std::min(numSpacePointVec.size()-1,tripletPair.second.size());
        numSpacePointVec[numSpacePoints]++;
    }
    std::cout << "====>> Found " << tripletMap.size() << " SpacePoints, numbers: ";
    for(const auto& count : numSpacePointVec) std::cout << count << " ";
    std::cout << std::endl;

    return;
}


// Useful for normalizing histograms
void SpacePointAnalysis::endJob(int numEvents)
{
    return;
}

DEFINE_ART_CLASS_TOOL(SpacePointAnalysis)
}
