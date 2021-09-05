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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

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

        tree->Branch("SimHitDeltaT0",      "std::vector<int>",   &fSimHitDeltaT0Vec);
        tree->Branch("SimHitDeltaT1",      "std::vector<int>",   &fSimHitDeltaT1Vec);
        tree->Branch("SimHitDeltaT2",      "std::vector<int>",   &fSimHitDeltaT2Vec);
        tree->Branch("SimDeltaT10",        "std::vector<int>",   &fSimDelta10Vec);
        tree->Branch("SimDeltaT11",        "std::vector<int>",   &fSimDelta21Vec);
        tree->Branch("HitDeltaT10",        "std::vector<int>",   &fHitDelta10Vec);
        tree->Branch("HitDeltaT11",        "std::vector<int>",   &fHitDelta21Vec);
        tree->Branch("MaxElectronDep0",    "std::vector<float>", &fBigElecDep0Vec);
        tree->Branch("MaxElectronDep1",    "std::vector<float>", &fBigElecDep1Vec);
        tree->Branch("MaxElectronDep2",    "std::vector<float>", &fBigElecDep2Vec);

        fTree = tree;
    }

    void fill()
    {
        if (fTree) fTree->Fill();
    }

    void clear()
    {
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

        fSimHitDeltaT0Vec.clear();
        fSimHitDeltaT1Vec.clear();
        fSimHitDeltaT2Vec.clear();
        fSimDelta10Vec.clear();
        fSimDelta21Vec.clear();
        fHitDelta10Vec.clear();
        fHitDelta21Vec.clear();
        fBigElecDep0Vec.clear();
        fBigElecDep1Vec.clear();
        fBigElecDep2Vec.clear();
    }

    // Define tuple vars, make public for direct access
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

    std::vector<int>   fSimHitDeltaT0Vec;
    std::vector<int>   fSimHitDeltaT1Vec;
    std::vector<int>   fSimHitDeltaT2Vec;
    std::vector<int>   fSimDelta10Vec;
    std::vector<int>   fSimDelta21Vec;
    std::vector<int>   fHitDelta10Vec;
    std::vector<int>   fHitDelta21Vec;
    std::vector<float> fBigElecDep0Vec;
    std::vector<float> fBigElecDep1Vec;
    std::vector<float> fBigElecDep2Vec;

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
    std::vector<art::InputTag>  fSpacePointLabelVec;
    art::InputTag               fBadChannelProducerLabel;
    std::vector<int>            fOffsetVec;              ///< Allow offsets for each plane

    // TTree variables
    mutable TTree*              fTree;

    mutable std::vector<int>    fTPCVec;
    mutable std::vector<int>    fCryoVec;
    mutable std::vector<int>    fPlaneVec;

    mutable HitSpacePointObj    fHitSpacePointObj;

    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
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
    fSpacePointLabelVec = pset.get< std::vector<art::InputTag>>("SpacePointLabelVec", std::vector<art::InputTag>() = {"cluster3d"});
    fOffsetVec          = pset.get<std::vector<int>           >("OffsetVec",          std::vector<int>()           = {0,0,0}      );

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

    clear();

    return;
}

void SpacePointAnalysis::clear() const
{
    fTPCVec.clear();
    fCryoVec.clear();
    fPlaneVec.clear();

    fHitSpacePointObj.clear();

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

    return;
}

void SpacePointAnalysis::processSpacePoints(const art::Event&                  event,
                                            const detinfo::DetectorClocksData& clockData) const
{
    // Diagnostics
    using Triplet    = std::tuple<const recob::Hit*, const recob::Hit*, const recob::Hit*>;
    using TripletMap = std::map<Triplet, std::vector<const recob::SpacePoint*>>;

    TripletMap tripletMap;

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
            float                       spQuality       = spacePointPtr->Chisq();
            float                       spCharge        = spacePointPtr->ErrXYZ()[1];
            float                       spAsymmetry     = spacePointPtr->ErrXYZ()[3];
            float                       smallestPH      = std::numeric_limits<float>::max();
            float                       largestPH       = 0.;
            int                         numHits         = 0;
            float                       averagePH       = 0.;
            float                       averagePT       = 0.;
            float                       largestDelT     = 0.;
            std::vector<float>          hitPeakTimeVec  = {-100.,-100.,-100.};
            std::vector<float>          bigElecDepVec   = {0.,0.,0.};
            std::vector<unsigned short> bigTDCVec       = {0,0,0};
            int                         numLongHits(0);
            int                         numIntersections(0);

            std::vector<const recob::Hit*> recobHitVec(3,nullptr);

            // Now we can use our maps to find out if the hits making up the SpacePoint are truly related...
            for(const auto& hitPtr : associatedHits)
            {
                float peakAmplitude = hitPtr->PeakAmplitude();
                float peakTime      = hitPtr->PeakTime();
                int   plane         = hitPtr->WireID().Plane;

                recobHitVec[plane] = hitPtr.get();
                numHits++;
                averagePH += peakAmplitude;
                averagePT += peakTime;
                smallestPH = std::min(peakAmplitude,smallestPH);
                largestPH  = std::max(peakAmplitude,largestPH);

                if (hitPtr->DegreesOfFreedom() < 2) numLongHits++;

                hitPeakTimeVec[plane] = clockData.TPCTick2TDC(peakTime);
            }
            Triplet hitTriplet(recobHitVec[0],recobHitVec[1],recobHitVec[2]);

            tripletMap[hitTriplet].emplace_back(spacePointPtr.get());

            averagePH /= float(numHits);
            averagePT /= float(numHits);

            for(const auto& hitPtr : associatedHits)
            {
                float delT = hitPtr->PeakTime() - averagePT;
                if (std::abs(delT) > std::abs(largestDelT)) largestDelT = delT;
            }

            // Fill for "all" cases
            fHitSpacePointObj.fSPQualityVec.push_back(spQuality);
            fHitSpacePointObj.fSPTotalChargeVec.push_back(spCharge);
            fHitSpacePointObj.fSPAsymmetryVec.push_back(spAsymmetry);
            fHitSpacePointObj.fSmallestPHVec.push_back(smallestPH);
            fHitSpacePointObj.fLargestPHVec.push_back(largestPH);
            fHitSpacePointObj.fAveragePHVec.push_back(averagePH);
            fHitSpacePointObj.fLargestDelTVec.push_back(largestDelT);
            fHitSpacePointObj.fNumLongHitsVec.emplace_back(numLongHits);
            fHitSpacePointObj.fNumIntersectSetVec.emplace_back(numIntersections);
            fHitSpacePointObj.fSimHitDeltaT0Vec.emplace_back(bigTDCVec[0] - hitPeakTimeVec[0]);
            fHitSpacePointObj.fSimHitDeltaT1Vec.emplace_back(bigTDCVec[1] - hitPeakTimeVec[1]);
            fHitSpacePointObj.fSimHitDeltaT2Vec.emplace_back(bigTDCVec[2] - hitPeakTimeVec[2]);
            fHitSpacePointObj.fSimDelta10Vec.emplace_back(bigTDCVec[1] - bigTDCVec[0]);
            fHitSpacePointObj.fSimDelta21Vec.emplace_back(bigTDCVec[2] - bigTDCVec[1]);
            fHitSpacePointObj.fHitDelta10Vec.emplace_back(hitPeakTimeVec[1] - hitPeakTimeVec[0]);
            fHitSpacePointObj.fHitDelta21Vec.emplace_back(hitPeakTimeVec[2] - hitPeakTimeVec[1]);
            fHitSpacePointObj.fBigElecDep0Vec.emplace_back(bigElecDepVec[0]);
            fHitSpacePointObj.fBigElecDep1Vec.emplace_back(bigElecDepVec[1]);
            fHitSpacePointObj.fBigElecDep2Vec.emplace_back(bigElecDepVec[2]);
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
