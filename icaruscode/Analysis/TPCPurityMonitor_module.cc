// TPCPurityMonitor_module.cc
// A basic "skeleton" to read in art::Event records from a file,
// access their information, and do something with them.

// See
// <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
// for a description of the ART classes used here.

// Almost everything you see in the code below may have to be changed
// by you to suit your task. The example task is to make histograms
// and n-tuples related to dE/dx of particle tracks in the detector.

// As you try to understand why things are done a certain way in this
// example ("What's all this stuff about 'auto const&'?"), it will help
// to read ADDITIONAL_NOTES.txt in the same directory as this file.

#ifndef TPCPurityMonitor_module
#define TPCPurityMonitor_module

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"

#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// One really does not like root but one has to use it
#include "TTree.h"

//purity info class
#include "sbnobj/Common/Analysis/TPCPurityInfo.hh"

// Eigen includes
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/Jacobi"

// C++ Includes
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

namespace TPCPurityMonitor
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class TPCPurityMonitor : public art::EDProducer
{
public:

    // Standard constructor and destructor for an ART module.
    explicit TPCPurityMonitor(fhicl::ParameterSet const& pset);
    virtual ~TPCPurityMonitor();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();
    void endJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
   // void beginRun(const art::Run& run);

    // The analysis routine, called once per event.
    void produce(art::Event& evt);

private:
    // Define a data structure to keep track of a hit's meta data
    using HitMetaPair            = std::pair<art::Ptr<recob::Hit>,const recob::TrackHitMeta*>;
    using HitMetaPairVec         = std::vector<HitMetaPair>;

    // Define the basic data structure we will use
    using StatusChargePair       = std::pair<bool,double>;
    using HitStatusChargePair    = std::pair<HitMetaPair,StatusChargePair>;
    using HitStatusChargePairVec = std::vector<HitStatusChargePair>;

    // We would also like to keep tracy of the trajectory points along the track
    using PointDirPair           = std::pair<geo::Point_t,geo::Vector_t>;
    using HitPointDirPairMap     = std::unordered_map<const recob::Hit*,PointDirPair>;

    // We also need to define a container for the output of the 2D PCA Analysis
    class PrincipalComponents2D
    {
    public:

        using EigenValues  = Eigen::Vector2d;
        using EigenVectors = Eigen::Matrix2d;

        PrincipalComponents2D() :
            fSVD_OK(false), fNumHitsUsed(0), fEigenValues(EigenValues::Zero()), fEigenVectors(EigenVectors::Zero()), fAvePosition(Eigen::Vector2d::Zero()) {}

        PrincipalComponents2D(bool ok, int nHits, const EigenValues& eigenValues, const EigenVectors& eigenVecs, const Eigen::Vector2d& avePos) :
            fSVD_OK(ok), fNumHitsUsed(nHits), fEigenValues(eigenValues), fEigenVectors(eigenVecs), fAvePosition(avePos) {}

        bool                   getSvdOK()                 const {return fSVD_OK;}
        int                    getNumHitsUsed()           const {return fNumHitsUsed;}
        const EigenValues&     getEigenValues()           const {return fEigenValues;}
        const EigenVectors&    getEigenVectors()          const {return fEigenVectors;}
        const Eigen::Vector2d& getAvePosition()           const {return fAvePosition;}
        void                   flipAxis(size_t axis)            { fEigenVectors.row(axis) = -fEigenVectors.row(axis);}

    private:

        bool            fSVD_OK;             ///< SVD Decomposition was successful
        int             fNumHitsUsed;        ///< Number of hits in the decomposition
        EigenValues     fEigenValues;        ///< Eigen values from SVD decomposition
        EigenVectors    fEigenVectors;       ///< The three principle axes
        Eigen::Vector2d fAvePosition;        ///< Average position of hits fed to PCA
    };

    // We also need to define a container for the output of the 2D PCA Analysis
    class PrincipalComponents3D
    {
    public:

        using EigenValues  = Eigen::Vector3d;
        using EigenVectors = Eigen::Matrix3d;

        PrincipalComponents3D() :
            fSVD_OK(false), fNumHitsUsed(0), fEigenValues(EigenValues::Zero()), fEigenVectors(EigenVectors::Zero()), fAvePosition(Eigen::Vector3d::Zero()) {}

        PrincipalComponents3D(bool ok, int nHits, const EigenValues& eigenValues, const EigenVectors& eigenVecs, const Eigen::Vector3d& avePos) :
            fSVD_OK(ok), fNumHitsUsed(nHits), fEigenValues(eigenValues), fEigenVectors(eigenVecs), fAvePosition(avePos) {}

        bool                   getSvdOK()                 const {return fSVD_OK;}
        int                    getNumHitsUsed()           const {return fNumHitsUsed;}
        const EigenValues&     getEigenValues()           const {return fEigenValues;}
        const EigenVectors&    getEigenVectors()          const {return fEigenVectors;}
        const Eigen::Vector3d& getAvePosition()           const {return fAvePosition;}
        void                   flipAxis(size_t axis)            { fEigenVectors.row(axis) = -fEigenVectors.row(axis);}

    private:

        bool            fSVD_OK;             ///< SVD Decomposition was successful
        int             fNumHitsUsed;        ///< Number of hits in the decomposition
        EigenValues     fEigenValues;        ///< Eigen values from SVD decomposition
        EigenVectors    fEigenVectors;       ///< The three principle axes
        Eigen::Vector3d fAvePosition;        ///< Average position of hits fed to PCA
    };

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // Compute the principle axes
    void GetPrincipalComponents2D(const HitStatusChargePairVec& hitPairVector, PrincipalComponents2D& pca)                      const;
    void GetPrincipalComponents3D(const HitStatusChargePairVec& hitPairVector, HitPointDirPairMap&, PrincipalComponents3D& pca) const;

    // Reject outliers
    void RejectOutliers(HitStatusChargePairVec& hitPairVector, const PrincipalComponents2D& pca) const;

    // The following typedefs will, obviously, be useful
    double length(const recob::Track* track);

    double projectedLength(const recob::Track* track);

    // The parameters we'll read from the .fcl file.
    std::vector<art::InputTag> fTrackLabelVec;      ///< labels for source of tracks
    unsigned                   fSelectedPlane;      ///< Select hits from this plane
    unsigned                   fMinNumHits;         ///< Minimum number of hits
    float                      fMinTickRange;       ///< Require track to span some number of ticks
    float                      fAssumedELifetime;   ///< Lifetime assumed for calculation
    float                      fMinRejectFraction;  ///< Reject this fraction of "low Q" hits
    float                      fMaxRejectFraction;  ///< Reject this fraction of "high Q" hits
    float                      fOutlierRejectFrac;  ///< Fraction of outliers to reject
    bool                       fUseHitIntegral;     ///< Setting to false swaps to "SummedADC"
    bool                       fWeightByChiSq;      ///< Weight the PCA by the hit chisquare
    bool                       fDiagnosticTuple;    ///< Output the diagnostic tuple

    float                      fSamplingRate;       ///< Recover the sampling rate from the clock data

    // Output tuple variables
    int                        fRunNumber;          ///< run number this event
    int                        fSubRunNumber;       ///< sub run number this event
    int                        fEventNumber;        ///< event number this event
    int                        fCryostat;           ///< Cryostat for given track
    int                        fTPC;                ///< TPC for given track
    int                        fTrackIdx;           ///< index of track
    int                        fWireRange;          ///< Last - first wire number
    int                        fWires;              ///< Number wires spanned
    int                        fTicks;              ///< Number ticks spanned
    double                     fAttenuation;        ///< Attenuation from calc
    double                     fError;              ///< Error from calc
    std::vector<double>        fTrackStartXVec;     ///< Starting x position of track
    std::vector<double>        fTrackStartYVec;     ///< Starting y position of track
    std::vector<double>        fTrackStartZVec;     ///< Starting z position of track
    std::vector<double>        fTrackDirXVec;       ///< Starting x direction of track
    std::vector<double>        fTrackDirYVec;       ///< Starting x direction of track
    std::vector<double>        fTrackDirZVec;       ///< Starting x direction of track
    std::vector<double>        fTrackEndXVec;       ///< Ending x position of track
    std::vector<double>        fTrackEndYVec;       ///< Ending y position of track
    std::vector<double>        fTrackEndZVec;       ///< Ending z position of track
    std::vector<double>        fTrackEndDirXVec;    ///< Ending x direction of track
    std::vector<double>        fTrackEndDirYVec;    ///< Ending x direction of track
    std::vector<double>        fTrackEndDirZVec;    ///< Ending x direction of track
    std::vector<double>        fPCAAxes2D;          ///< Axes for PCA
    std::vector<double>        fEigenValues2D;      ///< Eigen values 
    std::vector<double>        fMeanPosition2D;     ///< Mean position used for PCA
    std::vector<double>        fPCAAxes3D;          ///< Axes for PCA 3D
    std::vector<double>        fEigenValues3D;      ///< Eigen values 3D
    std::vector<double>        fMeanPosition3D;     ///< Mean position used for PCA
    std::vector<double>        fTickVec;            ///< vector of ticks
    std::vector<double>        fChargeVec;          ///< vector of hit charges
    std::vector<double>        fDeltaXVec;          ///< Keep track of hits path length from track fit
    std::vector<double>        fGoodnessOfFitVec;   ///< Goodness of the hit's fit
    std::vector<int>           fDegreesOfFreeVec;   ///< Degrees of freedom
    std::vector<int>           fSnippetLengthVec;   ///< Lenght from start/end of hit
    std::vector<bool>          fGoodHitVec;         ///< Hits were considered good

    TTree*                     fDiagnosticTree;     ///< Pointer to our tree

    int fNumEvents;

    // Other variables that will be shared between different methods.
    const geo::GeometryCore*                   fGeometry;       // pointer to Geometry service
}; // class TPCPurityMonitor

// This macro has to be defined for this module to be invoked from a
// .fcl file; see TPCPurityMonitor.fcl for more information.
DEFINE_ART_MODULE(TPCPurityMonitor)

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
TPCPurityMonitor::TPCPurityMonitor(fhicl::ParameterSet const& parameterSet)
    : EDProducer{parameterSet},
      fDiagnosticTree(nullptr),
      fNumEvents(0)

{
    fGeometry = lar::providerFrom<geo::Geometry>();

    // We're going to output purity objects
    produces<std::vector<anab::TPCPurityInfo>>("",art::Persistable::Yes);

    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
TPCPurityMonitor::~TPCPurityMonitor()
{}

//-----------------------------------------------------------------------
void TPCPurityMonitor::beginJob()
{
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    //double detectorLength = fGeometry->DetLength();

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    if (fDiagnosticTuple)
    {
        art::ServiceHandle<art::TFileService> tfs;
    
        fDiagnosticTree = tfs->make<TTree>("PurityMonitor","");

        fDiagnosticTree->Branch("run",         &fRunNumber,     "run/I");
        fDiagnosticTree->Branch("subrun",      &fSubRunNumber,  "subrun/I");
        fDiagnosticTree->Branch("event",       &fEventNumber,   "event/I");
        fDiagnosticTree->Branch("cryostat",    &fCryostat,      "cryostat/I");
        fDiagnosticTree->Branch("tpc",         &fTPC,           "tpc/I");
        fDiagnosticTree->Branch("trackidx",    &fTrackIdx,      "trackidx/I");
        fDiagnosticTree->Branch("wirerange",   &fWireRange,     "wirerange/I");
        fDiagnosticTree->Branch("nwires",      &fWires,         "nwires/I");
        fDiagnosticTree->Branch("nticks",      &fTicks,         "nticks/I");
        fDiagnosticTree->Branch("attenuation", &fAttenuation,   "attenuation/D");
        fDiagnosticTree->Branch("error",       &fError,         "error/D");

        fDiagnosticTree->Branch("trkstartx",   "std::vector<double>", &fTrackStartXVec);
        fDiagnosticTree->Branch("trkstarty",   "std::vector<double>", &fTrackStartYVec);
        fDiagnosticTree->Branch("trkstartz",   "std::vector<double>", &fTrackStartZVec);
        fDiagnosticTree->Branch("trkdirx",     "std::vector<double>", &fTrackDirXVec);
        fDiagnosticTree->Branch("trkdiry",     "std::vector<double>", &fTrackDirYVec);
        fDiagnosticTree->Branch("trkdirz",     "std::vector<double>", &fTrackDirZVec);
        fDiagnosticTree->Branch("trkendx",     "std::vector<double>", &fTrackEndXVec);
        fDiagnosticTree->Branch("trkendy",     "std::vector<double>", &fTrackEndYVec);
        fDiagnosticTree->Branch("trkendz",     "std::vector<double>", &fTrackEndZVec);
        fDiagnosticTree->Branch("trkenddirx",  "std::vector<double>", &fTrackEndDirXVec);
        fDiagnosticTree->Branch("trkenddiry",  "std::vector<double>", &fTrackEndDirYVec);
        fDiagnosticTree->Branch("trkenddirz",  "std::vector<double>", &fTrackEndDirZVec);
        fDiagnosticTree->Branch("pcavec2d",    "std::vector<double>", &fPCAAxes2D);
        fDiagnosticTree->Branch("eigenvec2d",  "std::vector<double>", &fEigenValues2D);
        fDiagnosticTree->Branch("meanpos2d",   "std::vector<double>", &fMeanPosition2D);
        fDiagnosticTree->Branch("pcavec3d",    "std::vector<double>", &fPCAAxes3D);
        fDiagnosticTree->Branch("eigenvec3d",  "std::vector<double>", &fEigenValues3D);
        fDiagnosticTree->Branch("meanpos3d",   "std::vector<double>", &fMeanPosition3D);
        fDiagnosticTree->Branch("tickvec",     "std::vector<double>", &fTickVec);
        fDiagnosticTree->Branch("chargevec",   "std::vector<double>", &fChargeVec);
        fDiagnosticTree->Branch("deltaxvec",   "std::vector<double>", &fDeltaXVec);
        fDiagnosticTree->Branch("goodnessvec", "std::vector<double>", &fGoodnessOfFitVec);
        fDiagnosticTree->Branch("freedomvec",  "std::vector<int>",    &fDegreesOfFreeVec);
        fDiagnosticTree->Branch("snippetvec",  "std::vector<int>",    &fSnippetLengthVec);
        fDiagnosticTree->Branch("goodhitvec",  "std::vector<bool>",   &fGoodHitVec);
    }


    // zero out the event counter
    fNumEvents = 0;
}

//-----------------------------------------------------------------------
//void TPCPurityMonitor::beginRun(const art::Run& /*run*/)
//{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
//}

//-----------------------------------------------------------------------
void TPCPurityMonitor::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fTrackLabelVec      = p.get< std::vector<art::InputTag> >("TrackLabel",             {""});
    fSelectedPlane      = p.get< unsigned                   >("SelectedPlane",             2);
    fMinNumHits         = p.get< unsigned                   >("MinNumHits",              100);
    fMinTickRange       = p.get< float                      >("MinTickRange",           150.);
    fAssumedELifetime   = p.get< float                      >("AssumedELifetime",    600000.);
    fMinRejectFraction  = p.get< float                      >("MinRejectFraction",      0.05);
    fMaxRejectFraction  = p.get< float                      >("MaxRejectFraction",      0.95);
    fOutlierRejectFrac  = p.get< float                      >("OutlierRejectFrac",      0.70);
    fUseHitIntegral     = p.get< bool                       >("UseHitIntegral",         true);
    fWeightByChiSq      = p.get< bool                       >("WeightByChiSq",         false);
    fDiagnosticTuple    = p.get< bool                       >("DiagnosticTuple",        true);

    return;
}

//-----------------------------------------------------------------------
void TPCPurityMonitor::produce(art::Event& event)
{
    //setup output vector#include "lardataobj/RecoBase/TrackHitMeta.h"
    std::unique_ptr< std::vector<anab::TPCPurityInfo> > outputPtrVector(new std::vector<anab::TPCPurityInfo>());

    fNumEvents++;
    
    // Recover the detector properties.
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);

    fSamplingRate = 1.e-3 * sampling_rate(clockData); // Note sampling rate is in ns, convert to us
    
    for(const auto& trackLabel : fTrackLabelVec)
    {
        // Make a pass through all hits to make contrasting plots
        art::Handle< std::vector<recob::Track>> trackHandle;
        event.getByLabel(trackLabel, trackHandle);
        
        if (!trackHandle.isValid()) continue;

        // I don't know another way to do this... but we need to build a map from hit to spacepoint
        // since we don't seem to have a way to go that direction with what we have for kalman fit
        // tracks. 
        using HitToSpacePointMap = std::unordered_map<const recob::Hit*,const recob::SpacePoint*>;
        
        HitToSpacePointMap hitToSpacePointMap;

        using PointCloud = std::vector<geo::Point_t>;

        PointCloud pointCloud;

        // Recover the collection of associations between tracks and hits and hits and spacepoints
        art::FindManyP<recob::Hit,recob::TrackHitMeta> trackHitAssns(trackHandle, event, trackLabel);

        // Loop over tracks and recover hits
        for(size_t trackIdx = 0; trackIdx < trackHandle->size(); trackIdx++)
        {
            art::Ptr<recob::Track> track(trackHandle,trackIdx);

            const std::vector<art::Ptr<recob::Hit>>&      trackHitsVec(trackHitAssns.at(track.key()));
            const std::vector<const recob::TrackHitMeta*> metaHitsVec = trackHitAssns.data(track.key());

            // Focus on selected hits:
            // 1) Pick out hits on a single plane given by fhicl parameter
            // 2) multiplicity == 1 which should give us clean gaussian shaped pulses
            using TPCToHitMetaPairVecMap = std::unordered_map<unsigned int,HitMetaPairVec>;

            TPCToHitMetaPairVecMap selectedHitMetaVecMap;

            for(size_t idx=0; idx<trackHitsVec.size(); idx++)
            {
                art::Ptr<recob::Hit> hit(trackHitsVec.at(idx));

                if (hit->WireID().Plane == fSelectedPlane && hit->Multiplicity() == 1) selectedHitMetaVecMap[hit->WireID().TPC].emplace_back(hit,metaHitsVec.at(idx));
            }

            if (selectedHitMetaVecMap.empty()) continue;

            // Currently we need to limit the analysis to a single TPC and we have tracks which may have been stitched across the cathode... 
            // For now, we search and find the TPC with the most hits
            TPCToHitMetaPairVecMap::iterator bestMapItr = selectedHitMetaVecMap.begin();

            for(TPCToHitMetaPairVecMap::iterator mapItr = selectedHitMetaVecMap.begin(); mapItr != selectedHitMetaVecMap.end(); mapItr++)
            {
                if (mapItr->second.size() > bestMapItr->second.size()) bestMapItr = mapItr;
            }

            HitMetaPairVec& selectedHitMetaVec = bestMapItr->second;

            // Need a minimum number of hits
            if (selectedHitMetaVec.size() < fMinNumHits) continue;

            // Sort hits by increasing time 
            std::sort(selectedHitMetaVec.begin(),selectedHitMetaVec.end(),[](const auto& left, const auto& right){return left.first->PeakTime() < right.first->PeakTime();});

            // Require track to have a minimum range in ticks
            if (selectedHitMetaVec.back().first->PeakTime() - selectedHitMetaVec.front().first->PeakTime() < fMinTickRange) continue;

            // At this point we should have a vector of art::Ptrs to hits on the selected plane
            // So we should be able to now transition to computing the attenuation
            // Start by forming a vector of pairs of the time (in ticks) and the ln of charge derated by an assumed lifetime
            HitStatusChargePairVec hitStatusChargePairVec;
            HitPointDirPairMap     hitPointDirPairMap;

            float  firstHitTime(selectedHitMetaVec.front().first->PeakTime());
            double maxDeltaX(1.5);   // Assume a "long hit" would be no more than 1.5 cm in length
            double wirePitch(0.3);

            for(const auto& hitMetaPair: selectedHitMetaVec)
            {
                unsigned int trkHitIndex = hitMetaPair.second->Index();
                double       deltaX      = 0.3;                         // Set this to 3 mm just in case no corresponding point
                double       cosTheta    = -100.;

                if (trkHitIndex != std::numeric_limits<unsigned int>::max() && track->HasValidPoint(trkHitIndex))
                {
                    geo::Point_t        hitPos  = track->LocationAtPoint(trkHitIndex);
                    geo::Vector_t       hitDir  = track->DirectionAtPoint(trkHitIndex);
                    const geo::WireGeo& wireGeo = fGeometry->Wire(hitMetaPair.first->WireID());
                    geo::Vector_t       wireDir(wireGeo.Direction()[0],wireGeo.Direction()[1],wireGeo.Direction()[2]);

                    pointCloud.emplace_back(hitPos);

                    cosTheta = std::abs(hitDir.Dot(wireDir));

                    if (cosTheta < 1.)
                    {
                        deltaX = std::min(wirePitch / (1. - cosTheta), maxDeltaX);
                    }
                    else deltaX = maxDeltaX;

                    double charge = fUseHitIntegral ? hitMetaPair.first->Integral() : hitMetaPair.first->SummedADC(); 

                    hitStatusChargePairVec.emplace_back(hitMetaPair,StatusChargePair(true,charge/deltaX));
                    hitPointDirPairMap[hitMetaPair.first.get()] = PointDirPair(hitPos,hitDir);
                }
            }

            size_t numOrig   = hitStatusChargePairVec.size();
            size_t lowCutIdx = fMinRejectFraction * numOrig;
            size_t hiCutIdx  = fMaxRejectFraction * numOrig;

            // Will require a minimum number of hits left over to proceed
            if (lowCutIdx + 10 >= hiCutIdx)
            {
                std::cout << "*****>>>> lowCutIdx >= hiCutIdx: " << lowCutIdx << ", " << hiCutIdx << std::endl;
                continue;
            }

//            HitStatusChargePairVec(hitStatusChargePairVec.begin() + lowCutIdx, hitStatusChargePairVec.begin() + hiCutIdx).swap(hitStatusChargePairVec);

            // Put back in time order
//            std::sort(hitStatusChargePairVec.begin(),hitStatusChargePairVec.end(),[](const auto& left, const auto& right){return left.first->PeakTime() < right.first->PeakTime();});

            // Tag the leading and trailing hits so as to not use them
            std::transform(hitStatusChargePairVec.begin(),hitStatusChargePairVec.begin()+lowCutIdx,hitStatusChargePairVec.begin(),      [](const auto& hitPair){return HitStatusChargePair(hitPair.first,StatusChargePair(false,hitPair.second.second));});
            std::transform(hitStatusChargePairVec.begin()+hiCutIdx,hitStatusChargePairVec.end(),hitStatusChargePairVec.begin()+hiCutIdx,[](const auto& hitPair){return HitStatusChargePair(hitPair.first,StatusChargePair(false,hitPair.second.second));});

            PrincipalComponents2D pca;

            GetPrincipalComponents2D(hitStatusChargePairVec, pca);

            // Reject the outliers
            RejectOutliers(hitStatusChargePairVec, pca);

            // Recompute the pca
            GetPrincipalComponents2D(hitStatusChargePairVec, pca);

            // If the PCA faild then we should bail out 
            if (!pca.getSvdOK()) continue;

            // Now get the 3D PCA so we can use this to help select on track straightness
            PrincipalComponents3D pca3D;

            GetPrincipalComponents3D(hitStatusChargePairVec, hitPointDirPairMap, pca3D);

            const PrincipalComponents2D::EigenVectors& eigenVectors = pca.getEigenVectors();

            double attenuation = eigenVectors.row(1)[1] / eigenVectors.row(1)[0];

            // Want to find the wire range (or should it be the number of wires?)
            unsigned maxWire(0);
            unsigned minWire(100000);

            std::set<unsigned> usedWiresSet;

            for(const auto& hitPair : hitStatusChargePairVec)
            {
                unsigned wire = hitPair.first.first->WireID().Wire;

                usedWiresSet.insert(wire);

                if (wire > maxWire) maxWire = wire;
                if (wire < minWire) minWire = wire;
            }

            geo::WireID wireID  = hitStatusChargePairVec.front().first.first->WireID();
			  
			anab::TPCPurityInfo purityInfo;

			purityInfo.Run         = event.run();
			purityInfo.Subrun      = event.subRun();
			purityInfo.Event       = event.event();
            purityInfo.Cryostat    = wireID.Cryostat;
			purityInfo.TPC         = wireID.TPC;
			purityInfo.Wires       = usedWiresSet.size(); //maxWire - minWire;
            purityInfo.Ticks       = hitStatusChargePairVec.back().first.first->PeakTime() - firstHitTime;
            purityInfo.Attenuation = -attenuation;
			purityInfo.FracError   = std::sqrt(pca.getEigenValues()[0] / pca.getEigenValues()[1]);

            outputPtrVector->emplace_back(purityInfo);

            if (fDiagnosticTuple)
            {
                fRunNumber    = event.run();
                fSubRunNumber = event.subRun();
                fEventNumber  = event.event();
                fCryostat     = wireID.Cryostat;
                fTPC          = wireID.TPC;
                fTrackIdx     = trackIdx; 
                fWireRange    = maxWire - minWire;
                fWires        = usedWiresSet.size();
                fTicks        = hitStatusChargePairVec.back().first.first->PeakTime() - firstHitTime;
                fAttenuation  = -attenuation;
                fError        = std::sqrt(pca.getEigenValues()[0] / pca.getEigenValues()[1]);

                const geo::Point_t& trackStartPos = track->LocationAtPoint(hitStatusChargePairVec.front().first.second->Index());
                const geo::Vector_t trackStartDir = track->DirectionAtPoint(hitStatusChargePairVec.front().first.second->Index());

                fTrackStartXVec.emplace_back(trackStartPos.X());
                fTrackStartYVec.emplace_back(trackStartPos.Y());
                fTrackStartZVec.emplace_back(trackStartPos.Z());
                fTrackDirXVec.emplace_back(trackStartDir.X());
                fTrackDirYVec.emplace_back(trackStartDir.Y());
                fTrackDirZVec.emplace_back(trackStartDir.Z());

                const geo::Point_t& trackEndPos = track->LocationAtPoint(hitStatusChargePairVec.back().first.second->Index());
                const geo::Vector_t trackEndDir = track->DirectionAtPoint(hitStatusChargePairVec.back().first.second->Index());

                fTrackEndXVec.emplace_back(trackEndPos.X());
                fTrackEndYVec.emplace_back(trackEndPos.Y());
                fTrackEndZVec.emplace_back(trackEndPos.Z());
                fTrackEndDirXVec.emplace_back(trackEndDir.X());
                fTrackEndDirYVec.emplace_back(trackEndDir.Y());
                fTrackEndDirZVec.emplace_back(trackEndDir.Z());

                // 2D PCA of time vs charge
                for(size_t rowIdx = 0; rowIdx < 2; rowIdx++)
                {
                    for(size_t colIdx = 0; colIdx < 2; colIdx++) fPCAAxes2D.emplace_back(eigenVectors.row(rowIdx)[colIdx]);
                }

                fEigenValues2D.emplace_back(pca.getEigenValues()[0]); 
                fEigenValues2D.emplace_back(pca.getEigenValues()[1]); 

                fMeanPosition2D.emplace_back(pca.getAvePosition()[0]);
                fMeanPosition2D.emplace_back(pca.getAvePosition()[1]);

                // 3D PCA of track trajectory points
                const PrincipalComponents3D::EigenVectors& eigenVectors3D = pca3D.getEigenVectors();

                for(size_t rowIdx = 0; rowIdx < 3; rowIdx++)
                {
                    for(size_t colIdx = 0; colIdx < 3; colIdx++) fPCAAxes3D.emplace_back(eigenVectors3D.row(rowIdx)[colIdx]);
                }

                fEigenValues3D.emplace_back(pca3D.getEigenValues()[0]); 
                fEigenValues3D.emplace_back(pca3D.getEigenValues()[1]); 
                fEigenValues3D.emplace_back(pca3D.getEigenValues()[2]); 

                fMeanPosition3D.emplace_back(pca3D.getAvePosition()[0]);
                fMeanPosition3D.emplace_back(pca3D.getAvePosition()[1]);
                fMeanPosition3D.emplace_back(pca3D.getAvePosition()[2]);

                for(const auto& hitPair : hitStatusChargePairVec)
                {
                    fTickVec.emplace_back(hitPair.first.first->PeakTime());
                    fChargeVec.emplace_back(hitPair.second.second);
                    fDeltaXVec.emplace_back(hitPair.first.second->Dx());
                    fGoodnessOfFitVec.emplace_back(hitPair.first.first->GoodnessOfFit());
                    fDegreesOfFreeVec.emplace_back(hitPair.first.first->DegreesOfFreedom());
                    fSnippetLengthVec.emplace_back(hitPair.first.first->EndTick() - hitPair.first.first->StartTick());
                    fGoodHitVec.emplace_back(hitPair.second.first);
                }

                fDiagnosticTree->Fill();

                fTrackStartXVec.clear();
                fTrackStartYVec.clear();
                fTrackStartZVec.clear();
                fTrackDirXVec.clear();
                fTrackDirYVec.clear();
                fTrackDirZVec.clear();
                fTrackEndXVec.clear();
                fTrackEndYVec.clear();
                fTrackEndZVec.clear();
                fTrackEndDirXVec.clear();
                fTrackEndDirYVec.clear();
                fTrackEndDirZVec.clear();
                fPCAAxes2D.clear(); 
                fEigenValues2D.clear(); 
                fMeanPosition2D.clear();
                fPCAAxes3D.clear(); 
                fEigenValues3D.clear(); 
                fMeanPosition3D.clear();
                fTickVec.clear(); 
                fChargeVec.clear();
                fDeltaXVec.clear();
                fGoodnessOfFitVec.clear();
                fDegreesOfFreeVec.clear();
                fSnippetLengthVec.clear();
                fGoodHitVec.clear();
            }
        }
    }
    
    //put info onto the event
    event.put(std::move(outputPtrVector));

    return;
}

void TPCPurityMonitor::endJob()
{
    return;
}

// Length of reconstructed track.
//----------------------------------------------------------------------------
double TPCPurityMonitor::length(const recob::Track* track)
{
    double   result(0.);
/*
    TVector3 disp(track->LocationAtPoint(0));
    TVector3 lastPoint(track->LocationAtPoint(0));
    TVector3 lastDir(0.,0.,0.);
    int      n(track->NumberTrajectoryPoints());

    for(int i = 1; i < n; ++i)
    {
        const TVector3& pos = track->LocationAtPoint(i);

        TVector3 trajDir = pos - lastPoint;

        if (trajDir.Mag2()) trajDir.SetMag(1.);

//        if (lastDir.Dot(trajDir) >= 0.)
//        {
            disp   -= pos;
            result += disp.Mag();
            disp    = pos;
//        }

        lastPoint = pos;
        lastDir   = trajDir;
    }
*/
    return result;
}

// Length of reconstructed track.
//----------------------------------------------------------------------------
double TPCPurityMonitor::projectedLength(const recob::Track* track)
{
    double   result(0.);
/*
    TVector3 lastPoint(track->LocationAtPoint(0));
    TVector3 lastDir(track->DirectionAtPoint(0));
    int      n(track->NumberTrajectoryPoints());

    for(int i = 1; i < n; ++i)
    {
        const TVector3& newPoint = track->LocationAtPoint(i);

        TVector3 lastToNewPoint = newPoint - lastPoint;
        double   arcLenToDoca   = lastDir.Dot(lastToNewPoint);

        result    += arcLenToDoca;
        lastPoint  = lastPoint + arcLenToDoca * lastDir;
        lastDir    = track->DirectionAtPoint(i);
    }
*/
    return result;
}


void TPCPurityMonitor::GetPrincipalComponents2D(const HitStatusChargePairVec& hitPairVector, PrincipalComponents2D& pca) const
{
    // Run through the HitPairList and get the mean position of all the hits
    Eigen::Vector2d meanPos(Eigen::Vector2d::Zero());
    double          meanWeightSum(0.);
    int             numPairsInt(0);

    float startTime = 0.; //hitPairVector.front().first->PeakTime();

    for (const auto& hitPair : hitPairVector) 
    {
        if (!hitPair.second.first) continue;

        const recob::Hit* hit = hitPair.first.first.get();

        // Weight the hit by the peak time difference significance
        double weight = fWeightByChiSq ? 1./hit->GoodnessOfFit() : 1.; 

        meanPos(0) += fSamplingRate * (hit->PeakTime() - startTime) * weight;
        meanPos(1) += std::log(hitPair.second.second) * weight;
        numPairsInt++;

        meanWeightSum += weight;
    }

    meanPos /= meanWeightSum;

    // Define elements of our covariance matrix
    double xi2(0.);
    double xiyi(0.);
    double yi2(0.0);
    double weightSum(0.);

    // Back through the hits to build the matrix
    for (const auto& hitPair : hitPairVector) 
    {
        if (!hitPair.second.first) continue;

        const recob::Hit* hit = hitPair.first.first.get();

        double weight = fWeightByChiSq ? 1./hit->GoodnessOfFit() : 1.;

        double x = (fSamplingRate * (hit->PeakTime() - startTime) - meanPos(0)) * weight;
        double y = (std::log(hitPair.second.second) - meanPos(1)) * weight;

        weightSum += weight * weight;

        xi2  += x * x;
        xiyi += x * y;
        yi2  += y * y;
    }

    // Using Eigen package
    Eigen::Matrix2d sig;

    sig << xi2, xiyi, xiyi, yi2;

    sig *= 1. / weightSum;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenMat(sig);

    if (eigenMat.info() == Eigen::ComputationInfo::Success) 
    {
        // Now copy outputPrincipalCo
        // The returned eigen values and vectors will be returned in an xyz system where x is the smallest spread,
        // y is the next smallest and z is the largest. Adopt that convention going forward
        PrincipalComponents2D::EigenValues  eigenVals = eigenMat.eigenvalues();
        PrincipalComponents2D::EigenVectors eigenVecs = eigenMat.eigenvectors().transpose();

        // Store away
        // NOTE: the major axis will be the second entry, the minor axis will be the first
        pca = PrincipalComponents2D(true, numPairsInt, eigenVals, eigenVecs, meanPos);
    }
    else 
    {
        mf::LogDebug("Cluster3D") << "PCA decompose failure, numPairs = " << numPairsInt << std::endl;
        pca = PrincipalComponents2D();
    }

    return;
}

void TPCPurityMonitor::GetPrincipalComponents3D(const HitStatusChargePairVec& hitPairVector, HitPointDirPairMap& hitPointDirPairMap, PrincipalComponents3D& pca) const
{
    // Run through the HitPairList and get the mean position of all the hits
    Eigen::Vector3d meanPos(Eigen::Vector3d::Zero());
    double          meanWeightSum(0.);
    int             numPairsInt(0);

    for (const auto& hitPair : hitPairVector) 
    {
        if (!hitPair.second.first) continue;

        const recob::Hit* hit = hitPair.first.first.get();

        geo::Point_t hitPos = hitPointDirPairMap[hit].first;

        // Weight the hit by the peak time difference significance
        double weight = fWeightByChiSq ? 1./hit->GoodnessOfFit() : 1.; 

        meanPos += Eigen::Vector3d(hitPos.X(),hitPos.Y(),hitPos.Z());

        numPairsInt++;

        meanWeightSum += weight;
    }

    meanPos /= meanWeightSum;

    // Define elements of our covariance matrix
    double xi2(0.);
    double xiyi(0.);
    double xizi(0.);
    double yi2(0.);
    double yizi(0.);
    double zi2(0.);
    double weightSum(0.);

    // Back through the hits to build the matrix
    for (const auto& hitPair : hitPairVector) 
    {
        if (!hitPair.second.first) continue;

        const recob::Hit* hit = hitPair.first.first.get();

        double weight = fWeightByChiSq ? 1./hit->GoodnessOfFit() : 1.;

        geo::Point_t    hitPos         = hitPointDirPairMap[hit].first;
        Eigen::Vector3d weightedHitPos = Eigen::Vector3d(hitPos.X(),hitPos.Y(),hitPos.Z()) - meanPos;

        weightSum += weight * weight;

        xi2  += weightedHitPos[0] * weightedHitPos[0];
        xiyi += weightedHitPos[0] * weightedHitPos[1];
        xizi += weightedHitPos[0] * weightedHitPos[2];
        yi2  += weightedHitPos[1] * weightedHitPos[1];
        yizi += weightedHitPos[1] * weightedHitPos[2];
        zi2  += weightedHitPos[2] * weightedHitPos[2];
    }

    // Using Eigen package
    Eigen::Matrix3d sig;

    sig << xi2, xiyi, xizi, xiyi, yi2, yizi, xizi, yizi, zi2;

    sig *= 1. / weightSum;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenMat(sig);

    if (eigenMat.info() == Eigen::ComputationInfo::Success) 
    {
        // Now copy outputPrincipalCo
        // The returned eigen values and vectors will be returned in an xyz system where x is the smallest spread,
        // y is the next smallest and z is the largest. Adopt that convention going forward
        PrincipalComponents3D::EigenValues  eigenVals = eigenMat.eigenvalues();
        PrincipalComponents3D::EigenVectors eigenVecs = eigenMat.eigenvectors().transpose();

        // Store away
        // NOTE: the major axis will be the second entry, the minor axis will be the first
        pca = PrincipalComponents3D(true, numPairsInt, eigenVals, eigenVecs, meanPos);
    }
    else 
    {
        mf::LogDebug("Cluster3D") << "PCA decompose failure, numPairs = " << numPairsInt << std::endl;
        pca = PrincipalComponents3D();
    }

    return;
}

void TPCPurityMonitor::RejectOutliers(HitStatusChargePairVec& hitPairVector, const PrincipalComponents2D& pca) const
{
    double                 slope  = pca.getEigenVectors().row(1)[1] / pca.getEigenVectors().row(1)[0];
    const Eigen::Vector2d& avePos = pca.getAvePosition(); 

    using HitPairDeltaLogChargePair    = std::pair<HitStatusChargePair*,double>;
    using HitPairDeltaLogChargePairVec = std::vector<HitPairDeltaLogChargePair>;

    HitPairDeltaLogChargePairVec hitPairDeltaLogChargePairVec;

    for(auto& hitPair : hitPairVector)
    {
        // We are only interested in the "good" hits here
        if (hitPair.second.first)
        {
            double predLogCharge  = (fSamplingRate * hitPair.first.first->PeakTime() - avePos[0]) * slope + avePos[1];
            double deltaLogCharge = std::log(hitPair.second.second) - predLogCharge;

            hitPairDeltaLogChargePairVec.emplace_back(&hitPair,deltaLogCharge);
        }
    }

    // Sort hits by their deviation from the prediction
    std::sort(hitPairDeltaLogChargePairVec.begin(),hitPairDeltaLogChargePairVec.end(),[](const auto& left, const auto& right){return left.second < right.second;});

    // Go through and tag those we are rejecting
    size_t loRejectIdx = 0.01 * hitPairDeltaLogChargePairVec.size();
    size_t hiRejectIdx = fOutlierRejectFrac * hitPairDeltaLogChargePairVec.size();

    const double outlierReject = 0.75;

    for(size_t idx = 0; idx < hitPairDeltaLogChargePairVec.size(); idx++)
    {
        if (idx < loRejectIdx || idx > hiRejectIdx)                    hitPairDeltaLogChargePairVec[idx].first->second.first = false;
        if (hitPairDeltaLogChargePairVec[idx].second < -outlierReject) hitPairDeltaLogChargePairVec[idx].first->second.first = false;
    }

    return;
}


} // namespace TPCPurityMonitor

#endif // TPCPurityMonitor_module
