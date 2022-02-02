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
    void beginRun(const art::Run& run);

    // The analysis routine, called once per event.
    void produce(art::Event& evt);

private:
    // Definie the basic data structure we will use
    using StatusChargePair       = std::pair<bool,double>;
    using HitStatusChargePair    = std::pair<const recob::Hit*,StatusChargePair>;
    using HitStatusChargePairVec = std::vector<HitStatusChargePair>;

    // We also need to define a container for the output of the PCA Analysis
    class PrincipalComponents
    {
    public:

        using EigenValues  = Eigen::Vector2d;
        using EigenVectors = Eigen::Matrix2d;

        PrincipalComponents() :
            fSVD_OK(false), fNumHitsUsed(0), fEigenValues(EigenValues::Zero()), fEigenVectors(EigenVectors::Zero()), fAvePosition(Eigen::Vector2d::Zero()) {}

    private:

        bool            fSVD_OK;             ///< SVD Decomposition was successful
        int             fNumHitsUsed;        ///< Number of hits in the decomposition
        EigenValues     fEigenValues;        ///< Eigen values from SVD decomposition
        EigenVectors    fEigenVectors;       ///< The three principle axes
        Eigen::Vector2d fAvePosition;        ///< Average position of hits fed to PCA

    public:

        PrincipalComponents(bool ok, int nHits, const EigenValues& eigenValues, const EigenVectors& eigenVecs, const Eigen::Vector2d& avePos) :
            fSVD_OK(ok), fNumHitsUsed(nHits), fEigenValues(eigenValues), fEigenVectors(eigenVecs), fAvePosition(avePos) {}

        bool                   getSvdOK()                 const {return fSVD_OK;}
        int                    getNumHitsUsed()           const {return fNumHitsUsed;}
        const EigenValues&     getEigenValues()           const {return fEigenValues;}
        const EigenVectors&    getEigenVectors()          const {return fEigenVectors;}
        const Eigen::Vector2d& getAvePosition()           const {return fAvePosition;}

        void                   flipAxis(size_t axis)            { fEigenVectors.row(axis) = -fEigenVectors.row(axis);}
    };

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // Compute the principle axes
    void GetPrincipalComponents(const HitStatusChargePairVec& hitPairVector, PrincipalComponents& pca) const;

    // Reject outliers
    void RejectOutliers(HitStatusChargePairVec& hitPairVector, const PrincipalComponents& pca) const;

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
    std::vector<double>        fPCAAxes;            ///< Axes for PCA
    std::vector<double>        fEigenValues;        ///< Eigen values 
    std::vector<double>        fMeanPosition;       ///< Mean position used for PCA
    std::vector<double>        fTickVec;            ///< vector of ticks
    std::vector<double>        fChargeVec;          ///< vector of hit charges
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

        fDiagnosticTree->Branch("pcavec",     "std::vector<double>", &fPCAAxes);
        fDiagnosticTree->Branch("eigenvec",   "std::vector<double>", &fEigenValues);
        fDiagnosticTree->Branch("meanpos",    "std::vector<double>", &fMeanPosition);
        fDiagnosticTree->Branch("tickvec",    "std::vector<double>", &fTickVec);
        fDiagnosticTree->Branch("chargevec",  "std::vector<double>", &fChargeVec);
        fDiagnosticTree->Branch("goodhitvec", "std::vector<bool>",   &fGoodHitVec);
    }


    // zero out the event counter
    fNumEvents = 0;
}

//-----------------------------------------------------------------------
void TPCPurityMonitor::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void TPCPurityMonitor::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fTrackLabelVec      = p.get< std::vector<art::InputTag> >("TrackLabel",    std::vector<art::InputTag>() = {""});
    fSelectedPlane      = p.get< unsigned                   >("SelectedPlane",                                   2);
    fMinNumHits         = p.get< unsigned                   >("MinNumHits",                                    100);
    fMinTickRange       = p.get< float                      >("MinTickRange",                                 150.);
    fAssumedELifetime   = p.get< float                      >("AssumedELifetime",                          600000.);
    fMinRejectFraction  = p.get< float                      >("MinRejectFraction",                            0.05);
    fMaxRejectFraction  = p.get< float                      >("MaxRejectFraction",                            0.95);
    fOutlierRejectFrac  = p.get< float                      >("OutlierRejectFrac",                            0.70);
    fUseHitIntegral     = p.get< bool                       >("UseHitIntegral",                               true);
    fWeightByChiSq      = p.get< bool                       >("WeightByChiSq",                               false);
    fDiagnosticTuple    = p.get< bool                       >("DiagnosticTuple",                              true);

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

        // Recover the collection of associations between tracks and hits
        art::FindManyP<recob::Hit,recob::TrackHitMeta> trackHitAssns(trackHandle, event, trackLabel);

        // Loop over tracks and recover hits
        for(size_t trackIdx = 0; trackIdx < trackHandle->size(); trackIdx++)
        {
            art::Ptr<recob::Track> track(trackHandle,trackIdx);

            const std::vector<art::Ptr<recob::Hit>>&      trackHitsVec(trackHitAssns.at(track.key()));
//            const std::vector<const recob::TrackHitMeta*> metaHitsVec = trackHitAssns.data(track.key());

            // Focus on selected hits:
            // 1) Pick out hits on a single plane given by fhicl parameter
            // 2) multiplicity == 1 which should give us clean gaussian shaped pulses
            std::vector<art::Ptr<recob::Hit>> selectedTrackHitsVec;

            for(auto& hit : trackHitsVec)
            {
                if (hit->WireID().Plane == fSelectedPlane && hit->Multiplicity() == 1) selectedTrackHitsVec.emplace_back(hit);
            }

            // Need a minimum number of hits
            if (selectedTrackHitsVec.size() < fMinNumHits) continue;

//            std::cout << "--> trackMeta index: " << metaHitsVec.front()->Index() << ", Dx: " << metaHitsVec.front()->Dx() << ", meta size: " << metaHitsVec.size() << ", track:" << trackHitsVec.size() << std::endl;

            // Sort hits by increasing time 
            std::sort(selectedTrackHitsVec.begin(),selectedTrackHitsVec.end(),[](const auto& left, const auto& right){return left->PeakTime() < right->PeakTime();});

            // Require track to have a minimum range in ticks
            if (selectedTrackHitsVec.back()->PeakTime() - selectedTrackHitsVec.front()->PeakTime() < fMinTickRange) continue;

            // At this point we should have a vector of art::Ptrs to hits on the selected plane
            // So we should be able to now transition to computing the attenuation
            // Start by forming a vector of pairs of the time (in ticks) and the ln of charge derated by an assumed lifetime
            HitStatusChargePairVec hitStatusChargePairVec;

            float firstHitTime(selectedTrackHitsVec.front()->PeakTime());

            for(const auto& hit : selectedTrackHitsVec)
            {
                float charge    = fUseHitIntegral ? hit->Integral() : hit->SummedADC(); 
//                float logCharge = charge * exp(fSamplingRate * (hit->PeakTime() - firstHitTime) / fAssumedELifetime);

                hitStatusChargePairVec.emplace_back(hit.get(),StatusChargePair(true,log(charge)));
            }

            // Drop the smallest and largest charges
//            std::sort(HitStatusChargePairVec.begin(),HitStatusChargePairVec.end(),[](const auto& left, const auto& right){return left.second.second < right.second.second;});

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

            PrincipalComponents pca;

            GetPrincipalComponents(hitStatusChargePairVec, pca);

            // Reject the outliers
            RejectOutliers(hitStatusChargePairVec, pca);

            // Recompute the pca
            GetPrincipalComponents(hitStatusChargePairVec, pca);

            const PrincipalComponents::EigenVectors& eigenVectors = pca.getEigenVectors();

            double attenuation = eigenVectors.row(1)[1] / eigenVectors.row(1)[0];

            // Want to find the wire range (or should it be the number of wires?)
            unsigned maxWire(0);
            unsigned minWire(100000);

            std::set<unsigned> usedWiresSet;

            for(const auto& hitPair : hitStatusChargePairVec)
            {
                unsigned wire = hitPair.first->WireID().Wire;

                usedWiresSet.insert(wire);

                if (wire > maxWire) maxWire = wire;
                if (wire < minWire) minWire = wire;
            }

            geo::WireID wireID  = hitStatusChargePairVec.front().first->WireID();
			  
			anab::TPCPurityInfo purityInfo;

			purityInfo.Run         = event.run();
			purityInfo.Subrun      = event.subRun();
			purityInfo.Event       = event.event();
            purityInfo.Cryostat    = wireID.Cryostat;
			purityInfo.TPC         = wireID.TPC;
			purityInfo.Wires       = usedWiresSet.size(); //maxWire - minWire;
            purityInfo.Ticks       = hitStatusChargePairVec.back().first->PeakTime() - firstHitTime;
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
                fTicks        = hitStatusChargePairVec.back().first->PeakTime() - firstHitTime;
                fAttenuation  = -attenuation;
                fError        = std::sqrt(pca.getEigenValues()[0] / pca.getEigenValues()[1]);

                for(size_t rowIdx = 0; rowIdx < 2; rowIdx++)
                {
                    for(size_t colIdx = 0; colIdx < 2; colIdx++) fPCAAxes.emplace_back(eigenVectors.row(rowIdx)[colIdx]);
                }

                fEigenValues.emplace_back(pca.getEigenValues()[0]); 
                fEigenValues.emplace_back(pca.getEigenValues()[1]); 

                fMeanPosition.emplace_back(pca.getAvePosition()[0]);
                fMeanPosition.emplace_back(pca.getAvePosition()[1]);

                for(const auto& hitPair : hitStatusChargePairVec)
                {
                    double charge = fUseHitIntegral ? hitPair.first->Integral() : hitPair.first->SummedADC();

                    fTickVec.emplace_back(hitPair.first->PeakTime());
                    fChargeVec.emplace_back(charge);
                    fGoodHitVec.emplace_back(hitPair.second.first);
                }

                fDiagnosticTree->Fill();

                fPCAAxes.clear(); 
                fEigenValues.clear(); 
                fMeanPosition.clear();
                fTickVec.clear(); 
                fChargeVec.clear(); 
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



void TPCPurityMonitor::GetPrincipalComponents(const HitStatusChargePairVec& hitPairVector, PrincipalComponents& pca) const
{
    // Run through the HitPairList and get the mean position of all the hits
    Eigen::Vector2d meanPos(Eigen::Vector2d::Zero());
    double          meanWeightSum(0.);
    int             numPairsInt(0);

    float startTime = hitPairVector.front().first->PeakTime();

    for (const auto& hitPair : hitPairVector) 
    {
        if (!hitPair.second.first) continue;

        const recob::Hit* hit = hitPair.first;

        // Weight the hit by the peak time difference significance
        double weight = fWeightByChiSq ? 1./hit->GoodnessOfFit() : 1.; 
        double charge = fUseHitIntegral ? hit->Integral() : hit->SummedADC();

        meanPos(0) += fSamplingRate * (hit->PeakTime() - startTime) * weight;
        meanPos(1) += std::log(charge) * weight;
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

        const recob::Hit* hit = hitPair.first;

        double weight = fWeightByChiSq ? 1./hit->GoodnessOfFit() : 1.;
        double charge = fUseHitIntegral ? hit->Integral() : hit->SummedADC();

        double x = (fSamplingRate * (hit->PeakTime() - startTime) - meanPos(0)) * weight;
        double y = (std::log(charge) - meanPos(1)) * weight;

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
        // Now copy output
        // The returned eigen values and vectors will be returned in an xyz system where x is the smallest spread,
        // y is the next smallest and z is the largest. Adopt that convention going forward
        PrincipalComponents::EigenValues  eigenVals = eigenMat.eigenvalues();
        PrincipalComponents::EigenVectors eigenVecs = eigenMat.eigenvectors().transpose();

        // Store away
        // NOTE: the major axis will be the second entry, the minor axis will be the first
        pca = PrincipalComponents(true, numPairsInt, eigenVals, eigenVecs, meanPos);
    }
    else 
    {
        mf::LogDebug("Cluster3D") << "PCA decompose failure, numPairs = " << numPairsInt << std::endl;
        pca = PrincipalComponents();
    }

    return;
}

void TPCPurityMonitor::RejectOutliers(HitStatusChargePairVec& hitPairVector, const PrincipalComponents& pca) const
{
    double                 slope  = pca.getEigenVectors().row(1)[1] / pca.getEigenVectors().row(1)[0];
    const Eigen::Vector2d& avePos = pca.getAvePosition(); 

    // We assume the input vector has been time ordered
    double firstHitTime = hitPairVector.front().first->PeakTime();

    for(auto& hitPair : hitPairVector)
    {
        // We are only interested in the "good" hits here
        if (hitPair.second.first)
        {
            double charge = fUseHitIntegral ? hitPair.first->Integral() : hitPair.first->SummedADC();

            double predLogCharge  = (fSamplingRate * (hitPair.first->PeakTime() - firstHitTime) - avePos[0]) * slope + avePos[1];
            double deltaLogCharge = std::log(charge) - predLogCharge;

            hitPair.second.second = deltaLogCharge;
        }
        else hitPair.second.second = 0.;  // These hits already rejected, don't double count them. 
    }

    // Sort hits by their deviation from the prediction
    std::sort(hitPairVector.begin(),hitPairVector.end(),[](const auto& left, const auto& right){return abs(left.second.second) < abs(right.second.second);});

    // Go through and tag those we are rejecting
    size_t rejectIdx = fOutlierRejectFrac * hitPairVector.size();

    for(size_t idx = rejectIdx; idx < hitPairVector.size(); idx++) hitPairVector[idx].second.first = false;

    // Put back in time order
    std::sort(hitPairVector.begin(), hitPairVector.end(), [](const auto& left, const auto& right){return left.first->PeakTime() < right.first->PeakTime();});

    return;
}


} // namespace TPCPurityMonitor

#endif // TPCPurityMonitor_module
