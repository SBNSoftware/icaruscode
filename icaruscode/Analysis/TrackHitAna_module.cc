// TrackHitAna_module.cc
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

#ifndef TrackHitAna_module
#define TrackHitAna_module

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

//#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "icaruscode/Analysis/tools/IHitHistogramTool.h"

// C++ Includes
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

#include <iostream>
#include <fstream>

namespace TrackHitAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class TrackHitAna : public art::EDAnalyzer
{
public:

    // Standard constructor and destructor for an ART module.
    explicit TrackHitAna(fhicl::ParameterSet const& pset);
    virtual ~TrackHitAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();
    void endJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event.
    void analyze (const art::Event& evt);

private:

    // Traverse PFParticle hierarchy
    int traversePFParticleHierarchy(art::Handle<std::vector<recob::PFParticle>>&,
                                    size_t,
                                    const art::FindManyP<recob::Track>&,
                                    const art::FindManyP<recob::Vertex>&,
                                    int&,
                                    int&) const;

    // The following typedefs will, obviously, be useful
    double length(const recob::Track* track);

    double projectedLength(const recob::Track* track);

    // The parameters we'll read from the .fcl file.
    art::InputTag fHitProducerLabel;
    art::InputTag fPFParticleProducerLabel;
    art::InputTag fTrackProducerLabel;
    art::InputTag fWireProducerLabel;

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fNumEvents;
    
    // Keep track of the hit histogramming tools here
    std::vector<std::unique_ptr<IHitHistogramTool>> fHitHistogramToolVec;

    std::vector<std::vector<double>> fChannelPedVec;

    // Other variables that will be shared between different methods.
    const geo::GeometryCore*           fGeometry;       // pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;
    const lariov::DetPedestalProvider& fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
}; // class TrackHitAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
TrackHitAna::TrackHitAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet),
      fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())

{
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
TrackHitAna::~TrackHitAna()
{}

//-----------------------------------------------------------------------
void TrackHitAna::beginJob()
{
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    //double detectorLength = fGeometry->DetLength();

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;

    // The arguments to 'make<whatever>' are the same as those passed
    // to the 'whatever' constructor, provided 'whatever' is a ROOT
    // class that TFileService recognizes.
    for (auto& hitHistTool : fHitHistogramToolVec) hitHistTool->initializeHists(tfs, "AllHits");

    // zero out the event counter
    fNumEvents = 0;
}

//-----------------------------------------------------------------------
void TrackHitAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void TrackHitAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fHitProducerLabel        = p.get< std::string >("HitModuleLabel",          "gauss");
    fPFParticleProducerLabel = p.get< std::string >("PFParticleProducerLabel", "cluster3d");
    fTrackProducerLabel      = p.get< std::string >("TrackProducerLabel",      "trackkalmanhit");
    fWireProducerLabel       = p.get< std::string >("WireProducerLabel",       "caldata");
    
    // Implement the tools for handling the responses
    const std::vector<fhicl::ParameterSet>& hitHistogramToolVec = p.get<std::vector<fhicl::ParameterSet>>("HitHistogramToolList");
    
    for(auto& hitHistogramTool : hitHistogramToolVec)
        fHitHistogramToolVec.push_back(art::make_tool<IHitHistogramTool>(hitHistogramTool));

    return;
}

//-----------------------------------------------------------------------
void TrackHitAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();

    fNumEvents++;
    
    // Make a pass through all hits to make contrasting plots
    art::Handle< std::vector<recob::Hit> > hitHandle;
    event.getByLabel(fHitProducerLabel, hitHandle);

    if (hitHandle.isValid())
    {
        IHitHistogramTool::HitPtrVec allHitVec;
        art::fill_ptr_vector(allHitVec, hitHandle);

        for(auto& hitHistTool : fHitHistogramToolVec) hitHistTool->fillHistograms(allHitVec);
    }

    return;
}

void TrackHitAna::endJob()
{
    // Make a call to normalize histograms if so desired
    for(auto& hitHistTool : fHitHistogramToolVec) hitHistTool->endJob(fNumEvents);

    return;
}

// Length of reconstructed track.
//----------------------------------------------------------------------------
double TrackHitAna::length(const recob::Track* track)
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
double TrackHitAna::projectedLength(const recob::Track* track)
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

int TrackHitAna::traversePFParticleHierarchy(art::Handle<std::vector<recob::PFParticle>>& pfParticleHandle,
                                             size_t                                       pfParticleIdx,
                                             const art::FindManyP<recob::Track>&          trackAssns,
                                             const art::FindManyP<recob::Vertex>&         vertexAssns,
                                             int&                                         nTracks,
                                             int&                                         nVertices) const
{
    // So far no daughters...
    int nDaughters(0);
/*
    // Get pointer to PFParticle
    art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfParticleIdx);

    // Recover tracks/vertices associated to this PFParticle
    std::vector<art::Ptr<recob::Track>>  pfPartTrackVec  = trackAssns.at(pfParticle.key());
    std::vector<art::Ptr<recob::Vertex>> pfPartVertexVec = vertexAssns.at(pfParticle.key());

    nTracks    += pfPartTrackVec.size();
    nVertices  += pfPartVertexVec.size();
    nDaughters += pfParticle->Daughters().size();

    for(auto& daughterIdx : pfParticle->Daughters())
    {
        nDaughters += traversePFParticleHierarchy(pfParticleHandle, daughterIdx, trackAssns, vertexAssns, nTracks, nVertices);
    }
*/
    return nDaughters;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see TrackHitAna.fcl for more information.
DEFINE_ART_MODULE(TrackHitAna)

} // namespace TrackHitAna

#endif // TrackHitAna_module
