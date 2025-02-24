#include "icaruscode/Analysis/tools/ITrackHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/WireReadout.h"

#include "lardataobj/RecoBase/Hit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <cmath>
#include <algorithm>

namespace BasicTrackAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       BasicTrackAnalysis
    // Module Type: producer
    // File:        BasicTrackAnalysis.h
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

class BasicTrackAnalysis : virtual public ITrackHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit BasicTrackAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~BasicTrackAnalysis();
    
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
    
    // Fcl parameters.
    std::vector<art::InputTag> fTrackProducerLabelVec; ///< tag for finding the tracks
    std::string                fLocalDirName;          ///< Directory name for histograms
    
    // Pointers to the histograms we'll create.
    TH1D*     fHitsByWire[3];
    TH2D*     fPulseHVsWidth[3];
    TProfile* fPulseHVsHitNo[3];
    
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
BasicTrackAnalysis::BasicTrackAnalysis(fhicl::ParameterSet const & pset)
{
    configure(pset);
    
    // Report.
    mf::LogInfo("BasicTrackAnalysis") << "BasicTrackAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
BasicTrackAnalysis::~BasicTrackAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void BasicTrackAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fTrackProducerLabelVec = pset.get<std::vector<art::InputTag>>("TrackProducerLabel",   std::vector<art::InputTag>()={""});
    fLocalDirName          = pset.get<std::string               >(  "LocalDirName",       std::string("wow"));
}

//----------------------------------------------------------------------------
/// Begin job method.
void BasicTrackAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());

    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const n_wires_0 = wireReadout.Nwires(geo::PlaneID{0, 0, 0});
    auto const n_wires_1 = wireReadout.Nwires(geo::PlaneID{0, 0, 1});
    auto const n_wires_2 = wireReadout.Nwires(geo::PlaneID{0, 0, 2});
    fHitsByWire[0]    = dir.make<TH1D>("HitsByWire0", ";Wire #", n_wires_0, 0., n_wires_0);
    fHitsByWire[1]    = dir.make<TH1D>("HitsByWire1", ";Wire #", n_wires_1, 0., n_wires_1);
    fHitsByWire[2]    = dir.make<TH1D>("HitsByWire2", ";Wire #", n_wires_2, 0., n_wires_2);
    
    fPulseHVsWidth[0] = dir.make<TH2D>("PHVsWidth0",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fPulseHVsWidth[1] = dir.make<TH2D>("PHVsWidth1",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fPulseHVsWidth[2] = dir.make<TH2D>("PHVsWidth2",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    
    fPulseHVsHitNo[0] = dir.make<TProfile>("PHVsNo0",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);
    fPulseHVsHitNo[1] = dir.make<TProfile>("PHVsNo1",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);
    fPulseHVsHitNo[2] = dir.make<TProfile>("PHVsNo2",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);

    return;
}
    
void BasicTrackAnalysis::fillHistograms(const art::Event& event) const
{
    // The game plan for this module is to look at recob::Tracks and objects associated to tracks
    // To do this we need a valid track collection for those we are hoping to look at
    for(const auto& trackLabel : fTrackProducerLabelVec)
    {
        art::Handle<std::vector<recob::Track> > trackHandle;
        event.getByLabel(trackLabel, trackHandle);
        
        if (trackHandle.isValid())
        {
            // Recover the collection of associations between tracks and hits
            art::FindManyP<recob::Hit> trackHitAssns(trackHandle, event, trackLabel);
            
            for(size_t trackIdx = 0; trackIdx < trackHandle->size(); trackIdx++)
            {
                art::Ptr<recob::Track> track(trackHandle,trackIdx);
                
            }
        }
    }

    return;
}
    
// Useful for normalizing histograms
void BasicTrackAnalysis::endJob(int numEvents)
{
    // Normalize wire profiles to be hits/event
    double normFactor(1./numEvents);
    
    for(size_t idx = 0; idx < 3; idx++) fHitsByWire[idx]->Scale(normFactor);
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(BasicTrackAnalysis)
}
