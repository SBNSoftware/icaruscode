
#include "icaruscode/Analysis/tools/MCTruth/IMCTruthMatching.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTracker.h"

#include <cmath>
#include <algorithm>

namespace truth
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       BackTrackerTruth
    // Module Type: art tool
    // File:        BackTrackerTruth.h
    //
    //              This provides MC truth information by interfacing to the
    //              LarSoft BackTracker service
    //
    // Configuration parameters:
    //
    // TruncMeanFraction     - the fraction of waveform bins to discard when
    //
    // Created by Tracy Usher (usher@slac.stanford.edu) on November 21, 2017
    //
    ////////////////////////////////////////////////////////////////////////

class BackTrackerTruth : virtual public IMCTruthMatching
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit BackTrackerTruth(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~BackTrackerTruth();
    
    /**
     *  @brief This rebuilds the internal maps, is a noop for this module since
     *         the BackTracker is a service and rebuilds its maps at begin run
     */
    void Rebuild(const art::Event& evt) {};

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset) override;

    /**
     *  @brief Get a reference to the ParticleList
     */
    const sim::ParticleList& ParticleList() const override;
    
    //    // Set the EveIdCalculator for the owned ParticleList
    //    void  SetEveIdCalculator(sim::EveIdCalculator *ec) { fParticleList.AdoptEveIdCalculator(ec); }
    
    // Return a pointer to the simb::MCParticle object corresponding to
    // the given TrackID
    const simb::MCParticle* TrackIDToParticle(int const& id)       const override;
    const simb::MCParticle* TrackIDToMotherParticle(int const& id) const override;
    
    // Get art::Ptr<> to simb::MCTruth and related information
    const art::Ptr<simb::MCTruth>&                TrackIDToMCTruth(int const& id)                        const override;
    const art::Ptr<simb::MCTruth>&                ParticleToMCTruth(const simb::MCParticle* p)           const override;
    std::vector<const simb::MCParticle*>          MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const override;
    const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector()                                        const override;
    
    // this method will return the Geant4 track IDs of
    // the particles contributing ionization electrons to the identified hit
    std::vector<sim::TrackIDE> HitToTrackID(recob::Hit const& hit)           const override;
    std::vector<sim::TrackIDE> HitToTrackID(art::Ptr<recob::Hit> const& hit) const override;
    
    // method to return a subset of allhits that are matched to a list of TrackIDs
    const std::vector<std::vector<art::Ptr<recob::Hit>>> TrackIDsToHits(std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                                std::vector<int> const& tkIDs) const override;
    
    // method to return the EveIDs of particles contributing ionization
    // electrons to the identified hit
    std::vector<sim::TrackIDE> HitToEveID(art::Ptr<recob::Hit> const& hit) const override;
    
    // method to return the XYZ position of the weighted average energy deposition for a given hit
    std::vector<double>  HitToXYZ(art::Ptr<recob::Hit> const& hit) const override;
    
    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    std::vector<double> SpacePointToXYZ(art::Ptr<recob::SpacePoint> const& spt,
                                                art::Event                  const& evt,
                                                std::string                 const& label) const override;
    
    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    std::vector<double> SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits) const override;
    
    // method to return the fraction of hits in a collection that come from the specified Geant4 track ids
    double HitCollectionPurity(std::set<int>                              trackIDs,
                                       std::vector< art::Ptr<recob::Hit> > const& hits) const override;
    
    // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
    // represented in a collection of hits
    double HitCollectionEfficiency(std::set<int>                              trackIDs,
                                           std::vector< art::Ptr<recob::Hit> > const& hits,
                                           std::vector< art::Ptr<recob::Hit> > const& allhits,
                                           geo::View_t                         const& view) const override;
    
    // method to return the fraction of charge in a collection that come from the specified Geant4 track ids
    double HitChargeCollectionPurity(std::set<int>                              trackIDs,
                                             std::vector< art::Ptr<recob::Hit> > const& hits) const override;
    
    // method to return the fraction of all charge in an event from a specific set of Geant4 track IDs that are
    // represented in a collection of hits
    double HitChargeCollectionEfficiency(std::set<int>                              trackIDs,
                                                 std::vector< art::Ptr<recob::Hit> > const& hits,
                                                 std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                 geo::View_t                         const& view) const override;
    
    // method to return all EveIDs corresponding to the current sim::ParticleList
    std::set<int> GetSetOfEveIDs() const override;
    
    // method to return all TrackIDs corresponding to the current sim::ParticleList
    std::set<int> GetSetOfTrackIDs() const override;
    
    // method to return all EveIDs corresponding to the given list of hits
    std::set<int> GetSetOfEveIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const override;
    
    // method to return all TrackIDs corresponding to the given list of hits
    std::set<int> GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const override;

private:
    
    // Fcl parameters.
    art::InputTag fTrackProducerLabel; ///< tag for finding the tracks
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*                      fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties*            fDetectorProperties;   ///< Detector properties service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
BackTrackerTruth::BackTrackerTruth(fhicl::ParameterSet const & pset)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("BackTrackerTruth") << "BackTrackerTruth configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
BackTrackerTruth::~BackTrackerTruth()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void BackTrackerTruth::reconfigure(fhicl::ParameterSet const & pset)
{
    fTrackProducerLabel = pset.get<art::InputTag>("TrackProducerLabel", "");
}
    
//----------------------------------------------------------------------
const sim::ParticleList& BackTrackerTruth::ParticleList() const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->ParticleList();
}

    
//----------------------------------------------------------------------
const simb::MCParticle* BackTrackerTruth::TrackIDToParticle(int const& id) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->TrackIDToParticle(id);
}
    
//----------------------------------------------------------------------
const simb::MCParticle* BackTrackerTruth::TrackIDToMotherParticle(int const& id) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->TrackIDToMotherParticle(id);
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& BackTrackerTruth::TrackIDToMCTruth(int const& id) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->TrackIDToMCTruth(id);
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& BackTrackerTruth::ParticleToMCTruth(const simb::MCParticle* p) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->ParticleToMCTruth(p);
}
    
//----------------------------------------------------------------------
std::vector<const simb::MCParticle*> BackTrackerTruth::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->MCTruthToParticles(mct);
}
    
//----------------------------------------------------------------------
const std::vector< art::Ptr<simb::MCTruth> >&  BackTrackerTruth::MCTruthVector() const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->MCTruthVector();
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> BackTrackerTruth::HitToTrackID(recob::Hit const& hit) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitToTrackID(hit);
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> BackTrackerTruth::HitToTrackID(art::Ptr<recob::Hit> const& hit) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitToTrackID(*hit);
}

//----------------------------------------------------------------------
const std::vector<std::vector<art::Ptr<recob::Hit>>> BackTrackerTruth::TrackIDsToHits(std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                                 std::vector<int> const& tkIDs) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->TrackIDsToHits(allhits,tkIDs);
}
    
//----------------------------------------------------------------------
// plist is assumed to have adopted the appropriate EveIdCalculator prior to
// having been passed to this method. It is likely that the EmEveIdCalculator is
// the one you always want to use
std::vector<sim::TrackIDE> BackTrackerTruth::HitToEveID(art::Ptr<recob::Hit> const& hit) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitToEveID(hit);
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfEveIDs() const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->GetSetOfEveIDs();
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfTrackIDs() const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->GetSetOfTrackIDs();
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfEveIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->GetSetOfEveIDs(hits);
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->GetSetOfTrackIDs(hits);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitCollectionPurity(std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitCollectionPurity(trackIDs, hits);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitChargeCollectionPurity(std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitChargeCollectionPurity(trackIDs, hits);
}
    
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitCollectionEfficiency(std::set<int> trackIDs,
                                            std::vector< art::Ptr<recob::Hit> > const& hits,
                                            std::vector< art::Ptr<recob::Hit> > const& allhits,
                                            geo::View_t const& view) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitCollectionEfficiency(trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitChargeCollectionEfficiency(std::set<int>                              trackIDs,
                                                  std::vector< art::Ptr<recob::Hit> > const& hits,
                                                  std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                  geo::View_t                         const& view) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitChargeCollectionEfficiency(trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::HitToXYZ(art::Ptr<recob::Hit> const& hit) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->HitToXYZ(hit);
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::SpacePointToXYZ(art::Ptr<recob::SpacePoint>         const& spt,
                                                 art::Event                          const& evt,
                                                 std::string                         const& label) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->SpacePointToXYZ(spt, evt, label);
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits) const
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    return backTracker->SpacePointHitsToXYZ(hits);
}

//----------------------------------------------------------------------------
    
DEFINE_ART_CLASS_TOOL(BackTrackerTruth)
}
