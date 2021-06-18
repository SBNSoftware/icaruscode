#include "icaruscode/Analysis/tools/MCTruth/IMCTruthMatching.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

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
    void Rebuild(const art::Event& evt) override {};

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
    const simb::MCParticle* TrackIDToParticle(int const id)       const override;
    const simb::MCParticle* TrackIDToMotherParticle(int const id) const override;
    
    // Get art::Ptr<> to simb::MCTruth and related information
    const art::Ptr<simb::MCTruth>&                TrackIDToMCTruth(int const id)                        const override;
    const art::Ptr<simb::MCTruth>&                ParticleToMCTruth(const simb::MCParticle* p)           const override;
    std::vector<const simb::MCParticle*>          MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const override;
    const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector()                                        const override;
    
    // this method will return the Geant4 track IDs of
    // the particles contributing ionization electrons to the identified hit
    std::vector<sim::TrackIDE> HitToTrackID(detinfo::DetectorClocksData const& clockData,
                                            recob::Hit const& hit)           const override;
    std::vector<sim::TrackIDE> HitToTrackID(detinfo::DetectorClocksData const& clockData,
                                            art::Ptr<recob::Hit> const& hit) const override;
    
    // method to return a subset of allhits that are matched to a list of TrackIDs
    std::vector<std::vector<art::Ptr<recob::Hit>>> TrackIDsToHits(detinfo::DetectorClocksData const& clockData,
                                                                  std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                  std::vector<int> const& tkIDs) const override;

    // method to return the EveIDs of particles contributing ionization
    // electrons to the identified hit
    std::vector<sim::TrackIDE> HitToEveID(detinfo::DetectorClocksData const& clockData,
                                          art::Ptr<recob::Hit>        const& hit) const override;

    // method to return the XYZ position of the weighted average energy deposition for a given hit
    std::vector<double> HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                 art::Ptr<recob::Hit>        const& hit) const override;

    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    std::vector<double> SpacePointToXYZ(detinfo::DetectorClocksData const&,
                                        art::Ptr<recob::SpacePoint> const& spt,
                                        art::Event                  const& evt,
                                        std::string                 const& label) const override;

    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    std::vector<double> SpacePointHitsToXYZ(detinfo::DetectorClocksData const& clockData,
                                            art::PtrVector<recob::Hit>  const& hits) const override;

    // method to return the fraction of hits in a collection that come from the specified Geant4 track ids
    double HitCollectionPurity(detinfo::DetectorClocksData         const& clockData,
                               std::set<int>                       const& trackIDs,
                               std::vector< art::Ptr<recob::Hit> > const& hits) const override;

    // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
    // represented in a collection of hits
    double HitCollectionEfficiency(detinfo::DetectorClocksData         const& clockData,
                                   std::set<int>                       const& trackIDs,
                                   std::vector< art::Ptr<recob::Hit> > const& hits,
                                   std::vector< art::Ptr<recob::Hit> > const& allhits,
                                   geo::View_t                         const  view) const override;
    
    // method to return the fraction of charge in a collection that come from the specified Geant4 track ids
    double HitChargeCollectionPurity(detinfo::DetectorClocksData         const& clockData,
                                     std::set<int>                       const& trackIDs,
                                     std::vector< art::Ptr<recob::Hit> > const& hits) const override;

    // method to return the fraction of all charge in an event from a specific set of Geant4 track IDs that are
    // represented in a collection of hits
    double HitChargeCollectionEfficiency(detinfo::DetectorClocksData         const& clockData,
                                         std::set<int>                       const& trackIDs,
                                         std::vector< art::Ptr<recob::Hit> > const& hits,
                                         std::vector< art::Ptr<recob::Hit> > const& allhits,
                                         geo::View_t                         const  view) const override;
    
    // method to return all EveIDs corresponding to the current sim::ParticleList
    std::set<int> GetSetOfEveIDs() const override;
    
    // method to return all TrackIDs corresponding to the current sim::ParticleList
    std::set<int> GetSetOfTrackIDs() const override;
    
    // method to return all EveIDs corresponding to the given list of hits
    std::set<int> GetSetOfEveIDs(detinfo::DetectorClocksData const& clockData,
                                 std::vector< art::Ptr<recob::Hit> > const& hits) const override;

    // method to return all TrackIDs corresponding to the given list of hits
    std::set<int> GetSetOfTrackIDs(detinfo::DetectorClocksData const& clockData,
                                   std::vector< art::Ptr<recob::Hit> > const& hits) const override;

private:
    
    // Fcl parameters.
    art::InputTag fTrackProducerLabel; ///< tag for finding the tracks
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*                      fGeometry;             ///< pointer to Geometry service
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
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;
    return partInventory->ParticleList();
}

    
//----------------------------------------------------------------------
const simb::MCParticle* BackTrackerTruth::TrackIDToParticle(int const id) const
{
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;
    return partInventory->TrackIdToParticle_P(id);
}
    
//----------------------------------------------------------------------
const simb::MCParticle* BackTrackerTruth::TrackIDToMotherParticle(int const id) const
{
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;
    return partInventory->TrackIdToMotherParticle_P(id);
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& BackTrackerTruth::TrackIDToMCTruth(int const id) const
{
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;
    return partInventory->TrackIdToMCTruth_P(id);
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& BackTrackerTruth::ParticleToMCTruth(const simb::MCParticle* p) const
{
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;
    return partInventory->ParticleToMCTruth_P(p);
}
    
//----------------------------------------------------------------------
std::vector<const simb::MCParticle*> BackTrackerTruth::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
{
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;
    return partInventory->MCTruthToParticles_Ps(mct);
}
    
//----------------------------------------------------------------------
const std::vector< art::Ptr<simb::MCTruth> >&  BackTrackerTruth::MCTruthVector() const
{
    art::ServiceHandle<cheat::ParticleInventoryService> partInventory;
    return partInventory->MCTruthVector_Ps();
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> BackTrackerTruth::HitToTrackID(detinfo::DetectorClocksData const& clockData,
                                                          recob::Hit const& hit) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->HitToTrackIDEs(clockData, hit);
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> BackTrackerTruth::HitToTrackID(detinfo::DetectorClocksData const& clockData,
                                                          art::Ptr<recob::Hit> const& hit) const
{
    return HitToTrackID(clockData, *hit);
}

//----------------------------------------------------------------------
std::vector<std::vector<art::Ptr<recob::Hit>>> BackTrackerTruth::TrackIDsToHits(detinfo::DetectorClocksData const& clockData,
                                                                                std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                                std::vector<int> const& tkIDs) const
{
    std::vector<std::vector<art::Ptr<recob::Hit>>> tkIDsToHitsVec;
    art::ServiceHandle<cheat::BackTrackerService>  backTracker;

    for(const auto& tkID : tkIDs)
    {
        std::vector<art::Ptr<recob::Hit>> hitVec = backTracker->TrackIdToHits_Ps(clockData, tkID, allhits);
        tkIDsToHitsVec.push_back(hitVec);
    }

    return tkIDsToHitsVec;
}
    
//----------------------------------------------------------------------
// plist is assumed to have adopted the appropriate EveIdCalculator prior to
// having been passed to this method. It is likely that the EmEveIdCalculator is
// the one you always want to use
std::vector<sim::TrackIDE> BackTrackerTruth::HitToEveID(detinfo::DetectorClocksData const& clockData,
                                                        art::Ptr<recob::Hit> const& hit) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->HitToEveTrackIDEs(clockData, hit);
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfEveIDs() const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->GetSetOfEveIds();
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfTrackIDs() const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->GetSetOfTrackIds();
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfEveIDs(detinfo::DetectorClocksData const& clockData,
                                               std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->GetSetOfEveIds(clockData, hits);
}
    
//----------------------------------------------------------------------
std::set<int> BackTrackerTruth::GetSetOfTrackIDs(detinfo::DetectorClocksData const& clockData,
                                                 std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->GetSetOfTrackIds(clockData, hits);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitCollectionPurity(detinfo::DetectorClocksData const& clockData,
                                             std::set<int> const& trackIDs,
                                             std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->HitCollectionPurity(clockData, trackIDs, hits);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitChargeCollectionPurity(detinfo::DetectorClocksData const& clockData,
                                                   std::set<int> const& trackIDs,
                                                   std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->HitChargeCollectionPurity(clockData, trackIDs, hits);
}
    
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitCollectionEfficiency(detinfo::DetectorClocksData         const& clockData,
                                                 std::set<int>                       const& trackIDs,
                                                 std::vector< art::Ptr<recob::Hit> > const& hits,
                                                 std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                 geo::View_t                         const  view) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->HitCollectionEfficiency(clockData, trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
double BackTrackerTruth::HitChargeCollectionEfficiency(detinfo::DetectorClocksData         const& clockData,
                                                       std::set<int>                       const& trackIDs,
                                                       std::vector< art::Ptr<recob::Hit> > const& hits,
                                                       std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                       geo::View_t                         const  view) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->HitChargeCollectionEfficiency(clockData, trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                               art::Ptr<recob::Hit> const& hit) const
{
    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->HitToXYZ(clockData, hit);
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::SpacePointToXYZ(detinfo::DetectorClocksData         const&,
                                                      art::Ptr<recob::SpacePoint>         const& spt,
                                                      art::Event                          const& evt,
                                                      std::string                         const& label) const
{
    std::vector<double> hitVec = {0.,0.,0.};

    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    //backTracker->SpacePointHitsToWeightedXYZ(spt, evt, label);

    return hitVec;
}
    
//----------------------------------------------------------------------
std::vector<double> BackTrackerTruth::SpacePointHitsToXYZ(detinfo::DetectorClocksData const& clockData,
                                                          art::PtrVector<recob::Hit>  const& hits) const
{
    std::vector<art::Ptr<recob::Hit>> hitVec;
    for(size_t idx=0; idx<hits.size(); idx++) hitVec.push_back(hits.at(idx));

    art::ServiceHandle<cheat::BackTrackerService> backTracker;
    return backTracker->SpacePointHitsToWeightedXYZ(clockData, hitVec);
}

//----------------------------------------------------------------------------
    
DEFINE_ART_CLASS_TOOL(BackTrackerTruth)
}
