////////////////////////////////////////////////////////////////////////
///
/// \file  IMCTruthMatching.h
/// \brief This provides an interface which defines truth matching functions
///        made available to downstream analysis code
///
/// \author  usher@slac.stanford.edu
///
////////////////////////////////////////////////////////////////////////
#ifndef IMCTRUTHMATCHING_H
#define IMCTRUTHMATCHING_H

#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nug4/ParticleNavigation/EveIdCalculator.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "lardataobj/RecoBase/Hit.h"

namespace detinfo {
  class DetectorClocksData;
}

namespace recob{
  class SpacePoint;
}

///code to link reconstructed objects back to the MC truth information
namespace truth
{
    
class IMCTruthMatching
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IMCTruthMatching() noexcept = default;

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;
    
    /**
     *  @brief This rebuilds the internal maps
     */
    virtual void Rebuild(const art::Event& evt) = 0;

    /**
     *  @brief Get a reference to the ParticleList
     */
    virtual const sim::ParticleList& ParticleList() const = 0;

//    // Set the EveIdCalculator for the owned ParticleList
//    void  SetEveIdCalculator(sim::EveIdCalculator *ec) { fParticleList.AdoptEveIdCalculator(ec); }

    // Return a pointer to the simb::MCParticle object corresponding to
    // the given TrackID
    virtual const simb::MCParticle* TrackIDToParticle(int id)       const = 0;
    virtual const simb::MCParticle* TrackIDToMotherParticle(int id) const = 0;

    // Get art::Ptr<> to simb::MCTruth and related information
    virtual const art::Ptr<simb::MCTruth>&                TrackIDToMCTruth(int id)                        const = 0;
    virtual const art::Ptr<simb::MCTruth>&                ParticleToMCTruth(const simb::MCParticle* p)           const = 0;
    virtual std::vector<const simb::MCParticle*>          MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const = 0;
    virtual const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector()                                        const = 0;

    // this method will return the Geant4 track IDs of 
    // the particles contributing ionization electrons to the identified hit
    virtual std::vector<sim::TrackIDE> HitToTrackID(detinfo::DetectorClocksData const& clockData,
                                                    recob::Hit const& hit)           const = 0;
    virtual std::vector<sim::TrackIDE> HitToTrackID(detinfo::DetectorClocksData const& clockData,
                                                    art::Ptr<recob::Hit> const& hit) const = 0;
    
    // method to return a subset of allhits that are matched to a list of TrackIDs
    virtual std::vector<std::vector<art::Ptr<recob::Hit>>> TrackIDsToHits(detinfo::DetectorClocksData const& clockData,
                                                                          std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                          std::vector<int> const& tkIDs) const = 0;

    // method to return the EveIDs of particles contributing ionization
    // electrons to the identified hit
    virtual std::vector<sim::TrackIDE> HitToEveID(detinfo::DetectorClocksData const& clockData,
                                                  art::Ptr<recob::Hit> const& hit) const = 0;

    // method to return the XYZ position of the weighted average energy deposition for a given hit
    virtual std::vector<double> HitToXYZ(detinfo::DetectorClocksData const& clockData,
                                          art::Ptr<recob::Hit> const& hit) const = 0;

    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    virtual std::vector<double> SpacePointToXYZ(detinfo::DetectorClocksData const& clockData,
                                                art::Ptr<recob::SpacePoint> const& spt,
                                                art::Event                  const& evt,
                                                std::string                 const& label) const = 0;

    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    virtual std::vector<double> SpacePointHitsToXYZ(detinfo::DetectorClocksData const& clockData,
                                                    art::PtrVector<recob::Hit>  const& hits) const = 0;

    // method to return the fraction of hits in a collection that come from the specified Geant4 track ids
    virtual double HitCollectionPurity(detinfo::DetectorClocksData         const& clockData,
                                       std::set<int>                       const& trackIDs,
                                       std::vector< art::Ptr<recob::Hit> > const& hits) const = 0;
    
    // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are 
    // represented in a collection of hits
    virtual double HitCollectionEfficiency(detinfo::DetectorClocksData         const& clockData,
                                           std::set<int>                       const& trackIDs,
                                           std::vector< art::Ptr<recob::Hit> > const& hits,
                                           std::vector< art::Ptr<recob::Hit> > const& allhits,
                                           geo::View_t                                view) const = 0;

    // method to return the fraction of charge in a collection that come from the specified Geant4 track ids
    virtual double HitChargeCollectionPurity(detinfo::DetectorClocksData         const& clockData,
                                             std::set<int>                       const& trackIDs,
                                             std::vector< art::Ptr<recob::Hit> > const& hits) const = 0;
    
    // method to return the fraction of all charge in an event from a specific set of Geant4 track IDs that are 
    // represented in a collection of hits
    virtual double HitChargeCollectionEfficiency(detinfo::DetectorClocksData         const& clockData,
                                                 std::set<int>                       const& trackIDs,
                                                 std::vector< art::Ptr<recob::Hit> > const& hits,
                                                 std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                 geo::View_t                         const  view) const = 0;
  
    // method to return all EveIDs corresponding to the current sim::ParticleList
    virtual std::set<int> GetSetOfEveIDs() const = 0;

    // method to return all TrackIDs corresponding to the current sim::ParticleList
    virtual std::set<int> GetSetOfTrackIDs() const = 0;

    // method to return all EveIDs corresponding to the given list of hits
    virtual std::set<int> GetSetOfEveIDs(detinfo::DetectorClocksData const& clockData,
                                         std::vector< art::Ptr<recob::Hit> > const& hits) const = 0;

    // method to return all TrackIDs corresponding to the given list of hits
    virtual std::set<int> GetSetOfTrackIDs(detinfo::DetectorClocksData const& clockData,
                                           std::vector< art::Ptr<recob::Hit> > const& hits) const = 0;
};
    
} // namespace
#endif // IMCTRUTHMATCHING_H
