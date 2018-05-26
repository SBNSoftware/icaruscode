/**
 * @file    MCTruthAssociations.h
 * @brief   This algorithm attempts to decode Track and Hit <--> MCParticle assocations
 * @author  Tracy Usher (usher@slac.stanford.edu)
 * @date    October 25, 2017
 * @see     galleryAnalysis.cpp
 * 
 */

#ifndef MCTruthAssociations_H
#define MCTruthAssociations_H

// LArSoft libraries
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
//#include "lardataobj/Simulation/SimChannel.h"   // need this for TrackIDE definition
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindOneP.h"

// nutools
#include "icaruscode/gallery/MCTruthBase/MCTruthParticleList.h"

// canvas libraries
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <vector>
#include <memory> // std::unique_ptr<>

namespace truth
{
using MCTruthTruthVec             = std::vector<art::Ptr<simb::MCTruth>>;
using MCParticleVec               = std::vector<art::Ptr<simb::MCParticle>>;
using MCTruthAssns                = art::FindOneP<simb::MCTruth>;
using HitParticleAssociations     = art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>;
using HitParticleAssociationsVec  = std::vector<const HitParticleAssociations*>;
using MCTruthParticleAssociations = art::Assns<simb::MCTruth, simb::MCParticle, void>;

// Definition of TrackIDE here to avoid pulling in SimChannel.h
/// Ionization energy from a Geant4 track
struct TrackIDE{
    int   trackID;      ///< Geant4 supplied trackID
    float energyFrac;   ///< fraction of hit energy from the particle with this trackID
    float energy;       ///< energy from the particle with this trackID [MeV]
    float numElectrons; ///< number of electrons from the particle detected on the wires
};

/**
 * @brief Obtains truth matching by using hit <--> MCParticle associations
 * 
 * Configuration
 * --------------
 * 
 */
class MCTruthAssociations
{
public:
  
    MCTruthAssociations(fhicl::ParameterSet const& config);
  
    void setup(const HitParticleAssociationsVec&,
               const MCParticleVec&,
               const MCTruthAssns&,
               const geo::GeometryCore&,
               const detinfo::DetectorProperties&);
    
    const MCTruthParticleList& getParticleList() const;

    // Return a pointer to the simb::MCParticle object corresponding to
    // the given TrackID
    const simb::MCParticle* TrackIDToParticle(int const& id)       const;
    const simb::MCParticle* TrackIDToMotherParticle(int const& id) const;
    
    // Get art::Ptr<> to simb::MCTruth and related information
    const art::Ptr<simb::MCTruth>&       TrackIDToMCTruth(int const& id)                        const;
    const art::Ptr<simb::MCTruth>&       ParticleToMCTruth(const simb::MCParticle* p)           const;
    std::vector<const simb::MCParticle*> MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const;
    const MCTruthTruthVec&               MCTruthVector()                                        const;
    
    // this method will return the Geant4 track IDs of
    // the particles contributing ionization electrons to the identified hit
    std::vector<TrackIDE> HitToTrackID(const recob::Hit*)           const;
    std::vector<TrackIDE> HitToTrackID(art::Ptr<recob::Hit> const&) const;

    // method to return a subset of allhits that are matched to a list of TrackIDs
    const std::vector<std::vector<art::Ptr<recob::Hit>>> TrackIDsToHits(std::vector<art::Ptr<recob::Hit>> const&,
                                                                        std::vector<int> const&) const;

    // method to return the EveIDs of particles contributing ionization
    // electrons to the identified hit
    std::vector<TrackIDE> HitToEveID(art::Ptr<recob::Hit> const& hit) const;

    // method to return the XYZ position of the weighted average energy deposition for a given hit
    std::vector<double>  HitToXYZ(art::Ptr<recob::Hit> const& hit) const;
    
    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    std::vector<double> SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits) const;

    // method to return the fraction of hits in a collection that come from the specified Geant4 track ids
    double HitCollectionPurity(std::set<int>,
                               std::vector< art::Ptr<recob::Hit> > const&) const;
    
    // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
    // represented in a collection of hits
    double HitCollectionEfficiency(std::set<int>,
                                   std::vector< art::Ptr<recob::Hit> > const&,
                                   std::vector< art::Ptr<recob::Hit> > const&,
                                   geo::View_t                         const&) const;
    
    // method to return the fraction of charge in a collection that come from the specified Geant4 track ids
    double HitChargeCollectionPurity(std::set<int>,
                                     std::vector< art::Ptr<recob::Hit> > const&) const;
    
    // method to return the fraction of all charge in an event from a specific set of Geant4 track IDs that are
    // represented in a collection of hits
    double HitChargeCollectionEfficiency(std::set<int>,
                                         std::vector< art::Ptr<recob::Hit> > const&,
                                         std::vector< art::Ptr<recob::Hit> > const&,
                                         geo::View_t                         const&) const;
    
    // method to return all EveIDs corresponding to the current sim::ParticleList
    std::set<int> GetSetOfEveIDs() const;
    
    // method to return all TrackIDs corresponding to the current sim::ParticleList
    std::set<int> GetSetOfTrackIDs() const;
    
    // method to return all EveIDs corresponding to the given list of hits
    std::set<int> GetSetOfEveIDs(std::vector< art::Ptr<recob::Hit> > const&) const;
    
    // method to return all TrackIDs corresponding to the given list of hits
    std::set<int> GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const&) const;

private:
    
    // Declare the containers for the basic maps
    using HitMatchDataPair  = std::pair<const recob::Hit*,const anab::BackTrackerHitMatchingData*>;
    using PartMatchDataPair = std::pair<const simb::MCParticle*,const anab::BackTrackerHitMatchingData*>;
    using HitToPartVecMap   = std::map<const recob::Hit*,std::set<PartMatchDataPair>>;
    using PartToHitVecMap   = std::map<const simb::MCParticle*, std::set<HitMatchDataPair>>;
    using MCTruthTrackIDMap = std::unordered_map<int, art::Ptr<simb::MCTruth>>;

    // Must allow for the case of multiple instances of hit <--> MCParticle associations
    // You ask "why do it this way? Can't these all be in a single set of containers?"
    // The answer is no because you want to avoid multiple counting
    struct HitPartAssnsStruct
    {
        // Declare containers as member variables
        HitToPartVecMap     fHitToPartVecMap;       ///< Mapping from hits to associated MCParticle/data pairs
        PartToHitVecMap     fPartToHitVecMap;       ///< Mapping from MCParticle to associated hit/data pairs
    };
    
    using HitPartAssnsList = std::list<HitPartAssnsStruct>;

    int    calculateEvdID(int) const;
    double length(const recob::Track*) const;
    double length(const simb::MCParticle& part, double dx,
                  TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
                  unsigned int tpc = 0, unsigned int cstat = 0) const;
    
    HitPartAssnsList                   fHitPartAssnsVec;       ///< Container for the (multiple) associations
    MCTruthParticleList                fParticleList;          ///< ParticleList to map track ID to
    MCTruthTruthVec                    fMCTruthVec;            ///< all the MCTruths for the event
    MCTruthTrackIDMap                  fTrackIDToMCTruthIndex; ///< map of track ids to MCTruthList entry

    float                              fMinHitEnergyFraction;  ///< minimum fraction of energy a track id has to
                                                               ///< contribute to a hit to be counted in
                                                               ///< purity and efficiency calculations
                                                               ///< based on hit collections

    geo::GeometryCore const*           fGeometry           = nullptr;
    const detinfo::DetectorProperties* fDetectorProperties = nullptr;   ///< Detector properties service
}; // class MCTruthAssociations
    
}  // End of namespace

#endif // MCTruthAssociations_H
