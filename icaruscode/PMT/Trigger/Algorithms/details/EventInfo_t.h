/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h
 * @brief  Class hosting selected information about the event.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 * @see    icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFO_T_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFO_T_H


// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time
#include "lardataalg/Utilities/quantities/energy.h" // gigaelectronvolt
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t

// C/C++ standard libraries
#include <vector>
#include <ostream>
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
namespace icarus::trigger::details {
  
  struct EventInfo_t;
  
  std::ostream& operator<< (std::ostream& out, EventInfo_t const& info);
  
} // namespace icarus::trigger::details


//------------------------------------------------------------------------------
/**
 * @brief Selected information about the event.
 *
 * This is intended as a record for transferring relevant event information
 * between object.
 */
struct icarus::trigger::details::EventInfo_t {
  
  /// @name Local type aliases
  /// @{
  using GeV = util::quantities::gigaelectronvolt;
  using simulation_time = detinfo::timescales::simulation_time;
  /// @}
  
  /// Constructor. As if nobody noticed.
  EventInfo_t() { fInteractions.fill(0U); }
  
  // --- BEGIN Query interface -----------------------------------------------
  /// @name Query interface
  /// @{

  /// Returns the number of weak charged current interactions in the event.
  unsigned int nWeakChargedCurrentInteractions() const
    { return fInteractions[itWCC]; }

  /// Returns the number of weak neutral current interactions in the event.
  unsigned int nWeakNeutralCurrentInteractions() const
    { return fInteractions[itWNC]; }

  /// Returns the number of weak current interactions in the event.
  unsigned int nWeakCurrentInteractions() const
    {
      return
        nWeakChargedCurrentInteractions() + nWeakNeutralCurrentInteractions();
    }

  /// Returns whether the event is generated as a neutrino CC interaction.
  bool isWeakChargedCurrent() const
    { return nWeakChargedCurrentInteractions() > 0U; }

  /// Returns whether the event is generated as a neutrino NC interaction.
  bool isWeakNeutralCurrent() const
    { return nWeakNeutralCurrentInteractions() > 0U; }

  /// Returns whether the event is generated as a neutrino interaction.
  bool isNeutrino() const
    { return isWeakChargedCurrent() || isWeakNeutralCurrent(); }
  
  /// Returns which neutrino flavor is present in an event
  bool isNu_mu() const { return nu_mu; }
  bool isNu_e() const { return nu_e; }

  /// Returns the neutrino PDG code
  int NeutrinoPDG() const { return fNeutrinoPDG; }

  /// Returns the total energy deposited in the detector during the event [GeV]
  GeV DepositedEnergy() const { return fEnergyDepTotal; }
  
  /// Returns the total energy deposited in the detector during beam [GeV]
  GeV DepositedEnergyInSpill() const { return fEnergyDepSpill; }

  /// Returns the energy deposited in the active volume during the event [GeV]
  GeV DepositedEnergyInActiveVolume() const { return fEnergyDepActive; }
  
  /// Returns the energy deposited in the active volume during the beam [GeV]
  GeV DepositedEnergyInSpillInActiveVolume() const
    { return fEnergyDepSpillActive; }

  /// Returns the neutrino energy [GeV]
  GeV NeutrinoEnergy() const { return fNeutrinoEnergy; }

  /// Returns the lepton energy [GeV]
  GeV LeptonEnergy() const { return fLeptonEnergy; }

  /// Returns the interaction type
  int InteractionType() const { return fInteractionType; }

  /// Returns the time of the first interaction, in simulation time scale [ns]
  simulation_time InteractionTime() const { return fInteractionTime; }
  
  /// Returns whether this type of event has a known vertex.
  bool hasVertex() const { return !fVertices.empty(); }
  
  /// Returns the number of known interaction vertices.
  unsigned int nVertices() const { return fVertices.size(); }
  
  /// Returns the list of a known interaction vertex.
  std::vector<geo::Point_t> const& Vertices() const { return fVertices; }
  
  /// Returns whether there is an interaction within the active volume.
  bool isInActiveVolume() const { return fInActiveVolume; }
  
  /// @}
  // --- END Query interface -------------------------------------------------


  // --- BEGIN Set interface -------------------------------------------------
  /// @name Set interface
  /// @{

  /// Marks this event as including _n_ more weak charged current interactions.
  void AddWeakChargedCurrentInteractions(unsigned int n = 1U)
    { fInteractions[itWCC] += n; }

  /// Marks this event as including _n_ more weak neutral current interactions.
  void AddWeakNeutralCurrentInteractions(unsigned int n = 1U)
    { fInteractions[itWNC] += n; }

  /// Marks the flavor of the neutrino in the first interaction.
  void SetNu_mu(bool numu) { nu_mu = numu; }
  void SetNu_e(bool nue) { nu_e = nue; }

  /// Marks the neutrino type of the first interaction in the event.
  void SetNeutrinoPDG(int NU) { fNeutrinoPDG = NU; }

  /// Sets the total deposited energy of the event [GeV]
  void SetDepositedEnergy(GeV e) { fEnergyDepTotal = e; }

  /// Sets the energy of the event deposited during beam gate [GeV]
  void SetDepositedEnergyInSpill(GeV e) { fEnergyDepSpill = e; }

  /// Sets the total deposited energy of the event in active volume [GeV]
  void SetDepositedEnergyInActiveVolume(GeV e) { fEnergyDepActive = e; }

  /// Sets the energy of the event deposited during beam gate in active volume
  /// [GeV]
  void SetDepositedEnergyInSpillInActiveVolume(GeV e)
    { fEnergyDepSpillActive = e; }

  /// Sets the neutrino energy.
  void SetNeutrinoEnergy(GeV eNu) { fNeutrinoEnergy = eNu; }

  /// Sets the lepton energy.
  void SetLeptonEnergy(GeV eL) { fLeptonEnergy = eL; }

  /// Sets the interaction type
  void SetInteractionType(int type) { fInteractionType = type; }

  /// Sets the time of the first interaction.
  void SetInteractionTime(simulation_time time) { fInteractionTime = time; }
  
  /// Set whether the event has relevant activity in the active volume.
  void SetInActiveVolume(bool active = true) { fInActiveVolume = active; }
  
  /// Adds a point to the list of interaction vertices in the event.
  void AddVertex(geo::Point_t const& vertex) { fVertices.push_back(vertex); }
  
  /// Inserts a point in the specified position of the list of interaction
  /// vertices in the event.
  void InsertVertex(geo::Point_t const& vertex, std::size_t beforeIndex)
    { fVertices.insert(next(begin(fVertices), beforeIndex), vertex); }
  
  /// @}
  // --- END Set interface ---------------------------------------------------

  /// Prints the content of the object into a stream.
  void dump(std::ostream& out) const;

    private:

  // --- BEGIN interaction type constants ------------------------------------
  
  static constexpr std::size_t itWCC { 0U }; ///< Charged weak current.
  static constexpr std::size_t itWNC { 1U }; ///< Neutral weak current.
  static constexpr std::size_t NInteractionTypes { 2U };

  // --- END interaction type constants --------------------------------------

  std::array<unsigned int, NInteractionTypes> fInteractions;

  GeV fEnergyDepTotal       { 0.0 }; ///< Total deposited energy.
  GeV fEnergyDepSpill       { 0.0 }; ///< Energy deposited in spill.
  GeV fEnergyDepActive      { 0.0 }; ///< Energy deposited in active volume.
  /// Energy deposited in active volume in spill.
  GeV fEnergyDepSpillActive { 0.0 };

  int fNeutrinoPDG { 0 };
  int fInteractionType { 0 };
  
  simulation_time fInteractionTime; ///< Time of the first interaction [ns]

  GeV fNeutrinoEnergy { 0.0 };
  GeV fLeptonEnergy { 0.0 };
  //GeV fNucleonEnergy { 0.0 };
  
  bool nu_mu { false };
  bool nu_e { false };
  
  /// Whether the event has activity inside the active volume.
  bool fInActiveVolume { false };
  
  std::vector<geo::Point_t> fVertices; ///< Position of all vertices.

}; // struct icarus::trigger::details::EventInfo_t

inline std::ostream& icarus::trigger::details::operator<<
  (std::ostream& out, EventInfo_t const& info)
  { info.dump(out); return out; }


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFO_T_H
