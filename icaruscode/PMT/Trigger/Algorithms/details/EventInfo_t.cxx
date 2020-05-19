/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.cxx
 * @brief  Class hosting selected information about the event (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 * @see    icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h
 */


// library header
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h"

// LArSoft libraries
#include "lardataalg/MCDumpers/MCDumperUtils.h" // sim::TruthInteractionTypeName

// C/C++ standard libraries
#include <ostream>


//------------------------------------------------------------------------------
void icarus::trigger::details::EventInfo_t::dump(std::ostream& out) const {
  out << "Event contains:";
  if (isNeutrino()) {
    if (nWeakChargedCurrentInteractions())
      out << " " << nWeakChargedCurrentInteractions() << " CC";
    if (nWeakNeutralCurrentInteractions())
      out << " " << nWeakNeutralCurrentInteractions() << " NC";
    if (isNu_mu()) out << " nu_mu";
    if (isNu_e()) out << " nu_e";
    out << "\nThe first neutrino has E=" << NeutrinoEnergy()
      << " and becomes a lepton with E=" << LeptonEnergy()
      << " with a " << sim::TruthInteractionTypeName(InteractionType())
      << " interaction at " << InteractionTime()
      ;
  }
  else {
    out << " no neutrino interaction";
  }
  out << "\nTotal deposited energy: " << DepositedEnergy()
    << ", of which in spill " << DepositedEnergyInSpill()
    << ", in active volume " << DepositedEnergyInActiveVolume()
    << ", in active volume and in spill "
      << DepositedEnergyInSpillInActiveVolume();
  if (fVertices.empty()) {
    out << "\nNo interaction vertex found.";
  }
  else {
    auto iVertex = fVertices.begin();
    auto const vend = fVertices.end();
    out
      << "\n" << fVertices.size() << " interaction vertices: " << *iVertex;
    while (++iVertex != vend) out << "; " << *iVertex;
    out << ".";
  }
  out << "\nThe event is" << (isInActiveVolume()? "": " NOT")
    << " marked as in the active volume of the detector.";
  out << "\n";
  out.flush();
} // icarus::trigger::details::EventInfo_t::dump()


//------------------------------------------------------------------------------
