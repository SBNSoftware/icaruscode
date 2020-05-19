/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/EventInfoTree.h
 * @brief  Class writing event information into a ROOT tree.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 * @see    icaruscode/PMT/Trigger/Algorithms/details/EventInfoTree.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFOTREE_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFOTREE_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/details/TreeHolder.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// ROOT libraries
#include "Rtypes.h"

// C/C++ standard libraries
#include <vector>


// -----------------------------------------------------------------------------
namespace icarus::trigger::details { struct EventInfoTree; }

/**
 * @brief Class managing the serialization of event information in a simple ROOT
 *        tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assignEvent()`, the branch addresses are assigned the values from the
 * event information. The tree is not `Fill()`-ed.
 *
 * The tree structure is:
 * `CC/i:NC/i:IntType/I:Time/D:NuE/D:OutLeptE/D:TotE/D:SpillE/D:InActive/O`,
 * with a single branch per element.
 *
 * Branches:
 *  * `CC` (unsigned integer): number of neutrino CC interactions in the event
 *  * `NC` (unsigned integer): number of neutrino NC interactions in the event
 *  * `IntType` (integer): code of interaction type (see `simb::int_type_`)
 *  * `Time` (double): time of the interaction in simulation time scale
 *  * `NuE` (double): energy of the generated initial state neutrino [GeV]
 *  * `OutLeptE` (double): energy of the generated final state lepton [GeV]
 *  * `TotE` (double): total deposited energy in the event [GeV]
 *  * `SpillE` (double): total deposited energy during the beam gate [GeV]
 *  * `InActive` (bool): whether an interaction happened in active volume'
 *      this requires an interaction vertex (e.g. cosmic rays are out)
 *  * `NVertices` (unsigned integer): number of interaction vertices in event
 *  * `Vertices_` (list of points): the location of all the interaction vertices
 *    in the event; it's a vector of GenVector 3D points (can access coordinates
 *    as `Vertices.X()` or `Vertices.fCoordinates.fX`)
 *
 */
struct icarus::trigger::details::EventInfoTree: public TreeHolder {

  /// Creates the required branches and assigns addresses to them.
  EventInfoTree(TTree& tree);

  /**
   * @brief Fills the information of the specified event.
   * @param info event information to fill the tree with
   * @param inPlots whether this event is plotted (as opposed to filtered out)
   */
  void assignEvent(EventInfo_t const& info);

  UInt_t fCC;
  UInt_t fNC;
  Int_t fIntType;
  Double_t fTime;
  Double_t fNuE;
  Double_t fOutLeptE;
  Double_t fTotE;
  Double_t fSpillE;
  Double_t fActiveE;
  Double_t fSpillActiveE;
  UInt_t fNVertices;
  std::vector<geo::Point_t> fVertices; // is this ROOT tree friendly?
  
  Bool_t fInActive;
  
}; // struct icarus::trigger::details::EventInfoTree


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFOTREE_H
