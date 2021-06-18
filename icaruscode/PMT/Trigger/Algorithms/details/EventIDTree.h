/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/EventIDTree.h
 * @brief  Class hosting selected information about the event.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTIDTREE_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTIDTREE_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/details/TreeHolder.h"

// framework libraries
#include "canvas/Persistency/Provenance/EventID.h"

// ROOT libraries
#include "Rtypes.h"


//------------------------------------------------------------------------------
namespace icarus::trigger::details { struct EventIDTree; }
/**
 * @brief Class managing the serialization of event ID in a simple ROOT tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assignEvent()`, the branch addresses are assigned the values from the
 * event ID. The tree is not `Fill()`-ed.
 *
 * The tree structure is: `Run/i:SubRun/i:Event/i`, with a single branch per
 * element.
 */
struct icarus::trigger::details::EventIDTree: public TreeHolder {

  /// Creates the required branches and assigns addresses to them.
  EventIDTree(TTree& tree);

  /// Fills the information of the specified event.
  void assignID(art::EventID const& id);

  UInt_t fRunNo;
  UInt_t fSubRunNo;
  UInt_t fEventNo;

}; // struct icarus::trigger::details::EventIDTree


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTIDTREE_H
