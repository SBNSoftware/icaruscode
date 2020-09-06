/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/EventIDTree.cxx
 * @brief  Class storing event information in a ROOT tree (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 * @see    icaruscode/PMT/Trigger/Algorithms/details/EventIDTree.h
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/details/EventIDTree.h"

// ROOT libraries
#include "TTree.h"


//------------------------------------------------------------------------------
icarus::trigger::details::EventIDTree::EventIDTree(TTree& tree)
  : TreeHolder(tree)
{

  this->tree().Branch("RunNo", &fRunNo);
  this->tree().Branch("SubRunNo", &fSubRunNo);
  this->tree().Branch("EventNo", &fEventNo);

} // icarus::trigger::details::EventIDTree::EventIDTree()


//------------------------------------------------------------------------------
void icarus::trigger::details::EventIDTree::assignID(art::EventID const& id) {

  fRunNo = id.run();
  fSubRunNo = id.subRun();
  fEventNo = id.event();

} // icarus::trigger::details::EventIDTree::assignID()


//------------------------------------------------------------------------------
