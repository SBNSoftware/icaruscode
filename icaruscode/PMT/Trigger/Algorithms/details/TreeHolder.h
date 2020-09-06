/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/TreeHolder.h
 * @brief  Class holding a ROOT tree, to be shared by other classes.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_TREEHOLDER_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_TREEHOLDER_H


// ROOT libraries
#include "TTree.h"


//------------------------------------------------------------------------------
namespace icarus::trigger::details { struct TreeHolder; }

/**
 * @brief Simple class holding a tree.
 * 
 * To be shared.
 */
struct icarus::trigger::details::TreeHolder {

  TreeHolder() = default;
  TreeHolder(TTree& tree): fTree(&tree) {}

  TTree& tree() { return *fTree; }
  TTree const& tree() const { return *fTree; }

    private:
  TTree* fTree = nullptr;

}; // struct TreeHolder


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_TREEHOLDER_H
