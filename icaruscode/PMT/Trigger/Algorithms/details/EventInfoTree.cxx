/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/EventInfoTree.cxx
 * @brief  Class writing event information into a ROOT tree (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 * @see    icaruscode/PMT/Trigger/Algorithms/details/EventInfoTree.h
 */


// library header
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfoTree.h"

// ROOT libraries
#include "TTree.h"


//------------------------------------------------------------------------------
//--- icarus::trigger::details::EventInfoTree
//------------------------------------------------------------------------------
icarus::trigger::details::EventInfoTree::EventInfoTree(TTree& tree)
  : TreeHolder(tree)
{

  this->tree().Branch("CC",       &fCC);
  this->tree().Branch("NC",       &fNC);
  this->tree().Branch("IntType",  &fIntType);
  this->tree().Branch("NuE",      &fNuE);
  this->tree().Branch("OutLeptE", &fOutLeptE);
  
  this->tree().Branch("TotE",         &fTotE);
  this->tree().Branch("SpillE",       &fSpillE);
  this->tree().Branch("ActiveE",      &fActiveE);
  this->tree().Branch("SpillActiveE", &fSpillActiveE);
  this->tree().Branch("InActive",     &fInActive);
  this->tree().Branch("Vertices",     &fVertices);
  
} // icarus::trigger::details::EventInfoTree::EventInfoTree()


//------------------------------------------------------------------------------
void icarus::trigger::details::EventInfoTree::assignEvent
  (EventInfo_t const& info)
{

  fCC       = info.nWeakChargedCurrentInteractions();
  fNC       = info.nWeakNeutralCurrentInteractions();
  fIntType  = info.InteractionType();
  fNuE      = static_cast<Double_t>(info.NeutrinoEnergy());
  fOutLeptE = static_cast<Double_t>(info.LeptonEnergy());

  fTotE         = static_cast<Double_t>(info.DepositedEnergy());
  fSpillE       = static_cast<Double_t>(info.DepositedEnergyInSpill());
  fActiveE      = static_cast<Double_t>(info.DepositedEnergyInActiveVolume());
  fSpillActiveE = static_cast<Double_t>(info.DepositedEnergyInSpillInActiveVolume());
  fInActive     = static_cast<Bool_t>(info.isInActiveVolume());
  fVertices     = info.Vertices();
  
} // icarus::trigger::details::EventInfoTree::assignEvent()


//------------------------------------------------------------------------------
