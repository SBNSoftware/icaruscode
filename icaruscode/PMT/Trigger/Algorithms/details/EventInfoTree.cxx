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
icarus::trigger::details::EventInfoTree::EventInfoTree
  (TTree& tree, bool fillGen /* = true */, bool fillEDep /* = true */)
  : TreeHolder(tree), fDoGen(fillGen), fDoEDep(fillEDep)
{
  if (fDoGen) {
    this->tree().Branch("CC",       &fCC);
    this->tree().Branch("NC",       &fNC);
    this->tree().Branch("IntType",  &fIntType);
    this->tree().Branch("Time",     &fTime);
    this->tree().Branch("NuE",      &fNuE);
    this->tree().Branch("OutLeptE", &fOutLeptE);
    this->tree().Branch("OutLeptAngle", &fOutLeptAngle);
//  this->tree().Branch("IsNu_mu", &fIsNu_mu);
//  this->tree().Branch("IsNu_e", &fIsNu_e);
  } // if generated
  
  if (fDoEDep) {
    this->tree().Branch("TotE",            &fTotE);
    this->tree().Branch("SpillE",          &fSpillE);
    this->tree().Branch("PreSpillE",       &fPreSpillE);
    this->tree().Branch("ActiveE",         &fActiveE);
    this->tree().Branch("SpillActiveE",    &fSpillActiveE);
    this->tree().Branch("PreSpillActiveE", &fPreSpillActiveE);
  } // if energy deposition
  
  if (fDoGen) {
    this->tree().Branch("InActive",        &fInActive);
    this->tree().Branch("Vertices",        &fVertices);
  } // if generated
  
} // icarus::trigger::details::EventInfoTree::EventInfoTree()


//------------------------------------------------------------------------------
void icarus::trigger::details::EventInfoTree::assignEvent
  (EventInfo_t const& info)
{
  
  if (fDoGen) {
    fCC       = info.nWeakChargedCurrentInteractions();
    fNC       = info.nWeakNeutralCurrentInteractions();
    fIntType  = info.InteractionType();
    fTime     = static_cast<Double_t>(info.InteractionTime());
    fNuE      = static_cast<Double_t>(info.NeutrinoEnergy());
    fOutLeptE = static_cast<Double_t>(info.LeptonEnergy());
    fOutLeptAngle = info.LeptonAngle();
 // fIsNu_mu = static_cast<Bool_t>(info.isNu_mu());
 // fIsNu_e = static_cast<Bool_t>(info.isNu_e());

  } // if generated
  
  if (fDoEDep) {
    fTotE         = static_cast<Double_t>(info.DepositedEnergy());
    fSpillE       = static_cast<Double_t>(info.DepositedEnergyInSpill());
    fPreSpillE    = static_cast<Double_t>(info.DepositedEnergyInPreSpill());
    fActiveE      = static_cast<Double_t>(info.DepositedEnergyInActiveVolume());
    fSpillActiveE
      = static_cast<Double_t>(info.DepositedEnergyInSpillInActiveVolume());
    fPreSpillActiveE
      = static_cast<Double_t>(info.DepositedEnergyInPreSpillInActiveVolume());
  } // if energy deposition
  
  if (fDoGen) {
    fInActive     = static_cast<Bool_t>(info.isInActiveVolume());
    fVertices     = info.Vertices();
  } // if generated
  
} // icarus::trigger::details::EventInfoTree::assignEvent()


//------------------------------------------------------------------------------
