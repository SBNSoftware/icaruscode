/**
 * @file   icaruscode/IcarusObj/SimEnergyDepositSummary.h
 * @brief  Object storing a summary of energy depositions in the detector.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * 
 */

#ifndef ICARUSCODE_ICARUSOBJ_SIMENERGYDEPOSITSUMMARY_H
#define ICARUSCODE_ICARUSOBJ_SIMENERGYDEPOSITSUMMARY_H

// -----------------------------------------------------------------------------
namespace icarus { struct SimEnergyDepositSummary; }

/// Data structure containing summary information about deposited energy.
struct icarus::SimEnergyDepositSummary {
  
  float Total       { 0.0 }; ///< Total deposited energy [GeV]
  float Spill       { 0.0 }; ///< Energy deposited in spill [GeV]
  float PreSpill    { 0.0 }; ///< Energy deposited in pre-spill [GeV]
  float Active      { 0.0 }; ///< Energy deposited in active volume [GeV]
  /// Energy deposited in active volume in spill [GeV]
  float SpillActive { 0.0 };
  /// Energy deposited in active volume in pre-spill window [GeV]
  float PreSpillActive { 0.0 };
  
}; // icarus::SimEnergyDepositSummary


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_ICARUSOBJ_SIMENERGYDEPOSITSUMMARY_H
