/**
 * @file   icaruscode/PMT/Trigger/Algorithms/FillTriggerGates.h
 * @brief  Utilities for I/O of trigger gate data.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 7, 2019
 * 
 * his library is currently header only.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FILLTRIGGERGATES_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FILLTRIGGERGATES_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h" // convenience

// // LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/Exception.h"

// C/C++ standard libraries
#include <vector>
#include <utility> // std::move()


namespace icarus::trigger {
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Creates a gate object out of trigger gate data products.
   * @tparam GateObject type of gate object to create
   * @tparam TriggerGateData the type of data stored
   * @param gates collection of data from all the gates
   * @param gateToWaveforms association between each gate and its waveforms
   * @return a STL vector of `GateObject` with copy of data from `gates`
   * 
   * Objects like `icarus::trigger::OpticalTriggerGate` are complex enough that
   * they are not saved directly into an _art_ event. Rather, they are diced
   * into pieces and the pieces are stored.
   * This function stitches the pieces and returns back an object like
   * `icarus::trigger::OpticalTriggerGate`.
   * 
   * The supported types for `GateObject` need to expose a
   * `icarus::trigger::OpticalTriggerGate`-like interface (including also e.g.
   * `icarus::trigger::SingleChannelOpticalTriggerGate` and 
   * `icarus::trigger::MultiChannelOpticalTriggerGate`).
   * They must accept:
   * * to be constructed with a `raw::OpDetWaveform` reference;
   * * to be added a `raw::OpDetWaveform` reference via `add()` call;
   * * to be assigned a `TriggerGateData` type.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * icarus::trigger::TriggerGateBuilder::TriggerGates readTriggerGates(
   *   art::Event const& event,
   *   icarus::trigger::ADCCounts_t const threshold,
   *   art::InputTag const& dataTag
   * ) {
   *   
   *   auto const& gates
   *     = *(event.getValidHandle<std::vector<TriggerGateData_t>>(dataTag));
   *   auto const& gateToWaveforms = *(
   *     event.getValidHandle<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
   *       (dataTag)
   *     );
   *   return {
   *     threshold,
   *     icarus::trigger::FillTriggerGates
   *       <icarus::trigger::SingleChannelOpticalTriggerGate>
   *       (gates, gateToWaveforms)
   *     };
   *   
   * } // readTriggerGates()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will return a `icarus::trigger::TriggerGateBuilder::TriggerGates` object
   * for the specified threshold, reading the information from `dateTag` data
   * products in an _art_ `event`.
   * 
   * 
   * Return value
   * -------------
   * 
   * The returned collection contains one element for each `TriggerGateData`
   * object in `gates`, in the same order. Each of these elements is of type
   * `GateObject`, contains a copy of the data of the corresponding gate,
   * and a list of optical waveforms (`raw::OpDetWaveform` objects) it is
   * associated to.
   * 
   * 
   * Requirements
   * -------------
   * 
   * The requirements bind the gates to their association to waveforms:
   * * each gate must be associated to at least one waveform
   * * the associations must be grouped so that all the association pairs
   *   pertaining a gate are contiguous
   *     * within each of these groups, which is made of at least one waveform,
   *       the waveforms must be ordered by increasing timestamp
   *     * the groups must be in the same order as their associated gates
   *   This constitutes the requirement of
   *   @ref LArSoftProxyDefinitionOneToManySeqAssn "one-to-many sequential association"
   *   with the addition that each element in `gates` must have at least one
   *   associated waveform.
   * 
   */
  template <typename GateObject, typename TriggerGateData>
  std::vector<GateObject> FillTriggerGates(
    std::vector<TriggerGateData> const& gates,
    art::Assns<TriggerGateData, raw::OpDetWaveform> const& gateToWaveforms
    );
  
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename GateObject, typename TriggerGateData>
std::vector<GateObject> icarus::trigger::FillTriggerGates(
  std::vector<TriggerGateData> const& gates,
  art::Assns<TriggerGateData, raw::OpDetWaveform> const& gateToWaveforms
) {
  using GateData_t = GateObject; ///< Type of object being returned,
  
  std::vector<GateData_t> allGates;
  
  auto iGateToWaveform = gateToWaveforms.begin();
  auto const gwend = gateToWaveforms.end();
  
  for (auto const& [ iGate, gate ]: util::enumerate(gates)) {
    
    //
    // find the relevant waveforms for this gate;
    // match by index of the element in the data product collection
    //
    while (iGateToWaveform != gwend) {
      if (iGateToWaveform->first.key() == iGate) break;
      ++iGateToWaveform;
    } // while
    auto const iFirstWaveform = iGateToWaveform;
    while (iGateToWaveform != gwend) {
      if (iGateToWaveform->first.key() != iGate) break;
      ++iGateToWaveform;
    } // while
    if (iFirstWaveform == iGateToWaveform) {
      throw cet::exception("FillTriggerGates")
        << " Could not find any waveform associated to trigger gate #" << iGate
        << "\n";
    }
    
    // NOTE we do not control that all waveforms come from the same channel
    GateData_t gateData(*(iFirstWaveform->second)); // add the first waveform
    auto iEndWaveform = iFirstWaveform;
    while (++iEndWaveform != iGateToWaveform) // add the other waveforms
      gateData.add(*(iEndWaveform->second));
    
    gateData = gate; // copy the gate data from the data product
    
    allGates.push_back(std::move(gateData));
    
  } // for gates
  
  return allGates;
} // icarus::trigger::TriggerDesignPlots::FillTriggerGates()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FILLTRIGGERGATES_H
