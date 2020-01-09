/**
 * @file   icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h
 * @brief  Class to create an object representing a beam gate.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_BEAMGATEMAKER_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_BEAMGATEMAKER_H

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // detinfo::timescales



//------------------------------------------------------------------------------
//
// declarations
//
namespace icarus::trigger {
  class BeamGateMaker;
} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Simple utility to generate gates around beam time.
 * 
 * Example: to create an object describing a beam gate of 10 &micro;s, run:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * using namespace util::quantities::time_literals;
 * 
 * icarus::trigger::BeamGateMaker makeBeamGate
 *   (*lar::providerFrom<detinfo::DetectorClocks>());
 * 
 * auto const beamGate = makeBeamGate(10_us);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * The time of the gate is in
 * @ref DetectorClocksOpticalElectronicsTime "optical detector electronics ticks".
 */
class icarus::trigger::BeamGateMaker {
  
  using time_interval = detinfo::timescales::time_interval;
  
  detinfo::DetectorTimings const fDetTimings; ///< Detector timing provider.
  
  
  /// Value used for default delay.
  static constexpr util::quantities::microsecond DefaultDelay { 0.0 };
  
  
    public:
  
  /// Constructor: uses a copy of the specified detector timing provider.
  BeamGateMaker(detinfo::DetectorTimings const& detTimings)
    : fDetTimings(detTimings)
    {}
  
  /// Constructor: uses the specified detector clocks service.
  BeamGateMaker(detinfo::DetectorClocks const& detClocks)
    : BeamGateMaker(detinfo::makeDetectorTimings(detClocks)) {}
  
  // @{
  /**
   * @brief Returns a gate object of the specified duration.
   * @tparam Gate type of gate object to be returned
   *              (default: `icarus::trigger::OpticalTriggerGate`)
   * @param length duration of the gate to be open
   * @param delay time from beam gate time when this gate opens (default: `0`)
   * @return a `Gate` object representing a gate opened at beam gate time
   * 
   * The returned object represents a gate opening at
   * @ref DetectorClocksBeamGateTime "beam gate time" (with an optional `delay`)
   * for the specified duration.
   * The time of the gate is in
   * @ref DetectorClocksOpticalElectronicsTime "optical ticks".
   */
  template <typename Gate = icarus::trigger::OpticalTriggerGate>
  Gate make(time_interval length, time_interval delay = DefaultDelay) const
    {
      icarus::trigger::OpticalTriggerGate beamGate;
      beamGate.gateLevels().openFor(
        fDetTimings.toOpticalTick(fDetTimings.BeamGateTime() + delay).value(),
        fDetTimings.toOpticalTicks(length).value()
        );
      return beamGate;
    }
  
  // alias
  template <typename Gate = icarus::trigger::OpticalTriggerGate>
  Gate operator()
    (time_interval length, time_interval delay = DefaultDelay) const
    { return make<Gate>(length, delay); }
  // @}
  
}; // class icarus::trigger::BeamGateMaker


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_BEAMGATEMAKER_H
