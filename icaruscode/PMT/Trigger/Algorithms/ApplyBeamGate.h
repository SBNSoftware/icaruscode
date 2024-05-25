/**
 * @file   icaruscode/PMT/Trigger/Algorithms/ApplyBeamGate.h
 * @brief  Helper to manage a beam gate.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 16, 2020
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_APPLYBEAMGATE_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_APPLYBEAMGATE_H


// ICARUS libraries
// #include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateStruct.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icarusalg/Utilities/mfLoggingClass.h"

// // LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // optical_tick, ...
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds

// // C/C++ standard libraries
// #include <utility> // std::pair<>
#include <ostream>
#include <vector>
#include <string>


// -----------------------------------------------------------------------------
namespace icarus::trigger {
  
  class ApplyBeamGateClass;
  
  /**
   * @brief Returns a new `ApplyBeamGateClass` object with the specified gate.
   * @param duration the desired duration of the gate
   * @param delay the delay of gate start relative to the global beam gate time
   * @param clockData the current detector clocks service provider data
   * @param logCategory message category the returned object will use
   * @return a new `ApplyBeamGateClass` object with the specified parameters
   */
  ApplyBeamGateClass makeApplyBeamGate(
    util::quantities::intervals::microseconds duration,
    util::quantities::intervals::microseconds delay,
    detinfo::DetectorClocksData const& clockData,
    std::string const& logCategory = "ApplyBeamGateClass"
    );
  
  /**
   * @brief Returns a new `ApplyBeamGateClass` object with the specified gate.
   * @param duration the desired duration of the gate
   * @param clockData the current detector clocks service provider data
   * @param logCategory message category the returned object will use
   * @return a new `ApplyBeamGateClass` object with the specified parameters
   * 
   * The gate starts at the nominal beam gate time as reported by `clockData`.
   */
  ApplyBeamGateClass makeApplyBeamGate(
    util::quantities::intervals::microseconds duration,
    detinfo::DetectorClocksData const& clockData,
    std::string const& logCategory = "ApplyBeamGateClass"
    );
  
  /**
   * @brief Returns a new `ApplyBeamGateClass` object with the specified gate.
   * @param start the start of the gate, in ticks on the optical time scale
   * @param duration the desired duration of the gate, in optical ticks
   * @param clockData the current detector clocks service provider data
   * @param logCategory message category the returned object will use
   * @return a new `ApplyBeamGateClass` object with the specified parameters
   * 
   * The gate starts at the nominal beam gate time as reported by `clockData`.
   */
  ApplyBeamGateClass makeApplyBeamGate(
    detinfo::timescales::optical_tick start,
    detinfo::timescales::optical_time_ticks duration,
    detinfo::DetectorClocksData const& clockData,
    std::string const& logCategory = "ApplyBeamGateClass"
    );
  
  std::ostream& operator<<
    (std::ostream& out, icarus::trigger::ApplyBeamGateClass const& gate);
  
} // namespace icarus::trigger

/**
 * @brief Helper applying a beam gate to any gate.
 * 
 * The gate starts from `detinfo::DetectorClocks::BeamGateTime()` (with an
 * optional delay) and has length specified on construction.
 * 
 * The assumption that the optical tick clock starts with the electronics time
 * is used.
 * 
 * Access to beam gate range may happen in a few ways: in simulation time scale
 * and in optical time ticks. In all cases, the returned object supports the
 * calls to `start()`, `stop()` and `duration()`.
 */
class icarus::trigger::ApplyBeamGateClass
  : protected icarus::ns::util::mfLoggingClass
{

    public:
  // @{
  /// Type aliases
  using microseconds = util::quantities::intervals::microseconds;
  using optical_tick = detinfo::timescales::optical_tick;
  using optical_time_ticks = detinfo::timescales::optical_time_ticks;
  // @}
  
  /// Constructor: gets the gate (in optical ticks) and its duration (in time).
  ApplyBeamGateClass(
    icarus::trigger::BeamGateStruct&& beamGate,
    std::string const& logCategory = "ApplyBeamGateClass"
    )
    : icarus::ns::util::mfLoggingClass(logCategory)
    , fGate(std::move(beamGate))
    {}
  
  /// Returns the beam gate as a `icarus::trigger::OpticalTriggerGate`.
  icarus::trigger::OpticalTriggerGate const& gate() const
    { return fGate.asGate(); }
  
  /// Returns the beam gate as the specified gate type `Gate`.
  template <typename Gate>
  Gate gateAs() const { return { gate() }; }
  
  /// Returns a copy of `gate` in AND with this beam gate.
  template <typename Gate>
  Gate apply(Gate gate) const
    { return icarus::trigger::mulGates(gate, this->gate()); }
  
  /// Returns a collection of copies of the specified `gates` each in AND with
  /// this beam gate.
  template
    <typename GateColl, typename GateObj = typename GateColl::value_type>
  std::vector<GateObj> applyToAll(GateColl const& gates) const;
  
  /// Returns the gate as a simulation time range.
  auto asSimulationTime() const { return fGate.asSimulationRange(); }
  
  /// Returns the range of the beam gate as start and stop tick.
  auto tickRange() const { return fGate.asOptTickRange(); }
  
  /// Returns the length of the gate (in time units).
  microseconds length() const { return fGate.duration(); }

  /// Returns the length of the gate (in time units).
  optical_time_ticks lengthTicks() const { return tickRange().duration(); }
  
  //@{
  /// Comparison operators.
  bool operator== (ApplyBeamGateClass const& other) const
    { return fGate == other.fGate; }
  bool operator!= (ApplyBeamGateClass const& other) const
    { return fGate != other.fGate; }
  //@}
  
    private:
  
  icarus::trigger::BeamGateStruct const fGate;
  
}; // class icarus::trigger::ApplyBeamGateClass


// -----------------------------------------------------------------------------
// --- inline implementation
// -----------------------------------------------------------------------------
inline auto icarus::trigger::makeApplyBeamGate(
  ::util::quantities::intervals::microseconds duration,
  ::util::quantities::intervals::microseconds delay,
  detinfo::DetectorClocksData const& clockData,
  std::string const& logCategory /* = "ApplyBeamGateClass" */
  ) -> ApplyBeamGateClass
{
  return {
    icarus::trigger::makeBeamGateStruct
      (detinfo::DetectorTimings{ clockData }, duration, delay),
    logCategory
    };
} // icarus::trigger::makeApplyBeamGate()


// -----------------------------------------------------------------------------
inline auto icarus::trigger::makeApplyBeamGate(
  ::util::quantities::intervals::microseconds duration,
  detinfo::DetectorClocksData const& clockData,
  std::string const& logCategory /* = "ApplyBeamGateClass" */
  ) -> ApplyBeamGateClass
  {
    using namespace util::quantities::time_literals;
    return makeApplyBeamGate(duration, 0.0_us, clockData, logCategory);
  }


// -----------------------------------------------------------------------------
inline auto icarus::trigger::makeApplyBeamGate(
  detinfo::timescales::optical_tick start,
  detinfo::timescales::optical_time_ticks duration,
  detinfo::DetectorClocksData const& clockData,
  std::string const& logCategory /* = "ApplyBeamGateClass" */
  ) -> ApplyBeamGateClass
{
  return {
    icarus::trigger::makeBeamGateStruct
      (detinfo::DetectorTimings{ clockData }, start, duration),
    logCategory
    };
} // icarus::trigger::makeApplyBeamGate()


// -----------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, icarus::trigger::ApplyBeamGateClass const& gate)
{
  out << gate.gate() << " (simulation time: " << gate.asSimulationTime()
    << ")";
  return out;
} // icarus::trigger::operator<< (icarus::trigger::ApplyBeamGateClass const&)


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
template
  <typename GateColl, typename GateObj /* = typename GateColl::value_type */>
std::vector<GateObj> icarus::trigger::ApplyBeamGateClass::applyToAll
  (GateColl const& gates) const
{
  std::vector<GateObj> res;
  std::transform(
    gates.begin(), gates.end(), std::back_inserter(res),
    [this](GateObj const& gate){ return apply(gate); }
    );
  return res;
} // icarus::trigger::ApplyBeamGateClass::applyToAll()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_APPLYBEAMGATE_H


