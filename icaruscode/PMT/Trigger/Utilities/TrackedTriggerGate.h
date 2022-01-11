/**
 * @file   icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h
 * @brief  A wrapper to trigger gate objects tracking the contributions.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 21, 2021
 * 
 * This library is header-only.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_TRACKEDTRIGGERGATE_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_TRACKEDTRIGGERGATE_H


// LArSoft libraries
#include "larcorealg/CoreUtils/span.h" // util::make_transformed_span()

// C/C++ standard libraries
#include <iosfwd>
#include <set>
#include <functional> // std::mem_fn()
#include <utility> // std::in_place_t, std::forward(), ...
#include <type_traits> // std::true_type...



// -----------------------------------------------------------------------------
namespace icarus::trigger {
  
  /// Trait returning whether `T` is some form of TrackedTriggerGate`.
  template <typename T>
  struct isTrackedTriggerGate;
  
  /// Whether `T` is some form of TrackedTriggerGate`.
  template <typename T>
  constexpr bool isTrackedTriggerGate_v = isTrackedTriggerGate<T>::value;
  
  template <typename Gate, typename TrackedType>
  class TrackedTriggerGate;
  
  template <typename Gate, typename TrackedType>
  std::ostream& operator<<
    (std::ostream& out, TrackedTriggerGate<Gate, TrackedType> const& gate);
  
  /// Returns an iterable with a reference to the gate of each entry of
  /// `trackingGates`.
  template <typename TrackingGateColl>
  auto gatesIn(TrackingGateColl& trackingGates);
  
} // namespace icarus::trigger

/**
  * @brief A wrapper to trigger gate objects tracking the input of operations.
  * @tparam Gate type of trigger gate object being wrapped
  * @tparam TrackedType type of the objects being tracked
  * 
  * This object includes its own `Gate` object, plus `tracking()`.
  * 
  * Functions taking as argument a generic sequence of `Gate` objects can be
  * passed the result of `gatesIn()`, which behaves like a sequence of `Gate`.
  * 
  * Supporting functions can add specific tracking operations when operating
  * of a `TrackedTriggerGate` object. The trait
  * `icarus::trigger::isTrackedTriggerGate_v` will return whether that is the
  * case.
  * 
  * Tracking is performed by accessing `tracking()` and `add()`'ing objects
  * to be tracked. Tracking is currently just a collection of objects of a
  * specific type (`TrackedType`). That type is expected to support strict
  * comparison. A copy of each tracking object is used and owned: if the object
  * is complex and externally owned, a pointer or handle should be stored
  * as tracking object instead of the object itself.
  * If an object is already present, it is not added again into the tracking.
  */
template <typename Gate, typename TrackedType>
class icarus::trigger::TrackedTriggerGate {
  
    public:
  
  using TriggerGate_t = Gate; ///< Gate type being wrapped.
  using Tracked_t = TrackedType; ///< Type for tracking.
  
  /// Tracked information. Interface is pretty minimal so far.
  class TrackingInfo {
    
    std::set<Tracked_t> fTracked; ///< All tracked objects.
    
      public:
    
    /// Add an object to the list of tracked objects, if it's not present yet.
    void add(TrackedType tracked);
    
    /// Adds all the objects `tracked` in the specified information object.
    void add(TrackingInfo const& tracked);
    
    /// Returns the number of objects currently tracked;
    std::size_t nTracked() const;
    
    /// Whether any tracked object is present.
    bool hasTracked() const;
    
    /// Returns an iterable of all tracked objects.
    auto getTracked() const;
    
  }; // class TrackingInfo
  
  /// Default constructor: default-constructed gate, no tracking.
  TrackedTriggerGate() = default;
  
  /// Constructor: copies the data of the specified gate (no tracking added).
  TrackedTriggerGate(TriggerGate_t gate): fGate(std::move(gate)) {}
  
  
  // @{
  
  /// Returns the tracking information.
  TrackingInfo const& tracking() const { return fTracking; }
  TrackingInfo& tracking() { return fTracking; }
  
  // @}
  
  
  // @{
  
  /// Returns the enclosed gate.
  TriggerGate_t const& gate() const { return fGate; }
  TriggerGate_t& gate() { return fGate; }
  
  // @}
  
  
  /// Returns the list of channels of the enclosed gate.
  decltype(auto) channels() const { return gate().channels(); }
  
  
    private:
  TriggerGate_t fGate; ///< Local copy of the gate information.
  TrackingInfo fTracking; ///< Tracking information.
  
}; // class icarus::trigger::TrackedTriggerGate


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
namespace icarus::trigger {
  
  template <typename T>
  struct isTrackedTriggerGate: std::false_type {};
  
  template <typename Gate, typename TrackedType>
  struct isTrackedTriggerGate<TrackedTriggerGate<Gate, TrackedType>>
    : std::true_type
  {};
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
template <typename TrackingGateColl>
auto icarus::trigger::gatesIn(TrackingGateColl& trackingGates) {
  
  using TrackingGate_t = typename TrackingGateColl::value_type;
  
  // C++20: this is definitely a template concept...
  static_assert(isTrackedTriggerGate_v<TrackingGate_t>,
    "icarus::trigger::gatesIn() requires a collection of TrackedTriggerGate objects");
  
  // constantness is driven by the one of type `TrackedTriggerGate`;
  // decltype(auto) return preserves referencehood
  auto getGate = [](auto& gate) -> decltype(auto) { return gate.gate(); };
  
#if 0
  // C++20: this is some
  return trackingGates | std::ranges::view::transform(getGate);
#endif // 0
  
  return util::make_transformed_span(trackingGates, getGate);
  
} // icarus::trigger::gatesIn()


// -----------------------------------------------------------------------------
// ---  icarus::trigger::TrackedTriggerGate<>
// -----------------------------------------------------------------------------
template <typename Gate, typename TrackedType>
void icarus::trigger::TrackedTriggerGate<Gate, TrackedType>::TrackingInfo::add
  (TrackedType tracked)
  { fTracked.insert(std::move(tracked)); }


// -----------------------------------------------------------------------------
template <typename Gate, typename TrackedType>
void icarus::trigger::TrackedTriggerGate<Gate, TrackedType>::TrackingInfo::add
  (TrackingInfo const& tracked)
  { fTracked.insert(begin(tracked.fTracked), end(tracked.fTracked)); }


// -----------------------------------------------------------------------------
template <typename Gate, typename TrackedType>
std::size_t
icarus::trigger::TrackedTriggerGate<Gate, TrackedType>::TrackingInfo::nTracked()
  const
  { return fTracked.size(); }


// -----------------------------------------------------------------------------
template <typename Gate, typename TrackedType>
bool
icarus::trigger::TrackedTriggerGate<Gate, TrackedType>::TrackingInfo::hasTracked()
  const
  { return !(fTracked.empty()); }


// -----------------------------------------------------------------------------
template <typename Gate, typename TrackedType>
auto icarus::trigger::TrackedTriggerGate<Gate, TrackedType>::TrackingInfo::getTracked()
  const
  { return fTracked; }


// -----------------------------------------------------------------------------
// ---  free function implementation
// -----------------------------------------------------------------------------
template <typename Gate, typename TrackedType>
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, TrackedTriggerGate<Gate, TrackedType> const& gate)
  { return out << gate.gate(); }


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_TRACKEDTRIGGERGATE_H
