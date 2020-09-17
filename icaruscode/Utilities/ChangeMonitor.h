/**
 * @file   icaruscode/Utilities/ChangeMonitor.h
 * @brief  Classes to detect the change of object values.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 16, 2020
 */

#ifndef ICARUSCODE_UTILITIES_CHANGEMONITOR_H
#define ICARUSCODE_UTILITIES_CHANGEMONITOR_H


// C/C++ standard libraries
#include <mutex>
#include <optional>
#include <utility> // std::move()
#include <functional> // std::equal_to<>


namespace icarus::ns::util {
  
  //----------------------------------------------------------------------------
  /**
   * @brief Helper to check if an object has changed.
   * @tparam T type of the object
   * @tparam Comp type of the comparison between `T` objects
   * 
   * This class can report if a value has changed from a previous check.
   * The usage pattern is:
   * 1. a `ChangeMonitor` is created
   * 2. the monitor object is given a value, which becomes the new reference
   * 3. arbitrary processing happens
   * 4. the monitor object is given another value, which becomes the new
   *   reference; if this value is different from the previous reference, the
   *   previous reference is returned
   * 5. steps 2 and 3 can be reiterated
   * 
   * Example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * // starts with no reference by default
   * icarus::ns::util::ChangeMonitor<int> monitor;
   * 
   * // first check just establishes the reference
   * int var = 0;
   * monitor(var); // this is also a `update()`, which returns no value
   * 
   * // reference is 0, new value is 1: a change is detected
   * if (monitor(1)) {
   *   std::cout << "Value has changed!" << std::endl;
   * }
   * 
   * var = 5; // this does not change the monitoring
   * // reference is now 1, new value is 1: no change is detected
   * if (monitor(1)) {
   *   std::cout << "Value has changed again!" << std::endl;
   * }
   * 
   * // more complex syntax (C++17) allows accessing the old reference value;
   * // reference is now 1, new value is 2: change is detected
   * if (auto prevVal = monitor(2); prevVal) {
   *   std::cout << "Value has changed from " << *prevVal << " to 2!"
   *     << std::endl;
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * A few observations:
   * * the object with default construction starts with no reference value;
   *   since changes are reported only when a reference value is available,
   *   the first `update()` will not report any change; a constructor is
   *   available to establish a reference at construction time;
   * * the monitor operates on _values_ but it does not monitor a variable the
   *   value might be stored in; therefore, even if the value of the variable
   *   `var` which established the reference is later changed, the monitor is
   *   oblivious to that change.
   * 
   * A copy of the "current" `T` object is kept in this object.
   * 
   * Requirements for `T`:
   * 
   * * must be _CopyConstructible_
   * * must be _EqualityComparable_ (if `Comp` is default)
   * 
   * @note This implementation is not thread-safe. For a thread-safe
   *       change monitor, see `icarus::ns::util::ThreadSafeChangeMonitor`.
   */
  template <typename T, typename Comp = std::equal_to<T>>
  struct ChangeMonitor {
    
    using Data_t = T; ///< Type of the object being monitored.
    using Comparison_t = Comp; ///< Type of object for reference comparison.
    
    /// Default constructor: starts with no reference value.
    ChangeMonitor(Comparison_t comp = Comparison_t{}): fComp(std::move(comp)) {}
    
    /// Constructor: starts with `ref` as the reference value.
    ChangeMonitor(Data_t const& ref, Comp comp = Comp{})
      : fRefObj(ref), fComp(std::move(comp)) {}
    
    
    /**
     * @brief Returns the old object if different from `newObj`.
     * @param currentObj the current object value
     * @return the old reference if different from `currentObj`, or no value
     * 
     * If there is no reference value, no value is returned.
     * Otherwise, the reference is updated and the value of the old reference
     * is returned.
     * After each call, the reference value will be equivalent to `currentObj`.
     * 
     * @note "No value returned" technically means that the `std::optional`
     *       object returned has no value, i.e. `has_value()` is `false`.
     */
    std::optional<Data_t> update(Data_t const& currentObj)
      {
        if (fRefObj && same(currentObj, fRefObj.value())) return {};
        auto lastObj = std::move(fRefObj);
        fRefObj.emplace(currentObj);
        return lastObj;
      }
    
    /// As `update()`.
    std::optional<Data_t> operator() (Data_t const& currentObj)
      { return update(currentObj); }
    
    /// Returns whether a reference value is present.
    bool hasReference() const { return fRefObj.has_value(); }
    
    /// Returns the reference value; undefined if `hasReference()` is `false`.
    Data_t const& reference() const { return fRefObj.value(); }
    
      private:
    
    /// The last object seen.
    std::optional<Data_t> fRefObj;
    
    Comparison_t fComp; ///< Comparison used for reference testing.
    
    /// Returns whether `A` and `B` represent the same value.
    bool same(Data_t const& A, Data_t const& B) const { return fComp(A, B); }
    
  }; // ChangeMonitor
  
  
  // Deduction guide: a single parameter is always a reference value.
  template <typename T>
  ChangeMonitor(T const&) -> ChangeMonitor<T>;
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Helper to check if an object has changed. Thread-safe.
   * @tparam T type of the object
   * @tparam Comp type of the comparison between `T` objects
   * 
   * This class operates like `ChangeMonitor`, but it is made thread-safe by
   * the use of a mutex.
   * 
   * @note This class is actually only _partially_ thread-safe:
   *       the member `reference()` is effectively not, since it returns a
   *       reference that can be then modified by another thread while accessed
   *       (read only) by another.
   */
  template <typename T, typename Comp = std::equal_to<T>>
  class ThreadSafeChangeMonitor: public ChangeMonitor<T, Comp> {
    
    using Base_t = ChangeMonitor<T, Comp>;
    using Data_t = typename Base_t::Data_t;
    
    mutable std::mutex fLock;
    
      public:
    
    // inherit constructors too
    using ChangeMonitor<T, Comp>::ChangeMonitor;
    
    
    /**
     * @brief Returns the old object if different from `newObj`.
     * @param currentObj the current object value
     * @return the old reference if different from `currentObj`, or no value
     * @see `ChangeMonitor::update()`
     * 
     * See `ChangeMonitor::update()` for details.
     */
    std::optional<Data_t> update(Data_t const& currentObj)
      { std::lock_guard lg { fLock }; return Base_t::update(currentObj); }
    
    /// As `update()`.
    std::optional<Data_t> operator() (Data_t const& currentObj)
      { return update(currentObj); }
    
    /// Returns whether a reference value is present.
    bool hasReference() const
      { std::lock_guard lg { fLock }; return Base_t::hasReference(); }
    
    /// Returns the reference value; undefined if `hasReference()` is `false`.
    /// @note While the method is thread-safe, it returns a reference to a
    ///       mutable object.
    Data_t const& reference() const
      { std::lock_guard lg { fLock }; return Base_t::reference(); }
    
  }; // ThreadSafeChangeMonitor
  
  
  // Deduction guide: a single parameter is always a reference value.
  template <typename T>
  ThreadSafeChangeMonitor(T const&) -> ThreadSafeChangeMonitor<T>;
  
  // ---------------------------------------------------------------------------
  
} // namespace icarus::ns::util


#endif // ICARUSCODE_UTILITIES_CHANGEMONITOR_H
