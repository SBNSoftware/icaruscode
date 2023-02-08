/**
 * @file   icaruscode/Utilities/intervals_utils.h
 * @brief  Utilities to make `util::Quantity` objects behave more like numbers.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 18, 2022
 * 
 * This library tries to provide utilities that a motivated developer may
 * utilize to keep its code `util::quantities::concepts::Quantity`-friendly.
 * 
 * The quantities library is a likely too ambitious project, and although has
 * proven useful at for some level convenient, for others it is a nightmare to
 * use, especially when mixed to standard variables.
 * These utilities are meant to make a syntax available that can be shared in
 * generic code whether the variables involved are quantities or plain values.
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_INTERVALS_UTILS_H
#define ICARUSCODE_UTILITIES_INTERVALS_UTILS_H

// ICARUS libraries
#include "icaruscode/Utilities/quantities_utils.h"

// LArSoft libraries
#include "lardataalg/Utilities/intervals.h"

// C/C++ standard libraries
#include <type_traits>


namespace util::quantites::concepts { // same namespace where `Quantity` lives
  
  namespace details {
    
    template <typename IP>
    using is_interval_or_point
      = std::disjunction<is_interval<IP>, is_point<IP>>;
    
    template <typename IP>
    constexpr bool is_interval_or_point_v = is_interval_or_point<IP>::value;
    
  } // namespace details
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Absolute value of a quantity.
   * @tparam IP the `Interval` or `Point` type
   * @param q interval or point to take absolute value of
   * @return an object of type `IP` with its value made non-negative
   * 
   * The generic way to use this is:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * void f(T value) {
   *   using std::abs;
   *   T = abs(value);
   *   // ...
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  template <
    typename IP,
    typename = std::enable_if_t<details::is_interval_or_point_v<IP>>
    >
  IP abs(IP q) { return q.abs(); }
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Scalar value of a interval or point.
   * @tparam IP the `Interval` or `Point` type
   * @param q interval or point to be returned the value of
   * @return a scalar value of the scalar type of `IP`
   * 
   * The generic way to use this is:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * float to_float(T v) {
   *   using util::value;
   *   return static_cast<float>(value(v));
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template<
    typename IP,
    typename = std::enable_if_t<details::is_interval_or_point_v<IP>>
    >
  typename IP::value_t value(IP q) { return q.value(); }
  
  
  // ---------------------------------------------------------------------------
  
} // namespace util


#endif // ICARUSCODE_UTILITIES_INTERVALS_UTILS_H
