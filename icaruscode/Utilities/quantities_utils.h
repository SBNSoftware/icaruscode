/**
 * @file   icaruscode/Utilities/quantities_utils.h
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

#ifndef ICARUSCODE_UTILITIES_QUANTITIES_UTILS_H
#define ICARUSCODE_UTILITIES_QUANTITIES_UTILS_H

// LArSoft libraries
#include "lardataalg/Utilities/quantities.h"

// C/C++ standard libraries
#include <type_traits>


namespace util::quantities::concepts { // same namespace where `Quantity` lives
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Absolute value of a quantity.
   * @tparam Q the `Quantity`-like type
   * @param q quantity object to take absolute value of
   * @return an object of type `T` with its value made non-negative
   * 
   * The generic way to use this is:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * void f(T value) {
   *   using std::abs;
   *   T const abs_value = abs(value);
   *   // ...
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template
    <typename Q, typename = std::enable_if_t<details::is_quantity<Q>::value>>
  Q abs(Q q) { return q.abs(); }
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Scalar value of a quantity.
   * @tparam Q the `Quantity`-like type
   * @param v quantity to be returned the value of
   * @return a scalar value of the scalar type of `Q`
   * @see `util::value()`
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
  template
    <typename Q, typename = std::enable_if_t<details::is_quantity<Q>::value>>
  typename Q::value_t value(Q q) { return q.value(); }
  
  
  // ---------------------------------------------------------------------------
  
} // namespace util::quantities::concepts


// -----------------------------------------------------------------------------
namespace util {
  
  /**
   * @brief Scalar value of a scalar quantity.
   * @tparam T the value type
   * @param v value to be returned
   * @return the value `v`
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
  template
    <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
  T value(T v) { return v; }
  
  /// The type returned by `util::value(T)`.
  template <typename T>
  using value_t = decltype(value(std::declval<T>()));
  
} // namespace util
// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_QUANTITIES_UTILS_H
