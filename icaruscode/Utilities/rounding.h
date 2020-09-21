/**
 * @file   icaruscode/Utilities/rounding.h
 * @brief  Utilities for numerical rounding.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 17, 2020
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_ROUNDING_H
#define ICARUSCODE_UTILITIES_ROUNDING_H


// C/C++ standard libraries
#include <cmath> // std::floor()


// ---------------------------------------------------------------------------
namespace icarus::ns::util {
  
  template <typename T, typename U>
  constexpr T rounddown(T const value, U const quantum, T const offset = T{});
  
  template <typename T, typename U>
  constexpr T roundup(T const value, U const quantum, T const offset = T{});
  
} // namespace icarus::ns::util


// ---------------------------------------------------------------------------
/**
 * @brief Returns the `value`, rounded down.
 * @tparam T type of the value to be rounded
 * @tparam U type of the quantization value
 * @param value the value to be rounded down
 * @param quantum the rounding unit
 * @param offset an offset to subtract for the rounding
 * @return `value` rounded down to multiples of `quantum`
 * 
 * The `value` is returned rounded _down_ into multiples of `quantum`.
 * Optionally, the rounding happens only for the amount of `value` above an
 * `offset`.
 * Examples:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::cout << icarus::ns::util::rounddown(23.0, 2.5)
 *   << "\n" << icarus::ns::util::rounddown(23.0, 2.5, 1.0)
 *   << std::endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will print `22.5` and `21`.
 * 
 * 
 * Type requirements
 * ------------------
 * 
 * * `T` must support addition and subtraction;
 * * `U` must support product to `double`: `P operator* (U, double)`, with
 *   the result `P` implicitly convertible to `T`;
 * * `T` and `U` must support multiplication (`P operator* (T, U)`), and the
 *   result `P` must be convertible back to `T`;
 * * `T` and `U` must support division (`R operator/ (T, U)`), and the result
 *   `R` must be implicitly convertible to an integral or floating point type;
 * * `T` and `U` must support addition (`S operator+ (T, U)`), and the result
 *   `S` must be convertible back to type `T`.
 * 
 */
template <typename T, typename U>
constexpr T icarus::ns::util::rounddown
  (T const value, U const quantum, T const offset /* = T{} */)
{
  return
    static_cast<T>(offset + quantum * std::floor((value - offset) / quantum));
} // icarus::ns::util::rounddown()


/**
 * @brief Returns the `value`, rounded up.
 * @tparam T type of the value to be rounded
 * @tparam U type of the quantization value
 * @param value the value to be rounded up
 * @param quantum the rounding unit
 * @param offset an offset to subtract for the rounding
 * @return `value` rounded up to multiples of `quantum`
 * 
 * The `value` is returned rounded _up_ into multiples of `quantum`.
 * Optionally, the rounding happens only for the amount of `value` above an
 * `offset`.
 * Examples:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::cout << icarus::ns::util::roundup(23.0, 2.5)
 *   << "\n" << icarus::ns::util::roundup(23.0, 2.5, 1.0)
 *   << std::endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will print `25` and `23.5`.
 * 
 * 
 * Type requirements
 * ------------------
 * 
 * * `T` must support comparison (`operator < (T, T)`;
 * * `T` must support addition and subtraction;
 * * `double` must be convertivle _into_ `T`;
 * * `T` and `U` must support multiplication (`P operator* (T, U)`), and the
 *   result `P` must be convertible back to type `T`;
 * * `T` and `U` must support division (`R operator/ (T, U)`), and the result
 *   `R` must be implicitly convertible to an integral or floating point type;
 * * `T` and `U` must support addition (`S operator+ (T, U)`), and the result
 *   `S` must be convertible back to type `T`.
 * 
 */
template <typename T, typename U>
constexpr T icarus::ns::util::roundup
  (T const value, U const quantum, T const offset /* = T{} */)
{
  T const rounded = rounddown(value, quantum, offset);
  return (rounded < value)? static_cast<T>(rounded + quantum): rounded;
} // icarus::ns::util::roundup()


// ---------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_ROUNDING_H

