/**
 * @file   icaruscode/Utilities/Expected.h
 * @brief  Class providing a value or an error.
 * @author petrillo@slac,stanford.edu
 * @date   March 13, 2024
 * 
 * This library is header-only.
 */

#ifndef ICARUSCODE_UTILITIES_EXPECTED_H
#define ICARUSCODE_UTILITIES_EXPECTED_H

// C/C++ standard libraries
#include <variant>
#include <utility> // std::move(), std::forward()


// -----------------------------------------------------------------------------
namespace util { template <typename T, typename E = int> class Expected; }

/**
 * @brief Class holding a value or an error.
 * @tparam T type of the value
 * @tparam E (default: `int`) type of the error
 * 
 * This class is loosely based on C++-23 `std::expected`. When C++23 is
 * available, this class should be removed.
 * The main difference is the default value of the error type `E`, which is not
 * present at all in the standard (so the upgrade path would be to define
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * namespace util {
 *   template <typename T, typename E = int>
 *   using Expected = std::expected<T, E>;
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
template <typename T, typename E /* = int */>
class util::Expected {
  
  std::variant<T, E> fValue; ///< The stored value or error.
  
    public:
  
  /// Constructor: starts with a default-constructed error value.
  constexpr Expected(E e = E{}): fValue{ std::move(e) } {}
  
  /// Constructor: starts with a copy of the value `t`.
  constexpr Expected(T t): fValue{ std::move(t) } {}
  
  Expected<T, E>& operator= (E e) { fValue = std::move(e); return *this; }
  Expected<T, E>& operator= (T t) { fValue = std::move(t); return *this; }
  
  /// Constructs the expected _value_ in place.
  template <typename... Args>
  T& emplace(Args&&... args)
    { return fValue.emplace(std::forward<Args>(args)...); }
  
  /// Returns whether a value is stored (as opposed to an error).
  bool has_value() const { return std::holds_alternative<T>(fValue); }
  
  /// Returns the error (undefined behaviour when storing a value).
  E error() const { return std::get<E>(fValue); }
  
  /// Returns the value (undefined behaviour when storing an error).
  T const& value() const { return std::get<T>(fValue); }
  
  /// Returns a pointer to the value (`nullptr` when storing an error).
  T const* valuePtr() const { return has_value()? &value(): nullptr; }
  
  /// Implicit conversion to value (undefined behaviour if storing an error).
  T const& operator*() const { return value(); }
  
  /// Implicit conversion to value (undefined behaviour if storing an error).
  T const* operator-> () const { return valuePtr(); }
  
  /// Returns whether a value is stored (as opposed to an error).
  explicit operator bool() const { return has_value(); }
  
}; // util::Expected


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_UTILITIES_EXPECTED_H
