/**
 * @file   icaruscode/Utilities/quantities/datasize.h
 * @brief  Dimensioned variables representing data size.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 2, 2018
 * @see    icaruscode/Utilities/quantities.h
 * 
 * Set of basic quantities related to data size.
 * This is mostly a proof of concept for custom prefixes.
 * 
 * This is a header-only library.
 * 
 */

#ifndef ICARUSCODE_UTILITIES_QUANTITIES_DATASIZE_H
#define ICARUSCODE_UTILITIES_QUANTITIES_DATASIZE_H

// LArSoft libraries
#include "icaruscode/Utilities/quantities.h"

// C/C++ standard libraries
#include <string_view>
#include <ratio>


//------------------------------------------------------------------------------
namespace util::quantities {
  
  namespace units {
    
    using namespace std::string_view_literals; // for operator""sv()
    
    struct Byte: public concepts::UnitBase {
      static constexpr auto symbol = "B"sv;
      static constexpr auto name   = "byte"sv;
    };
    
  } // namespace units
  
  
  namespace prefixes {
    
    /// Factor 1'024 (2^10).
    using kibi = std::ratio<(1LL << 10)>;
    
    /// Factor 1'048'576 (2^20).
    using mebi = std::ratio<(1LL << 20)>;
    
    /// Factor 1'073'741'824 (2^30).
    using gibi = std::ratio<(1LL << 30)>;
    
    /// Factor 2^40.
    using tebi = std::ratio<(1LL << 40)>;
    
    /// Factor 2^50.
    using pebi = std::ratio<(1LL << 50)>;
    
    /// Factor 2^60.
    using exbi = std::ratio<(1LL << 60)>;
    
  } // namespace prefixes
  
  
  // -- BEGIN Data size --------------------------------------------------------
  /**
   * @name Data size
   * 
   * These time quantities are tied to `util::quantities::units::Byte`.
   * 
   * * most general template, `scaled_byte`, allowing to choose both the scale
   *     of the unit (e.g. `util::quantities::prefixes::mebi` for mebibyte) and
   *     the type of the numerical representation
   * * generic templates (e.g. `byte_as`), allowing to choose which numerical
   *     representation to use
   * * unsigned integral number (e.g. `byte`), 64-bit on amd64 architecture,
   *     ready for use
   * 
   */
  /// @{
  
  
  /// The most generic `units::Byte`-based quantity.
  template <typename R, typename T = unsigned long long int>
  using scaled_byte = concepts::scaled_quantity<units::Byte, R, T>;
  
  //
  // byte
  //
  /// Type of data size stored in bytes.
  template <typename T = unsigned long long int>
  using byte_as = scaled_byte<std::ratio<1>, T>;
  
  /// Alias for common language habits.
  template <typename T = unsigned long long int>
  using bytes_as = byte_as<T>;
  
  /// Type of data size stored in bytes, in `long long` precision.
  using byte = byte_as<>;
  
  /// Alias for common language habits.
  using bytes = byte;
  
  //
  // kibibyte
  //
  /// Type of data size stored in kibibytes.
  template <typename T = unsigned long long int>
  using kibibyte_as = concepts::rescale<byte_as<T>, prefixes::kibi>;
  
  /// Alias for common language habits.
  template <typename T = unsigned long long int>
  using kibibytes_as = kibibyte_as<T>;

  /// Type of data size stored in kibibytes, in `long long` precision.
  using kibibyte = kibibyte_as<>;
  
  /// Alias for common language habits.
  using kibibytes = kibibyte;
  
  //
  // mebibyte
  //
  /// Type of data size stored in mebibytes.
  template <typename T = unsigned long long int>
  using mebibyte_as = concepts::rescale<byte_as<T>, prefixes::mebi>;
  
  /// Alias for common language habits.
  template <typename T = unsigned long long int>
  using mebibytes_as = mebibyte_as<T>;

  /// Type of data size stored in mebibytes, in `long long` precision.
  using mebibyte = mebibyte_as<>;
  
  /// Alias for common language habits.
  using mebibytes = mebibyte;
  
  //
  // gibibyte
  //
  /// Type of data size stored in gibibytes.
  template <typename T = unsigned long long int>
  using gibibyte_as = concepts::rescale<byte_as<T>, prefixes::gibi>;
  
  /// Alias for common language habits.
  template <typename T = unsigned long long int>
  using gibibytes_as = gibibyte_as<T>;

  /// Type of data size stored in pebibytes, in `long long` precision.
  using gibibyte = gibibyte_as<>;
  
  /// Alias for common language habits.
  using gibibytes = gibibyte;
  
  //
  // tebibyte
  //
  /// Type of data size stored in tebibytes.
  template <typename T = unsigned long long int>
  using tebibyte_as = concepts::rescale<byte_as<T>, prefixes::tebi>;
  
  /// Alias for common language habits.
  template <typename T = unsigned long long int>
  using tebibytes_as = tebibyte_as<T>;

  /// Type of data size stored in tebibytes, in `long long` precision.
  using tebibyte = tebibyte_as<>;
  
  /// Alias for common language habits.
  using tebibytes = tebibyte;
  
  //
  // pebibyte
  //
  /// Type of data size stored in pebibytes.
  template <typename T = unsigned long long int>
  using pebibyte_as = concepts::rescale<byte_as<T>, prefixes::pebi>;
  
  /// Alias for common language habits.
  template <typename T = unsigned long long int>
  using pebibytes_as = pebibyte_as<T>;

  /// Type of data size stored in pebibytes, in `long long` precision.
  using pebibyte = pebibyte_as<>;
  
  /// Alias for common language habits.
  using pebibytes = pebibyte;
  
  //
  // exbibyte
  //
  /// Type of data size stored in exbibytes.
  template <typename T = unsigned long long int>
  using exbibyte_as = concepts::rescale<byte_as<T>, prefixes::exbi>;
  
  /// Alias for common language habits.
  template <typename T = unsigned long long int>
  using exbibytes_as = exbibyte_as<T>;

  /// Type of data size stored in exbibytes, in `long long` precision.
  using exbibyte = exbibyte_as<>;
  
  /// Alias for common language habits.
  using exbibytes = exbibyte;
  
  
  /**
   * @brief Literal constants for data size quantities.
   * 
   * These functions allow a simplified syntax for specifying a data size
   * quantity.
   * In order to use these, their namespace must be used:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using namespace util::quantities::datasize_literals;
   * 
   * // definition of `util::quantities::byte` constant:
   * constexpr auto s_B = 12_B;
   * 
   * // assignment (likely to a quantity) of
   * // `util::quantities::kibibyte{512}`
   * s_B = 512_kiB;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  namespace datasize_literals {
    
    // @{
    /// Literal second value.
    constexpr byte operator""_B (long double v)
      { return byte{ static_cast<unsigned long long int>(v) }; }
    constexpr byte operator""_B (unsigned long long int v)
      { return byte{ v }; }
    // @}
    
    // @{
    /// Literal kibibyte value.
    constexpr kibibyte operator""_kiB (long double v)
      { return kibibyte{ static_cast<unsigned long long int>(v) }; }
    constexpr kibibyte operator""_kiB (unsigned long long int v)
      { return kibibyte{ v }; }
    // @}
    
    // @{
    /// Literal mebibyte value.
    constexpr mebibyte operator""_MiB (long double v)
      { return mebibyte{ static_cast<unsigned long long int>(v) }; }
    constexpr mebibyte operator""_MiB (unsigned long long int v)
      { return mebibyte{ v }; }
    // @}
    
    // @{
    /// Literal gibibyte value.
    constexpr gibibyte operator""_GiB (long double v)
      { return gibibyte{ static_cast<unsigned long long int>(v) }; }
    constexpr gibibyte operator""_GiB (unsigned long long int v)
      { return gibibyte{ v }; }
    // @}
    
    // @{
    /// Literal tebibyte value.
    constexpr tebibyte operator""_TiB (long double v)
      { return tebibyte{ static_cast<unsigned long long int>(v) }; }
    constexpr tebibyte operator""_TiB (unsigned long long int v)
      { return tebibyte{ v }; }
    // @}
    
    // @{
    /// Literal pebibyte value.
    constexpr pebibyte operator""_PiB (long double v)
      { return pebibyte{ static_cast<unsigned long long int>(v) }; }
    constexpr pebibyte operator""_PiB (unsigned long long int v)
      { return pebibyte{ v }; }
    // @}
    
    // @{
    /// Literal exbibyte value.
    constexpr exbibyte operator""_EiB (long double v)
      { return exbibyte{ static_cast<unsigned long long int>(v) }; }
    constexpr exbibyte operator""_EiB (unsigned long long int v)
      { return exbibyte{ v }; }
    // @}
    
    
  } // datasize_literals
  
  
  /// @}
  // -- END Data size ----------------------------------------------------------
  
} // namespace util::quantities

//------------------------------------------------------------------------------
//--- Template specializations
//------------------------------------------------------------------------------
namespace util::quantities::concepts {
  
  /// Prefix for 1024 (2^10).
  template <>
  struct Prefix<prefixes::kibi> {
    
      /// Returns the symbol of the prefix.
      static constexpr auto symbol() { return "ki"sv; }
      
      /// Returns the full name of the prefix.
      static constexpr auto name() { return "kibi"sv; }
    
  }; // struct Prefix<prefixes::kibi>
  
  
  /// Prefix for 1048576 (2^20).
  template <>
  struct Prefix<prefixes::mebi> {
    
      /// Returns the symbol of the prefix.
      static constexpr auto symbol() { return "Mi"sv; }
      
      /// Returns the full name of the prefix.
      static constexpr auto name() { return "mebi"sv; }
    
  }; // struct Prefix<prefixes::mebi>
  
  
  /// Prefix for 1073741824 (2^30).
  template <>
  struct Prefix<prefixes::gibi> {
    
      /// Returns the symbol of the prefix.
      static constexpr auto symbol() { return "Gi"sv; }
      
      /// Returns the full name of the prefix.
      static constexpr auto name() { return "gibi"sv; }
    
  }; // struct Prefix<prefixes::gibi>
  
  
  /// Prefix for 2^40.
  template <>
  struct Prefix<prefixes::tebi> {
    
      /// Returns the symbol of the prefix.
      static constexpr auto symbol() { return "Ti"sv; }
      
      /// Returns the full name of the prefix.
      static constexpr auto name() { return "tebi"sv; }
    
  }; // struct Prefix<prefixes::tebi>
  
  
  /// Prefix for 2^50.
  template <>
  struct Prefix<prefixes::pebi> {
    
      /// Returns the symbol of the prefix.
      static constexpr auto symbol() { return "Pi"sv; }
      
      /// Returns the full name of the prefix.
      static constexpr auto name() { return "pebi"sv; }
    
  }; // struct Prefix<prefixes::pebi>
  
  
  /// Prefix for 2^60.
  template <>
  struct Prefix<prefixes::exbi> {
    
      /// Returns the symbol of the prefix.
      static constexpr auto symbol() { return "Ei"sv; }
      
      /// Returns the full name of the prefix.
      static constexpr auto name() { return "exbi"sv; }
    
  }; // struct Prefix<prefixes::exbi>
  
  
} // namespace util::quantities::concepts


//------------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_QUANTITIES_DATASIZE_H
