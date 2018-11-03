/**
 * @file   icaruscode/Utilities/quantities/electronics.h
 * @brief  Dimensioned variables related to electronics.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 2, 2018
 * @see    icaruscode/Utilities/quantities.h
 * 
 * Set of basic quantities related to electronics. Currently, quantities are
 * defined based on the following units:
 * * tick
 * 
 * This is a header-only library.
 * 
 */

#ifndef ICARUSCODE_UTILITIES_QUANTITIES_ELECTRONICS_H
#define ICARUSCODE_UTILITIES_QUANTITIES_ELECTRONICS_H

// LArSoft libraries
#include "icaruscode/Utilities/quantities.h"

// C/C++ standard libraries
#include <string_view>
#include <ratio>
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
namespace util::quantities {
  
  namespace units {
    
    using namespace std::string_view_literals; // for operator""sv()
    
    struct Tick: public concepts::UnitBase {
      static constexpr auto symbol = "#"sv;
      static constexpr auto name   = "tick"sv;
    };
    
    struct Counts: public concepts::UnitBase {
      static constexpr auto symbol = "#"sv;
      static constexpr auto name   = "counts"sv;
    };
    
  } // namespace units
  
  
  // -- BEGIN Ticks ------------------------------------------------------------
  /**
   * @name Ticks
   * 
   * These tick quantities are tied to `util::quantities::units::Tick`.
   * A few options are provided:
   * 
   * * generic template (`tick_as`), allowing to choose which numerical
   *     representation to use
   * * unsigned integer (`tick`), based on `std::size_t`, ready for use
   * 
   * For this unit in particular, additional options are provided to accommodate
   * the custom of using the unit in plural form: `ticks_as` and `ticks`
   * are exactly equivalent to the singular-named counterparts.
   */
  /// @{
  
  /// Tick number, represented by the specified type `T`.
  template <typename T = std::size_t>
  using tick_as = concepts::scaled_quantity<units::Tick, std::ratio<1>, T>;
  
  /// Alias for common language habits.
  template <typename T = std::size_t>
  using ticks_as = tick_as<T>;
  
  /// Tick number, represented by `std::size_t`.
  using tick = tick_as<>;
  
  /// Alias for common language habits.
  using ticks = tick;
  
  
  /// @}
  // -- END Ticks --------------------------------------------------------------
  
  
  // -- BEGIN ADC counts -------------------------------------------------------
  /**
   * @name ADC counts
   * 
   * These ADC count quantities are tied to `util::quantities::units::Counts`.
   * A few options are provided:
   * 
   * * generic template (`counts_as`), allowing to choose which numerical
   *     representation to use
   * * unsigned integer (`counts`), based on `signed short int`, ready for use
   * 
   */
  /// @{
  
  /// Number of ADC counts, represented by the specified type `T`.
  template <typename T = signed short int>
  using counts_as = concepts::scaled_quantity<units::Counts, std::ratio<1>, T>;
  
  /// Number of ADC counts, represented by `signed short int`.
  using counts = counts_as<>;
  
  // -- END ADC counts ---------------------------------------------------------
  
  
  /**
   * @brief Literal constants for electronics quantities.
   * 
   * These functions allow a simplified syntax for specifying a tick quantity.
   * In order to use these, their namespace must be used:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using namespace util::quantities::electronics_literals;
   * 
   * // definition of `util::quantities::tick` constant:
   * constexpr auto i = 56_tick;
   * 
   * // definition of `util::quantities::counts` constant:
   * constexpr auto q = 675_ADC;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  namespace electronics_literals {
    
    // @{
    /// Literal tick value.
    constexpr tick operator""_tick (long double v)
      { return tick{ static_cast<std::size_t>(v) }; }
    constexpr tick operator""_tick (unsigned long long int v)
      { return tick{ static_cast<std::size_t>(v) }; }
    // @}
    
    // @{
    /// Literal ADC count value.
    constexpr counts operator""_ADC (long double v)
      { return counts{ static_cast<signed short int>(v) }; }
    constexpr counts operator""_ADC (unsigned long long int v)
      { return counts{ static_cast<signed short int>(v) }; }
    // @}
    
    
  } // electronics_literals
  
  
  /// @}
  
} // namespace util::quantities

//------------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_QUANTITIES_ELECTRONICS_H
