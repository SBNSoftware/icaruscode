/**
 * @file   icaruscode/Utilities/quantities/electromagnetism.h
 * @brief  Dimensioned variables representing electromagnetic quantities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 2, 2018
 * @see    icaruscode/Utilities/quantities.h
 * 
 * Set of basic quantities related to electromagnetism. Currently, quantities
 * are defined based on the following units:
 * * coulomb (fC, pC, nC, uC, mC, C)
 * 
 * This is a header-only library.
 * 
 * @todo Also belong here: ampere, volt, farad, ohm...
 * 
 */

#ifndef ICARUSCODE_UTILITIES_QUANTITIES_ELECTROMAGNETISM_H
#define ICARUSCODE_UTILITIES_QUANTITIES_ELECTROMAGNETISM_H

// LArSoft libraries
#include "icaruscode/Utilities/quantities.h"

// C/C++ standard libraries
#include <string_view>
#include <ratio>


//------------------------------------------------------------------------------
namespace util::quantities {
  
  namespace units {
    
    using namespace std::string_view_literals; // for operator""sv()
    
    struct Coulomb: public concepts::UnitBase {
      static constexpr auto symbol = "C"sv;
      static constexpr auto name   = "coulomb"sv;
    };
    
    struct Volt: public concepts::UnitBase {
      static constexpr auto symbol = "V"sv;
      static constexpr auto name   = "volt"sv;
    };
    
  } // namespace units
  
  
  // -- BEGIN Charge -----------------------------------------------------------
  /**
   * @name Charge quantities
   * 
   * These time quantities are tied to `util::quantities::units::Coulomb`.
   * A few options are provided:
   * 
   * * most general template, `scaled_coulomb`, allowing to choose both the
   *     scale of the unit (e.g. `std::pico` for picocoulomb) and the type of
   *     the numerical representation
   * * generic templates (e.g. `coulomb_as`), allowing to choose which numerical
   *     representation to use
   * * double precision (e.g. `coulomb`), ready for use
   * 
   */
  /// @{
  
  
  /// The most generic `units::Coulomb`-based quantity.
  template <typename R, typename T = double>
  using scaled_coulomb = concepts::scaled_quantity<units::Coulomb, R, T>;
  
  //
  // coulomb
  //
  /// Type of charge stored in coulomb.
  template <typename T = double>
  using coulomb_as = scaled_coulomb<std::ratio<1>, T>;
  
  /// Type of time stored in coulombs, in double precision.
  using coulomb = coulomb_as<>;
  
  //
  // millicoulomb
  //
  /// Type of charge stored in millicoulomb.
  template <typename T = double>
  using millicoulomb_as = concepts::rescale<coulomb_as<T>, std::milli>;
  
  /// Type of time stored in millicoulomb, in double precision.
  using millicoulomb = millicoulomb_as<>;
  
  //
  // microcoulomb
  //
  /// Type of charge stored in microcoulomb.
  template <typename T = double>
  using microcoulomb_as = concepts::rescale<coulomb_as<T>, std::micro>;
  
  /// Type of time stored in microcoulomb, in double precision.
  using microcoulomb = microcoulomb_as<>;
  
  //
  // nanocoulomb
  //
  /// Type of charge stored in nanocoulomb.
  template <typename T = double>
  using nanocoulomb_as = concepts::rescale<coulomb_as<T>, std::nano>;
  
  /// Type of time stored in nanocoulomb, in double precision.
  using nanocoulomb = nanocoulomb_as<>;
  
  //
  // picocoulomb
  //
  /// Type of charge stored in picocoulomb.
  template <typename T = double>
  using picocoulomb_as = concepts::rescale<coulomb_as<T>, std::pico>;
  
  /// Type of time stored in picocoulomb, in double precision.
  using picocoulomb = picocoulomb_as<>;
  
  //
  // femtocoulomb
  //
  /// Type of charge stored in femtocoulomb.
  template <typename T = double>
  using femtocoulomb_as = concepts::rescale<coulomb_as<T>, std::femto>;
  
  /// Type of time stored in femtocoulomb, in double precision.
  using femtocoulomb = femtocoulomb_as<>;
  
  
  // -- END Charge -------------------------------------------------------------
  
  
  // -- BEGIN Electric potential -----------------------------------------------
  /**
   * @name Electric potential quantities
   * 
   * These time quantities are tied to `util::quantities::units::Volt`.
   * A few options are provided:
   * 
   * * most general template, `scaled_volt`, allowing to choose both the
   *     scale of the unit (e.g. `std::kilo` for kilovolt) and the type of
   *     the numerical representation
   * * generic templates (e.g. `volt_as`), allowing to choose which numerical
   *     representation to use
   * * double precision (e.g. `volt`), ready for use
   * 
   */
  /// @{
  
  
  /// The most generic `units::Volt`-based quantity.
  template <typename R, typename T = double>
  using scaled_volt = concepts::scaled_quantity<units::Volt, R, T>;
  
  //
  // volt
  //
  /// Type of potential stored in volt.
  template <typename T = double>
  using volt_as = scaled_volt<std::ratio<1>, T>;
  
  /// Type of time stored in volts, in double precision.
  using volt = volt_as<>;
  
  //
  // millivolt
  //
  /// Type of potential stored in millivolt.
  template <typename T = double>
  using millivolt_as = concepts::rescale<volt_as<T>, std::milli>;
  
  /// Type of time stored in millivolt, in double precision.
  using millivolt = millivolt_as<>;
  
  //
  // microvolt
  //
  /// Type of potential stored in microvolt.
  template <typename T = double>
  using microvolt_as = concepts::rescale<volt_as<T>, std::micro>;
  
  /// Type of time stored in microvolt, in double precision.
  using microvolt = microvolt_as<>;
  
  //
  // kilovolt
  //
  /// Type of potential stored in kilovolt.
  template <typename T = double>
  using kilovolt_as = concepts::rescale<volt_as<T>, std::kilo>;
  
  /// Type of time stored in kilovolt, in double precision.
  using kilovolt = kilovolt_as<>;
  
  //
  // megavolt
  //
  /// Type of potential stored in megavolt.
  template <typename T = double>
  using megavolt_as = concepts::rescale<volt_as<T>, std::mega>;
  
  /// Type of time stored in megavolt, in double precision.
  using megavolt = megavolt_as<>;
  
  //
  // femtovolt
  //
  /// Type of potential stored in gigavolt.
  template <typename T = double>
  using gigavolt_as = concepts::rescale<volt_as<T>, std::giga>;
  
  /// Type of time stored in gigavolt, in double precision.
  using gigavolt = gigavolt_as<>;
  
  /// @}
  // -- END Electric potential -------------------------------------------------
  
  
  /**
   * @brief Literal constants for time quantities.
   * 
   * These functions allow a simplified syntax for specifying a time quantity.
   * In order to use these, their namespace must be used:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using namespace util::quantities::time_literals;
   * 
   * // definition of `util::quantities::picocoulomb` constant:
   * constexpr auto Q_pC = 230_pC;
   * 
   * // assignment (likely to a quantity) of
   * // `util::quantities::femtocoulomb{500.0}`
   * Q = 500_fC;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  namespace electromagnetism_literals {
    
    // @{
    /// Literal coulomb value.
    constexpr coulomb operator""_C (long double v)
      { return coulomb{ static_cast<double>(v) }; }
    constexpr coulomb operator""_C (unsigned long long int v)
      { return coulomb{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal millicoulomb value.
    constexpr millicoulomb operator""_mC (long double v)
      { return millicoulomb{ static_cast<double>(v) }; }
    constexpr millicoulomb operator""_mC (unsigned long long int v)
      { return millicoulomb{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal microcoulomb value.
    constexpr microcoulomb operator""_uC (long double v)
      { return microcoulomb{ static_cast<double>(v) }; }
    constexpr microcoulomb operator""_uC (unsigned long long int v)
      { return microcoulomb{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal nanocoulomb value.
    constexpr nanocoulomb operator""_nC (long double v)
      { return nanocoulomb{ static_cast<double>(v) }; }
    constexpr nanocoulomb operator""_nC (unsigned long long int v)
      { return nanocoulomb{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal picocoulomb value.
    constexpr picocoulomb operator""_pC (long double v)
      { return picocoulomb{ static_cast<double>(v) }; }
    constexpr picocoulomb operator""_pC (unsigned long long int v)
      { return picocoulomb{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal femtocoulomb value.
    constexpr femtocoulomb operator""_fC (long double v)
      { return femtocoulomb{ static_cast<double>(v) }; }
    constexpr femtocoulomb operator""_fC (unsigned long long int v)
      { return femtocoulomb{ static_cast<double>(v) }; }
    // @}
    
    
    // @{
    /// Literal volt value.
    constexpr volt operator""_V (long double v)
      { return volt{ static_cast<double>(v) }; }
    constexpr volt operator""_V (unsigned long long int v)
      { return volt{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal millivolt value.
    constexpr millivolt operator""_mV (long double v)
      { return millivolt{ static_cast<double>(v) }; }
    constexpr millivolt operator""_mV (unsigned long long int v)
      { return millivolt{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal microvolt value.
    constexpr microvolt operator""_uV (long double v)
      { return microvolt{ static_cast<double>(v) }; }
    constexpr microvolt operator""_uV (unsigned long long int v)
      { return microvolt{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal kilovolt value.
    constexpr kilovolt operator""_kV (long double v)
      { return kilovolt{ static_cast<double>(v) }; }
    constexpr kilovolt operator""_kV (unsigned long long int v)
      { return kilovolt{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal megavolt value.
    constexpr megavolt operator""_MV (long double v)
      { return megavolt{ static_cast<double>(v) }; }
    constexpr megavolt operator""_MV (unsigned long long int v)
      { return megavolt{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal gigavolt value.
    constexpr gigavolt operator""_GV (long double v)
      { return gigavolt{ static_cast<double>(v) }; }
    constexpr gigavolt operator""_GV (unsigned long long int v)
      { return gigavolt{ static_cast<double>(v) }; }
    // @}
    
  } // electromagnetism_literals
  
  
} // namespace util::quantities

//------------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_QUANTITIES_ELECTROMAGNETISM_H
