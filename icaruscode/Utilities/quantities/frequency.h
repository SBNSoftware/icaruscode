/**
 * @file   icaruscode/Utilities/quantities/frequency.h
 * @brief  Dimensioned variables representing frequency quantities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 30, 2018
 * @see    icaruscode/Utilities/quantities.h,
 *         icaruscode/Utilities/quantities/spacetime.h
 * 
 * Set of quantities related to frequency (inverse of time). Currently,
 * quantities are defined based on the following units:
 * * hertz (Hz, kHz, MHz, GHz)
 * 
 * Also, special operations with units from`util::quantities::unit::Second` are
 * supported:
 * * _t_ x _f_ = _s_
 * * _s_ / _t_ = _f_
 * * _s_ / _f_ = _t_
 * where _t_ is a time, _f_ a frequency and _s_ a pure number.
 * 
 * This is a header-only library.
 * 
 */

#ifndef ICARUSCODE_UTILITIES_QUANTITIES_FREQUENCY_H
#define ICARUSCODE_UTILITIES_QUANTITIES_FREQUENCY_H

// LArSoft libraries
#include "icaruscode/Utilities/quantities/spacetime.h" // ...::units::Second
#include "icaruscode/Utilities/quantities.h"


// C/C++ standard libraries
#include <string_view>
#include <ratio>


//------------------------------------------------------------------------------
namespace util::quantities {
  
  namespace units {
    
    using namespace std::string_view_literals; // for operator""sv()
    
    struct Hertz: public concepts::UnitBase {
      static constexpr auto symbol = "Hz"sv;
      static constexpr auto name   = "hertz"sv;
    };
    
  } // namespace units
  
  
  // -- BEGIN Frequency --------------------------------------------------------
  /**
   * @name Frequency quantities
   * 
   * These frequency quantities are tied to `util::quantities::units::Hertz`.
   * A few options are provided:
   * 
   * * most general template, `scaled_hertz`, allowing to choose both the scale
   *     of the unit (e.g. `std::kilo` for kilohertz) and the type of the
   *     numerical representation
   * * generic template (e.g. `hertz_as`), allowing to choose which numerical
   *     representation to use
   * * double precision (e.g. `hertz`), ready for use
   * 
   */
  /// @{
  
  /// The most generic `units::Hertz`-based quantity.
  template <typename R, typename T = double>
  using scaled_hertz = concepts::scaled_quantity<units::Hertz, R, T>;
  
  //
  // hertz
  //
  /// Type of frequency stored in hertz.
  template <typename T = double>
  using hertz_as = scaled_hertz<std::ratio<1>, T>;
  
  /// Type of frequency stored in hertz, in double precision.
  using hertz = hertz_as<>;
  
  //
  // kilohertz
  //
  /// Type of frequency stored in kilohertz.
  template <typename T = double>
  using kilohertz_as = concepts::rescale<hertz_as<T>, std::kilo>;
  
  /// Type of frequency stored in kilohertz, in double precision.
  using kilohertz = kilohertz_as<>;
  
  //
  // megahertz
  //
  /// Type of frequency stored in megahertz.
  template <typename T = double>
  using megahertz_as = concepts::rescale<hertz_as<T>, std::mega>;
  
  /// Type of frequency stored in megahertz, in double precision.
  using megahertz = megahertz_as<>;
  
  //
  // gigahertz
  //
  /// Type of frequency stored in gigahertz.
  template <typename T = double>
  using gigahertz_as = concepts::rescale<hertz_as<T>, std::giga>;
  
  /// Type of frequency stored in gigahertz, in double precision.
  using gigahertz = gigahertz_as<>;
  
  
  /**
   * @brief Literal constants for frequency quantities.
   * 
   * These functions allow a simplified syntax for specifying a frequency
   * quantity. In order to use these, their namespace must be used:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using namespace util::quantities::frequency_literals;
   * 
   * // definition of `util::quantities::hertz` constant:
   * constexpr auto f_Hz = 12_Hz;
   * 
   * // assignment (likely to a quantity) of `util::quantities::megahertz{50.0}`
   * f_Hz = 50_MHz;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  namespace frequency_literals {
    
    // @{
    /// Literal hertz value.
    constexpr hertz operator""_Hz (long double v)
      { return hertz{ static_cast<double>(v) }; }
    constexpr hertz operator""_Hz (unsigned long long int v)
      { return hertz{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal kilohertz value.
    constexpr kilohertz operator""_kHz (long double v)
      { return kilohertz{ static_cast<double>(v) }; }
    constexpr kilohertz operator""_kHz (unsigned long long int v)
      { return kilohertz{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal megahertz value.
    constexpr megahertz operator""_MHz (long double v)
      { return megahertz{ static_cast<double>(v) }; }
    constexpr megahertz operator""_MHz (unsigned long long int v)
      { return megahertz{ static_cast<double>(v) }; }
    // @}
    
    // @{
    /// Literal gigahertz value.
    constexpr gigahertz operator""_GHz (long double v)
      { return gigahertz{ static_cast<double>(v) }; }
    constexpr gigahertz operator""_GHz (unsigned long long int v)
      { return gigahertz{ static_cast<double>(v) }; }
    // @}
    
  } // namespace frequency_literals
  
  
  /// @}
  // -- END Frequency ----------------------------------------------------------
  
  
  // -- BEGIN Special operations between Second and Hertz units ----------------
  /**
   * @name Special operations between Second and Hertz units.
   * 
   * The following operations are supported:
   * * _t_ x _f_ = _s_
   * * _s_ / _t_ = _f_ // TODO
   * * _s_ / _f_ = _t_ // TODO
   * where _t_ is a time, _f_ a frequency and _s_ a pure number.
   * 
   */
  /// @{
  
  namespace concepts {
    //@{
    /**
     * @brief Returns the product (as scalar) of a time and a frequency.
     * @tparam TR type of time unit scale (e.g. `std::micro`)
     * @tparam TT type of time value representation
     * @tparam FR type of frequency unit scale (e.g. `std::mega`)
     * @tparam FT type of frequency value representation
     * @param t time quantity, based on `util::quantities::units::Second`
     * @param f frequency quantity, based on `util::quantities::units::Hertz`
     * @return the product of the two (using the C++ type of `TT * FT`)
     * 
     */
    template <typename TR, typename TT, typename FR, typename FT>
    constexpr auto operator*(
      scaled_quantity<util::quantities::units::Second, TR, TT> t,
      scaled_quantity<util::quantities::units::Hertz, FR, FT> f
      )
      -> decltype(std::declval<TT>() * std::declval<FT>())
      ;
    template <typename FR, typename FT, typename TR, typename TT>
    constexpr auto operator*(
      scaled_quantity<util::quantities::units::Hertz, FR, FT> f,
      scaled_quantity<util::quantities::units::Second, TR, TT> t
      )
      { return t * f; }
    //@}
    
    
    /**
     * @brief Returns a frequency as the inverse of a time.
     * @tparam T type of pure number
     * @tparam TR type of time unit scale (e.g. `std::micro`)
     * @tparam TT type of time value representation
     * @param v scalar value to be divided
     * @param t time quantity, based on `util::quantities::unit::Second`
     * @return a frequency _f_ so that `f * t` equals `v`
     * 
     * The scale of the frequency unit is the inverse of the time one (e.g.,
     * a division by `util::quantities::millisecond` gives
     * `util::quantities::kilohertz`).
     */
    template <typename T, typename TR, typename TT>
    constexpr auto operator/(
      T v,
      scaled_quantity<util::quantities::units::Second, TR, TT> t
      )
      -> std::enable_if_t<
        std::is_fundamental_v<T> && std::is_arithmetic_v<T>,
        scaled_quantity<
          util::quantities::units::Hertz,
          details::invert_t<TR>,
          decltype(std::declval<T>() / std::declval<TT>())
        >
        >
      ;
    
    
    /**
     * @brief Returns a time as the inverse of a frequency.
     * @tparam T type of pure number
     * @tparam FR type of frequency unit scale (e.g. `std::mega`)
     * @tparam FT type of frequency value representation
     * @param v scalar value to be divided
     * @param t frequency quantity, based on `util::quantities::unit::Hertz`
     * @return a time _t_ so that `t * f` equals `v`
     * 
     * The scale of the time unit is the inverse of the frequency one (e.g.,
     * a division by `util::quantities::kilohertz` gives
     * `util::quantities::millisecond`).
     */
    template <typename T, typename FR, typename FT>
    constexpr auto operator/(
      T v,
      scaled_quantity<util::quantities::units::Hertz, FR, FT> f
      )
      -> std::enable_if_t<
        std::is_fundamental_v<T> && std::is_arithmetic_v<T>,
        scaled_quantity<
          util::quantities::units::Second,
          details::invert_t<FR>,
          decltype(std::declval<T>() / std::declval<FT>())
        >
        >
      ;
    
    
  } // namespace concepts
  
  /// @}
  // -- END Special operations between Second and Hertz units ------------------

  
} // namespace util::quantities

//------------------------------------------------------------------------------
//--- template implementation
template <typename TR, typename TT, typename FR, typename FT>
constexpr auto util::quantities::concepts::operator*(
  scaled_quantity<util::quantities::units::Second, TR, TT> t,
  scaled_quantity<util::quantities::units::Hertz, FR, FT> f
  )
  -> decltype(std::declval<TT>() * std::declval<FT>())
{
  return details::applyRatioToValue<simplified_ratio_multiply<TR, FR> >
    (t.value() * f.value());
} // util::quantities::operator*(Second, Hertz)


//------------------------------------------------------------------------------
template <typename T, typename TR, typename TT>
constexpr auto util::quantities::concepts::operator/(
  T v,
  scaled_quantity<util::quantities::units::Second, TR, TT> t
  )
  -> std::enable_if_t<
    std::is_fundamental_v<T> && std::is_arithmetic_v<T>,
    scaled_quantity<
      util::quantities::units::Hertz,
      details::invert_t<TR>,
      decltype(std::declval<T>() / std::declval<TT>())
    >
    >
{
  return scaled_quantity<
    util::quantities::units::Hertz,
    details::invert_t<TR>,
    decltype(std::declval<T>() / std::declval<TT>())
    >
    { v / t.value() };
} // util::quantities::operator/(Second)


//------------------------------------------------------------------------------
template <typename T, typename FR, typename FT>
constexpr auto util::quantities::concepts::operator/(
  T v,
  scaled_quantity<util::quantities::units::Hertz, FR, FT> f
  )
  -> std::enable_if_t<
    std::is_fundamental_v<T> && std::is_arithmetic_v<T>,
    scaled_quantity<
      util::quantities::units::Second,
      details::invert_t<FR>,
      decltype(std::declval<T>() / std::declval<FT>())
    >
    >
{
  return scaled_quantity<
    util::quantities::units::Second,
    details::invert_t<FR>,
    decltype(std::declval<T>() / std::declval<FT>())
    >
    { v / f.value() };
} // util::quantities::operator/(Hertz)


//------------------------------------------------------------------------------


//------------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_QUANTITIES_FREQUENCY_H
