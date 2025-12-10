/**
 * @file   icaruscode/PMT/Trigger/Utilities/DiscrThresholdParser.h
 * @brief  Utility to parse a value (threshold) in voltage or ADC counts.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 9, 2025
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_DISCRTHRESHOLDPARSER_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_DISCRTHRESHOLDPARSER_H


// ICARUS and LArSoft libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Algorithms/ADCsettings.h"
#include "lardataalg/Utilities/quantities/electromagnetism.h" // Volt
#include "lardataalg/Utilities/quantities.h"

// C++ standard libraries
#include <cstdint> // std::size_t
#include <cctype> // std::isblank()
#include <charconv> // std::from_char()
#include <functional> // std::bind()
#include <stdexcept> // std::out_of_range, std::logic_error, ...
#include <string_view>
#include <system_error> // std::errc
#include <type_traits> // std::enable_if_t()
#include <utility> // std::move()
#include <vector>


//------------------------------------------------------------------------------
namespace icarus::trigger{ class DiscrThresholdParser; }
/**
 * @brief Converter of voltage level specifications into ADC counts.
 * 
 * This converter seamlessly supports specifications in ADC counts (pure
 * numbers) or voltage (`util::quantities::units::Volt` quantity objects).
 * 
 * The interface allows for the conversion of a single value or of a sequence.
 * 
 * The type values are converted to is `icarus::trigger::ADCCounts_t`, which is
 * an integral type. Rounding (as well as conversions) is performed via the
 * `icarus::trigger::ADCsettings` helper.
 * 
 */
class icarus::trigger::DiscrThresholdParser {
  
  /// Type trait: `T` is a string (actually: binds/converts to a string_view).
  template <typename T>
  static constexpr bool is_string_v
    = std::is_convertible_v<T, std::string_view>;
  
  /// Conversion utility.
  icarus::ADCsettings<double> const fADCsettings;
  
    public:
  
  /// Constructor: assigns the conversion parameters (or uses CAEN V1730 defaults).
  DiscrThresholdParser(icarus::ADCsettings<double> ADCsettings = {})
    : fADCsettings{ std::move(ADCsettings) }
    {}
  
  // --- BEGIN ---  Parsing functions  -----------------------------------------
  /// @name Parsing functions
  /// @{
  
  /**
   * @brief Parses a string and returns the corresponding ADC count value.
   * @tparam String a string type for the argument
   * @param spec single specification of the level, as ADC or voltage
   * @return the value of `spec` converted in `ADCCounts_t`
   * @see @ref icarus_trigger_DiscrThresholdParser_parse_Coll "parse(Coll const&)"
   * 
   * Converts `spec` into a ADC value (quantity object of type
   * `icarus::trigger::ADCCounts_t`).
   * The value is rounded using `icarus::trigger::ADCsettings::roundADC()`.
   * 
   * 
   */
  ADCCounts_t parse(std::string_view spec) const { return parseImpl(spec); }
  
  /**
   * @brief  Parses a string and returns the corresponding ADC count value.
   * @anchor icarus_trigger_DiscrThresholdParser_parse_Coll
   * @tparam Coll type of a collection of string type objects
   * @param  specs collection of specifications to convert
   * @return collection of values from the sequence, converted in `ADCCounts_t`
   * @see `parse(std::string_view)`
   */
  template <typename Coll>
  std::enable_if_t<
    !is_string_v<Coll> && is_string_v<typename Coll::value_type>,
    std::vector<ADCCounts_t>
  >
  parse(Coll const& specs) const;
  
  /**
   * @brief  Parses a string and returns the corresponding ADC count value.
   * @tparam BIter type of iterator to the strings to convert
   * @tparam EIter type of end-of-sequence iterator
   * @param  sbegin iterator pointing to the first specification to convert
   * @param  send iterator pointing past the last specification to convert
   * @return collection of values from the sequence, converted in `ADCCounts_t`
   */
  template <typename BIter, typename EIter>
  std::vector<ADCCounts_t> parse(BIter sbegin, EIter send) const;
  
  /// @}
  // ---- END ----  Parsing functions  -----------------------------------------
  
  
  /**
   * @brief  Converts `sv` to a numerical value.
   * @tparam T (default: `double`) type of real or integral value to read
   * @param  sv the string to be converted to a number
   * @return the value from `sv`
   * @throw  std::out_of_range the value is beyond the range of `T`
   * @throw  std::invalid_argument `sv` does not represent a number
   * @throw  std::logic_error on other (unexpected) errors
   * 
   * The argument `sv` is stripped of heading and trailing spaces.
   * After that, the full remaining string must be converted into a numerical
   * value.
   */
  template <typename T = double>
  static double toNumber(std::string_view sv);
  
  /**
   * @brief Returns whether `spec` is a voltage (`util::quantities::units::Volt`).
   * @param spec candidate voltage specification string
   * @return whether `spec` is a voltage
   * 
   * This function is less strict than an actual conversion (i.e. for a `spec`
   * that returns `true`, `makeQuantity<volt>(spec)` could still fail, mostly
   * due to missing spaces or invalid prefixes).
   */
  static bool isVoltSpec(std::string_view spec);
  
    private:
  
  /// Implementation of `parse(string_view)`.
  /// @see `parse(std::string_view)`
  ADCCounts_t parseImpl(std::string_view spec) const;
  
  
}; // icarus::trigger::DiscrThresholdParser



//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
template <typename Coll>
auto icarus::trigger::DiscrThresholdParser::parse(Coll const& specs) const
  -> std::enable_if_t<
    !is_string_v<Coll> && is_string_v<typename Coll::value_type>,
    std::vector<ADCCounts_t>
  >
{
  return parse(begin(specs), end(specs));
}


//------------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::DiscrThresholdParser::parse
  (BIter sbegin, EIter send) const -> std::vector<ADCCounts_t>
{
  std::vector<icarus::trigger::ADCCounts_t> ADCs;
  std::transform(sbegin, send, back_inserter(ADCs),
    std::bind(&DiscrThresholdParser::parseImpl, this, std::placeholders::_1));
  return ADCs;
} // icarus::trigger::DiscrThresholdParser::parse(BIter, EIter)


//------------------------------------------------------------------------------
template <typename T /* = double */ >
double icarus::trigger::DiscrThresholdParser::toNumber
  (std::string_view sv)
{
  
  /*
   * std::from_char() is... difficult:
   *  * it does not want leading spaces
   *  * it does not want leading plus sign (minus is fine)
   *  * it may stop before the end of the string
   */
  auto const lstrip = [](std::string_view sv) -> std::string_view
    {
      while (!sv.empty()) {
        if (!std::isblank(sv.front())) return sv;
        sv.remove_prefix(1);
      }
      return {};
    };
  
  // remove all the stuff std::from_char() does not like at the front
  sv = lstrip(sv);
  if (!sv.empty() && (sv.front() == '+'))
    { sv.remove_prefix(1); sv = lstrip(sv); }
  
  T val;
  auto const svendptr = sv.data() + sv.length();
  std::from_chars_result const res = std::from_chars(sv.data(), svendptr, val);
  switch (res.ec) {
    case std::errc::invalid_argument:
      throw std::invalid_argument
        { "Invalid double argument '" + std::string{ sv } + "'" };
    case std::errc::result_out_of_range:
      throw std::out_of_range
        { "Argument out of double range: '" + std::string{ sv } + "'" };
    default:
      if (res.ec != std::errc{}) { // can't stay in `case:`
        throw std::logic_error{
          "Unexpected error converting '" + std::string{ sv } + "': "
          + std::make_error_condition(res.ec).message()
          };
      }
      if (res.ptr != svendptr) {
        throw std::invalid_argument
          { "Extra characters after value: '" + std::string{ sv } + "'" };
      }
      return val;
  } // switch
  
} // icarus::trigger::DiscrThresholdParser::toNumber(()


//------------------------------------------------------------------------------
//---  Inline implementation
//------------------------------------------------------------------------------
inline bool icarus::trigger::DiscrThresholdParser::isVoltSpec
  (std::string_view sv)
{
  std::string_view const& unitStr = util::quantities::units::Volt::symbol;
  
  // strip trailing blanks
  while (!sv.empty()) {
    if (!std::isblank(sv.back())) break;
    sv.remove_suffix(1);
  }
  
  // C++20: sv.ends_with(unitStr)!!
  return sv.substr(sv.length() - unitStr.length()) == unitStr;
} // icarus::trigger::DiscrThresholdParser::isVoltSpec()


//------------------------------------------------------------------------------
inline auto icarus::trigger::DiscrThresholdParser::parseImpl
  (std::string_view sv) const -> ADCCounts_t
{
  using util::quantities::volt;
  
  short int const ADC = isVoltSpec(sv)
    ? fADCsettings.to_ADC(util::quantities::makeQuantity<volt>(sv))  // volt
    : fADCsettings.roundADC(toNumber(sv))                     // pure number
    ;
  
  return icarus::trigger::ADCCounts_t{ ADC };
} // icarus::trigger::DiscrThresholdParser::parseImpl()


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_DISCRTHRESHOLDPARSER_H
