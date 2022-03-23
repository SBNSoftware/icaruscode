/**
 * @file   sbnobj/Common/Trigger/BeamBits.h
 * @brief  Definitions of the trigger bits for SBN.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 19, 2021
 * 
 * This is a header-only library.
 * 
 * @todo This header should be moved into the proper place
 *       (proposal: `sbnobj/Common/Trigger` in `sbnobj`)
 */

#ifndef SBNOBJ_COMMON_TRIGGER_BEAMBITS_H
#define SBNOBJ_COMMON_TRIGGER_BEAMBITS_H


// C/C++ standard libraries
#include <stdexcept> // std::runtime_error
#include <string>
#include <vector>
#include <initializer_list>
#include <type_traits> // std::underlying_type_t


// -----------------------------------------------------------------------------
/**
 * @brief Set of values to identify a beam.
 * 
 * The constants are used to identify the gate during which a trigger is
 * acquired, in the `raw::Trigger::TriggerBits()` bit mask.
 */
namespace sbn {
  
  /**
   * @brief Simple utilities to deal with bit enumerators.
   * 
   * A few simple concepts:
   * * a bit is identified by its position (integral number: `0`, `1`, `2`...)
   * * representation of a bit is via a enumerator (`enum class`)
   * * bits may be associated with a name (`name()`)
   * * to test bits, `mask()`, `value()` and `hasBitSet()` help with the
   *   conversions between `enum class` and the underlying integral type
   * 
   */
  namespace bits {
    
    // --- BEGIN -- Generic bit functions --------------------------------------
    /// @name Generic bit functions
    /// @{
    
    /// Type for bit masks.
    /// @note This is a glorified integral type.
    template <typename EnumType>
    struct mask_t {
      using bits_t = EnumType; ///< Enumeration type of the bits.
      using maskbits_t = std::underlying_type_t<EnumType>; ///< Bit data type.
      maskbits_t bits { 0 };
      operator maskbits_t() const { return bits; }
    }; // mask_t
    
    /// Returns the value of specified `bit` (conversion like `enum` to `int`).
    template <typename EnumType>
    constexpr auto value(EnumType bit);
    
    /// Returns a mask with all specified bits set.
    template <typename EnumType, typename... OtherBits>
    constexpr mask_t<EnumType> mask(EnumType bit, OtherBits... otherBits);
    
    /// Returns whether the specified `bit` is set in `bitMask`.
    template <typename EnumType>
    constexpr bool hasBitSet(mask_t<EnumType> bitMask, EnumType bit);
    
    /// Returns the name of the specified `bit`. Delegates to `bitName()`.
    template <typename EnumType>
    std::string name(EnumType bit);
    
    /// Returns a list of the names of all the bits set in `mask`.
    /// Mask is interpreted as made of only bits of type `EnumType`.
    template <typename EnumType>
    std::vector<std::string> names(mask_t<EnumType> mask);
    
    /// @}
    // --- END ---- Generic bit functions --------------------------------------
    
    
    // --- BEGIN -- Beam bits --------------------------------------------------
    /// @name Beam bits
    /// @{
    
    /// Type of beam or beam gate or other trigger source.
    enum class triggerSource: unsigned int {
      Unknown,     ///< Type of beam unknown.
      BNB,         ///< Type of beam: BNB.
      NuMI,        ///< Type of beam: NuMI.
      OffbeamBNB,  ///< Type of Offbeam: BNB
      OffbeamNuMI, ///< Type of Offbeam: NuMI
      // ==> add here if more are needed <==
      NBits    ///< Number of bits currently supported.
    }; // triggerSource
    
    using triggerSourceMask = mask_t<triggerSource>;
    
    /// Returns a mnemonic short name of the beam type.
    std::string bitName(triggerSource bit);
    
    /// @}
    // --- END ---- Beam bits --------------------------------------------------

  } // namespace bits
  
  using bits::triggerSource; // import symbol
  using bits::triggerSourceMask; // import symbol
  
} // namespace sbn


// -----------------------------------------------------------------------------
// ---  inline and template implementation
// -----------------------------------------------------------------------------
template <typename EnumType>
constexpr auto sbn::bits::value(EnumType bit)
  { return static_cast<std::underlying_type_t<EnumType>>(bit); }


// -----------------------------------------------------------------------------
template <typename EnumType, typename... OtherBits>
constexpr auto sbn::bits::mask(EnumType bit, OtherBits... otherBits)
  -> mask_t<EnumType>
{
  unsigned int m { 1U << value(bit) };
  if constexpr(sizeof...(OtherBits) > 0U) m |= mask(otherBits...);
  return { m };
} // sbn::mask()


// -----------------------------------------------------------------------------
template <typename EnumType>
constexpr bool sbn::bits::hasBitSet(mask_t<EnumType> bitMask, EnumType bit)
  { return bitMask & mask(bit); }


// -----------------------------------------------------------------------------
template <typename EnumType>
std::string sbn::bits::name(EnumType bit)
  { return bitName(bit); }
  
  
// -----------------------------------------------------------------------------
template <typename EnumType>
std::vector<std::string> sbn::bits::names(mask_t<EnumType> mask) {
  static_assert(value(EnumType::NBits) >= 0);
  
  using namespace std::string_literals;
  
  constexpr std::size_t MaxBits = sizeof(mask) * 8U;
  constexpr std::size_t NSupportedBits = value(EnumType::NBits);
  
  std::vector<std::string> names;
  for (std::size_t bit = 0U; bit < MaxBits; ++bit) {
    auto const typedBit = static_cast<EnumType>(bit);
    if (!hasBitSet(mask, typedBit)) continue;
    names.push_back((bit > NSupportedBits)
      ? "<unsupported ["s + std::to_string(bit) + "]>"s: name(typedBit)
      );
  } // for
  return names;
  
} // sbn::bits::names()


// -----------------------------------------------------------------------------
/// @todo Move into an implementation file once this header in the final location.
inline std::string sbn::bits::bitName(triggerSource bit) {
  
  using namespace std::string_literals;
  switch (bit) {
    case triggerSource::Unknown:     return "unknown"s;
    case triggerSource::BNB:         return "BNB"s;
    case triggerSource::NuMI:        return "NuMI"s;
    case triggerSource::OffbeamBNB:  return "OffbeamBNB"s;
    case triggerSource::OffbeamNuMI: return "OffbeamNuMI"s;
    case triggerSource::NBits:       return "<invalid>"s;
  } // switch
  throw std::runtime_error("sbn::bits::bitName(triggerSource{ "s
    + std::to_string(value(bit)) + " }): unknown bit"s);
} // sbn::bitName()

// -----------------------------------------------------------------------------

#endif // SBNOBJ_COMMON_TRIGGER_BEAMBITS_H
