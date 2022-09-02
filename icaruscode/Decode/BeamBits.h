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

    /// Converts a integral type into a mask.
    template <typename EnumType>
    mask_t<EnumType> makeMask
      (typename mask_t<EnumType>::maskbits_t bits) noexcept;
    
    
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
      Calib,     ///< Type of Offbeam: Calibration
      // ==> add here if more are needed <==
      NBits    ///< Number of bits currently supported.
    }; // triggerSource

    /// Type of mask with `triggerSource` bits.
    using triggerSourceMask = mask_t<triggerSource>;

    
    /// Location or locations generating a trigger.
    enum class triggerLocation: unsigned int {
      CryoEast,    ///< A trigger happened in the east cryostat.
      CryoWest,    ///< A trigger happened in the west cryostat.
      TPCEE,       ///< A trigger happened in the east cryostat, east TPC.
      TPCEW,       ///< A trigger happened in the east cryostat, west TPC.
      TPCWE,       ///< A trigger happened in the west cryostat, east TPC.
      TPCWW,       ///< A trigger happened in the west cryostat, west TPC.
      // ==> add here if more are needed <==
      NBits    ///< Number of bits currently supported.
    }; // triggerLocation

    /// Type of mask with `triggerLocation` bits.
    using triggerLocationMask = mask_t<triggerLocation>;

    /// Type representing the type(s) of this trigger.
    enum class triggerType: unsigned int {
      Majority,    ///< A minimum number of close-by PMT pairs above threshold was reached.
      MinimumBias, ///< Data collected at gate opening with no further requirement imposed.
      // ==> add here if more are needed <==
      NBits    ///< Number of bits currently supported.
    }; // triggerType

    /// Type of mask with `triggerType` bits.
    using triggerTypeMask = mask_t<triggerType>;
    
    /// Trigger window mode
    enum class triggerWindowMode: unsigned int {
      Separated,    ///< Separated, non-overlapping contigous window
      Overlapping,  ///< Overlaping windows
      //==> add here if more are needed <==
      NBits     ///< Number of Bits currently supported
    };

    

    /// Enabled gates in the trigger configuration. See register 0X00050008 in docdb SBN-doc-23778-v1
    enum class gateSelection: unsigned int {
      GateBNB,                     ///<Enable receiving BNB early warning signal (gatedBES) to open BNB gates
      DriftGateBNB,                ///<Enable BNB early-early warning signal ($1D) for light out-of-time in BNB gates
      GateNuMI,                    ///<Enable NuMI early warning signal (MIBS$74) to open NuMI gates
      DriftGateNuMI,               ///<Enable receiving NuMI early-early warning signal ($AD) for light out-of-time in NuMI gates
      GateOffbeamBNB,              ///<Enable Offbeam gate for BNB
      DriftGateOffbeamBNB,         ///<Enable Offbeam drift gate BNB (for light out-of-time in offbeam gates)
      GateOffbeamNuMI,             ///<Enable Offbeam gate for NuMI
      DriftGateOffbeamNuMI,        ///<Enable Offbeam drift gate NuMI (for light out-of-time in offbeam gates)
      GateCalibration,             ///<Enable Calibration gate
      DriftGateCalibration,        ///<Enable Calibration drift gate (for light out-of-time in calibration gates)
      MinbiasGateBNB,              ///<Enable MinBias triggers for the BNB stream
      MinbiasGateNuMI,             ///<Enable MinBias triggers for the NuMI stream
      MinbiasGateOffbeamBNB,       ///<Enabke MinBias triggers for the Offbeam BNB stream
      MinbiasGateOffbeamNuMI,      ///<Enable MinBias triggers for the Offbeam NuMI stream
      MinbiasGateCalibration,      ///<Enable MinBias triggers for the Calibration stream
      MinbiasDriftGateBNB,         ///<Enable light out-of-time for MinBias triggers in BNB stream
      MinbiasDriftGateNuMI,        ///<Enable light out-of-time for MinBias triggers in NuMI stream
      MinbiasDriftGateOffbeamBNB,  ///<Enable light out-of-time for MinBias triggers in Offbeam BNB stream
      MinbiasDriftGateOffbeamNuMI, ///<Enable light out-of-time for MinBias triggers in Offbeam NuMI stream
      MinbiasDriftGateCalibration, ///<Enable light out-of-time for MinBias triggers in Calibration stream
      // ==> add here if more are needed <==
      NBits
    }; // gateSelection

    using gateSelectionMask = mask_t<gateSelection>;

    /// Returns a mnemonic short name of the beam type.
    std::string bitName(triggerSource bit);
    /// Returns a mnemonic short name of the trigger location.
    std::string bitName(triggerLocation bit);
    /// Returns a mnemonic short name of the trigger type.
    std::string bitName(triggerType bit);
    /// Returns a mnemonic short name for the trigger window mode.
    std::string bitName(triggerWindowMode bit);
    /// Returns a mnemonic short name for the trigger window mode.
    std::string bitName(gateSelection bit);

    /// @}
    // --- END ---- Beam bits --------------------------------------------------

  } // namespace bits
  
  using bits::triggerSource; // import symbol
  using bits::triggerSourceMask;
  using bits::triggerType;
  using bits::triggerTypeMask;
  using bits::triggerLocation;
  using bits::triggerLocationMask;
  using bits::triggerWindowMode;
  using bits::gateSelection;
  
} // namespace sbn


// -----------------------------------------------------------------------------
// ---  inline and template implementation
// -----------------------------------------------------------------------------
template <typename EnumType>
auto sbn::bits::makeMask  (typename mask_t<EnumType>::maskbits_t bits) noexcept
  -> mask_t<EnumType>
  { return { bits }; }


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
    case triggerSource::Calib:       return "Calib"s;
    case triggerSource::NBits:       return "<invalid>"s;
  } // switch
  throw std::runtime_error("sbn::bits::bitName(triggerSource{ "s
    + std::to_string(value(bit)) + " }): unknown bit"s);
} // sbn::bitName()

// -----------------------------------------------------------------------------

inline std::string sbn::bits::bitName(triggerLocation bit) {

  using namespace std::string_literals;
  switch (bit) {
    case triggerLocation::CryoEast: return "Cryo E"s;
    case triggerLocation::CryoWest: return "Cryo W"s;
    case triggerLocation::TPCEE:    return "TPC EE"s;
    case triggerLocation::TPCEW:    return "TPC EW"s;
    case triggerLocation::TPCWE:    return "TPC WE"s;
    case triggerLocation::TPCWW:    return "TPC WW"s;
    case triggerLocation::NBits:    return "<invalid>"s;
  } // switch
  throw std::runtime_error("sbn::bits::bitName(triggerLocation{ "s
    + std::to_string(value(bit)) + " }): unknown bit"s);
} // sbn::bitName(triggerLocation)

// -----------------------------------------------------------------------------
inline std::string sbn::bits::bitName(triggerType bit) {

  using namespace std::string_literals;
  switch (bit) {
    case triggerType::Majority:    return "majority"s;
    case triggerType::MinimumBias: return "minimum bias"s;
    case triggerType::NBits:       return "<invalid>"s;
  } // switch
  throw std::runtime_error("sbn::bits::bitName(triggerType{ "s
    + std::to_string(value(bit)) + " }): unknown bit"s);
} // sbn::bitName(triggerType)

inline std::string sbn::bits::bitName(triggerWindowMode bit) {

  using namespace std::string_literals;
  switch (bit) {
    case sbn::bits::triggerWindowMode::Separated:    return "Separated Window"s;
    case sbn::bits::triggerWindowMode::Overlapping:  return "Overlapping Window"s;
    case sbn::bits::triggerWindowMode::NBits:    return "<invalid>"s;
  } // switch
  throw std::runtime_error("sbn::bits::bitName(triggerWindowMode{ "s
    + std::to_string(value(bit)) + " }): unknown bit"s);
} // triggerWindowMode


inline std::string sbn::bits::bitName(gateSelection bit) {

  using namespace std::string_literals;
  switch (bit) {
    case gateSelection::GateBNB:                     return "GateBNB"s;
    case gateSelection::DriftGateBNB:                return "DriftGateBNB"s;
    case gateSelection::GateNuMI:                    return "GateNuMI"s;
    case gateSelection::DriftGateNuMI:               return "DriftGateNuMI"s;
    case gateSelection::GateOffbeamBNB:              return "GateOffbeamBNB"s;
    case gateSelection::DriftGateOffbeamBNB:         return "DriftGateOffbeamBNB"s;
    case gateSelection::GateOffbeamNuMI:             return "GateOffbeamNuMI"s;
    case gateSelection::DriftGateOffbeamNuMI:        return "DriftGateOffbeamNuMI"s;
    case gateSelection::GateCalibration:             return "GateCalibration"s;
    case gateSelection::DriftGateCalibration:        return "DriftGateCalibration"s;
    case gateSelection::MinbiasGateBNB:              return "MinbiasGateBNB"s;
    case gateSelection::MinbiasGateNuMI:             return "MinbiasGateNuMI"s;
    case gateSelection::MinbiasGateOffbeamBNB:       return "MinbiasGateOffbeamBNB"s;
    case gateSelection::MinbiasGateOffbeamNuMI:      return "MinbiasGateOffbeamNuMI"s;
    case gateSelection::MinbiasGateCalibration:      return "MinbiasGateCalibration"s;
    case gateSelection::MinbiasDriftGateBNB:         return "MinbiasDriftGateBNB"s;
    case gateSelection::MinbiasDriftGateNuMI:        return "MinbiasDriftGateNuMI"s;
    case gateSelection::MinbiasDriftGateOffbeamBNB:  return "MinbiasDriftGateOffbeamBNB"s;
    case gateSelection::MinbiasDriftGateOffbeamNuMI: return "MinbiasDriftGateOffbeamNuMI"s;
    case gateSelection::MinbiasDriftGateCalibration: return "MinbiasDriftGateCalibration"s;
    case gateSelection::NBits:                       return "NBits"s;
  } // switch
  throw std::runtime_error("sbn::bits::bitName(gateSelection{ "s
    + std::to_string(value(bit)) + " }): unknown bit"s);
} // sbn::bitName()


namespace icarus::trigger {

  using triggerLocation = sbn::triggerLocation;
  using triggerSource   = sbn::triggerSource;

  static constexpr std::size_t kEast             = sbn::bits::value<triggerLocation>(triggerLocation::CryoEast);
  static constexpr std::size_t kWest             = sbn::bits::value<triggerLocation>(triggerLocation::CryoWest);
  static constexpr std::size_t kNTriggerLocation = sbn::bits::value<triggerLocation>(triggerLocation::NBits);

  static constexpr std::size_t kBNB            = sbn::bits::value<triggerSource>(triggerSource::BNB);
  static constexpr std::size_t kNuMI           = sbn::bits::value<triggerSource>(triggerSource::NuMI);
  static constexpr std::size_t kOffBeamBNB     = sbn::bits::value<triggerSource>(triggerSource::OffbeamBNB);
  static constexpr std::size_t kOffBeamNuMI    = sbn::bits::value<triggerSource>(triggerSource::OffbeamNuMI);
  static constexpr std::size_t kCalibration    = sbn::bits::value<triggerSource>(triggerSource::Calib);
  static constexpr std::size_t kNTriggerSource = sbn::bits::value<triggerSource>(triggerSource::NBits);

}


// -----------------------------------------------------------------------------

#endif // SBNOBJ_COMMON_TRIGGER_BEAMBITS_H
