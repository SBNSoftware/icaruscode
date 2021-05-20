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
#include <string>
#include <initializer_list>


// -----------------------------------------------------------------------------
/**
 * @brief Set of values to identify a beam.
 * 
 * The constants are used to identify the gate during which a trigger is
 * acquired, in the `raw::Trigger::TriggerBits()` bit mask.
 */
namespace sbn {
  
  /// Type of beam or beam gate or other trigger source.
  enum class triggerSource: unsigned int {
    Unknown, ///< Type of beam unknown.
    BNB,     ///< Type of beam: BNB.
    NuMI,    ///< Type of beam: NuMI.
    // ==> add here if more are needed <==
    NBits    ///< Number of bits currently supported.
  }; // triggerSource
  
  /// Returns the value of the specified bit.
  constexpr unsigned int value(triggerSource bit);
  
  /// Returns a mask with all specified bits set.
  template <typename... OtherBits>
  constexpr unsigned int mask(triggerSource bit, OtherBits... otherBits);
  //@}
  
  /// Returns a mnemonic short name of the beam type.
  std::string bitName(triggerSource bit);
  
} // namespace sbn


// -----------------------------------------------------------------------------
inline constexpr unsigned int sbn::value(triggerSource bit)
  { return static_cast<int>(bit); }


// -----------------------------------------------------------------------------
template <typename... OtherBits>
constexpr unsigned int sbn::mask(triggerSource bit, OtherBits... otherBits) {
  unsigned int m { 1U << value(bit) };
  if constexpr(sizeof...(OtherBits) > 0U) m |= mask(otherBits...);
  return m;
} // sbn::mask()


// -----------------------------------------------------------------------------
/// @todo Move into an implementation file once this header in the final location.
inline std::string sbn::bitName(triggerSource bit) {
  
  using namespace std::string_literals;
  switch (bit) {
    case triggerSource::Unknown: return "unknown"s;
    case triggerSource::BNB:     return "BNB"s;
    case triggerSource::NuMI:    return "NuMI"s;
    case triggerSource::NBits:   return "<invalid>"s;
  } // switch
  
} // sbn::bitName()


// -----------------------------------------------------------------------------

#endif // SBNOBJ_COMMON_TRIGGER_BEAMBITS_H
