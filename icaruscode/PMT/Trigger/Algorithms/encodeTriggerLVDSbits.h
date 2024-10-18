/**
 * @file   icaruscode/PMT/Trigger/Algorithms/encodeTriggerLVDSbits.h
 * @brief  Functions to encode ICARUS trigger LVDS bits.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 28, 2024
 * @see    icaruscode/PMT/Trigger/Algorithms/encodeTriggerLVDSbits.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ENCODETRIGGERLVDSBITS_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ENCODETRIGGERLVDSBITS_H


// C/C++ standard libraries
#include <array>
#include <cstdint> // std::uint64_t


// -----------------------------------------------------------------------------
namespace icarus::trigger {
  
  class LVDSbitMaps; // forward declaration
  
  /**
   * @brief Encodes all the trigger LVDS bits for the specified `cryostat`.
   * @param PMTpairMap interface to the PMT mapping database
   * @param cryostat `0` for east cryostat, `1` for west cryostat
   * @param connector01word the raw value of connectors 0 and 1
   * @param connector23word the raw value of connectors 2 and 3
   * @param compact (default: `false`) remove alignment bits between connectors
   * @return an array with the encoded LVDS bits.
   * 
   * This function encodes the bits from the trigger hardware FPGA input
   * connectors into a defined format as required by `sbn::ExtraTriggerInfo`.
   * 
   * The connector words are encoded like in ICARUS hardware, with connectors 0
   * and 2 as most significant words and connector 1 and 3 as least significant
   * words. For each 32-bit connector word, the bits 0-23 report the LVDS state
   * for PMT pair discrimination as mapped in the database, while adder bits are
   * somewhere between bits 26 and 29, as also mapped in the database
   * (here adders are not encoded).
   * 
   * The returned encoding has the bits in the (LArSoft) PMT channel ID order,
   * the most significant bit including PMT channel 0. For LVDS bits
   * representing a PMT pair, the position is assigned considering the smallest
   * of the two channel ID. The bits are still structured in two separate words,
   * in the form `0x00AABBCC'00DDEEFF`, with `A`'s the lowest PMT channel IDs,
   * unless `compact` flag is set, in which case the word is encoded as
   * `0x0000AABB'CCDDEEFF`.
   * 
   * This function requires access to the PMT mapping database.
   * 
   */
  std::array<std::uint64_t, 2U> encodeTriggerLVDSbits(
    icarus::trigger::LVDSbitMaps const& PMTpairMap,
    short int cryostat,
    std::uint64_t connector01word, std::uint64_t connector23word,
    bool compact = false
    );
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ENCODETRIGGERLVDSBITS_H
