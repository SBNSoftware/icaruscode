/**
 * @file   icaruscode/PMT/Trigger/Algorithms/encodeTriggerLVDSbits.cxx
 * @brief  Functions to encode ICARUS trigger LVDS bits.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 28, 2024
 * @see    icaruscode/PMT/Trigger/Algorithms/encodeTriggerLVDSbits.h
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/encodeTriggerLVDSbits.h"

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/LVDSbitMaps.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"

// C/C++ standard libraries
#include <cassert>


// -----------------------------------------------------------------------------
std::array<std::uint64_t, 2U> icarus::trigger::encodeTriggerLVDSbits(
  LVDSbitMaps const& PMTpairMap, short int cryostat,
  std::uint64_t connector01word, std::uint64_t connector23word,
  bool compact /* = false */
) {
  
  std::array<std::uint64_t, 2U> outputWords;
  outputWords.fill(0);
  
  /*
    * The (first) two `LVDSstatus` 64-bit words are initialized according to
    * the prescription in the class documentation: LVDS bits from PMT with
    * lower channel ID end in lower (least significant) bits in `LVDSstatus`.
    * 
    * The algorithm will fill the bits of the output words in sequence,
    * each time picking the value from the appropriate bits in the input
    * connector words. The correct bit is determined by a precooked map.
    */
  
  // connector words have the first connector in the most significant 32 bits:
  std::array<std::uint32_t, 4U> const connectorBits {
    static_cast<std::uint32_t>(connector01word >> 32ULL & 0x00FF'FFFFULL),
    static_cast<std::uint32_t>(connector01word          & 0x00FF'FFFFULL),
    static_cast<std::uint32_t>(connector23word >> 32ULL & 0x00FF'FFFFULL),
    static_cast<std::uint32_t>(connector23word          & 0x00FF'FFFFULL)
  };
    
  auto const mask
    = [](auto bitNo){ return 1ULL << static_cast<std::uint64_t>(bitNo); };
  
  for (auto const PMTwall: util::counter(2U)) {
    
    std::uint64_t& outputWord = outputWords[PMTwall];
    
    for (auto const bit: util::counter<PMTpairBitID::StatusBit_t>(64)) {
      
      PMTpairBitID const bitID{ (unsigned) cryostat, PMTwall, bit };
      
      LVDSHWbitID const source = PMTpairMap.bitSource(bitID).source;
      
      if (!source) continue; // bit not mapped to anything
      
      assert(source.cryostat == cryostat);
      bool const bitValue
        = connectorBits.at(source.connector) & mask(source.bit);
      
      if (bitValue) {
        PMTpairBitID::StatusBit_t const destBit
          = (compact && (bit >= 32))? bit - 8: bit;
        outputWord |= mask(destBit);
      }
      
    } // bit
    
  } // PMT wall
  
  return outputWords;
} // icarus::trigger::encodeTriggerLVDSbits()


// -----------------------------------------------------------------------------
