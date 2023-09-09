/**
 * @file   icaruscode/Decode/ChannelMapping/RunPeriods.h
 * @brief  Utilities to define run periods.
 * @date   August 8, 2023
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_RUNPERIODS_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_RUNPERIODS_H


// C/C++ standard libraries
#include <array>

// -----------------------------------------------------------------------------
namespace icarusDB {
  
  /// List of run periods.
  enum class RunPeriod: unsigned int {
    Runs0to2,       ///< Runs from the very beginning to the first part of Summer 2023.
    Run2shutdownB1, ///< A few runs in Summer 2023, break point 1 (August 23-29).
    Runs3andOn,     ///< Runs from last part of Summer 2023 on.
    
    // leave the following one as last, of course:
    NPeriods ///< Number of run periods.
  };
  
  /// Returns an array with all run periods. Also in `RunPeriods::All`.
  constexpr std::array<RunPeriod, static_cast<unsigned int>(RunPeriod::NPeriods)>
    allPeriods() noexcept;
  
  /// Returns the number of supported periods.
  constexpr std::size_t NPeriods() noexcept
    { return static_cast<std::size_t>(RunPeriod::NPeriods); }
  
  struct RunPeriods;
  
} // namespace icarusDB


// -----------------------------------------------------------------------------
inline constexpr std::array<icarusDB::RunPeriod, icarusDB::NPeriods()>
icarusDB::allPeriods() noexcept {
  std::array<RunPeriod, NPeriods()> periods{};
  for (unsigned int p = 0; p < NPeriods(); ++p)
    periods[p] = static_cast<RunPeriod>(p);
  return periods;
}


// -----------------------------------------------------------------------------
/**
 * @brief Definition of supported run periods.
 * 
 * This static class currently provides only the transformation from a run
 * number to a run period (`withRun()`).
 */
struct icarusDB::RunPeriods {
  
  static constexpr std::size_t NPeriods = icarusDB::NPeriods();
  
  /// The list of all supported periods.
  static constexpr std::array<RunPeriod, NPeriods> All = allPeriods();
  
  
  /**
   * @brief Returns the period of the specified `run`.
   * @return the period `run` belongs to, or `NPeriods` if not supported
   */
  static constexpr RunPeriod withRun(int run)
    {
      
      if (run < 10369) return RunPeriod::Runs0to2;
      // new PMT mapping (swapped some digitizer channels)
      if (run < 10441) return RunPeriod::Run2shutdownB1;
      // new PMT mapping (swapped more digitizer channels)
      if (run < 0xBADCAFE) return RunPeriod::Runs3andOn; // large value; reminds me of Fermilab
      
      return RunPeriod::NPeriods;
      
    }
  
}; // icarusDB::RunPeriods


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_RUNPERIODS_H
