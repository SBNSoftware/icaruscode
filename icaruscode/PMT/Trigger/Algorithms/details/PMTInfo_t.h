/**
 * @file    icaruscode/PMT/Trigger/Algorithms/details/PMTInfo_t.h
 * @brief   Helper class to store discriminated PMT information.
 * @authors Ryan Howell (rhowell3@ur.rochester.edu)
 *          Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    January 28, 2021
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_PMTINFO_T_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_PMTINFO_T_H

// ICARUS libraries
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t, ...

// C++ standard libraries
#include <vector>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::trigger::details { struct PMTInfo_t; }

// -----------------------------------------------------------------------------
/// Helper data structure to store PMT activity information in the event.
struct icarus::trigger::details::PMTInfo_t {
  
  using ChannelID_t = raw::Channel_t; ///< Type to represent a channel ID.
  
  using ADCcount_t = raw::ADC_Count_t; ///< Type to represent threshold.
  
  
  // --- BEGIN -- Construction -------------------------------------------------
  
  PMTInfo_t() = default; // no trigger
  PMTInfo_t(ADCcount_t threshold, std::vector<ChannelID_t> activeChannels);
  
  // --- END -- Construction ---------------------------------------------------
  
  
  // --- BEGIN -- Access to PMT information ------------------------------------
  /// @name Access to PMT information
  /// @{
  
  /// The threshold this PMT data is extracted with.
  ADCcount_t threshold() const { return fThreshold; }
  
  /// Returns the list of channels with activity above threshold.
  std::vector<ChannelID_t> const& activeChannels() const
    { return fActiveChannels; }
  
  /// @}
  // --- END -- Access to PMT information --------------------------------------
  
    private:
  
  ADCcount_t fThreshold { 0 }; ///< Discrimination threshold.
  
  /// Channels whose activity is above threshold.
  std::vector<ChannelID_t> fActiveChannels;
  
}; // icarus::trigger::details::PMTInfo_t


// -----------------------------------------------------------------------------
// ---  Inline implementation
// -----------------------------------------------------------------------------
icarus::trigger::details::PMTInfo_t::PMTInfo_t
  (ADCcount_t threshold, std::vector<ChannelID_t> activeChannels)
  : fThreshold(threshold)
  , fActiveChannels(std::move(activeChannels))
{}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_DETAILS_PMTINFO_T_H
