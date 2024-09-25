/**
 * @file   icaruscode/IcarusObj/PMTBeamSignal.h
 * @brief  Holds the event-by-event RWM or EW times
 * @author Matteo Vicenzi (mvicenzi@bnl.gov)
 * @date   March 14 2024
 */

#ifndef ICARUSCODE_ICARUSOBJ_PMTBEAMSIGNAL_H
#define ICARUSCODE_ICARUSOBJ_PMTBEAMSIGNAL_H

// C/C++ standard libraries
#include <limits>
#include <string>
#include <cstddef>

namespace icarus::timing
{

  /// Special value to denote no special channel information
  static constexpr auto NoChannel = std::numeric_limits<unsigned int>::max();
  /// Special value to denote no time channel informatio
  static constexpr double NoTime = 0.;
  // Special value to denote no sample information
  static constexpr std::size_t NoSample = 0;

  struct PMTBeamSignal
  {

    /// Special channel this time was extracted from
    unsigned int specialChannel = NoChannel;

    /// Board on which the special channel is on
    std::string digitizerLabel = "";

    /// Crate this time applies to
    std::string crate = "";

    /// Sample within the waveform where the reference signal is found
    std::size_t sample = NoSample;

    /// Start time in electronics time [us]
    double startTimeAbs = NoTime;

    /// Start time relative to trigger time [us]
    double startTime = NoTime;

    PMTBeamSignal(unsigned int ch, std::string b, std::string c,
                  std::size_t s, double t, double tt) : specialChannel(ch), digitizerLabel(b), crate(c), sample(s),
                                                        startTimeAbs(t), startTime(tt) {};

    PMTBeamSignal() {};

    /// Returns whether the time is valid.
    bool isValid() const { return ((sample != NoSample) && (startTime != NoTime)); }
  };

} // namespace icarus::timing

#endif // ICARUSCODE_ICARUSOBJ_PMTBEAMSIGNAL_H
