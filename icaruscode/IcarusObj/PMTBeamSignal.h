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

  /// Special value to denote no special channel information.
  static constexpr auto NoChannel = std::numeric_limits<unsigned int>::max();
  /// Special value to denote no time channel information.
  static constexpr double NoTime = std::numeric_limits<double>::max();
  // Special value to denote no sample information.
  static constexpr std::size_t NoSample = 0;

  /**
   * @brief Beam time as seen by a PMT readout board.
   *
   * This could either be an early warning (EW) or a resistive wall monitor (RWM) time.
   * These signals are delivered via fibers and digitized in special PMT channels.
   *
   * Both the time in @ref DetectorClocksElectronicsTime "electronics time scale"
   * and time time relative to the hardware trigger are included.
   *
   * The information in this object may be missing: its validity should
   * always be checked in advance with `isValid()`.
   */

  struct PMTBeamSignal
  {

    /// Special channel this time was extracted from.
    /// These are defined in `CAEN_V1730_setup_icarus.fcl`.
    unsigned int specialChannel = NoChannel;

    /// Board on which the special channel is on (e.g: WW-TOP-A).
    /// Should match the same format as `icarusDB::PMTChannelInfo_t::digitizerLabel`.
    std::string digitizerLabel = "";

    /// Crate this time applies to (e.g.: WW-TOP).
    /// Corresponds to the first part of `digitizerLabel`.
    std::string crate = "";

    /// Sample within the waveform where the reference signal is found.
    std::size_t sample = NoSample;

    /// Start time in electronics time [us].
    double startTimeAbs = NoTime;

    /// Start time relative to trigger time [us].
    double startTime = NoTime;

    /// Returns whether the time is valid.
    bool isValid() const { return (sample != NoSample); }
  };

} // namespace icarus::timing

#endif // ICARUSCODE_ICARUSOBJ_PMTBEAMSIGNAL_H
