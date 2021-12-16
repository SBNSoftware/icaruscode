/**
 * @file   icaruscode/PMT/Algorithms/PMTcoverageInfoUtils.h
 * @brief  Writes a collection of sbn::PMTcoverageInfo from PMT waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 22, 2021
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_PMTCOVERAGEINFOUTILS_H
#define ICARUSCODE_PMT_ALGORITHMS_PMTCOVERAGEINFOUTILS_H


// SBN libraries
#include "icaruscode/IcarusObj/PMTcoverageInfo.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <optional>


// -----------------------------------------------------------------------------
// forward declarations
namespace detinfo { class DetectorTimings; }


// -----------------------------------------------------------------------------
namespace sbn {
  
  // --- BEGIN -- Creation of sbn::PMTcoverageInfo from raw::OpDetWaveform -----
  /**
   * @name Creation of `sbn::PMTcoverageInfo` from `raw::OpDetWaveform`
   * 
   * The creation of summary objects `sbn::PMTcoverageInfo` from optical
   * detector waveforms is possible with two options:
   * 
   * * one shot: call to convert a single waveform (`makePMTcoverageInfo()`);
   * * bulk: converter object reused for multiple conversions
   *   (`sbn::PMTcoverageInfoMaker`).
   * 
   * For usage examples, see their respective documentation.
   */
  /// @{
  class PMTcoverageInfoMaker;
  
  /**
   * @brief Creates a `sbn::PMTcoverageInfo` out of a `raw::OpDetWaveform`.
   * @param waveform the input waveform
   * @param detTimings timing service provider
   * @return a `sbn::PMTcoverageInfo` object with the summary information
   * 
   * Returns a new summary object extracted from the `waveform`.
   * 
   * The timing information is used to determine whether the global trigger
   * and beam times are included.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings
   *   (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event));
   * 
   * sbn::PMTcoverageInfo info = sbn::makePMTcoverageInfo(waveform, detTimings);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  PMTcoverageInfo makePMTcoverageInfo(
    raw::OpDetWaveform const& waveform,
    detinfo::DetectorTimings const& detTimings
    );
  
  
  /**
   * @brief Creates a `sbn::PMTcoverageInfo` out of a `raw::OpDetWaveform`.
   * @param waveform the input waveform
   * @param opDetTickPeriod period of the optical detector digitizer
   * @return a `sbn::PMTcoverageInfo` object with the summary information
   * 
   * Returns a new summary object extracted from the `waveform`.
   * 
   * Information requiring timing is not saved.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * util::quantities::intervals::nanoseconds const opDetPeriod { 2.0 };
   * 
   * sbn::PMTcoverageInfo info = sbn::makePMTcoverageInfo(waveform, opDetPeriod);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  PMTcoverageInfo makePMTcoverageInfo(
    raw::OpDetWaveform const& waveform,
    util::quantities::intervals::microseconds opDetTickPeriod
    );
  
  /// @}
  // --- END ---- Creation of sbn::PMTcoverageInfo from raw::OpDetWaveform -----
  
  
} // namespace sbn


// -----------------------------------------------------------------------------
/**
 * @brief Converter from `raw::OpDetWaveform` into `sbn::PMTcoverageInfo`.
 * 
 * An object of this class is initialized once with some timings (e.g. once per
 * event), used to `make()` multiple `sbn::PMTcoverageInfo` and then discarded:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * 
 * detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings
 *   (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event));
 * 
 * std::vector<sbn::PMTcoverageInfo> PMTinfo;
 * for (raw::OpDetWaveform const& waveform: waveforms)
 *   PMTinfo.push_back(makePMTcoverageInfo(waveform));
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * It supports a scenario where times (trigger and beam gate) are not available
 * or not relevant, in which case only the duration in time of a optical
 * detector waveform tick is needed.
 */
class sbn::PMTcoverageInfoMaker {
  
    public:
  
  using microseconds = util::quantities::intervals::microseconds;
  
  /// Constructor: allows creation of `sbn::PMTcoverageInfo` with full
  /// information.
  PMTcoverageInfoMaker(detinfo::DetectorTimings const& detTimings);
  
  
  /// Constructor: allows creation of `sbn::PMTcoverageInfo` with no
  /// trigger/beam time information.
  PMTcoverageInfoMaker(microseconds opDetTickPeriod);
  
  //@{
  /// Creates a `sbn::PMTcoverageInfo` out of the specified `waveform`.
  sbn::PMTcoverageInfo make(raw::OpDetWaveform const& waveform) const;
  sbn::PMTcoverageInfo operator() (raw::OpDetWaveform const& waveform) const
    { return make(waveform); }
  //@}
  
    private:
  
  using electronics_time = detinfo::timescales::electronics_time;
  
  microseconds fOpDetTickPeriod; ///< The duration of a optical detector tick.
  
  std::optional<electronics_time> fTriggerTime; ///< Cached trigger time.
  std::optional<electronics_time> fBeamGateTime; ///< Cached beam gate time.
  
}; // sbn::PMTcoverageInfoMaker


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_PMTCOVERAGEINFOUTILS_H
