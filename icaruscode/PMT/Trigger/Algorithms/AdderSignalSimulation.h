/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulation.h
 * @brief  Provides a adder simulation manager.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulation.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALSIMULATION_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALSIMULATION_H


// ICARUS/SBN libraries
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelSimulator.h"
#include "icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulationTypes.h" // ADCsettings_t
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelMaps.h"
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelID.h"
#include "icaruscode/PMT/Trigger/Algorithms/WaveformWithBaseline.h"
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "icarusalg/Utilities/GroupByIndex.h"
#include "icarusalg/Utilities/TimeInterval.h"

// LArSoft/framework libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// // C/C++ standard libraries
#include <set>
#include <stdexcept> // std::runtime_error (protocol only)
#include <string>
#include <utility> // std::pair
#include <vector>


// -----------------------------------------------------------------------------
// forward declarations
class TDirectory;

// -----------------------------------------------------------------------------
namespace icarus::trigger { class AdderSignalSimulation; }
/**
 * @brief Algorithm coordinating the simulation of adder channels.
 * 
 * This algorithm produces a waveform for each adder board and each PMT trigger
 * primitive in the event (see "Timing" below).
 * 
 * The input is PMT optical waveforms: this algorithm composes them into added
 * analogue signals, reshapes and finally "calibrates" them:
 * 
 * 1. Waveforms are assigned to an adder channel according to the specified
 *    mapping, and then grouped by time.
 * 2. For each group of PMT waveforms, a simulation interval is defined on the
 *    time interval overlapping all non-missing PMT waveforms; if a channel is
 *    listed as missing (see below) it is used if present but not required.
 *    In the coincidence of an appropriate mapping (all PMT channels to an adder
 *    coming from the same readout board) and readout scheme (by readout board
 *    trigger), this step is trivial as all PMT waveforms completely overlap.
 *    This is the case of ICARUS set up.
 * 3. The PMT signals are baseline-subtracted, converted into voltage and summed
 *    into a single waveform (not necessarily in this order of steps).
 *    The signal is inverted: since PMT waveforms are expected to be with
 *    "negative polarity", the PMT sum will have positive polarity.
 * 4. A reshaping is applied on each PMT sum waveform if configured.
 * 5. A calibration is applied on each reshaped PMT sum if configured.
 * 6. The resulting waveforms are converted back to ADC counts, added an offset
 *    (from `AdderBaseline`), their samples converted to integer (by truncation:
 *    see `AdderChannelSimulator::AdderSignal_t::makeOpDetWaveform()`) and
 *    stored into LArSoft-standard `raw::OpDetWaveform` objects (which are based
 *    on signed ADC count type).
 * 
 * 
 * Mapping
 * --------
 * 
 * Each adder board and its 15 contributing PMT channels are connected by a
 * channel map. This algorithm requires such map to be provided by the caller
 * at `setup()` (or construction) time. A mapping is mandatory for simulation.
 * The algorithm `icarus::trigger::AdderChannelMapBuilder` can build one.
 * 
 * 
 * Timing
 * -------
 * 
 * A "PMT trigger primitive" is a signal sent to the PMT readout boards to
 * record (freeze and send) the digitized PMT data. On each primitive, all the
 * PMT in a readout board are recorded and stamped with the same time (which
 * can be recalibrated later on during reconstruction).
 * 
 * The algorithm identifies groups of PMT waveforms overlapping in time, then
 * assigns them to adder boards for independent simulation.
 * The adders are simulated in the time of overlap among all the contributing
 * PMT waveforms. The details of this procedure are described in
 * `findOverlappingTimeIntervals()`.
 * 
 * @note It is important for this procedure that the input waveform time is
 *       as consistent as possible among PMT channels. PMT channels on the same
 *       readout board are originally assigned the same consistent time.
 *       The standard mapping of PMT channels into adders in ICARUS Run3 (2024)
 *       and later routes all the PMT from one readout board into a single adder
 *       board, which is enough in principle to preserve for such consistency;
 *       but per-channel corrections have the chance of breaking that.
 *       The correctness of the algorithm implementation for odd waveform sizes
 *       that may result from misaligned PMT waveform times has not been
 *       verified (as of `icaruscode` `v10_06_01`). Here, "odd" is literally
 *       meant with the mathematical concept on the number of samples in the
 *       waveform, and the Fourier transform library can either take advantage
 *       or require that number to be even or a multiple of some power of 2.
 * 
 * 
 * Missing channels
 * -----------------
 * 
 * The algorithm allows for the specification of a list of missing PMT channels.
 * A missing channel is in general a channel whose output should be ignored.
 * In practice, depending on the type of input and of channel issue, the channel
 * may be completely absent from the input data product, be present with only
 * noise, normal or high (typical of detector data), or be present with physical
 * activity (in simulation where these channels are not masked).
 * 
 * The correct treatment of these channels depends on the input and this
 * algorithm is not spending time in figuring it out.
 * The current implementation is that waveforms from a "missing" channels, as
 * flagged by the configuration parameter `MissingChannels`, are used if present
 * but not required (as opposed to throwing an exception when an expected and
 * good channel of an adder is not found).
 * 
 * 
 * On calibration
 * ---------------
 * 
 * The calibration is bridging the gap between the first principle, approximate
 * simulation and the observed data.
 * 
 * While the interface allows for any calibration algorithm, only one is
 * implemented in `icaruscode` `v010_06_01`:
 * `icarus::trigger::AmplitudeAdderCalibration`.
 * This calibration just rescales the whole waveform so that the peak amplitude
 * reflects the one in data.
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * Enabling the writing of plots into the ROOT directory (by passing a non-null
 * pointer to one in the `simulate()` call) is not compatible with
 * multithreading.
 * 
 * 
 * Input
 * ======
 * 
 * Event data
 * -----------
 * 
 * The only physical input required is a single waveform for each recorded
 * optical detector activity; each waveform comes from a single channel, and
 * there may be multiple waveforms on the same channel at different times.
 * The time stamp is expected to be on the
 * @ref DetectorClocksElectronicsTime "electronics time scale"
 * and therefore expressed in microseconds. For each waveform, a precomputed
 * baseline is also required. This input is collected in the
 * `WaveformWithBaseline` object, which hosts a reference (pointer) to the
 * original data, which must be valid and available during all the simulation.
 * 
 * 
 * Additional information/algorithms
 * ----------------------------------
 * 
 * * Adder-to-PMT channel mapping (e.g. from `AdderChannelMapBuilder`).
 * * Digitizer settings (default values should be appropriate for CAEN V1730).
 * * Reshaping algorithm (essential although technically optional).
 * * Calibration algorithm (important but technically optional).
 * 
 * 
 * Output
 * =======
 * 
 * The call to `simulate()` returns a collection of waveforms representing
 * the simulated adder board response.
 * The waveforms are "reproduced" from the already digitized PMT waveforms, in
 * contrast with the hardware which performs an analogue sum.
 * The channel number is unique on each of them, and it is following a
 * convention documented in `icarus::trigger::AdderChannelMaps`.
 * Note that the baseline is expected to be `0` by construction.
 * 
 * 
 * Plots
 * ------
 * 
 * If a `TDirectory` is passed to `simulate()`, intermediate waveforms are
 * written into that directory as ROOT `TGraph` objects.
 * The plots are organised in one level of subdirectories: each subdirectory
 * includes the waveforms of a single channel.
 * See the "Plots" section of `simulate()` documentation for details.
 * 
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * The part of the configuration that is static is performed via FHiCL at
 * construction. The rest of the configuration ("dynamic", may in principle
 * change run by run or even event by event) is performed via the constructor
 * or `setup()` calls.
 * 
 * The FHiCL configuration (`Config`) includes:
 * 
 * * `MissingChannels` (list of integers): list of the PMT channel numbers
 *   (as reported in the `WaveformTag` input collection) flagged as "missing".
 *   Note that in the hardware bad channels may be included in the adder input.
 * * `SplitToAdders` (real, default: 5%): the fraction of PMT signal sent to the
 *   adder boards. It is assumed that the PMT signal in input received
 *   the complementary fraction of signal.
 * * `AdderBaseline` (real, default: `0`): offset in ADC counts to be added to
 *   the adder waveforms before rounding to integral sample.
 * * `LogCategory` (string, default: `"SimulateAdderSignal"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * 
 */
class icarus::trigger::AdderSignalSimulation
  : private icarus::ns::util::mfLoggingClass
{
  
  // imported types
  using nanoseconds = util::quantities::intervals::nanoseconds;
  using electronics_time = detinfo::timescales::electronics_time;
  using TimeInterval_t = icarus::ns::util::TimeInterval<electronics_time>;
  
  /// Map of waveform information records per channel.
  using WaveformsByChannel_t
    = icarus::ns::util::GroupByIndex<WaveformWithBaseline>;
  
  /// Special value to indicate no channel number information.
  static constexpr AdderChannelID NoChannel{};
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  std::vector<raw::Channel_t> const fMissingChannels; ///< Channels to skip.
  
  double const fSplitToAdders; ///< Fraction of PMT signal diverted to adders.
  
  float const fAdderBaseline; ///< Offset added to all samples in adder output.
  
  adder::types::ADCsettings_t const fADCsettings; ///< Digitizer settings.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------

  /// Mapping of adder to PMT channels.
  AdderChannelMap const* fAdderChannels = nullptr;

  /// Reshaping algorithm.
  FilteringCircuit const* fReshapingFilter = nullptr;
  
  /// Adder calibration algorithm.
  AdderCalibrationDatabase::RunCalibration const* fCalibration = nullptr;
  
  // --- END Algorithms --------------------------------------------------------
  
  /// Returns the sampling period of the digitizer.
  nanoseconds opticalTick() const { return fADCsettings.samplingTime; }
  
  
  // --- BEGIN ---  Algorithm implementation functions  ------------------------
  /// @name Algorithm implementation functions
  /// @{
  
  /// Returns a map channel -> waveforms on that channel.
  /// If `waveformInfo` changes, the returned map becomes invalid.
  static WaveformsByChannel_t groupByChannel
    (std::vector<WaveformWithBaseline> const& waveformInfo);
  
  /**
   * @brief Returns the waveforms contributing to a adder channel.
   * @param waveformInfoByChannel all available waveforms, grouped by channel
   * @param windowChannels list of PMT channels included in the adder channel
   * @param adderChannel (default: `NoChannel`) adder channel identifier
   * @return the list of contributing waveforms, and the non-missing channels
   * 
   * All waveforms are provided as a `WaveformWithBaseline` object (that is,
   * actual waveform plus its baseline information).
   * 
   * The input `waveformInfoByChannel` is a map type that takes a PMT channel
   * number as key and returns a list of waveforms for that channel.
   * Waveforms from channels not in `windowChannels` are ignored.
   * 
   * The algorithm returns the list of waveforms associated to `windowChannels`,
   * still grouped by channel, but in a straight sequence (so the index of the
   * element is not the PMT channel any more, and it is in fact carrying no
   * information).
   * In addition, it returns which of the channels in the window are not marked
   * as "missing", and therefore expected to have meaningful data.
   * However, also waveforms from "missing" channels are included in the
   * returned list if present.
   */
  std::pair<std::vector<std::vector<WaveformWithBaseline const*>>, std::set<raw::Channel_t>>
  collectAdderWindowWaveforms(
    WaveformsByChannel_t const& waveformInfoByChannel,
    icarus::trigger::AdderChannelInfo_t::PMTchannelList_t const& windowChannels,
    AdderChannelID adderChannel = NoChannel
    ) const;
  
  /**
   * @brief Returns the time intervals where all channels have information.
   * @param waveformsPerChannel input waveforms, grouped by channel
   * @param requiredChannels list of channels required to be present
   * @return sequence of pairs: the time intervals and the contributing waveforms
   * 
   * Each of the returned time intervals is entirely covered by _all_ the
   * required channels (and possibly some more). A channel "covers" an interval
   * if it has a waveform that spans at least the whole interval, without gaps;
   * the waveform may contain light signals, only electronics noise or even be
   * constant at `0` mV, and it will still "cover" the entire time interval where
   * it is defined.
   * 
   * Each entry of the input list `waveformsPerChannel` is a list of all
   * waveforms (with baseline) on the same channel. Neither the order nor the
   * index meaning for the channels are constrained.
   * Empty entries (i.e. empty lists) will be just ignored.
   * However, within each channel the waveforms are assumed to be sorted in
   * increasing time.
   * 
   * For each returned item, the time interval is included together with all the
   * input waveforms that cover it completely. By construction, in each time
   * interval there will be only one waveform from each of the channels.
   * The required channels, which define the overlap interval, will always
   * contribute one waveform each. A non-required channel will contribute a
   * waveform only if that waveform has full coverage of the interval.
   * 
   * It may happen that the same waveform contributes to different intervals.
   */
  std::vector<std::pair<TimeInterval_t, std::vector<WaveformWithBaseline const*>>>
  findOverlappingTimeIntervals(
    std::vector<std::vector<WaveformWithBaseline const*>> const& waveformsPerChannel,
    std::set<raw::Channel_t> const& requiredChannels
    ) const;
  
  /// Returns whether `channel` is in the configured missing channel list.
  bool isMissingChannel(raw::Channel_t channel) const;
  
  //@{
  /// Returns the time interval covered by the `waveform`.
  TimeInterval_t waveformInterval(raw::OpDetWaveform const& waveform) const;
  TimeInterval_t waveformInterval(WaveformWithBaseline const& waveform) const;
  //@}
  
  /// @}
  // ---- END ----  Algorithm implementation functions  ------------------------
  
    public:
  
  /// FHiCL configuration of the module.
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<raw::Channel_t> MissingChannels{
      Name{ "MissingChannels" },
      Comment("do not include the specified channels [default: all included]"),
      std::vector<raw::Channel_t>{}
      };
    
    fhicl::Atom<double> SplitToAdders{
      Name{ "SplitToAdders" },
      Comment{ "fraction of PMT signal diverted to the adder boards" },
      0.05
      };
    
    fhicl::Atom<float> AdderBaseline{
      Name{ "AdderBaseline" },
      Comment{ "baseline (offset) of the output adder signal [ADC#]" },
      0.0f
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "tag of the module output to console via message facility" },
      "SimulateAdderSignal"
      };
    
  }; // Config
  
  
  /**
   * @brief Constructor: algorithm configuration via FHiCL.
   * @param config FHiCL configuration
   * @param ADCsettings (default: use CAEN 1730B's) settings of digitizer ADC
   * @param adderChannels mapping of adder channels to their contributing PMTs
   * @param reshapingFilter filtering algorithm
   * @param calibration calibration algorithm
   * 
   * Configuration is read once and for all.
   * If the algorithms are specified, they will be integrated (equivalent to a
   * call to `setup()`), otherwise they can be `setup()` later.
   */
  AdderSignalSimulation(
    Config const& config,
    adder::types::ADCsettings_t ADCsettings = adder::types::ADCsettings_t{},
    icarus::trigger::AdderChannelMap const* adderChannels = nullptr,
    FilteringCircuit const* reshapingFilter = nullptr,
    AdderCalibrationDatabase::RunCalibration const* calibration = nullptr
    );
  
  /**
   * @brief Set up (i.e. dynamic configuration) of all algorithms.
   * @param adderChannels mapping of adder channels to their contributing PMTs
   * @param reshapingFilter filtering algorithm
   * @param calibration calibration algorithm
   * 
   * The specified algorithms will be used from now on.
   * Algorithm objects are owned by the caller and they are expected to stay
   * valid for the duration of this object lifetime, or up to the next
   * `setup()` call.
   * If an algorithm is `nullptr`, the corresponding functionality is not
   * executed: no reshaping if `reshapingFilter` is `nullptr` (PMT signal sum
   * is used directly), no calibration if `calibrationDB` is `nullptr`
   * (reshaped signals are used as they come out of the simulation).
   * The `adderChannel` map is, instead, mandatory.
   */
  void setup(
    icarus::trigger::AdderChannelMap const* adderChannels,
    FilteringCircuit const* reshapingFilter,
    AdderCalibrationDatabase::RunCalibration const* calibration
    );
  
  
  // --- BEGIN ---  Configuration query  ---------------------------------------
  /// @name Configuration query
  /// @{
  
  /// Returns the current reshaping filter algorithm (it may be `nullptr`).
  FilteringCircuit const* reshapingFilter() const { return fReshapingFilter; }
  
  /// Returns the current calibration algorithm (it may be `nullptr`).
  AdderCalibrationDatabase::RunCalibration const* calibration() const
    { return fCalibration; }
  
  /**
   * @brief Dumps the current configuration.
   * @param indent indentation on the start of each line
   * @param firstIndent indentation on the start of the first line only
   * @return a human-readable string describing the algorithm configuration
   * 
   * The string starts straight into the first line, and does not terminate the
   * last one (it is safe to assume the line where this dump ends needs to be
   * explicitly terminated).
   */
  std::string configurationDump
    (std::string const& indent, std::string const& firstIndent) const;
  
  /// Dumps the current configuration (first line indent same as others).
  std::string configurationDump(std::string const& indent) const
    { return configurationDump(indent, indent); }
  
  /// @}
  // ---- END ----  Configuration query  ---------------------------------------
  
  /**
   * @brief Simulates adder signals from the whole specified dataset.
   * @param waveformInfo collection of input PMT digitized waveforms
   * @param plotDir (default: `nullptr`) directory where to put debug plots
   * @return a list of simulated adder waveforms
   * @throw std::runtime_error in case of misconfiguration or incomplete setup
   * 
   * The main input, `waveformInfo`, consists of an unsorted list of waveforms
   * from contributing PMT, each one associated to its precomputed pedestal
   * value (baseline), both expressed in ADC counts on the digitizer scale
   * (the baseline may be still a real number as result of a computation).
   * The `WaveformWithBaseline` objects only hold pointers to the actual data
   * (including the baseline), so the original data needs to be available for
   * the duration of the execution of `simulate()`.
   * 
   * If `plotDir` is specified, plots are produced and put into it.
   *
   * The order of the simulated waveforms in the output is not prescribed.
   *
   * The simulation is outlined in the class documentation.
   *
   * Plots
   * ------
   * 
   * One subdirectory is created for each adder channel.
   * Subdirectories have the name pattern: `<parent name>_CH<adder channel ID>`,
   * where `<parent name>` is the ROOT name of the parent directory (the one
   * passed to `simulate()`) and the `<adder channel ID>` is in hexadecimal
   * (e.g. `0XA002`).
   * 
   * Within this directory, plots are saved as described in 
   * `icarus::trigger::AdderChannelSimulator`.
   * 
   */
  std::vector<raw::OpDetWaveform> simulate(
    std::vector<WaveformWithBaseline> const& waveformInfo,
    TDirectory* plotDir = nullptr
    ) const;
  
}; // icarus::trigger::AdderSignalSimulation


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALSIMULATION_H
