/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AdderChannelSimulator.h
 * @brief  Simulation of adder board output.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AdderChannelSimulator.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELSIMULATOR_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELSIMULATOR_H

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelID.h" // AdderChannelID
#include "icaruscode/PMT/Trigger/Algorithms/AdderCalibrationDatabase.h"
#include "icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.h"
#include "icaruscode/PMT/Trigger/Algorithms/WaveformWithBaseline.h"
#include "icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulationTypes.h"
#include "icarusalg/Utilities/TimeInterval.h"
#include "icarusalg/Utilities/mfLoggingClass.h"

// // LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataalg/Utilities/quantities.h" // util::quantities::concepts
#include "lardataalg/Utilities/quantities/electromagnetism.h" // volt...
#include "lardataobj/RawData/OpDetWaveform.h"

// C/C++ standard libraries
#include <cstdint> // std::size_t
#include <complex>
#include <limits>
#include <memory>
#include <ratio> // std::milli
#include <string>
#include <type_traits> // std::enable_if_t
#include <utility> // std::move(), std::pair
#include <valarray>
#include <vector>


// -----------------------------------------------------------------------------
// forward declarations
class TObject;
class TDirectory;
class TGraph;
namespace icarus::trigger { class WaveformFilterAlg; } // Pimpl'ing this one

// -----------------------------------------------------------------------------
namespace icarus::trigger { class AdderChannelSimulator; }
/**
 * @brief Algorithm class simulating the output of an ICARUS adder board.
 * 
 * The algorithm takes as input the PMT waveforms input of the adder board
 * and simulates the output of the board.
 * 
 * The steps of the simulation are:
 * 
 * 1. determination of the time interval to operate on
 * 2. analogue sum of the PMT waveforms in that interval
 * 3. signal reshaping via a configured filter response
 * 4. calibration
 * 
 * The output signal is expressed in the same type of sampling as the input
 * waveforms, but it is represented in potential (millivolt) rather than ADC
 * counts.
 * 
 * Plots
 * ------
 * 
 * If a plot directory is configured, it will be filled with plots from the
 * adder channel being simulated:
 *  * adder-related waveform (as a ROOT `TGraph`):
 *      * `PMTsum`: PMT sum waveform;
 *      * `ReshapedPMTsum`: waveform out of the reshaping of PMT sum;
 *      * `SimulatedAdder`: waveform after calibration of reshaped PMT sum;
 * 
 *    If a step is skipped, the waveform plot is still produced and it will look
 *    the same as the previous one (for example, if no calibration was
 *    configured, `SimulatedAdder` plot will show the same waveform as
 *    `ReshapedPMTsum`).
 *    The naming of the graph objects in the directory (described in detail in
 *    `makeWaveformPlot()`) is `<waveform type>_<channel directory>_T<timestamp>`
 *    with `<waveform type>` the name of the plot from the list above,
 *    `<channel directory>` the ROOT name of the configured ROOT directory,
 *    and `<timestamp>` a time for the group of waveforms (the start of the
 *    simulation interval for all the waveforms; for `SimulatedAdder` that time
 *    also includes any time calibration), in format `0000.000` and in
 *    electronics time (microseconds).
 *    Waveforms are plotted in full resolution, and each is taking 8 times the
 *    memory compared to an input PMT waveform (`double` vs. 16-bit integer,
 *    and storing the time points too), meaning that the size of the output
 *    will result _larger_ than the raw data the simulation starts with.
 * 
 */
class icarus::trigger::AdderChannelSimulator
  : private icarus::ns::util::mfLoggingClass
{
  
    public:
  
  // imported types
  using AdderChannelID = icarus::trigger::AdderChannelID;
  
  using Voltage_t = adder::types::Voltage_t;
  using ADCsettings_t = adder::types::ADCsettings_t;
  using WaveformSamples_t = adder::types::WaveformSamples_t;
  
  // --- BEGIN ---  local aliases  ---------------------------------------------
  using electronics_time = detinfo::timescales::electronics_time;
  using WaveformWithBaseline = icarus::trigger::WaveformWithBaseline;
  using nanoseconds = util::quantities::intervals::nanoseconds;
  // ---- END ----  local aliases  ---------------------------------------------
  
  
  // --- BEGIN ---  data type definitions  -------------------------------------
  
  /// Time interval in electronics time.
  using TimeInterval_t = icarus::ns::util::TimeInterval<electronics_time>;
  // ---- END ----  data type definitions  -------------------------------------
  
  
  // --- BEGIN ---  interface data structures  ---------------------------------
  
  /// The output of the adder signal.
  struct AdderSignal_t {
    AdderChannelID channel; ///< Adder channel number.
    WaveformSamples_t samples; ///< Values of the samples in millivolt.
    raw::TimeStamp_t startTime; ///< Timestamp of the first sample.
    unsigned int nClipped = 0; ///< Count of clipped contributing waveforms.
    
    /**
     * @brief Returns a `raw::OpDetWaveform` object out of this record content.
     * @param ADCsettings settings used in volt-to-ADC conversion
     * @param baseline offset added to the waveform [ADC counts]
     * @return a `raw::OpDetWaveform` object out of this record content
     * 
     * Conversion to integral is made via plain C++ assignment (which performs
     * a truncation like `static_cast<short int>()`).
     */
    raw::OpDetWaveform makeOpDetWaveform
      (ADCsettings_t const& ADCsettings, Voltage_t baseline = 0) const;
    
  }; // AdderSignal_t
  
  // ---- END ----  interface data structures  ---------------------------------
  
  
  /**
   * @brief Constructor: configures the algorithm for this specific channel.
   * @param channel number of the adder board channel (used for tracking only)
   * @param ADCsettings relevant settings of the waveform ADC
   * @param splitFraction how much of the PMT signal is diverted to the adders
   * @param reshapingFilter filtering object for signal reshaping
   * @param calibration object for adder signal calibration
   * @param plotDir (optional) ROOT directory where to store debugging plots
   * @param logCategory (default: `"AdderChannelSimulator"`) logging stream name
   * 
   * If the reshaper algorithm pointer is null, no reshaping will be performed.
   * If the calibration database pointer is null, no calibration will be
   * applied.
   */
  AdderChannelSimulator(
    AdderChannelID channel,
    ADCsettings_t ADCsettings,
    double splitFraction,
    FilteringCircuit const* reshapingFilter,
    AdderCalibrationDatabase::RunCalibration const* calibration,
    TDirectory* plotDir = nullptr,
    std::string logCategory = "AdderChannelSimulator"
    );
  
  ~AdderChannelSimulator();
  
  /**
   * @brief Builds an adder waveform from all channels in the specified time window.
   * @param timeInterval the interval to constrain the waveform into
   * @param waveforms the list of waveforms (plus baseline) to include
   * @return the added waveform and the list of the contributing waveforms
   * 
   * The returned waveform covers the entire `timeInterval` and includes all
   * the waveforms in the input that contribute to that interval.
   * If a PMT channel has no waveform information in the interval or part of it,
   * it will not contribute (that includes also no noise contribution).
   * 
   * This function performs the following steps:
   *  1. addition (`addWaveformsInInterval()`);
   *  2. reshaping (`reshapeWaveform()`)
   *  3. amplitude calibration (`calibrateReshapedPMTsum()`)
   * 
   * See the documentation of those methods for additional details.
   * The returned list includes all and only the waveforms which have at least
   * one sample in the requested time interval.
   */
  std::pair<AdderSignal_t, std::vector<WaveformWithBaseline const*>> simulate(
    TimeInterval_t const& timeInterval,
    std::vector<WaveformWithBaseline const*> const& waveforms
    ) const;
  
  
  /// Returns the ADC settings in use.
  ADCsettings_t const& ADCsettings() const { return fADCsettings; }
  
    private:
  
  // privately import types
  using WaveformAmplitudes_t = adder::types::WaveformAmplitudes_t;
  
  
  // --- BEGIN ---  Configuration  ---------------------------------------------
  
  AdderChannelID const fChannel; ///< Adder channel identifier.
  
  ADCsettings_t const fADCsettings; ///< ADC configuration.
  
  double const fPMTsplitterScale; ///< Scale factor from PMT digitized signal to sum.
  
  TDirectory* const fPlotDir; ///< Directory where to store plots.
  
  // ---- END ----  Configuration  ---------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  FilteringCircuit const* fReshapingFilter; ///< Reshaping algorithm.
  
  /// Adder calibration algorithm.
  AdderCalibrationDatabase::RunCalibration const* fCalibration;
  
  /// Waveform filtering algorithm.
  std::unique_ptr<WaveformFilterAlg> const fFilterAlg;
  // That was totally cheating: the access to `fFilterAlg` is all but `const`
  // (the object has state that changes when filtering), but the _pointer_
  // to the algorithm is constant so we can use it in const-qualified methods
  // (like `WaveformFilterAlg* const`, unlike `WaveformFilterAlg const*`).
  
  // --- END Algorithms --------------------------------------------------------
  
  
  // --- BEGIN ---  Main algorithm functions  ----------------------------------
  /// @name Main algorithm functions
  ///@{
  /// Describes the start and end of a waveform in low-level form.
  struct WaveformExtent_t {
    raw::TimeStamp_t startTime = std::numeric_limits<raw::TimeStamp_t>::min();
    std::size_t nSamples = 0;
  }; // WaveformExtent_t
  
  
  /**
   * @brief Returns the sum of PMT waveforms limited to the specified interval.
   * @param waveforms the waveforms (and baselines) to add
   * @param timeInterval interval of definition of the final waveform
   * @param[out] pNClipped storage for the number of contributing channels clipped
   * @return the samples of the PMT sum waveform and the contributing waveforms
   * 
   * Each waveform has its baseline subtracted (negative "polarity" is assumed).
   * For each waveform, the intersection with `timeInterval` is added to the
   * result with no relative time shift and no signal scale factor.
   * Samples are added directly, assuming that the value of sample #0 of a
   * waveform starting at time _t_ has the same value `waveform[0]` in the time
   * interval [ _t_, &Delta; _t_ ], where &Delta; _t_ is the duration of a ADC
   * sampling tick (e.g. 2 ns).
   * For example, if the requested interval starts at 3.2 ns and an input
   * waveform starts with time 0.5 ns, the second sample of the waveform
   * (covering [2.5;4.5] ns for 2 ns sampling) will be added to the first sample
   * of the sum.
   * 
   * The returned waveform exactly covers the `timeInterval`.
   * Waveforms that do not cover the full `timeInterval` are considered to be
   * uniformly null in the interval they do not cover.
   * 
   * The function also returns the list of waveforms which contributed to at
   * least one sample of the result.
   * 
   * If `pNClipped` address is not null, the number of channels contributing
   * with clipped samples is stored in that memory address.
   * A clipped sample is a sample that reaches the limits of the ADC range,
   * a waveform is considered clipped if it contributed has at least one clipped
   * sample within `timeInterval`, and a channel is considered clipped if it
   * contributes with at least one clipped waveform.
   * For this purpose, the input waveforms are considered to have "negative
   * polarity", and therefore the only clipped limit value is `0`.
   * 
   * 
   * Precision
   * ----------
   * 
   * Computation is performed with `Voltage_t` precision.
   * Assuming the standard ICARUS adder, with 14 bit dynamic range and 15
   * contributing channels, 18 bits will be needed to preserve the full range.
   * A single precision type, with its 23 bit mantissa, guarantees no loss.
   * 
   * In the current implementation, the samples are immediately converted into
   * `Voltage_t` numbers.
   */
  std::pair<WaveformSamples_t, std::vector<WaveformWithBaseline const*>>
  addWaveformsInInterval(
    std::vector<WaveformWithBaseline const*> const& waveforms,
    TimeInterval_t const& timeInterval,
    unsigned int* pNClipped = nullptr
    ) const;
  
  /**
   * @brief Applies reshaping to the specified waveform.
   * @param PMTsum samples of the waveform to be reshaped
   * @return reshaped waveform samples in the same time interval as `PMTsum`
   * 
   * The specified waveform is applied a transfer function in
   * frequency domain. The result is converted back to time domain.
   * 
   * The conversion between frequency and time domains is performed via
   * Fourier transform on the full waveform. As a consequence,
   * the lowest frequency and the number of frequencies included in the
   * transform depends on the length of the waveform.
   */
  WaveformSamples_t reshapeWaveform(WaveformSamples_t PMTsum) const;
  
  /**
   * @brief Applies a calibration to the specified waveform.
   * @param samples samples of the waveform to be calibrated
   * @param channel adder channel to use the calibration of
   * @return calibrated waveform samples
   * 
   * The calibration of the specified adder `channel` is applied to the input
   * waveform.
   * The calibration is a single factor, uniformly applied to all samples of
   * the waveform.
   * 
   */
  WaveformSamples_t calibrateReshapedPMTsum
    (WaveformSamples_t PMTsum, AdderChannelID channel) const;
  
  /// @}
  // ---- END ----  Main algorithm functions  ----------------------------------
  
  
  // --- BEGIN ---  Helper functions  ------------------------------------------
  /// @name Helper functions
  /// @{

  /**
   * @brief Creates and writes a ROOT graph from the specified waveform samples.
   * @tparam Samples type of container for the samples
   * @param startTime time of the start of the first sample
   * @param samples values of the samples
   * @param typeName name (tag) of the type of waveform
   * @param descr description of the type of waveform
   * @param write (default: `true`) immediately write the graph into the directory
   * @return a pointer to the created graph
   * 
   * A new graph is created with a name based on `typeName`, plus the name
   * of the current directory and a tag of the timestamp.
   * The title is also made out of the title of the current directory,
   * in addition to the `descr` text.
   * 
   * The time axis is constructed based on the configured tick duration and to
   * the `startTime`, which is on the electronics scale and presumably in
   * microseconds.
   * 
   * If the title does not contain any axis labels, some are made up assuming
   * the ordinate axis shows a voltage in millivolt.
   * 
   * The graph is written in the plot directory (if no plot directory was
   * configured, this function immediately exits).
   */
  template <typename Samples>
  std::unique_ptr<TGraph> makeWaveformPlot(
    electronics_time startTime, Samples const& samples,
    std::string const& typeName = "Waveform", std::string const& descr = "",
    bool write = true
    ) const;
  
  
  /// Write the specified ROOT object into the configured plot directory.
  /// @return the value returned by `TObject::Write()`, `0` if no directory
  int writePlot(TObject* obj) const;

  /// @}
  // ---- END ----  Helper functions  ------------------------------------------
  
    
  /// Returns the extend of the waveform.
  static WaveformExtent_t makeWaveformExtent
    (WaveformWithBaseline const& waveform);
  
  /// Returns the extend of the waveform.
  static WaveformExtent_t makeWaveformExtent(raw::OpDetWaveform const& waveform);
  

}; // icarus::trigger::AdderChannelSimulator


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERCHANNELSIMULATOR_H
