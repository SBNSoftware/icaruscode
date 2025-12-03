/**
 * @file   icaruscode/PMT/Trigger/Algorithms/AmplitudeAdderCalibration.h
 * @brief  Calibration for adder board output based on signal amplitude.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 5, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/AmplitudeAdderCalibration.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALCALIBRATION_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALCALIBRATION_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/AdderCalibrationDatabase.h"

// LAroft libraries
#include "lardataalg/Utilities/quantities_fhicl.h" // for voltage parameters
#include "lardataalg/Utilities/quantities/electromagnetism.h" // volt, millivolt
#include "lardataalg/Utilities/intervals_fhicl.h" // for time parameters

// framework libraries
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <algorithm> // std::minmax_element
#include <array>
#include <cmath> // std::sqrt()
#include <memory> // std::unique_ptr
#include <string>
#include <vector>


// -----------------------------------------------------------------------------
namespace icarus::trigger { class AmplitudeAdderCalibration; }
/**
 * @brief Adder calibration database with a single period and data from FHiCL.
 * 
 * This simple adder configuration database supports a single set of linear
 * calibration curves, one per channel.
 * The curve is parametrized by waveform amplitude, expressed in volt from the
 * baseline.
 * The amplitude is simply computed as highest sample value of the waveform.
 * 
 */
class icarus::trigger::AmplitudeAdderCalibration
  : public icarus::trigger::AdderCalibrationDatabase
{
  
    public:
  
  // import type
  using millivolt = util::quantities::millivolt;
  
  
  // --- BEGIN ---  Algorithm configuration records  ---------------------------
  /// @name Algorithm configuration records
  /// @{
  
  /**
   * @brief Linear calibration for a specific channel.
   *
   * The calibration "curve" describes the calibrated peak amplitude as function
   * of the original one.
   * The returned calibration `factor()` describes the ratio between the two,
   * while `computeAmplitude()` and `computeUncertainty()` calibrate a specific
   * amplitude.
   */
  struct ChannelCalibration_t {
    
    /// Type for calibration curve parameters.
    using par_matrix_t = std::array<double, 2U>;
    
    /// Type for calibration curve parameter covariance matrix.
    using cov_matrix_t = std::array<par_matrix_t, 2U> ;
    
    AdderChannelID channel; ///< Channel the calibration applies to.
    
    /// Parameters: [0] offset [mV]  [1] slope
    par_matrix_t parameters = { 0.0, 0.0 };
    
    /// Covariance matrix: [0] offset [mV]  [1] slope
    cov_matrix_t covariance = { 0.0, 0.0, 0.0, 0.0 };
    
    /// Time interval to add to the reshaped waveform.
    microseconds timeOffset { 0.0 };
    
    /// Uncertainty on the time interval to add to the reshaped waveform.
    microseconds timeOffsetErr { 0.0 };
    
    /// Returns the slope of the calibration line.
    double slope() const { return parameters[1]; }
    
    /// Returns the offset of the calibration line at 0 amplitude.
    millivolt offset() const { return millivolt{ parameters[0] }; }
    
    /// Returns the uncertainty on the `slope()`.
    double slopeErr() const { return std::sqrt(slopeVar()); }
    
    /// Returns the uncertainty on the `offset()`.
    millivolt offsetErr() const { return millivolt{ std::sqrt(offsetVar()) }; }
    
    /// Returns the variance on the `offset()` [mV^2]
    double offsetVar() const { return covariance[0][0]; }
    
    /// Returns the variance on the `slope()`.
    double slopeVar() const { return covariance[1][1]; }
    
    /// Returns the covariance between `offset()` and `slope()` [mV]
    double offsetSlopeCovar() const { return covariance[0][1]; }
    
    /// Constructor: all set to zero.
    ChannelCalibration_t();
    
    /// Constructor: sets all values.
    ChannelCalibration_t(
      AdderChannelID channel,
      double slope, millivolt offset,
      double slopeErr = 0.0, millivolt offsetErr = millivolt{ 0.0 },
      millivolt covariance = millivolt{ 0.0 },
      microseconds timeOffset = microseconds{ 0.0 },
      microseconds timeOffsetErr = microseconds{ 0.0 }
      );
    
    /// Returns the calibration factor for a waveform of given `amplitude`.
    double factor(millivolt amplitude) const;
    
    /// Returns the relative calibration factor uncertainty for a waveform of given `amplitude`.
    double factorUncertainty(millivolt amplitude) const;
    
    /// Returns the calibration factor uncertainty for a waveform of given `amplitude`.
    double factorRelUncertainty(millivolt amplitude) const;
    
    /// Returns the calibrated value of `amplitude`.
    millivolt computeAmplitude(millivolt amplitude) const;
    
    /// Returns the variance on the calibrated value of `amplitude` [mV&sup2;].
    double computeVariance(millivolt amplitude) const;
    
    /// Returns the uncertainty on the calibrated value of `amplitude`.
    millivolt computeUncertainty(millivolt amplitude) const;
    
    /// Returns a `waveform` calibrated for the specified `amplitude`.
    template <typename ValArray>
    auto applyToWaveform(millivolt amplitude, ValArray const& waveform) const
      { return waveform * factor(amplitude); }
    
  }; // ChannelCalibration_t
  
  /// Type to store the calibration constants from all adder channels.
  using CalibrationConstants_t = std::vector<ChannelCalibration_t>;
  
  /// @}
  // ---- END ----  Algorithm configuration records  ---------------------------
  
  // --- BEGIN ---  Algorithm configuration (FHiCL)  ---------------------------
  /// @name  Algorithm configuration (FHiCL)
  /// @{
  
  /// Adder channel calibration FHiCL configuration.
  struct ChannelCalibrationConfig {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<raw::Channel_t> channel {
      Name{ "channel" },
      Comment{ "adder channel number" }
      };
    
    fhicl::Atom<double> slope {
      Name{ "slope" },
      Comment{ "linear calibration slope coefficient" }
      };
    
    fhicl::Atom<millivolt> offset {
      Name{ "offset" },
      Comment{ "linear calibration offset" }
      };
    
    fhicl::Atom<double> slopeErr {
      Name{ "slopeErr" },
      Comment{ "uncertainty on the linear calibration slope coefficient" },
      0.0
      };
    
    fhicl::Atom<millivolt> offsetErr {
      Name{ "offsetErr" },
      Comment{ "uncertainty on the linear calibration offset" },
      millivolt{ 0.0 }
      };
    
    fhicl::Atom<millivolt> slopeOffsetCov {
      Name{ "slopeOffsetCov" },
      Comment{ "covariance of the linear calibration parameters" },
      millivolt{ 0.0 }
      };
    
    fhicl::Atom<microseconds> timeOffset {
      Name{ "timeOffset" },
      Comment{ "delay to be added to waveform time" },
      microseconds{ 0.0 }
      };
    
    fhicl::Atom<microseconds> timeOffsetErr {
      Name{ "timeOffsetErr" },
      Comment{ "uncertainty on the waveform time correction" },
      microseconds{ 0.0 }
      };
    
  }; // ChannelCalibrationConfig
  
  /// Adder calibration database FHiCL configuration.
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<fhicl::TableAs<ChannelCalibration_t, ChannelCalibrationConfig>>
    CalibrationCoefficients {
      Name{ "CalibrationCoefficients" },
      Comment{ "all calibration coefficients, one line per channel" }
      };
    
  }; // Config
  
  using Parameters = fhicl::Table<Config>;
  
  
  /// @}
  // ---- END ----  Algorithm configuration (FHiCL)  ---------------------------
  
  
  /**
   * @brief Implementation of the calibration of waveforms in a specific run.
   * @see RunCalibration
   * 
   * See `AmplitudeAdderCalibration` for details on the specific implementation
   * of this calibration, and `AdderCalibrationDatabase::RunCalibration` for
   * the documentation of the interface.
   * 
   */
  class RunCalibration: public AdderCalibrationDatabase::RunCalibration {
    
    using Base_t = AdderCalibrationDatabase::RunCalibration;
    
    CalibrationConstants_t fCalibration; ///< All calibration constants.
    
    /**
     * @brief Implementation of the `calibrationFactor()` method.
     * @param channel the number of adder channel the waveform is from
     * @param waveform the sequence of samples of the waveform
     * @return a calibration factor
     * 
     * Computes the amplitude of the waveform and returns the corresponding
     * calibration factor from the configured linear curve.
     */
    virtual double doCalibrationFactor
      (AdderChannelID channel, WaveformSamples_t const& waveform) const override;
    
    /**
     * @brief Implementation of the `timeOffset()` method.
     * @param channel the number of adder channel the waveform is from
     * @param waveform the sequence of samples of the waveform (ignored)
     * @return a calibration factor
     * 
     * Returns the time offset configured for the specified channel.
     * The content of the waveform is ignored.
     */
    virtual microseconds doTimeOffset
      (AdderChannelID channel, WaveformSamples_t const& waveform) const override;
    
    /// Returns the calibration record of the specified channel.
    /// @throws UnknownChannelError if channel is not found
    ChannelCalibration_t const& fetchChannel(AdderChannelID channel) const;
  
      public:
    
    /**
     * @brief Constructor: stores the calibration curve.
     * @param run the run this calibration is for
     * @param calibration all the calibration constants
     * @param logCategory name of the logging stream to emit output to
     */
    RunCalibration(
      int run, CalibrationConstants_t calibration,
      std::string logCategory = "AdderCalibrationDatabase"
      );
    
    /// Returns the calibration constants of `channel`, `nullptr` if not known.
    ChannelCalibration_t const* findChannel(AdderChannelID channel) const;
    
  }; // RunCalibration
  
  /**
   * @brief Constructor: reads from FHiCL all the `calibration` constants.
   * @param config FHiCL configuration
   * 
   * This constructor can be triggered by passing a `fhicl::ParameterSet`.
   */
  AmplitudeAdderCalibration(Parameters const& params)
    : AmplitudeAdderCalibration{ params() }
    {}
  
  /**
   * @brief Constructor: reads from FHiCL all the `calibration` constants.
   * @param config FHiCL configuration
   */
  AmplitudeAdderCalibration(Config const& config);
  
  /**
   * @brief Constructor: copies all the `calibration` constants.
   * @param calibration all calibration coefficients, one item per channel
   */
  AmplitudeAdderCalibration(CalibrationConstants_t calibration);
  
  
  /**
   * @brief Extracts the amplitude of the specified waveform.
   * @tparam Samples type of list of sample values in the waveform
   * @tparam Baseline type of baseline value
   * @param samples list of sample values in the waveform
   * @param baseline (default: 0) baseline value
   * @return the amplitude of the waveform peak
   * 
   * This simple algorithm returns the maximum deviation from the baseline.
   * The returned value is always positive.
   */
  template <typename Samples, typename Baseline = typename Samples::value_type>
  static auto computePeakAmplitude
    (Samples const& samples, Baseline baseline = Baseline{});
  
  
    protected:
  
  CalibrationConstants_t fCalibration; ///< All calibration constants.
  
  /// Returns the calibration factor for the specified waveform and channel.
  virtual std::unique_ptr<AdderCalibrationDatabase::RunCalibration>
  doCalibrationForRun
    (unsigned int run) const override;
  
  /// Implements configuration dumping. Prints all the configuration.
  virtual void doDumpConfig
    (std::ostream& out, details::Indenter& nextline) const override;
  
}; // icarus::trigger::AmplitudeAdderCalibration


namespace icarus::trigger {
  
  /// Conversion of channel calibration parameter for FHiCL.
  AmplitudeAdderCalibration::ChannelCalibration_t convert
    (AmplitudeAdderCalibration::ChannelCalibrationConfig const& config);
  
} // namespace icarus::trigger


// --- BEGIN --- Sorting for AmplitudeAdderCalibration::ChannelCalibration_t ---
namespace icarus::trigger::details {
  struct ChannelCalibration_t_sorter {
    using ChannelCalibration_t = AmplitudeAdderCalibration::ChannelCalibration_t;
    constexpr bool operator()
      (ChannelCalibration_t const& a, ChannelCalibration_t const& b) const noexcept
      { return a.channel < b.channel; }
    constexpr bool operator()
      (AdderChannelID a, ChannelCalibration_t const& b) const noexcept
      { return a < b.channel; }
    constexpr bool operator()
      (ChannelCalibration_t const& a, AdderChannelID b) const noexcept
      { return a.channel < b; }
    constexpr bool operator()
      (AdderChannelID a, AdderChannelID b) const noexcept
      { return a < b; }
  }; // ChannelCalibration_t_sorter
  
} // namespace icarus::trigger
// ---- END ---- Sorting for AmplitudeAdderCalibration::ChannelCalibration_t ---



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
// ---  icarus::trigger::AmplitudeAdderCalibration
//------------------------------------------------------------------------------
template
  <typename Samples, typename Baseline /* = typename Samples::value_type */>
auto icarus::trigger::AmplitudeAdderCalibration::computePeakAmplitude
  (Samples const& samples, Baseline baseline /* = Baseline{} */)
{
  using std::begin, std::end;
  auto const b = begin(samples), e = end(samples);
  auto const [ iMin, iMax ] = std::minmax_element(b, e);
  // hoping that C++ converts that `0` into the appropriate type
  return (b == e)? 0: std::max(baseline - *iMin, *iMax - baseline);
} // icarus::trigger::AmplitudeAdderCalibration::computePeakAmplitude()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_ADDERSIGNALCALIBRATION_H
