/**
 * @file   icaruscode/PMT/Data/WaveformBaseline.h
 * @brief  A baseline for a waveform.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 11, 2020
 * @see    icaruscode/PMT/Data/WaveformBaseline.cxx
 */

#ifndef ICARUSCODE_PMT_DATA_WAVEFORMBASELINE_H
#define ICARUSCODE_PMT_DATA_WAVEFORMBASELINE_H


// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <cmath> // std::round()


//------------------------------------------------------------------------------
namespace icarus {
  
  struct WaveformBaseline;
  
  /// Prints the value of the baseline into a stream.
  std::ostream& operator<<
    (std::ostream& out, icarus::WaveformBaseline const& baseline);
  
} // namespace icarus

/**
 * @brief Class containing a waveform baseline value.
 * 
 * The baseline is stored as a floating point value, not to lose precision.
 * 
 * This class is a data product wrapper for a simple value, with some usability
 * candies attached.
 * 
 * 
 * Example of simple usage
 * ------------------------
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * 
 * icarus::WaveformBaseline const baseline { 1.2f };
 * 
 * std::cout << "Baseline: " << baseline << " ADC" << std::endl;
 * 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will print `Baseline: 1.2 ADC`.
 * 
 * 
 * Example of usage with quantity values
 * --------------------------------------
 * 
 * The following is a more complex example showing baseline subtraction.
 * We assume we are using the type `util::quantities::counts` as ADC count type
 * (defined in `lardataalg/Utilities/quantities/electronics.h`) for the
 * waveform.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * using ADCCount_t = util::quantities::counts;
 * std::vector<ADCCount_t> subtractBaseline
 *   (std::vector<ADCCount_t> const& data, icarus::WaveformBaseline const& baseline)
 * {
 *   
 *   icarus::waveform_operations::NegativePolarityOperations<float> const
 *     waveOps { baseline() };
 *   
 *   std::vector<ADCCount_t> subtracted;
 *   subtracted.reserve(data.size());
 *   for (auto sample: data) {
 *     subtracted.emplace_back
 *       (std::round(waveOps.subtractBaseline(sample.value())));
 *   }
 *   
 *   return subtracted;
 * } // subtractBaseline()
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * The `subtracted` waveform is positive and baseline-subtracted.
 * We use ICARUS utilities to manage the polarity of the waveform
 * (hard-coded negative). Note that the subtraction is less than trivial because
 * of the integral type waveform on top of the floating point baseline.
 * The samples are converted into floating point for the subtraction, then
 * reconverted (rounded) back.
 */
struct icarus::WaveformBaseline {
  
  using Baseline_t = float; ///< Type of baseline value.
  
  Baseline_t fBaseline {}; ///< The current value of the baseline.
  
  
  // --- BEGIN -- Constructors -------------------------------------------------
  
  /// Constructor: default baseline (`0`).
  WaveformBaseline() = default;
  
  /// Constructor: sets the baseline.
  WaveformBaseline(Baseline_t baseline): fBaseline(baseline) {}
  
  // --- END -- Constructors ---------------------------------------------------
  
  
  // --- BEGIN -- Access to the baseline ---------------------------------------
  /**
   * @name Access to the baseline
   * 
   * In addition to the direct method (`baseline()`) a few candies are offered:
   * a function-like operator for baseline access.
   * 
   */
  /// @{
  /// Returns the current baseline value.
  Baseline_t baseline() const { return fBaseline; }
  
  /// Returns the current baseline value.
  Baseline_t operator() () const { return baseline(); }
  
  // --- END -- Access to the baseline -----------------------------------------
  
}; // icarus::trigger::WaveformBaseline


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_DATA_WAVEFORMBASELINE_H
