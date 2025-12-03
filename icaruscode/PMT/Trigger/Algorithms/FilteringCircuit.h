/**
 * @file   icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.h
 * @brief  Filtering circuit response interface.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.cxx
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FILTERINGCIRCUIT_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FILTERINGCIRCUIT_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/details/Indenter.h"

// C/C++ standard libraries
#include <cstdint> // std::size_t
#include <complex>
#include <iosfwd> // std::ostream
#include <utility> // std::move()
#include <valarray>


// -----------------------------------------------------------------------------
namespace icarus::trigger { class FilteringCircuit; }
/// Filtering circuit simulation interface.
class icarus::trigger::FilteringCircuit {
  
    public:
  
  /// Base data type used `Frequency_t` and `Voltage_t`.
  using Data_t = float;
  
  /// Type representing voltage, used internally for computation.
  using Voltage_t = Data_t;
  
  /// Representation of amplitudes in frequency domain for internal computation.
  using WaveformAmplitudes_t = std::valarray<std::complex<Voltage_t>>;

  /// Type for representing a frequency.
  using Frequency_t = std::complex<Data_t>;
  
  /// Representation of amplitudes in frequency domain for internal computation.
  using Frequencies_t = std::valarray<Frequency_t>;
  
  
  virtual ~FilteringCircuit() = default;
  
  /**
   * @brief Applies the transformation to the specified amplitudes.
   * @param period base period of the amplitudes [s]
   * @param amplitudes the vector of amplitudes to transform
   * @return the new amplitudes
   * @see `makeFrequencies()`
   * 
   * The amplitudes are assumed to be each relative to the frequency in the
   * `makeFrequencies(period)` list.
   */
  WaveformAmplitudes_t apply
    (Data_t period, WaveformAmplitudes_t&& amplitudes) const;
  
  /**
   * @brief Applies the transformation to the specified amplitudes.
   * @param s frequency corresponding to each amplitude [Hz]
   * @param amplitudes the vector of amplitudes to transform
   * @return the new amplitudes
   */
  WaveformAmplitudes_t apply
    (Frequencies_t const& s, WaveformAmplitudes_t&& amplitudes) const
    { return doApply(s, std::move(amplitudes)); }
  
  /**
   * @brief Applies the transformation in place to the specified amplitudes.
   * @param s frequency corresponding to each amplitude [Hz]
   * @param[in,out] amplitudes the vector of amplitudes to transform in place
   */
  void apply(Frequencies_t const& s, WaveformAmplitudes_t& amplitudes) const
    { amplitudes = doApply(s, std::move(amplitudes)); }
  
  
  // --- BEGIN --- Configuration dump ----------------------------------------
  /// @name Configuration dumping
  /// @{
  /**
   * @brief Dumps the configuration to a stream.
   * @tparam Stream type of stream to insert to
   * @param out the stream to dump the information into
   * @param indent string used to indent any new line
   * @param firstIndent string used to indent the starting line
   * 
   * The dump is in human-readable form and is intended for logging.
   * The output ends with a new line (and no indentation inserted).
   * 
   * The output can be multi-line, and it is indented using the `indent`
   * string, with the exception of the first line which uses `firstIndent`.
   */
  template <typename Stream>
  void dumpConfig
    (Stream& out, std::string indent, std::string firstIndent) const;
  
  /**
   * @brief Dumps the configuration to a stream.
   * @tparam Stream type of stream to insert to
   * @param out the stream to dump the information into
   * @param indent string used to indent all lines (including the first)
   * @see dumpConfig(Stream&, std::string, std::string)
   */
  template <typename Stream>
  void dumpConfig(Stream& out, std::string const& indent = "") const
    { dumpConfig(out, indent, indent); }
  
  /**
   * @brief Dumps the configuration to a stream.
   * @param indent string used to indent any new line
   * @param firstIndent string used to indent the starting line
   * @return an opaque object managing the dump
   * @see `dumpConfig()`
   * 
   * The output is from `dumpConfig()`. This function allows direct insertion:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << filter.dumpConfig("   ");
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Differently from `dumpConfig()`, this ends does not end insertion with
   * a new line.
   */
  auto printConfig(std::string indent, std::string firstIndent) const;

  /**
   * @brief Dumps the configuration to a stream.
   * @param indent string used to indent all lines (including the first)
   * @see printConfig(std::string, std::string)
   */
  auto printConfig(std::string const& indent = "") const;

  /// @}
  // ---- END ---- Configuration dump ----------------------------------------
  
  /**
   * @brief Returns a vector of frequencies with base sample period.
   * @param period sampling period (the highest frequency is twice this)
   * @param n number of samples
   * @return a vector of `n` / 2 + 1 frequencies (_s_)
   * 
   * Given a sequence of `n` samples with sampling time `period` (e.g. 2 ns),
   * the relevant frequencies _f_ in a spectrum analysis are from `0`
   * (continuum) to 1/(2*`period`), with steps of 1/(`n`*`period`)
   * (inverse of the duration of the full waveform).
   * This function returns a vector of values _s_ (defined as 2 &pi; i _f_),
   * from `0` to (&pi; i)/`period`.
   * 
   * For example, with 5000 samples of 2 ns each, and a waveform duration of
   * 10 microseconds, there are 2501 frequencies, from 0 to 250 MHz in steps of
   * 100 kHz.
   * 
   * The frequencies are bare values returned in units of inverse of `period`
   * (e.g. gigahertz if `period` is in nanoseconds). It is recommended to
   * use `period` in seconds (then frequencies will be returned in hertz).
   */
  static Frequencies_t makeFrequencies(Data_t period, std::size_t n);
  
    protected:
  
  struct Dumper;
  
  friend std::ostream& operator<< (std::ostream&, Dumper);
  
  /**
   * @brief Applies the transformation to the specified amplitudes.
   * @param s (complex) frequencies in the spectrum (e.g. 2&pi;i _f_)
   * @param amplitudes (complex) amplitude for each frequency
   * @return the filtered amplitudes
   * 
   * A way to obtain the frequencies `s` is via `makeFrequencies()`.
   */
  virtual WaveformAmplitudes_t doApply
    (Frequencies_t const& s, WaveformAmplitudes_t amplitudes) const
    = 0;
  
  /**
   * @brief Implementation of `dumpConfig`.
   * @param out stream to insert output to
   * @param nextLine an `Indenter` object to help with new lines
   * 
   * Implementation guidelines:
   *  * Use `out << nextLine` at the beginning of each output line
   *    (first included)
   *  * When passing the indenter around, pass it by reference.
   *  * Do not end the last line (it will be done by `dumpConfig()`).
   * 
   * See `SallenKeyFilter::doDumpConfig()` for an implementation example.
   */
  virtual void doDumpConfig
    (std::ostream& out, details::Indenter& nextLine) const;
  
}; // icarus::trigger::FilteringCircuit


namespace icarus::trigger {
  
  /// Dumps a filter configuration into the `out` stream.
  std::ostream& operator<< (std::ostream& out, FilteringCircuit::Dumper dumper);
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
// ---  Inline implementation
// -----------------------------------------------------------------------------
inline auto icarus::trigger::FilteringCircuit::apply
  (Data_t period, WaveformAmplitudes_t&& amplitudes) const -> WaveformAmplitudes_t
{
  return apply
    (makeFrequencies(period, 2 * (amplitudes.size() - 1)), std::move(amplitudes));
}


// -----------------------------------------------------------------------------
class icarus::trigger::FilteringCircuit::Dumper {
  FilteringCircuit const* fFilter = nullptr;
  details::Indenter fIndenter;
    public:
  
  Dumper(FilteringCircuit const& filter, details::Indenter indenter)
    : fFilter{ &filter }, fIndenter{ std::move(indenter) } {}

  void into(std::ostream& out) { fFilter->doDumpConfig(out, fIndenter); }
  
};


// -----------------------------------------------------------------------------
inline std::ostream& icarus::trigger::operator<<
  (std::ostream& out, FilteringCircuit::Dumper dumper)
{
  dumper.into(out);
  return out;
}


// -----------------------------------------------------------------------------
inline auto icarus::trigger::FilteringCircuit::printConfig
  (std::string indent, std::string firstIndent) const
{
  return Dumper
    (*this, details::Indenter{ std::move(indent), std::move(firstIndent) });
}

inline auto icarus::trigger::FilteringCircuit::printConfig
  (std::string const& indent /* = "" */) const
  { return printConfig(indent, indent); }


// -----------------------------------------------------------------------------
template <typename Stream>
void icarus::trigger::FilteringCircuit::dumpConfig(
  Stream& out, std::string indent, std::string firstIndent
) const {
  out << printConfig(std::move(indent), std::move(firstIndent)) << std::endl;
} // icarus::trigger::AdderSignalSimulation::FilteringCircuit::dumpConfig()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_FILTERINGCIRCUIT_H
