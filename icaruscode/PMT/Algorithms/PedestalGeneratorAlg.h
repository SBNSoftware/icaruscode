/**
 * @file   icaruscode/PMT/Algorithms/PedestalGeneratorAlg.h
 * @brief  Interface for a pedestal generating algorithm.
 * @date   November 19, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_PEDESTALGENERATORALG_H
#define ICARUSCODE_PMT_ALGORITHMS_PEDESTALGENERATORALG_H

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t

// C/C++ standard libraries
#include <ostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <limits> // std::numeric_limits
#include <cstdint> // std::uint64_t
#include <cstddef> // std::size_t
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename ADCT> class PedestalGeneratorAlg; 
  template <typename ADCT>
  std::ostream& operator<<
    (std::ostream& out, PedestalGeneratorAlg<ADCT> const& gen);
} // namespace icarus::opdet

// -----------------------------------------------------------------------------
/**
 * @brief Interface for a PMT electronics pedestal generator algorithm.
 * @param ADCT type of the ADC count this interface deals with.
 * 
 * A pedestal generator is an algorithm able to fill digitization samples with
 * a pedestal. In the ideal case, the pedestal is constant. In any realistic
 * case, the pedestal is at least affected by some electronics noise.
 * The interface allows to specify the channel that the pedestal will pertain to
 * and the absolute time at which the pedestal is produced. Neither is an
 * essential parameter.
 * 
 * Because it's such a common requirement, the interface also explicitly allows
 * a noise generator algorithm to be specified; typically the noise algorithm
 * object will become owned by the pedestal generator.
 * 
 * The configuration of the generator is expected to be performed on
 * construction, and it is therefore out of the interface.
 * Configuration and all additional setup elements must be passed through the
 * constructor; that includes the noise generator algorithm.
 * 
 * Currently the interface expects to know ahead how many noise samples are
 * going to be extracted by each call.
 * 
 * 
 * Implementer notes
 * ==================
 * 
 * The interface methods are deliberately marked non-`const`, with the
 * expectation that random engines will be used and that the implementers will
 * not be forced to use mutable state: as a baseline, the implementations of
 * this interface should be treated as thread-unsafe.
 * 
  * The interface currently demands the definition of two functions, one to fill
 * (overwrite) a buffer, the other to add to it. The rationale is that the
 * implementer can choose the ones which is more efficient to implement (also
 * depending on the expected most common use), and it is expected that most
 * implementations will have either of the two call the other.
 * 
 * The other overrideable function is a configuration dump that is convenient
 * for the log files but does not affect the noise generation.
 * 
 * Although the constructor is in principle not bound to any particular
 * signature from the point of view of C++ inheritance, because these objects
 * are going to be created via a factory written in C++, it is still demanded
 * (by the factory interface or by the callers of the factory) that the
 * constructor signature be standardized. To make these objects work with _art_
 * tool infrastructure, they need to be wrapped into a class which exposes a
 * constructor whose sole argument can be initialized from a
 * `fhicl::ParameterSet`.
 * 
 */
template <typename ADCT>
class icarus::opdet::PedestalGeneratorAlg {
  
    public:
  
  using ADCcount_t = ADCT; ///< Type of the sample value.
  
  using Timestamp_t = std::uint64_t; ///< Type of timestamp.
  
  
  /// Special pedestal value to express no knowledge of the pedestal.
  static constexpr ADCcount_t NoPedestalLevel
    = std::numeric_limits<ADCcount_t>::max();
  
  
  virtual ~PedestalGeneratorAlg() = default;
  
  
  /**
   * @brief Returns the pedestal level at a certain `channel` and `time`.
   * @param channel channel the pedestal belongs to
   * @param time absolute time at which to query the pedestal level [UTC, ns]
   * @return the average level at the specific channel and time
   * 
   * If the tool does not support this concept, this function will return the
   * value `NoPedestalLevel`.
   */
  ADCcount_t pedestalLevel(raw::Channel_t channel, Timestamp_t time) const;

  
  // --- BEGIN -- Pedestal addition --------------------------------------------
  /// @name Pedestal addition
  /// @{
  
  /**
   * @brief Adds pedestal to `n` samples starting at `begin`.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param begin pointer to the first sample to be added with pedestal
   * @param n number of samples to add pedestal to
   * @return the number of samples actually added pedestal
   * 
   * No check is performed on the validity of the destination buffer.
   * It is guaranteed that no more than `n` samples are changed, in the
   * contiguous sequence after `begin`.
   * 
   * The return value should be `n` unless the generator was unable to fully
   * fulfill the request.
   */
  std::size_t add(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, std::size_t n
    );
  
  /**
   * @brief Adds pedestal to the samples from `begin` up to `end` (excluded).
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param begin pointer to the first sample to be added pedestal
   * @param end pointer past the last sample to be added pedestal
   * @return pointer to the first sample not added pedestal
   * @see `fill(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`
   * 
   * The same considerations apply as expressed in
   * `add(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`,
   * with the difference that the expected return value is `end` unless the
   * generator was unable to fully fulfill the request.
   */
  ADCcount_t* add(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, ADCcount_t* end
    );
  
  /**
   * @brief Adds pedestal to all the `samples`.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param samples container of the samples to be added pedestal to
   * @return number of samples actually added pedestal
   * @see `fill(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`
   * 
   * The same general considerations apply as expressed in
   * `add(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`.
   * 
   * In addition, the buffer (`samples`) is now automatically guaranteed to be
   * valid and _all_ its samples are added pedestal. The size of the buffer is
   * never resized.
   * The expected return value is `samples.size()` unless the generator was
   * unable to fully fulfill the request.
   */
  std::size_t add(
    raw::Channel_t channel, Timestamp_t time,
    std::vector<ADCcount_t>& samples
    );
  
  /// @}
  // --- END ---- Pedestal addition --------------------------------------------
  
  
  // --- BEGIN -- Pedestal overwrite -------------------------------------------
  /// @name Pedestal overwrite
  /// @{
  
  /**
   * @brief Overwrites `n` samples starting at `begin` with pedestal.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param begin pointer to the first sample to be overwritten with pedestal
   * @param n number of samples to fill with pedestal
   * @return the number of samples actually overwritten with pedestal
   * 
   * No check is performed on the validity of the destination buffer.
   * Old values are typically completely overwritten. It is guaranteed that no
   * more than `n` samples are overwritten, in the contiguous sequence after
   * `begin`.
   * 
   * The return value should be `n` unless the generator was unable to fully
   * fulfill the request.
   */
  std::size_t fill(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, std::size_t n
    );
  
  /**
   * @brief Overwrites samples from `begin` to `end` with pedestal.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param begin pointer to the first sample to be overwritten with pedestal
   * @param end pointer past the last sample to be overwritten with pedestal
   * @return pointer to the first sample not overwritten with pedestal
   * @see `fill(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`
   * 
   * The same considerations apply as expressed in
   * `fill(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`,
   * with the difference that the expected return value is `end` unless the
   * generator was unable to fully fulfill the request.
   */
  ADCcount_t* fill(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, ADCcount_t* end
    );
  
  /**
   * @brief Overwrites all the `samples` with pedestal.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param samples container of the samples to be overwritten
   * @return number of samples actually overwritten with pedestal
   * @see `fill(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`
   * 
   * The same general considerations apply as expressed in
   * `fill(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`.
   * 
   * In addition, the buffer (`samples`) is now automatically guaranteed to be
   * valid and _all_ its samples are overwritten with pedestal. The size of the
   * buffer is never resized.
   * The expected return value is `samples.size()` unless the generator was
   * unable to fully fulfill the request.
   */
  std::size_t fill(
    raw::Channel_t channel, Timestamp_t time,
    std::vector<ADCcount_t>& samples
    );
  
  
  /**
   * @brief Adds to the end of `samples` additional `n` pedestal samples.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param samples container of the samples to be expanded
   * @param n number of pedestal samples to be added
   * @return number of samples actually overwritten with pedestal
   * @see `fill(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`
   * 
   * The same general considerations apply as expressed in
   * `fill(raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t)`.
   * 
   * In addition, the buffer (`samples`) is automatically expanded with `n`
   * more elements, which are assigned pedestal values.
   * The expected return value is `n` unless the generator was unable to fully
   * fulfill the request. In that case, the sample is resized to include only
   * up to the samples actually added with pedestal (the actually allocated
   * memory, `capacity()`, may exceed that though, and the buffer does not
   * `shrink_to_fit()`).
   */
  std::size_t append(
    raw::Channel_t channel, Timestamp_t time,
    std::vector<ADCcount_t>& samples, std::size_t n
    );
  
  /// @}
  // --- END ---- Pedestal overwrite -------------------------------------------
  
  
  // --- BEGIN -- Dump configuration on screen ---------------------------------
  /// @name Dump configuration on screen
  /// @{
  
  //@{
  /**
   * @brief Prints on stream the configuration parameters.
   * @param out the stream to write into
   * @param indent indentation string, prepended to all lines except first
   * @param firstIndent indentation string prepended to the first line
   */
  void dump(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const;
  void dump(std::ostream&& out, std::string const& indent = "") const;
  //@}


  //@{
  /**
   * @brief Returns the configuration parameter description.
   * @param indent indentation string, prepended to all lines except first
   * @param firstIndent indentation string prepended to the first line
   * @return a string with the parameters of this shape
   */
  std::string toString
    (std::string const& indent, std::string const& firstIndent) const;
  std::string toString(std::string const& indent = "") const;
  //@}
  
  /// @}
  // --- END ---- Dump configuration on screen ---------------------------------
  
  
    protected:
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /**
   * @brief Returns the pedestal level at a certain `channel` and `time`.
   * @param channel channel the pedestal belongs to
   * @param time absolute time at which to query the pedestal level [UTC, ns]
   * @return the average level at the specific channel and time
   * 
   * If the tool does not support this concept, this function will return the
   * value `NoPedestalLevel`.
   */
  virtual ADCcount_t doPedestalLevel
    (raw::Channel_t channel, Timestamp_t time) const = 0;
  
  
  /**
   * @brief Adds pedestal to `n` samples starting at `begin`.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param begin pointer to the first sample to be added with pedestal
   * @param n number of samples to add pedestal to
   * @return the number of samples actually added pedestal
   * 
   * No check is performed on the validity of the destination buffer.
   * It is guaranteed that no more than `n` samples are changed, in the
   * contiguous sequence after `begin`.
   * 
   * The return value should be `n` unless the generator was unable to fully
   * fulfill the request.
   */
  virtual std::size_t doAdd(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, std::size_t n
    ) = 0;
  
  /**
   * @brief Overwrites `n` samples starting at `begin` with pedestal.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param begin pointer to the first sample to be overwritten with pedestal
   * @param n number of samples to fill with pedestal
   * @return the number of samples actually overwritten with pedestal
   * 
   * No check is performed on the validity of the destination buffer.
   * Old values are typically completely overwritten. It is guaranteed that no
   * more than `n` samples are overwritten, in the contiguous sequence after
   * `begin`.
   * 
   * The return value should be `n` unless the generator was unable to fully
   * fulfill the request.
   */
  virtual std::size_t doFill(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, std::size_t n
    ) = 0;
  
  /**
   * @brief Prints into the stream the parameters of this algorithm.
   * @param out the C++ output stream to write into
   * @param indent indentation string, prepended to all lines except first
   * @param firstIndent indentation string prepended to the first line
   * 
   * Sends a complete description of the configuration of the algorithm to the
   * `out` stream. The output starts on the current line of the stream, which
   * is added a `firstIndent` indentation string. All following new lines are
   * started with an indentation string `indent`. The last line is not ended
   * with a newline character.
   * 
   * The default implementation prints nothing (also no indentation).
   * Custom implementations may choose to add line breaks in addition to the
   * ones described above.
   */
  virtual void doDump(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const {}
  
  // --- END ---- Virtual interface --------------------------------------------
  
  
  // --- BEGIN -- Convenience implementation functions -------------------------
  
  /**
   * @brief Implementation helper: implements `doFill()` based on `doAdd()`.
   * 
   * This implementation is trivial: the buffer is first zeroed, and then
   * pedestal is added.
   * 
   * The definition of `doFill()` becomes:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::size_t doFill(
   *   raw::Channel_t channel, Timestamp_t time,
   *   ADCcount_t* begin, std::size_t n
   *   ) override
   *   { return _fillByAdding(channel, time, begin, n); }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  std::size_t _fillByAdding(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, std::size_t n
    );
  
  /**
   * @brief Implementation helper: implements `doAdd()` based on `doFill()`.
   * 
   * This implementation is trivial: the pedestal is first generated in a new
   * buffer, and then it's added to the existing one.
   * 
   * The definition of `doAdd()` becomes:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::size_t doAdd(
   *   raw::Channel_t channel, Timestamp_t time,
   *   ADCcount_t* begin, std::size_t n
   *   ) override
   *   { return _addByFilling(channel, time, begin, n); }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  std::size_t _addByFilling(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, std::size_t n
    );
  
  // --- END ---- Convenience implementation functions -------------------------
  
}; // icarus::opdet::PedestalGeneratorAlg


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename ADCT>
auto icarus::opdet::PedestalGeneratorAlg<ADCT>::pedestalLevel
  (raw::Channel_t channel, Timestamp_t time) const -> ADCcount_t
{
  return doPedestalLevel(channel, time);
}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::PedestalGeneratorAlg<ADCT>::add(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  return doAdd(channel, time, begin, n);
}


// -----------------------------------------------------------------------------
template <typename ADCT>
auto icarus::opdet::PedestalGeneratorAlg<ADCT>::add(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, ADCcount_t* end
) -> ADCcount_t* {
  std::size_t const n = add(channel, time, begin, end - begin);
  return begin + n;
}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::PedestalGeneratorAlg<ADCT>::add(
  raw::Channel_t channel, Timestamp_t time,
  std::vector<ADCcount_t>& samples
) {
  return add(channel, time, samples.data(), samples.size());
}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::PedestalGeneratorAlg<ADCT>::fill(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  return doFill(channel, time, begin, n);
}


// -----------------------------------------------------------------------------
template <typename ADCT>
auto icarus::opdet::PedestalGeneratorAlg<ADCT>::fill(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, ADCcount_t* end
) -> ADCcount_t* {
  std::size_t const n = fill(channel, time, begin, end - begin);
  return begin + n;
}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::PedestalGeneratorAlg<ADCT>::fill(
  raw::Channel_t channel, Timestamp_t time,
  std::vector<ADCcount_t>& samples
) {
  return fill(channel, time, samples.data(), samples.size());
}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::PedestalGeneratorAlg<ADCT>::append(
  raw::Channel_t channel, Timestamp_t time,
  std::vector<ADCcount_t>& samples, std::size_t n
) {
  std::size_t const oldSize = samples.size();
  samples.resize(oldSize + n);
  std::size_t const nAdded = fill(channel, time, samples.data() + oldSize, n);
  samples.resize(oldSize + nAdded); // no-op if nAdded is equal to n
  return nAdded;
}


// -----------------------------------------------------------------------------
template <typename ADCT>
void icarus::opdet::PedestalGeneratorAlg<ADCT>::dump(
  std::ostream& out, std::string const& indent, std::string const& firstIndent
) const {
  doDump(out, indent, firstIndent);
}


// -----------------------------------------------------------------------------
template <typename ADCT>
void icarus::opdet::PedestalGeneratorAlg<ADCT>::dump
  (std::ostream&& out, std::string const& indent /* = "" */) const
  { dump(out, indent, indent); }


// -----------------------------------------------------------------------------
template <typename ADCT>
std::string icarus::opdet::PedestalGeneratorAlg<ADCT>::toString
  (std::string const& indent, std::string const& firstIndent) const
{
  std::ostringstream sstr;
  dump(sstr, indent, firstIndent);
  return std::move(sstr).str();
} // icarus::opdet::PedestalGeneratorAlg<>::toString()


// -----------------------------------------------------------------------------
template <typename ADCT>
std::string icarus::opdet::PedestalGeneratorAlg<ADCT>::toString
  (std::string const& indent /* = "" */) const
  { return toString(indent, indent); }


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::PedestalGeneratorAlg<ADCT>::_fillByAdding(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  std::fill_n(begin, n, ADCcount_t{ 0 });
  return doAdd(channel, time, begin, n);
}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::PedestalGeneratorAlg<ADCT>::_addByFilling(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  std::vector<ADCcount_t> noise(n);
  std::size_t const nSamples = fill(channel, time, noise);
  assert(nSamples <= n);
  
  auto it = noise.cbegin();
  auto const nend = it + nSamples;
  while (it != nend) *(begin++) += *it++;
  
  return nSamples;
} // icarus::opdet::PedestalGeneratorAlg::_addByFilling()


// -----------------------------------------------------------------------------
// --- free function implementation
// -----------------------------------------------------------------------------
template <typename ADCT>
std::ostream& icarus::opdet::operator<<
  (std::ostream& out, icarus::opdet::PedestalGeneratorAlg<ADCT> const& gen)
  { out << gen.toString(); return out; }


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_PEDESTALGENERATORALG_H


