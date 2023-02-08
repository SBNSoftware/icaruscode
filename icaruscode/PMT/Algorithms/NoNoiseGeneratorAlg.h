/**
 * @file   icaruscode/PMT/Algorithms/NoNoiseGeneratorAlg.h
 * @brief  A (faster?) Gaussian noise generation algorithm.
 * @date   November 17, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/PMT/Algorithms/PMTnoiseGenerator.h
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_NONOISEGENERATORALG_H
#define ICARUSCODE_PMT_ALGORITHMS_NONOISEGENERATORALG_H

// ICARUS libraries
#include "icaruscode/PMT/Algorithms/NoiseGeneratorAlg.h"
#include "icaruscode/Utilities/quantities_utils.h" // util::...::value()

// C/C++ standard libraries
#include <algorithm>
#include <vector>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename ADCT = double> class NoNoiseGeneratorAlg;
}
/**
 * @brief No-noise generator algorithm for PMT electronics simulation.
 * @param ADCT type of the ADC count this interface treats with.
 * 
 * This noise generator always emits no noise (`0`).
 * 
 * The CLHEP random engine passed in the configuration is ignored.
 */
template <typename ADCT /* = double */>
class icarus::opdet::NoNoiseGeneratorAlg
  : public icarus::opdet::NoiseGeneratorAlg<ADCT>
{
  using Base_t = icarus::opdet::NoiseGeneratorAlg<ADCT>; ///< Base class alias.
  
    public:
  
  using ADCcount_t = typename Base_t::ADCcount_t;
  using Timestamp_t = typename Base_t::Timestamp_t;
  
  /// Underlying fundamental type of `ADCcount_t`, e.g. `float`.
  using ADCvalue_t = util::value_t<typename Base_t::ADCcount_t>;
  
  
  // --- BEGIN -- Configuration ------------------------------------------------
  /// Algorithm configuration parameters.
  struct Params_t {
    // empty
  };
  
  
  /// Configuration parameters.
  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    // empty
    
  }; // struct Config
  
  
  // --- END ---- Configuration ------------------------------------------------
  
  
  // ---------------------------------------------------------------------------
  /// Default constructor.
  NoNoiseGeneratorAlg() = default;
  
  /**
   * @brief Constructor: acquires configuration parameters as data structure.
   * @param params configuration parameters
   * @param engine _ignored_
   */
  NoNoiseGeneratorAlg(Params_t params, CLHEP::HepRandomEngine& engine);
  
  
  /**
   * @brief Constructor: acquires configuration parameters as data structure.
   * @param params configuration parameters in FHiCL format
   * @param engine _ignored_
   */
  NoNoiseGeneratorAlg(Config const& config, CLHEP::HepRandomEngine& engine);
  
  
    private:
  
  Params_t const fParams; ///< Configuration parameters.
  
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /**
   * @brief Adds noise to `n` samples starting at `begin`.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param begin pointer to the first sample to be added with noise
   * @param n number of samples to add noise to
   * @return the number of samples actually added noise
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
    ) override;
  
  /**
   * @brief Overwrites `n` samples starting at `begin` with noise.
   * @param channel ID of the readout channel the samples are from
   * @param time the absolute time of the first sample being filled [UTC, ns]
   * @param begin pointer to the first sample to be overwritten with noise
   * @param n number of samples to fill with noise
   * @return the number of samples actually overwritten with noise
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
    ) override;
  
  /**
   * @brief Prints into the stream the parameters of this algorithm.
   * @param out the C++ output stream to write into
   * @param indent indentation string, prepended to all lines except first
   * @param firstIndent indentation string prepended to the first line
   */
  virtual void doDump(
    std::ostream& out,
    std::string const& indent, std::string const& firstIndent
    ) const override;
  
  // --- END ---- Virtual interface --------------------------------------------
  
  
  /// Extracts the configuration parameters from a FHiCL configuration.
  static Params_t convert(Config const& config);
  
  
}; // icarus::opdet::NoNoiseGeneratorAlg



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename ADCT>
icarus::opdet::NoNoiseGeneratorAlg<ADCT>::NoNoiseGeneratorAlg
  (Params_t params, CLHEP::HepRandomEngine&)
  : fParams{ std::move(params) }
{}


// -----------------------------------------------------------------------------
template <typename ADCT>
icarus::opdet::NoNoiseGeneratorAlg<ADCT>::NoNoiseGeneratorAlg
  (Config const& config, CLHEP::HepRandomEngine& engine)
  : NoNoiseGeneratorAlg{ convert(config), engine }
{}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::NoNoiseGeneratorAlg<ADCT>::doAdd
  (raw::Channel_t, Timestamp_t, ADCcount_t*, std::size_t n)
{
  return n;
} // icarus::opdet::NoNoiseGeneratorAlg<ADCT>::doAdd()


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::NoNoiseGeneratorAlg<ADCT>::doFill
  (raw::Channel_t, Timestamp_t, ADCcount_t* begin, std::size_t n)
{
  std::fill_n(begin, n, ADCcount_t{ 0 });
  return n;
} // icarus::opdet::NoNoiseGeneratorAlg<ADCT>::doFill()


// -----------------------------------------------------------------------------
template <typename ADCT>
void icarus::opdet::NoNoiseGeneratorAlg<ADCT>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  out << firstIndent << "No-noise generator algorithm (NoNoiseGeneratorAlg).";
} // icarus::opdet::NoNoiseGeneratorAlg<ADCT>::doDump()


// -----------------------------------------------------------------------------
template <typename ADCT>
auto icarus::opdet::NoNoiseGeneratorAlg<ADCT>::convert
  (Config const& config) -> Params_t
{
  return {};
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_NONOISEGENERATORALG_H


