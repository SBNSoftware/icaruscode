/**
 * @file   icaruscode/PMT/Algorithms/FastGaussianNoiseGeneratorAlg.h
 * @brief  A (faster?) Gaussian noise generation algorithm.
 * @date   November 17, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/PMT/Algorithms/PMTnoiseGenerator.h
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_FASTGAUSSIANNOISEGENERATORALG_H
#define ICARUSCODE_PMT_ALGORITHMS_FASTGAUSSIANNOISEGENERATORALG_H

// ICARUS libraries
#include "icaruscode/PMT/Algorithms/NoiseGeneratorAlg.h"
#include "icaruscode/Utilities/quantities_utils.h" // util::...::value()
#include "icarusalg/Utilities/FastAndPoorGauss.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C/C++ standard libraries
#include <vector>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename ADCT = double> class FastGaussianNoiseGeneratorAlg;
}
/**
 * @brief Gaussian noise generator algorithm for PMT electronics simulation.
 * @param ADCT type of the ADC count this interface treats with.
 * 
 * This noise generator produces a random Gaussian-distributed noise with
 * the specified RMS and no bias.
 * 
 * It requires a CLHEP random engine as input, and directly uses a custom
 * adapter.
 * 
 * Generating for any type `ADCT` other that the CLHEP-native `double` requires
 * a conversion and slows down the generation.
 * 
 * Compared to GaussianNoiseGeneratorAlg(), we use a somehow faster random
 * generator; to squeeze the CPU cycles, we avoid the CLHEP interface as much as
 * possible; the random number from the engine is immediately converted
 * to single precision, and the rest of the math happens in there as well.
 * No virtual interfaces nor indirection is involved within this function
 * (except for CLHEP random engine and the call to `add()`/`fill()`). We
 * generate a normal variable _z_ (standard deviation 1, mean 0) and we just
 * scale it to the desired standard deviation, not bothering to add the mean
 * offset of 0.
 * 
 * Note that unless the random engine is multi-thread safe, this function
 * won't gain anything from multi-threading.
 * 
 * @note Despite the name, the generator limits to generate one sample per call.
 *       It may be possible, if this proved to be a limitation, to extend
 *       `util::FastAndPoorGauss` to generate an array of numbers, in the hope
 *       that the compiler lets vectorization kick in.
 */
template <typename ADCT /* = double */>
class icarus::opdet::FastGaussianNoiseGeneratorAlg
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
    ADCcount_t RMS; ///< Gaussian sigma parameter for the noise.
  };
  
  
  /// Configuration parameters.
  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<ADCvalue_t> RMS{
      Name{ "RMS" },
      Comment{ "RMS of the Gaussian noise" }
      // mandatory
      };
    
  }; // struct Config
  
  
  // --- END ---- Configuration ------------------------------------------------
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Constructor: acquires configuration parameters as data structure.
   * @param params configuration parameters
   * @param engine random engines
   */
  FastGaussianNoiseGeneratorAlg(Params_t params, CLHEP::HepRandomEngine& engine);
  
  
  /**
   * @brief Constructor: acquires configuration parameters as data structure.
   * @param params configuration parameters in FHiCL format
   * @param engine random engines
   */
  FastGaussianNoiseGeneratorAlg
    (Config const& config, CLHEP::HepRandomEngine& engine);
  
  
    private:
  
  using GaussAdapter_t = util::FastAndPoorGauss<32768U, ADCvalue_t>;
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  Params_t const fParams; ///< All configuration parameters.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  
  /// Random engine used by this algorithm.
  CLHEP::HepRandomEngine& fRandomEngine;
  
  
  /// Random adapter.
  static GaussAdapter_t const FastGauss;
  
  
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
  
  
}; // icarus::opdet::FastGaussianNoiseGeneratorAlg



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename ADCT>
typename icarus::opdet::FastGaussianNoiseGeneratorAlg<ADCT>::GaussAdapter_t
const icarus::opdet::FastGaussianNoiseGeneratorAlg<ADCT>::FastGauss;


// -----------------------------------------------------------------------------
template <typename ADCT>
icarus::opdet::FastGaussianNoiseGeneratorAlg<ADCT>::FastGaussianNoiseGeneratorAlg
  (Params_t params, CLHEP::HepRandomEngine& engine)
  : fParams{ std::move(params) }
  , fRandomEngine{ engine }
{}


// -----------------------------------------------------------------------------
template <typename ADCT>
icarus::opdet::FastGaussianNoiseGeneratorAlg<ADCT>::FastGaussianNoiseGeneratorAlg
  (Config const& config, CLHEP::HepRandomEngine& engine)
  : FastGaussianNoiseGeneratorAlg{ convert(config), engine }
{}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::FastGaussianNoiseGeneratorAlg<ADCT>::doAdd(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  ADCcount_t* const end = begin + n;
  while (begin != end) {
    *(begin++)
      += static_cast<ADCcount_t>(fParams.RMS*FastGauss(fRandomEngine.flat()));
  }
  return n;
} // icarus::opdet::FastGaussianNoiseGeneratorAlg<>::doAdd()


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::FastGaussianNoiseGeneratorAlg<ADCT>::doFill(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  ADCcount_t* const end = begin + n;
  while (begin != end) {
    *(begin++)
      = static_cast<ADCcount_t>(fParams.RMS*FastGauss(fRandomEngine.flat()));
  }
  return n;
} // icarus::opdet::FastGaussianNoiseGeneratorAlg<>::doFill()


// -----------------------------------------------------------------------------
template <typename ADCT>
void icarus::opdet::FastGaussianNoiseGeneratorAlg<ADCT>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  out << firstIndent
    << "Gaussian noise generator algorithm (FastGaussianNoiseGeneratorAlg with "
    << GaussAdapter_t::NPoints << " precomputed values):"
    << "\n" << indent << "  RMS: " << fParams.RMS;
} // icarus::opdet::FastGaussianNoiseGeneratorAlg<>::doDump()


// -----------------------------------------------------------------------------
template <typename ADCT>
auto icarus::opdet::FastGaussianNoiseGeneratorAlg<ADCT>::convert
  (Config const& config) -> Params_t
{
  return {
    ADCcount_t{ config.RMS() } // RMS
    };
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_FASTGAUSSIANNOISEGENERATORALG_H


