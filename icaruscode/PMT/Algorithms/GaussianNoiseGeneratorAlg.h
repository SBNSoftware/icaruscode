/**
 * @file   icaruscode/PMT/Algorithms/GaussianNoiseGeneratorAlg.h
 * @brief  A Gaussian noise generation algorithm.
 * @date   November 17, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/PMT/Algorithms/PMTnoiseGenerator.h
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_GAUSSIANNOISEGENERATORALG_H
#define ICARUSCODE_PMT_ALGORITHMS_GAUSSIANNOISEGENERATORALG_H

// LArSoft/ICARUS libraries
#include "icaruscode/PMT/Algorithms/NoiseGeneratorAlg.h"
#include "icaruscode/Utilities/quantities_utils.h" // util::...::value()
#include "lardataalg/Utilities/quantities_fhicl.h" // count_f from FHiCL

// framework libraries
#include "fhiclcpp/types/Atom.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine
#include "CLHEP/Random/RandGaussQ.h"

// C/C++ standard libraries
#include <vector>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename ADCT = double> class GaussianNoiseGeneratorAlg;
}
/**
 * @brief Gaussian noise generator algorithm for PMT electronics simulation.
 * @param ADCT type of the ADC count this interface treats with.
 * 
 * This noise generator produces a random Gaussian-distributed noise with
 * the specified RMS and no bias.
 * 
 * It requires a CLHEP random engine as input, and directly uses
 * `CLHEP::RandGaussQ` adaptor on it.
 * 
 * Generating for any type `ADCT` other that the CLHEP-native `double` requires
 * a conversion and slows down the generation.
 * 
 */
template <typename ADCT /* = double */>
class icarus::opdet::GaussianNoiseGeneratorAlg
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
  GaussianNoiseGeneratorAlg(Params_t params, CLHEP::HepRandomEngine& engine);
  
  
  /**
   * @brief Constructor: acquires configuration parameters as data structure.
   * @param params configuration parameters in FHiCL format
   * @param engine random engines
   */
  GaussianNoiseGeneratorAlg
    (Config const& config, CLHEP::HepRandomEngine& engine);

  
  
    protected:
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  Params_t const fParams; ///< All configuration parameters.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  
  /// Random engine used by this algorithm.
  CLHEP::HepRandomEngine& fRandomEngine;
  
  CLHEP::RandGaussQ fGausRandom; ///< Gaussian random extractor adapter.
  
  
  
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
  
  
}; // icarus::opdet::GaussianNoiseGeneratorAlg



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename ADCT>
icarus::opdet::GaussianNoiseGeneratorAlg<ADCT>::GaussianNoiseGeneratorAlg
  (Params_t params, CLHEP::HepRandomEngine& engine)
  : fParams{ std::move(params) }
  , fRandomEngine{ engine }
  , fGausRandom{ fRandomEngine, 0.0, static_cast<double>(value(fParams.RMS)) }
{}


// -----------------------------------------------------------------------------
template <typename ADCT>
icarus::opdet::GaussianNoiseGeneratorAlg<ADCT>::GaussianNoiseGeneratorAlg
  (Config const& config, CLHEP::HepRandomEngine& engine)
  : GaussianNoiseGeneratorAlg{ convert(config), engine }
{}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::GaussianNoiseGeneratorAlg<ADCT>::doAdd(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  std::vector<double> noise(n);
  fGausRandom.fireArray(static_cast<int>(n), noise.data());
  for (double const sample: noise)
    *(begin++) += static_cast<ADCcount_t>(sample);
  return n;
} // icarus::opdet::GaussianNoiseGeneratorAlg<>::doAdd()


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::GaussianNoiseGeneratorAlg<ADCT>::doFill(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  // if destination type is the same as generator natively uses,
  // we avoid conversions
  if constexpr(std::is_same_v<ADCcount_t, double>) {
    fGausRandom.fireArray(static_cast<int>(n), begin);
    return n;
  }
  else {
    std::vector<double> noise(n);
    fGausRandom.fireArray(static_cast<int>(n), noise.data());
    for (double const sample: noise)
      *(begin++) = static_cast<ADCcount_t>(sample);
    return n;
  }
} // icarus::opdet::GaussianNoiseGeneratorAlg<>::doFill()


// -----------------------------------------------------------------------------
template <typename ADCT>
void icarus::opdet::GaussianNoiseGeneratorAlg<ADCT>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  out << firstIndent
    << "Gaussian noise generator algorithm (GaussianNoiseGeneratorAlg):"
    << "\n" << indent << "  RMS: " << fParams.RMS;
} // icarus::opdet::GaussianNoiseGeneratorAlg<>::doDump()


// -----------------------------------------------------------------------------
template <typename ADCT>
auto icarus::opdet::GaussianNoiseGeneratorAlg<ADCT>::convert
  (Config const& config) -> Params_t
{
  return {
    ADCcount_t{ config.RMS() } // RMS
    };
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_GAUSSIANNOISEGENERATORALG_H


