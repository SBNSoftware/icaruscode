/**
 * @file   icaruscode/PMT/Algorithms/ConstantPedestalGeneratorAlg.h
 * @brief  A pedestal generation algorithm using a single constant baseline.
 * @date   November 17, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/PMT/Algorithms/PMTnoiseGenerator.h
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_CONSTANTPEDESTALGENERATORALG_H
#define ICARUSCODE_PMT_ALGORITHMS_CONSTANTPEDESTALGENERATORALG_H

// LArSoft/ICARUS libraries
#include "icaruscode/PMT/Algorithms/PedestalGeneratorAlg.h"
#include "icaruscode/PMT/Algorithms/NoiseGeneratorAlg.h"
#include "icaruscode/Utilities/quantities_utils.h" // util::...::value()
#include "lardataalg/Utilities/quantities_fhicl.h" // counts_f from FHiCL

// framework libraries
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <algorithm> // std::fill_n()
#include <memory> // std::unique_ptr<>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename ADCT = double> class ConstantPedestalGeneratorAlg;
}
/**
 * @brief Constant pedestal generator algorithm for PMT electronics simulation.
 * @param ADCT type of the ADC count this interface treats with.
 * 
 * This generator produces a constant pedestal, with the addition of the noise
 * from the configured noise generator.
 * 
 */
template <typename ADCT /* = double */>
class icarus::opdet::ConstantPedestalGeneratorAlg
  : public icarus::opdet::PedestalGeneratorAlg<ADCT>
{
  ///< Base class alias.
  using Base_t = icarus::opdet::PedestalGeneratorAlg<ADCT>;
  
    public:
  
  using ADCcount_t = typename Base_t::ADCcount_t;
  using Timestamp_t = typename Base_t::Timestamp_t;
  
  /// Underlying fundamental type of `ADCcount_t`, e.g. `float`.
  using ADCvalue_t = util::value_t<typename Base_t::ADCcount_t>;
  
  /// Type of the noise generator algorithm interface.
  using NoiseGenerator_t = icarus::opdet::NoiseGeneratorAlg<ADCcount_t>;
  
  
  // --- BEGIN -- Configuration ------------------------------------------------
  /// Algorithm configuration parameters.
  struct Params_t {
    ADCcount_t level; ///< Pedestal level.
  };
  
  
  /// Configuration parameters.
  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<ADCvalue_t> Pedestal{
      Name{ "Pedestal" },
      Comment{ "level of the pedestal" }
      // mandatory
      };
    
  }; // struct Config
  
  
  // --- END ---- Configuration ------------------------------------------------
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Constructor: acquires configuration parameters as data structure.
   * @param params configuration parameters
   * @param noiseGen noise generation algorithm (will be acquired and owned)
   */
  ConstantPedestalGeneratorAlg
    (Params_t params, std::unique_ptr<NoiseGenerator_t>&& noiseGen);
  
  
  /**
   * @brief Constructor: acquires configuration parameters as data structure.
   * @param params configuration parameters in FHiCL format
   * @param noiseGen noise generation algorithm (will be acquired and owned)
   */
  ConstantPedestalGeneratorAlg
    (Config const& config, std::unique_ptr<NoiseGenerator_t>&& noiseGen);

  
  
    protected:
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  Params_t const fParams; ///< All configuration parameters.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  
  std::unique_ptr<NoiseGenerator_t> fNoiseGen; ///< Noise generator algorithm.
  
  
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  virtual ADCcount_t doPedestalLevel
    (raw::Channel_t channel, Timestamp_t time) const override;
  
  
  /**
   * @brief Adds pedestal (with noise) to `n` samples starting at `begin`.
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
   * 
   * The pedestal is added _before_ the noise is added (that may be relevant
   * if the noise level depends on the value of the sample).
   */
  virtual std::size_t doAdd(
    raw::Channel_t channel, Timestamp_t time,
    ADCcount_t* begin, std::size_t n
    ) override;
  
  /**
   * @brief Overwrites `n` samples starting at `begin` with pedestal plus noise.
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
   * 
   * The pedestal is set _before_ the noise is added (that may be relevant
   * if the noise level depends on the value of the sample).
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
  
  
}; // icarus::opdet::ConstantPedestalGeneratorAlg



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename ADCT>
icarus::opdet::ConstantPedestalGeneratorAlg<ADCT>::ConstantPedestalGeneratorAlg
  (Params_t params, std::unique_ptr<NoiseGenerator_t>&& noiseGen)
  : fParams{ std::move(params) }
  , fNoiseGen{ std::move(noiseGen) }
{}


// -----------------------------------------------------------------------------
template <typename ADCT>
icarus::opdet::ConstantPedestalGeneratorAlg<ADCT>::ConstantPedestalGeneratorAlg
  (Config const& config, std::unique_ptr<NoiseGenerator_t>&& noiseGen)
  : ConstantPedestalGeneratorAlg{ convert(config), std::move(noiseGen) }
{}


// -----------------------------------------------------------------------------
template <typename ADCT>
auto icarus::opdet::ConstantPedestalGeneratorAlg<ADCT>::doPedestalLevel
  (raw::Channel_t channel, Timestamp_t time) const -> ADCcount_t
{
  return fParams.level;
}


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::ConstantPedestalGeneratorAlg<ADCT>::doAdd(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  if (fNoiseGen) n = fNoiseGen->add(channel, time, begin, n);
  ADCcount_t* const end = begin + n;
  while (begin != end) *(begin++) += fParams.level;
  return n;
} // icarus::opdet::ConstantPedestalGeneratorAlg<>::doAdd()


// -----------------------------------------------------------------------------
template <typename ADCT>
std::size_t icarus::opdet::ConstantPedestalGeneratorAlg<ADCT>::doFill(
  raw::Channel_t channel, Timestamp_t time,
  ADCcount_t* begin, std::size_t n
) {
  std::fill_n(begin, n, fParams.level);
  return fNoiseGen? fNoiseGen->add(channel, time, begin, n): n;
} // icarus::opdet::ConstantPedestalGeneratorAlg<>::doFill()


// -----------------------------------------------------------------------------
template <typename ADCT>
void icarus::opdet::ConstantPedestalGeneratorAlg<ADCT>::doDump(
  std::ostream& out,
  std::string const& indent, std::string const& firstIndent
) const {
  out << firstIndent
    << "Constant pedestal generator algorithm (ConstantPedestalGeneratorAlg):"
    << "\n" << indent << "  level: " << fParams.level
    << "\n" << indent << "  noise: ";
  if (fNoiseGen) fNoiseGen->dump(out, indent + "  ", "");
  else out << "none";
} // icarus::opdet::ConstantPedestalGeneratorAlg<>::doDump()


// -----------------------------------------------------------------------------
template <typename ADCT>
auto icarus::opdet::ConstantPedestalGeneratorAlg<ADCT>::convert
  (Config const& config) -> Params_t
{
  return {
    ADCcount_t{ config.Pedestal() } // level
    };
}


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_CONSTANTPEDESTALGENERATORALG_H


