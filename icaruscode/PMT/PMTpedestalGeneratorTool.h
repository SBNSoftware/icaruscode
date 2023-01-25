/**
 * @file   icaruscode/PMT/PMTpedestalGeneratorTool.h
 * @brief  Tool to create a PMT pedestal generator.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 19, 2022
 * @see    `icaruscode/PMT/Algorithms/PMTpedestalGenerator.h`
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_PMT_PMTPEDESTALGENERATORTOOL_H
#define ICARUSCODE_PMT_PMTPEDESTALGENERATORTOOL_H


// ICARUS libraries
#include "icaruscode/PMT/Algorithms/PedestalGeneratorAlg.h"
#include "icaruscode/PMT/Algorithms/NoiseGeneratorAlg.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electronics.h" // counts_f

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C/C++ standard libraries
#include <memory> // std::unique_ptr()
#include <utility> // std::move()
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::opdet { struct PMTpedestalGeneratorTool; }
/**
 * @brief Creates a `PMTpedestalGenerator` algorithm object.
 * 
 * Implementations of this tool create an algorithm to emit PMT readout
 * electronics pedestal, as an object derived from
 * `icarus::opdet::PMTpedestalGenerator<util::quantities::counts_f>`
 * (the template parameter is the same used by the single photoelectron
 * response, but it's not programmatically tied to it to avoid entangling
 * dependencies).
 * 
 * The tool object is just a wrapper factory that creates the actual generator,
 * which is owned and can be returned for use.
 * 
 * The signature of the tool constructor should be compatible with:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * PMTpedestalGeneratorTool(fhicl::ParameterSet const&)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * in order to be compatible with `art::make_tool()`.
 * 
 * The only purpose of the tool itself is to create and return the pedestal
 * generation algorithm. This may be done by calling `makeGenerator()` either
 * with the noise generator to be used or with the random engine to use to
 * create a new noise generator (the latter functionality is not guaranteed for
 * all tool implementations).
 * 
 */
struct icarus::opdet::PMTpedestalGeneratorTool {
  
  using ADCcount_t = util::quantities::counts_f; ///< Type of ADC count used.
  
  /// Base type of algorithm returned.
  using Generator_t = icarus::opdet::PedestalGeneratorAlg<ADCcount_t>;
  
  /// Type of noise generation algorithm passed to the pedestal generator.
  using NoiseGenerator_t = icarus::opdet::NoiseGeneratorAlg<ADCcount_t>;
  
  
  virtual ~PMTpedestalGeneratorTool() = default;
  
  /**
   * @brief Returns an instance of the pedestal generator.
   * @param noiseGen noise generation algorithm (will be acquired and owned)
   * 
   * This function is guaranteed to return a valid pulse function only once.
   * The return value on following calls is undefined, but it is considered
   * an error to return a null pointer.
   */
  std::unique_ptr<Generator_t> makeGenerator
    (std::unique_ptr<NoiseGenerator_t>&& noiseGen);
  
  /**
   * @brief Returns an instance of the pedestal generator.
   * @param engine the random engine used for the noise generation algorithm
   * 
   * This function will attempt to create the proper noise generation algorithm
   * from the configuration, using the specified random `engine`, and then
   * proceed to the creation of the pedestal generation algorithm utilizing
   * the noise generator just created.
   * 
   * It is possible that the creation of a certain generator with this procedure
   * is not possible, in which case an exception should be expected.
   * This function is guaranteed to return a valid pulse function only once.
   * The return value on following calls is undefined, but it is considered
   * an error to return a null pointer.
   */
  std::unique_ptr<Generator_t> makeGenerator
    (CLHEP::HepRandomEngine& engine);
  
  
    private:
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /**
   * @brief Returns an instance of the pedestal generator
   * @param noiseGen noise generation algorithm (will be acquired and owned)
   * 
   * Only the first call is guaranteed to return a non-null pointer.
   * 
   */
  virtual std::unique_ptr<Generator_t> doMakeGenerator
    (std::unique_ptr<NoiseGenerator_t>&& noiseGen) = 0;
  
  
  /**
   * @brief Returns an instance of the pedestal generator
   * @param engine the random engine used for the noise generation algorithm
   * 
   * Only the first call is guaranteed to return a non-null pointer.
   */
  virtual std::unique_ptr<Generator_t> doMakeGenerator
    (CLHEP::HepRandomEngine& engine) = 0;


  
  // --- END -- Virtual interface ----------------------------------------------
  
  
}; // icarus::opdet::PMTpedestalGeneratorTool


//------------------------------------------------------------------------------
//--- icarus::opdet::PMTpedestalGeneratorTool implementation
//------------------------------------------------------------------------------
inline auto icarus::opdet::PMTpedestalGeneratorTool::makeGenerator
  (std::unique_ptr<NoiseGenerator_t>&& noiseGen) -> std::unique_ptr<Generator_t>
{
  std::unique_ptr<Generator_t> ptr = doMakeGenerator(std::move(noiseGen));
  assert(ptr);
  return ptr;
} // icarus::opdet::PMTpedestalGeneratorTool::makeGenerator()


//------------------------------------------------------------------------------
inline auto icarus::opdet::PMTpedestalGeneratorTool::makeGenerator
  (CLHEP::HepRandomEngine& engine) -> std::unique_ptr<Generator_t>
{
  std::unique_ptr<Generator_t> ptr = doMakeGenerator(engine);
  assert(ptr);
  return ptr;
} // icarus::opdet::PMTpedestalGeneratorTool::makeGenerator()


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_PMTPEDESTALGENERATORTOOL_H
