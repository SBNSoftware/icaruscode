/**
 * @file   icaruscode/PMT/PMTnoiseGeneratorTool.h
 * @brief  Tool to create a PMT noise generator.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 17, 2022
 * @see    `icaruscode/PMT/Algorithms/PMTnoiseGenerator.h`
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_PMT_PMTNOISEGENERATORTOOL_H
#define ICARUSCODE_PMT_PMTNOISEGENERATORTOOL_H


// ICARUS libraries
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
namespace icarus::opdet { struct PMTnoiseGeneratorTool; }
/**
 * @brief Creates a `PMTnoiseGenerator` algorithm object.
 * 
 * Implementations of this tool create an algorithm to emit PMT readout
 * electronics noise, as an object derived from
 * `icarus::opdet::PMTnoiseGenerator<util::quantities::counts_f>`
 * (the template parameter is the same used by the single photoelectron
 * response, but it's not programmatically tied to it to avoid entangling
 * dependencies).
 * 
 * The tool object is just a wrapper factory that creates the actual generator,
 * which is owned and can be returned for use.
 * 
 * The signature of the tool constructor should be compatible with:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * PMTnoiseGeneratorTool(fhicl::ParameterSet const&)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * in order to be compatible with `art::make_tool()`.
 * 
 * 
 */
struct icarus::opdet::PMTnoiseGeneratorTool {
  
  using ADCcount_t = util::quantities::counts_f; ///< Type of ADC count used.
  
  /// Type of function returned.
  using Generator_t = icarus::opdet::NoiseGeneratorAlg<ADCcount_t>;
  
  
  virtual ~PMTnoiseGeneratorTool() = default;
  
  /**
   * @brief Returns an instance of the noise generator.
   * @param engine the random engine to be used by the generator.
   * 
   * This function is guaranteed to return a valid generator only once.
   * The return value on following calls is undefined, but it is considered
   * an error to return a null pointer.
   */
  std::unique_ptr<Generator_t> makeGenerator(CLHEP::HepRandomEngine& engine);
  
  
    private:
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /**
   * @brief Returns an instance of the noise generator.
   * @param engine the random engine to be used by the generator.
   * 
   * Only the first call is guaranteed to return a non-null pointer.
   * 
   */
  virtual std::unique_ptr<Generator_t> doMakeGenerator
    (CLHEP::HepRandomEngine& engine) = 0;

  
  // --- END -- Virtual interface ----------------------------------------------
  
  
}; // icarus::opdet::PMTnoiseGeneratorTool


//------------------------------------------------------------------------------
//--- icarus::opdet::PMTnoiseGeneratorTool implementation
//------------------------------------------------------------------------------
inline auto icarus::opdet::PMTnoiseGeneratorTool::makeGenerator
  (CLHEP::HepRandomEngine& engine) -> std::unique_ptr<Generator_t>
{
  std::unique_ptr<Generator_t> ptr = doMakeGenerator(engine);
  assert(ptr);
  return ptr;
} // icarus::opdet::PMTnoiseGeneratorTool::makeGenerator()


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_PMTNOISEGENERATORTOOL_H
