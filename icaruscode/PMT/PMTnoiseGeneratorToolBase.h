/**
 * @file   icaruscode/PMT/PMTnoiseGeneratorToolBase.h
 * @brief  Standard base for implementation of tools of PMT noise generators.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 17, 2022
 * @see    `icaruscode/PMT/Algorithms/PMTgausNoiseGenerator.h`
 */

#ifndef ICARUSCODE_PMT_PMTNOISEGENERATORTOOLBASE_H
#define ICARUSCODE_PMT_PMTNOISEGENERATORTOOLBASE_H


// ICARUS libraries
#include "icaruscode/PMT/PMTnoiseGeneratorTool.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities_fhicl.h" // quantities from FHiCL

// framework libraries
#include "art/Utilities/ToolConfigTable.h" // for the derived classes
#include "art/Utilities/ToolMacros.h" // for the derived classes
#include "canvas/Utilities/Exception.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C/C++ standard libraries
#include <memory> // std::unique_ptr()
#include <string>


//------------------------------------------------------------------------------
namespace icarus::opdet {
  template <template <typename> class ProvidedGenerator>
  struct PMTnoiseGeneratorToolBase;
}
/**
 * @brief Base class for a standard tool wrapper for noise generator algorithms.
 * @see `icarus::opdet::PMTnoiseGeneratorTool`
 * 
 * This is a base class aimed to make the addition of tools for PMT noise
 * generators faster by implementing most of the boilerplate code, at the price
 * of a loss of flexibility.
 * 
 * 
 * Usage pattern
 * --------------
 * 
 * To have a tool serving the algorithm `GenAlg`, a new class must be defined
 * and declared a tool:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * class icarus::opdet::PMTGenAlgTool
 *   : public icarus::opdet::PMTnoiseGeneratorToolBase<GenAlg>
 * {
 *   using Base_t = icarus::opdet::PMTnoiseGeneratorToolBase<GenAlg>;
 *   
 *   using Base_t::Base_t;
 *   
 * }; // icarus::opdet::PMTGenAlgTool
 * 
 * DEFINE_ART_CLASS_TOOL(icarus::opdet::PMTGenAlgTool)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * To use this tool, given the necessity of specifying more than just the FHiCL
 * configuration (in particular, the random generator engine is needed), a
 * two-step procedure is necessary:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::unique_ptr<icarus::opdet::PMTnoiseGeneratorTool> genTool
 *   = art::make_tool<icarus::opdet::PMTnoiseGeneratorTool>(config);
 * std::unique_ptr<icarus::opdet::NoiseGeneratorAlg<util::quantities::counts_f>> genAlg
 *   = genTool->makeGenerator(engine);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * where `config` is better served as a `fhicl::DelegatedParameter` object.
 * The tool can be destroyed after the algorithm is created, as the ownership
 * of the algorithm is given to the caller.
 * Of course, to initialize a member variable in the initialization list of a
 * constructor the two steps may be merged, provided that the `engine` is
 * already available:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * : fGenAlg{
 *     art::make_tool<icarus::opdet::PMTnoiseGeneratorTool>(params().GenAlgo())
 *       ->makeGenerator(engine)
 *     }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 
 * Configuration
 * --------------
 * 
 * The configuration is specific to the featured algorithm.
 * Run `lar --print-description` on the tool or read the algorithm `Config` data
 * structure for a short explanation of the meaning of the parameters, or check
 * the documentation of the generator algorithms.
 * 
 * 
 * 
 * Requirements
 * =============
 * 
 * Generator algorithm:
 *  * `Config` data structure with FHiCL configuration;
 *  * construction from a configuration table and a CLHEP random engine
 * 
 * 
 * Implementations:
 *  * (optional) override the name of the tool (`toolName()`)
 * 
 */
template <template <typename ADCT> class ProvidedGenerator>
struct icarus::opdet::PMTnoiseGeneratorToolBase
  : public icarus::opdet::PMTnoiseGeneratorTool
{
  
  /// Type of generator actually being provided.
  using ProvidedGenerator_t = ProvidedGenerator<ADCcount_t>;
  
  /// Tool parameter configuration.
  using Parameters = art::ToolConfigTable<typename ProvidedGenerator_t::Config>;
  
  /// Constructor: sets the configuration.
  PMTnoiseGeneratorToolBase(Parameters const& params)
    : fConfigParams(std::move(params))
    {}
  
  
    protected:
  
  /// Returns the name of this tool for error messages.
  virtual std::string toolName() const
    { return fConfigParams.get_PSet().template get<std::string>("tool_type"); }
  
  
    private:
  
  Parameters const fConfigParams; ///< Configuration parameters.
  
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /// Creates a generator with the stored configuration and the specified random
  /// `engine`.
  virtual std::unique_ptr<Generator_t> doMakeGenerator
    (CLHEP::HepRandomEngine& engine) override
    { return std::make_unique<ProvidedGenerator_t>(fConfigParams(), engine); }
  
  // --- END -- Virtual interface ----------------------------------------------
  
}; // icarus::opdet::PMTgausNoiseGeneratorTool


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_PMTNOISEGENERATORTOOLBASE_H
