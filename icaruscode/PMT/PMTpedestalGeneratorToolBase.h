/**
 * @file   icaruscode/PMT/PMTpedestalGeneratorToolBase.h
 * @brief  Standard base for implementation of tools of PMT pedestal generators.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 19, 2022
 * @see    `icaruscode/PMT/Algorithms/PMTgausPedestalGenerator.h`
 */

#ifndef ICARUSCODE_PMT_PMTNOISEGENERATORTOOLBASE_H
#define ICARUSCODE_PMT_PMTNOISEGENERATORTOOLBASE_H


// ICARUS libraries
#include "icaruscode/PMT/PMTpedestalGeneratorTool.h"
#include "icaruscode/PMT/PMTnoiseGeneratorTool.h"
#include "icaruscode/PMT/Algorithms/NoiseGeneratorAlg.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electronics.h" // counts_f
#include "lardataalg/Utilities/quantities_fhicl.h" // quantities from FHiCL

// framework libraries
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolConfigTable.h" // for the derived classes
#include "art/Utilities/ToolMacros.h" // for the derived classes
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/TableFragment.h"
#include "fhiclcpp/ParameterSet.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C/C++ standard libraries
#include <memory> // std::unique_ptr
#include <utility> // std::move()
#include <string>


//------------------------------------------------------------------------------
namespace icarus::opdet {
  template <template <typename> class ProvidedGenerator>
  struct PMTpedestalGeneratorToolBase;
}
/**
 * @brief Base for a standard tool wrapper for pedestal generator algorithms.
 * @see `icarus::opdet::PMTpedestalGeneratorTool`
 * 
 * This is a base class aimed to make the addition of tools for PMT pedestal
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
 *   : public icarus::opdet::PMTpedestalGeneratorToolBase<GenAlg>
 * {
 *   using Base_t = icarus::opdet::PMTpedestalGeneratorToolBase<GenAlg>;
 *   
 *   using Base_t::Base_t;
 *   
 * }; // icarus::opdet::PMTGenAlgTool
 * 
 * DEFINE_ART_CLASS_TOOL(icarus::opdet::PMTGenAlgTool)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * To use this tool, given the necessity of specifying more than just the FHiCL
 * configuration (in particular, a random number engine for the noise generator
 * algorithm is needed), a two-step procedure is necessary:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::unique_ptr<icarus::opdet::PMTpedestalGeneratorTool> genTool
 *   = art::make_tool<icarus::opdet::PMTpedestalGeneratorTool>(config);
 * std::unique_ptr<icarus::opdet::PedestalGeneratorAlg<util::quantities::counts_f>> genAlg
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
 *     art::make_tool<icarus::opdet::PMTpedestalGeneratorTool>(params().GenAlgo())
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
 * It is expected, though, that the configuration contains a `NoiseGenerator`
 * table for the configuration of the noise generation algorithm and tool.
 * If that mandatory table is empty, no noise generator algorithm will be
 * provided (which the pedestal generator may or may not be happy with).
 * 
 * 
 * Requirements
 * =============
 * 
 * Generator algorithm:
 *  * `Config` data structure with FHiCL configuration;
 *  * construction from a configuration table and a `std::unique_ptr` to a
 *    noise generator algorithm.
 * 
 * Implementations:
 *  * (optional) override the name of the tool (`toolName()`)
 * 
 */
template <template <typename ADCT> class ProvidedGenerator>
struct icarus::opdet::PMTpedestalGeneratorToolBase
  : public icarus::opdet::PMTpedestalGeneratorTool
{
  
  /// Type of generator actually being provided.
  using ProvidedGenerator_t = ProvidedGenerator<ADCcount_t>;
  
  /// FHiCL configuration of the provided generator.
  using ProvidedGeneratorConfig_t = typename ProvidedGenerator_t::Config;
  
  /// Configuration class for pedestal generation tools.
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::TableFragment<ProvidedGeneratorConfig_t> PedestalConfig;
    fhicl::DelegatedParameter NoiseGenerator{
      Name{ "NoiseGenerator" },
      Comment{
        "Configuration of the noise generation tool"
        " (from PMTnoiseGenerationTool)"
        }
      };
    
  }; // Config
  
  /// Tool parameter configuration.
  using Parameters = art::ToolConfigTable<Config>;
  
  /// Constructor: sets the configuration.
  PMTpedestalGeneratorToolBase(Parameters const& params)
    : fConfigParams(std::move(params))
    {}
  
  
    protected:
  
  /// Returns the name of this tool for error messages.
  virtual std::string toolName() const
    { return fConfigParams.get_PSet().template get<std::string>("tool_type"); }
  
  
    private:
  
  Parameters const fConfigParams; ///< Configuration parameters.
  
  
  // --- BEGIN -- Virtual interface --------------------------------------------
  
  /// Creates a generator with the stored configuration and the specified noise
  /// generator.
  virtual std::unique_ptr<Generator_t> doMakeGenerator
    (std::unique_ptr<NoiseGenerator_t>&& noiseGen) override;
  
  /// Creates a generator with the stored configuration and the specified random
  /// `engine`.
  virtual std::unique_ptr<Generator_t> doMakeGenerator
    (CLHEP::HepRandomEngine& engine) override;
  
  // --- END -- Virtual interface ----------------------------------------------
  
}; // icarus::opdet::PMTgausPedestalGeneratorTool


//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
template <template <typename ADCT> class ProvidedGenerator>
auto
icarus::opdet::PMTpedestalGeneratorToolBase<ProvidedGenerator>::doMakeGenerator
  (std::unique_ptr<NoiseGenerator_t>&& noiseGen) -> std::unique_ptr<Generator_t>
{
  return std::make_unique<ProvidedGenerator_t>
    (fConfigParams().PedestalConfig(), std::move(noiseGen));
} // icarus::opdet::PMTpedestalGeneratorToolBase::doMakeGenerator(noiseGen)


//------------------------------------------------------------------------------
template <template <typename ADCT> class ProvidedGenerator>
auto
icarus::opdet::PMTpedestalGeneratorToolBase<ProvidedGenerator>::doMakeGenerator
  (CLHEP::HepRandomEngine& engine) -> std::unique_ptr<Generator_t>
{
  using NoiseGenerator_t
    = icarus::opdet::NoiseGeneratorAlg<util::quantities::counts_f>;
  
  std::unique_ptr<NoiseGenerator_t> noiseGen;
  fhicl::ParameterSet const& noiseConfig
    = fConfigParams().NoiseGenerator.template get<fhicl::ParameterSet>();
  if (!noiseConfig.is_empty()) {
    noiseGen = art::make_tool<icarus::opdet::PMTnoiseGeneratorTool>(noiseConfig)
      ->makeGenerator(engine);
  }
  return makeGenerator(std::move(noiseGen));
} // icarus::opdet::PMTpedestalGeneratorToolBase::doMakeGenerator(engine)


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_PMTNOISEGENERATORTOOLBASE_H
