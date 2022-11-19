/**
 * @file   icaruscode/PMT/PMTnoNoiseGeneratorTool_tool.cc
 * @brief  Toolization of `icarus::opdet::NoNoiseGeneratorAlg<counts_f>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/PMTgausNoiseGenerator.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::PMTnoiseGeneratorTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/PMTnoiseGeneratorToolBase.h"
#include "icaruscode/PMT/Algorithms/NoNoiseGeneratorAlg.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"


//------------------------------------------------------------------------------
namespace icarus::opdet { struct PMTnoNoiseGeneratorTool; }
/**
 * @brief Creates a `icarus::opdet::NoNoiseGeneratorAlg` algorithm.
 * @see `icarus::opdet::NoNoiseGeneratorAlg`
 * 
 * This tool creates a `icarus::opdet::NoNoiseGeneratorAlg<counts_f>`
 * "noise" generator algorithm.
 * 
 * The usage pattern of this algorithm is the standard one documented in
 * `icarus::opdet::PMTnoiseGeneratorToolBase`.
 * 
 * 
 * Configuration
 * --------------
 * 
 * Run `lar --print-description PMTgausNoiseGeneratorTool` (or read
 * `icarus::opdet::NoNoiseGeneratorAlg::Config` data structure) for a
 * short explanation of the meaning of the parameters, or check the
 * documentation of `icarus::opdet::NoNoiseGeneratorAlg`.
 * 
 */
class icarus::opdet::PMTnoNoiseGeneratorTool
  : public icarus::opdet::PMTnoiseGeneratorToolBase<NoNoiseGeneratorAlg>
{
  using Base_t
    = icarus::opdet::PMTnoiseGeneratorToolBase<NoNoiseGeneratorAlg>;
  
  using Base_t::Base_t;
  
}; // icarus::opdet::PMTnoNoiseGeneratorTool


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::PMTnoNoiseGeneratorTool)


//------------------------------------------------------------------------------

