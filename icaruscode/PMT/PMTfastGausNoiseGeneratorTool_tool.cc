/**
 * @file   icaruscode/PMT/PMTfastGausNoiseGeneratorTool_tool.cc
 * @brief  Toolization of `icarus::opdet::FastGaussianNoiseGeneratorAlg<counts_f>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/GastGaussianNoiseGenerator.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::PMTnoiseGeneratorTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/PMTnoiseGeneratorToolBase.h"
#include "icaruscode/PMT/Algorithms/FastGaussianNoiseGeneratorAlg.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"


//------------------------------------------------------------------------------
namespace icarus::opdet { struct PMTfastGausNoiseGeneratorTool; }
/**
 * @brief Creates a `icarus::opdet::FastGaussianNoiseGeneratorAlg` algorithm.
 * @see `icarus::opdet::FastGaussianNoiseGeneratorAlg`
 * 
 * This tool creates a `icarus::opdet::FastGaussianNoiseGeneratorAlg<counts_f>`
 * noise generator algorithm.
 * 
 * The usage pattern of this algorithm is the standard one documented in
 * `icarus::opdet::PMTnoiseGeneratorToolBase`.
 * 
 * 
 * Configuration
 * --------------
 * 
 * Run `lar --print-description PMTgausNoiseGeneratorTool` (or read
 * `icarus::opdet::FastGaussianNoiseGeneratorAlg::Config` data structure) for a
 * short explanation of the meaning of the parameters, or check the
 * documentation of `icarus::opdet::FastGaussianNoiseGeneratorAlg`.
 * 
 */
class icarus::opdet::PMTfastGausNoiseGeneratorTool
  : public icarus::opdet::PMTnoiseGeneratorToolBase<FastGaussianNoiseGeneratorAlg>
{
  using Base_t
    = icarus::opdet::PMTnoiseGeneratorToolBase<FastGaussianNoiseGeneratorAlg>;
  
  using Base_t::Base_t;
  
}; // icarus::opdet::PMTfastGausNoiseGeneratorTool


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::PMTfastGausNoiseGeneratorTool)


//------------------------------------------------------------------------------

