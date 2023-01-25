/**
 * @file   icaruscode/PMT/PMTgausNoiseGeneratorTool_tool.cc
 * @brief  Toolization of `icarus::opdet::GaussianNoiseGeneratorAlg<counts_f>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/GaussianNoiseGenerator.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::PMTnoiseGeneratorTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/PMTnoiseGeneratorToolBase.h"
#include "icaruscode/PMT/Algorithms/GaussianNoiseGeneratorAlg.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"


//------------------------------------------------------------------------------
namespace icarus::opdet { struct PMTgausNoiseGeneratorTool; }
/**
 * @brief Creates a `icarus::opdet::GaussianNoiseGeneratorAlg` algorithm.
 * @see `icarus::opdet::GaussianNoiseGeneratorAlg`
 * 
 * This tool creates a `icarus::opdet::GaussianNoiseGeneratorAlg<counts_f>`
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
 * `icarus::opdet::GaussianNoiseGeneratorAlg::Config` data structure) for a
 * short explanation of the meaning of the parameters, or check the
 * documentation of `icarus::opdet::GaussianNoiseGeneratorAlg`.
 * 
 */
struct icarus::opdet::PMTgausNoiseGeneratorTool
  : public icarus::opdet::PMTnoiseGeneratorToolBase<GaussianNoiseGeneratorAlg>
{
  using Base_t
    = icarus::opdet::PMTnoiseGeneratorToolBase<GaussianNoiseGeneratorAlg>;
  
  using Base_t::Base_t;
  
}; // icarus::opdet::PMTgausNoiseGeneratorTool


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::PMTgausNoiseGeneratorTool)


//------------------------------------------------------------------------------

