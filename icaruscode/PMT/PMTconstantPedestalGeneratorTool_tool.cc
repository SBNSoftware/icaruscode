/**
 * @file   icaruscode/PMT/PMTconstantPedestalGeneratorTool_tool.cc
 * @brief  Toolization of `icarus::opdet::ConstantPedestalGeneratorAlg<counts_f>`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 17, 2020
 * @see    `icaruscode/PMT/Algorithms/PMTconstantPedestalGenerator.h`
 * 
 * This is an implementation of tool interface
 * `icarus::opdet::PMTpedestalGeneratorTool`.
 */


// ICARUS libraries
#include "icaruscode/PMT/PMTpedestalGeneratorToolBase.h"
#include "icaruscode/PMT/Algorithms/ConstantPedestalGeneratorAlg.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"


//------------------------------------------------------------------------------
namespace icarus::opdet { struct PMTconstantPedestalGeneratorTool; }
/**
 * @brief Creates a `icarus::opdet::ConstantPedestalGeneratorAlg` algorithm.
 * @see `icarus::opdet::ConstantPedestalGeneratorAlg`
 * 
 * This tool creates a `icarus::opdet::ConstantPedestalGeneratorAlg<counts_f>`
 * pedestal generator algorithm, armed with the noise generator specified in
 * the configuration.
 * 
 * The usage pattern of this algorithm is the standard one documented in
 * `icarus::opdet::PMTpedestalGeneratorToolBase`.
 * 
 * 
 * Configuration
 * --------------
 * 
 * Run `lar --print-description PMTconstantPedestalGeneratorTool` for a
 * short explanation of the meaning of the parameters; or read
 * `icarus::opdet::ConstantPedestalGeneratorAlg::Config` data structure, or
 * check the documentation of `icarus::opdet::ConstantPedestalGeneratorAlg`,
 * plus the one in `icarus::opdet::PMTpedestalGeneratorToolBase` which explain
 * the configuration of the noise generator tool within.
 * 
 */
struct icarus::opdet::PMTconstantPedestalGeneratorTool
  : public icarus::opdet::PMTpedestalGeneratorToolBase<ConstantPedestalGeneratorAlg>
{
  using Base_t
    = icarus::opdet::PMTpedestalGeneratorToolBase<ConstantPedestalGeneratorAlg>;
  
  using Base_t::Base_t;
  
}; // icarus::opdet::PMTconstantPedestalGeneratorTool


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::opdet::PMTconstantPedestalGeneratorTool)


//------------------------------------------------------------------------------

