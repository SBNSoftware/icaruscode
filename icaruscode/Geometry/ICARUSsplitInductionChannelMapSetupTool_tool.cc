/**
 * @file   icaruscode/Geometry/ICARUSsplitInductionChannelMapSetupTool_tool.cc
 * @brief  Tool to create ICARUS channel mapper with split induction wires.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 7, 2019
 * @see    `larcore/Geometry/ChannelMapSetupTool.h`
 * 
 * This library is header-only.
 */

// ICARUS libraries
#include "icaruscode/Geometry/ICARUSChannelMapAlg.h"

// LArSoft libraries
#include "larcore/Geometry/ChannelMapSetupTool.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"


// -----------------------------------------------------------------------------
namespace icarus {
  class ICARUSsplitInductionChannelMapSetupTool;
} // namespace icarus


/**
 * @brief Interface for a tool creating the standard ICARUS channel mapper.
 * 
 * This tool creates a `icarus::ICARUSChannelMapAlg` mapper.
 * 
 */
class icarus::ICARUSsplitInductionChannelMapSetupTool
  : public geo::ChannelMapSetupTool
{
  
  /// Type of channel mapping algorithm being created.
  using Mapper_t = icarus::ICARUSChannelMapAlg;
  
    public:
  
  using Parameters = art::ToolConfigTable<Mapper_t::Config>;
  
  
  /// Constructor: passes all parameters to the channel mapping algorithm.
  ICARUSsplitInductionChannelMapSetupTool(Parameters const& config);
  
  
    protected:
  
  std::unique_ptr<Mapper_t> fChannelMap;
  
  
  // --- BEGIN -- Virtual interface implementation ---------------------------
  /// @name Virtual interface
  /// @{
  
  /// Returns a pointer to the channel mapping.
  virtual std::unique_ptr<geo::ChannelMapAlg> doChannelMap() override
    { return std::move(fChannelMap); }
  
  /// @}
  // --- END -- Virtual interface implementation -----------------------------
  
  
  /// Creates and returns the channel mapping algorithm.
  std::unique_ptr<Mapper_t> makeChannelMap
    (Mapper_t::Config const& config) const;
  
}; // class icarus::ICARUSsplitInductionChannelMapSetupTool


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarus::ICARUSsplitInductionChannelMapSetupTool::ICARUSsplitInductionChannelMapSetupTool
  (Parameters const& config)
  : fChannelMap(makeChannelMap(config()))
{}


// -----------------------------------------------------------------------------
auto icarus::ICARUSsplitInductionChannelMapSetupTool::makeChannelMap
  (Mapper_t::Config const& config) const -> std::unique_ptr<Mapper_t>
{
  return std::make_unique<Mapper_t>(config);
} // icarus::ICARUSsplitInductionChannelMapSetupTool::makeChannelMap()



// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::ICARUSsplitInductionChannelMapSetupTool)


// -----------------------------------------------------------------------------
