/**
 * @file   icaruscode/Geometry/ICARUSsingleInductionChannelMapSetupTool_tool.cc
 * @brief  Tool to create ICARUS channel mapper with split induction wires.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 7, 2019
 * @see    `larcore/Geometry/ChannelMapSetupTool.h`
 * 
 * This library is header-only.
 */

// ICARUS libraries
#include "icaruscode/Geometry/ChannelMapIcarusAlg.h"

// LArSoft libraries
#include "larcore/Geometry/ChannelMapSetupTool.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"


// -----------------------------------------------------------------------------
namespace icarus {
  class ICARUSsingleInductionChannelMapSetupTool;
} // namespace icarus


/**
 * @brief Interface for a tool creating the standard ICARUS channel mapper.
 * 
 * This tool creates a `geo::ChannelMapIcarusAlg` mapper.
 * 
 */
class icarus::ICARUSsingleInductionChannelMapSetupTool
  : public geo::ChannelMapSetupTool
{
  
  /// Type of channel mapping algorithm being created.
  using Mapper_t = geo::ChannelMapIcarusAlg;
  
    public:
  
  /// Constructor: passes all parameters to the channel mapping algorithm.
  ICARUSsingleInductionChannelMapSetupTool(fhicl::ParameterSet const& pset);
  
  
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
    (fhicl::ParameterSet const& pset) const;
  
}; // class icarus::ICARUSsingleInductionChannelMapSetupTool


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarus::ICARUSsingleInductionChannelMapSetupTool::ICARUSsingleInductionChannelMapSetupTool
  (fhicl::ParameterSet const& pset)
  : fChannelMap(makeChannelMap(pset))
{}


// -----------------------------------------------------------------------------
auto icarus::ICARUSsingleInductionChannelMapSetupTool::makeChannelMap
  (fhicl::ParameterSet const& pset) const -> std::unique_ptr<Mapper_t>
{
  return std::make_unique<Mapper_t>(pset);
} // icarus::ICARUSsingleInductionChannelMapSetupTool::makeChannelMap()



// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::ICARUSsingleInductionChannelMapSetupTool)


// -----------------------------------------------------------------------------
