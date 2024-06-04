/**
 * @file   icaruscode/Geometry/ICARUSsplitInductionWireReadoutSetupTool_tool.cc
 * @brief  Tool to create ICARUS channel mapper with split induction wires.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 7, 2019
 * @see    `larcore/Geometry/WireReadoutSetupTool.h`
 *
 * This library is header-only.
 */

// ICARUS libraries
#include "icarusalg/Geometry/ICARUSWireReadoutGeom.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadoutSetupTool.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"


// -----------------------------------------------------------------------------
namespace icarus {
  class ICARUSsplitInductionWireReadoutSetupTool;
} // namespace icarus


/**
 * @brief Interface for a tool creating the standard ICARUS channel mapper.
 *
 * This tool creates a `icarus::ICARUSWireReadoutGeom` mapper.
 *
 */
class icarus::ICARUSsplitInductionWireReadoutSetupTool
  : public geo::WireReadoutSetupTool
{

  /// Type of channel mapping algorithm being created.
  using Mapper_t = icarus::ICARUSWireReadoutGeom;

    public:

  using Parameters = art::ToolConfigTable<Mapper_t::Config>;


  /// Constructor: passes all parameters to the channel mapping algorithm.
  ICARUSsplitInductionWireReadoutSetupTool(Parameters const& config);


    protected:

  Mapper_t::Config fParams;


  // --- BEGIN -- Virtual interface implementation ---------------------------
  /// @name Virtual interface
  /// @{

  /// Returns a pointer to the channel mapping.
  std::unique_ptr<geo::WireReadoutGeom> doWireReadout(
    geo::GeometryCore const* geom,
    std::unique_ptr<geo::WireReadoutSorter> sorter) override;

  /// @}
  // --- END -- Virtual interface implementation -----------------------------


}; // class icarus::ICARUSsplitInductionWireReadoutSetupTool


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarus::ICARUSsplitInductionWireReadoutSetupTool::ICARUSsplitInductionWireReadoutSetupTool
  (Parameters const& config)
    : fParams{config()}
{}


// -----------------------------------------------------------------------------
std::unique_ptr<geo::WireReadoutGeom>
icarus::ICARUSsplitInductionWireReadoutSetupTool::doWireReadout
  (geo::GeometryCore const* geom, std::unique_ptr<geo::WireReadoutSorter> sorter)
{
  return std::make_unique<Mapper_t>(fParams, geom, std::move(sorter));
}



// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(icarus::ICARUSsplitInductionWireReadoutSetupTool)


// -----------------------------------------------------------------------------
