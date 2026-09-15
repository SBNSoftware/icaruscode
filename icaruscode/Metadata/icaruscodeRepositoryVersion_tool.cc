/**
 * @file   RepositoryVersion_icaruscode_tool.cc
 * @brief  _art_ tool reporting the version of `icaruscode`-related packages.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 18, 2025
 * 
 */

// SBN libraries
#include "icaruscode/Metadata/RepositoryVersion_icaruscode.h"
#include "icarusalg/Metadata/RepositoryVersion_icarusalg.h"
#include "sbncode/Metadata/RepositoryVersionReportUtils.h"
#include "sbncode/Metadata/RepositoryVersionReporter.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard libraries
#include <string>
#include <vector>


// -----------------------------------------------------------------------------
namespace sbn { struct icaruscodeRepositoryVersion; }
/**
 * @brief Implements the `sbn::RepositoryVersionReporter` interface for
 * `icaruscode`.
 * 
 * It collects information from the following repositories: `icaruscode`,
 * `icarusalg` and `sbncode` (and its dependencies).
 * 
 */
struct sbn::icaruscodeRepositoryVersion: public sbn::RepositoryVersionReporter {
  
  struct Config {};
  
  using Parameters = art::ToolConfigTable<Config>;
  
  icaruscodeRepositoryVersion(Parameters const&);
  
}; // sbn::icaruscodeRepositoryVersion()


// -----------------------------------------------------------------------------
// ---  implementation
// -----------------------------------------------------------------------------
sbn::icaruscodeRepositoryVersion::icaruscodeRepositoryVersion(Parameters const&)
{
  
  // this will consult with `sbncode` plugin and report its dependencies too:
  addVersionFromRepository(packageVersions, "sbncode");
  
  packageVersions.items.emplace_back("icarusalg", ::RepositoryVersion_icarusalg);
  packageVersions.items.emplace_back("icaruscode", ::RepositoryVersion_icaruscode);
  
  packageVersions.finish();
  
}


// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(sbn::icaruscodeRepositoryVersion)


// -----------------------------------------------------------------------------
