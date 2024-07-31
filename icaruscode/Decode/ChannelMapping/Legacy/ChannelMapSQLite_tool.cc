/**
 * @file icaruscode/Decode/ChannelMapping/Legacy/ChannelMapSQLite_tool.cc
 * @brief _art_ tool wrapping ICARUS channel mapping SQLite database helper.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see icaruscode/Decode/ChannelMapping/ChannelMapSQLite.cxx
 * 
 * This tool is a bit unusual in _art_ standards in that the name of the tool
 * class, `ChannelMapSQLiteTool`, does not match the name of the file,
 * `ChannelMapSQLite_tool.cc` (normally it would be
 * `ChannelMapSQLiteTool_tool.cc`).
 * The choice of this source file name is for backward compatibility,
 * since before this refactoring a tool called `ChannelMapSQLite` existed.
 * The choice to change the name to _the class_ is that in the refactoring
 * `icarusDB::ChannelMapSQLite` still exists but it's not an _art_ tool any more
 * but rather a plain algorithm.
 * 
 * When choosing a tool, the configuration parameter `tool_type` is used to
 * put together the name expected for the library, and by default than name is
 * based on the name of the source file without suffix, i.e. `ChannelMapSQLite`:
 * so, `tool_type: ChannelMapSQLite` will be enough to find the library.
 * Once in the library, _art_ uses as the tool object the one specified by the
 * `DEFINE_ART_CLASS_TOOL()` macro, which is `icarusDB::ChannelMapSQLiteTool`.
 * So everything works.
 */

#include "icaruscode/Decode/ChannelMapping/ChannelMapSQLite.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/TableFragment.h"


namespace icarusDB { struct ChannelMapSQLiteTool; }

/// Toolification of `icarusDB::ChannelMapSQLite`.
struct icarusDB::ChannelMapSQLiteTool: public icarusDB::ChannelMapSQLite {
  
  struct Config {
    
    fhicl::TableFragment<icarusDB::ChannelMapSQLite::Config> ChannelMapSQLiteConfig;
    
  };
  
  using Parameters = art::ToolConfigTable<Config>;
  
  ChannelMapSQLiteTool(Parameters const& params)
    : icarusDB::ChannelMapSQLite{ params().ChannelMapSQLiteConfig() }
    {}
  
}; // icarusDB::ChannelMapSQLiteTool


DEFINE_ART_CLASS_TOOL(icarusDB::ChannelMapSQLiteTool)


