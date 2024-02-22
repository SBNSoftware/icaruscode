/**
 * @file icaruscode/Decode/ChannelMapping/ChannelMapPostGres_tool.cc
 * @brief _art_ tool wrapping ICARUS channel mapping PostGres database helper.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see icaruscode/Decode/ChannelMapping/ChannelMapPostGres.cxx
 * 
 * This tool is a bit unusual in _art_ standards in that the name of the tool
 * class, `ChannelMapPostGresTool`, does not match the name of the file,
 * `ChannelMapPostGres_tool.cc` (normally it would be
 * `ChannelMapPostGresTool_tool.cc`).
 * The choice of this source file name is for backward compatibility,
 * since before this refactoring a tool called `ChannelMapPostGres` existed.
 * The choice to change the name to _the class_ is that in the refactoring
 * `icarusDB::ChannelMapPostGres` still exists but it's not an _art_ tool any
 * more but rather a plain algorithm.
 * 
 * When choosing a tool, the configuration parameter `tool_type` is used to
 * put together the name expected for the library, and by default than name is
 * based on the name of the source file without suffix, i.e.
 * `ChannelMapPostGres`: so, `tool_type: ChannelMapPostGres` will be enough to
 * find the library. Once in the library, _art_ uses as the tool object the one
 * specified by the `DEFINE_ART_CLASS_TOOL()` macro, which is
 * `icarusDB::ChannelMapPostGresTool`. So everything works.
 */

#include "icaruscode/Decode/ChannelMapping/ChannelMapPostGres.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/TableFragment.h"


namespace icarusDB { class ChannelMapPostGresTool; }

/// Toolification of `icarusDB::ChannelMapPostGres`.
struct icarusDB::ChannelMapPostGresTool: public icarusDB::ChannelMapPostGres {
  
  struct Config {
    
    fhicl::TableFragment<icarusDB::ChannelMapPostGres::Config> ChannelMapPostGresConfig;
    
  };
  
  using Parameters = art::ToolConfigTable<Config>;
  
  ChannelMapPostGresTool(Parameters const& params)
    : icarusDB::ChannelMapPostGres{ params().ChannelMapPostGresConfig() }
    {}
  
}; // icarusDB::ChannelMapPostGresTool


DEFINE_ART_CLASS_TOOL(icarusDB::ChannelMapPostGresTool)
