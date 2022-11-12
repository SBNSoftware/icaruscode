/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WindowChannelMap.cxx
 * @brief  Applies sliding window trigger patterns.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/WindowChannelMap.h
 * 
 * Linkage to the implementation file is required when using `dump()` methods.
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/WindowChannelMap.h"

// ICARUS libraries
#include "icarusalg/Utilities/IntegerRanges.h"

// C/C++ standard libraries
#include <ostream>


//------------------------------------------------------------------------------
void icarus::trigger::WindowChannelMap::WindowInfo_t::dump
  (std::ostream& out, std::string const& indent /* = "" */) const
{
  
  out << indent << "window #" << topology.index
    << " at " << composition.center << " cm";
  if (composition.hasCryostat()) out << " in " << composition.cryoid;
  else               out << " (cryostat undefined)";
  out << " includes " << composition.channels.size() << " channels";
  if (!composition.channels.empty())
    out << " (" << icarus::makeIntegerRanges(composition.channels) << ")";
  if (topology.hasOppositeWindow())
    out << " opposite to [#" << topology.opposite << "]";
  if (topology.hasUpstreamWindow())
    out << " downstream of [#" << topology.upstream << "]";
  if (topology.hasDownstreamWindow())
    out << " upstream of [#" << topology.downstream << "]";
} // icarus::trigger::WindowChannelMap::WindowInfo_t::dump()


//------------------------------------------------------------------------------
void icarus::trigger::WindowChannelMap::dump
  (std::ostream& out, std::string const& indent, std::string const& firstIndent)
  const
{
  out << firstIndent << "Map has " << nWindows() << " windows:";
  for (WindowInfo_t const& info: fWindows) {
    out << "\n  "; // additional indentation
    info.dump(out, indent);
  } // for
  out << '\n';
} // icarus::trigger::WindowChannelMap::dump()


//------------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, icarus::trigger::WindowChannelMap const& wi)
  { wi.dump(out); return out; }


//------------------------------------------------------------------------------
