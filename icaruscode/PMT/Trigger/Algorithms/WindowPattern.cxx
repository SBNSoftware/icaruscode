/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WindowPattern.cxx
 * @brief  Defines a (sliding) window trigger pattern (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 21, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/WindowPattern.h
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/WindowPattern.h"


//------------------------------------------------------------------------------
std::string icarus::trigger::WindowPattern::tag() const {
  using namespace std::string_literals;
  
  std::string s;
  
  s += "M"s + std::to_string(minInMainWindow);
  
  if (minInOppositeWindow > 0U)
    s += "O"s + std::to_string(minInOppositeWindow);
  
  if ((minInDownstreamWindow > 0U) || requireDownstreamWindow) {
    s += "D"s + std::to_string(minInDownstreamWindow);
    if (requireDownstreamWindow) s+= "req"s;
  } // if downstream
  
  if ((minInUpstreamWindow > 0U) || requireUpstreamWindow) {
    s += "U"s + std::to_string(minInUpstreamWindow);
    if (requireUpstreamWindow) s+= "req"s;
  } // if upstream
  
  return s;
} // icarus::trigger::WindowPattern::description()


//------------------------------------------------------------------------------
std::string icarus::trigger::WindowPattern::description() const {
  using namespace std::string_literals;
  
  std::string s = "required:";
  
  s += " "s + std::to_string(minInMainWindow);
  
  if (minInOppositeWindow > 0U)
    s += " + " + std::to_string(minInOppositeWindow) + " (opposite)"s;
  
  if ((minInDownstreamWindow > 0U) || requireDownstreamWindow) {
    if (minInDownstreamWindow) s += " + " + std::to_string(minInDownstreamWindow);
    s += " (downstream";
    if (requireDownstreamWindow) s+= ", mandatory)"s;
    s += ")"s;
  } // if downstream
  
  if ((minInUpstreamWindow > 0U) || requireUpstreamWindow) {
    if (minInUpstreamWindow) s += " + " + std::to_string(minInUpstreamWindow);
    s += " (upstream";
    if (requireUpstreamWindow) s+= ", mandatory)"s;
    s += ")"s;
  } // if upstream
  
  return s;
} // icarus::trigger::WindowPattern::description()


//------------------------------------------------------------------------------
