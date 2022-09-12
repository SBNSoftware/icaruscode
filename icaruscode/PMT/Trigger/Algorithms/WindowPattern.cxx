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
bool icarus::trigger::WindowPattern::isMainRequirementRelevant() const {
  
  /*
   * If the main window requirement is no larger than half the requested sum
   * (rounded up) then the requirement is not relevant, although in the special
   * cases where `S` is odd and `M` is its half rounded up (e.g. `M3S5`), the
   * main window requirement may be used to decide to which side of the TPC
   * the trigger should be assigned.
   */
  
  unsigned int const maxIrrelevant
    = minSumInOppositeWindows - (minSumInOppositeWindows / 2);
  
   return minInMainWindow > maxIrrelevant;
  
} // icarus::trigger::WindowPattern::isMainRequirementRelevant()
  
  
//------------------------------------------------------------------------------
bool icarus::trigger::WindowPattern::isSumRequirementRelevant() const {
  
  return minSumInOppositeWindows > (minInMainWindow + minInOppositeWindow);
  
} // icarus::trigger::WindowPattern::isSumRequirementRelevant()
  
  
//------------------------------------------------------------------------------
std::string icarus::trigger::WindowPattern::tag() const {
  using namespace std::string_literals;
  
  std::string s;
  
  // `M` is added last, because it may be needed even if irrelevant
  if (minInOppositeWindow > 0U)
    s += "O"s + std::to_string(minInOppositeWindow);
  
  if (isSumRequirementRelevant())
    s += "S"s + std::to_string(minSumInOppositeWindows);
  
  if ((minInDownstreamWindow > 0U) || requireDownstreamWindow) {
    s += "D"s + std::to_string(minInDownstreamWindow);
    if (requireDownstreamWindow) s+= "req"s;
  } // if downstream
  
  if ((minInUpstreamWindow > 0U) || requireUpstreamWindow) {
    s += "U"s + std::to_string(minInUpstreamWindow);
    if (requireUpstreamWindow) s+= "req"s;
  } // if upstream
  
  if (isMainRequirementRelevant() || s.empty())
    s = "M"s + std::to_string(minInMainWindow) + s;
  
  return s;
} // icarus::trigger::WindowPattern::description()


//------------------------------------------------------------------------------
std::string icarus::trigger::WindowPattern::description() const {
  using namespace std::string_literals;
  
  bool const useMain = isMainRequirementRelevant();
  bool const useSum = isSumRequirementRelevant();
  
  std::string s = "required:";
  
  if (useMain || !useSum)
    s += " "s + std::to_string(minInMainWindow);
  else
    s += " "s + std::to_string(minSumInOppositeWindows) + " (main+opposite)"s;
  
  if (minInOppositeWindow > 0U)
    s += " + "s + std::to_string(minInOppositeWindow) + " (opposite)"s;
  
  if (useMain && useSum)
    s += " (and "s + std::to_string(minSumInOppositeWindows) + " main+opposite)"s;
  
  if ((minInDownstreamWindow > 0U) || requireDownstreamWindow) {
    if (minInDownstreamWindow)
      s += " + "s + std::to_string(minInDownstreamWindow);
    s += " (downstream"s;
    if (requireDownstreamWindow) s+= ", mandatory)"s;
    s += ")"s;
  } // if downstream
  
  if ((minInUpstreamWindow > 0U) || requireUpstreamWindow) {
    if (minInUpstreamWindow) s += " + "s + std::to_string(minInUpstreamWindow);
    s += " (upstream"s;
    if (requireUpstreamWindow) s+= ", mandatory)"s;
    s += ")"s;
  } // if upstream
  
  return s;
} // icarus::trigger::WindowPattern::description()


//------------------------------------------------------------------------------
