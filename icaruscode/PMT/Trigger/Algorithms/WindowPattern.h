/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WindowPattern.h
 * @brief  Defines a (sliding) window trigger pattern.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 21, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/WindowPattern.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWPATTERN_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWPATTERN_H


// C/C++ standard libraries
#include <vector>
#include <string>
#include <limits>


// -----------------------------------------------------------------------------
namespace icarus::trigger {
  
  struct WindowPattern;
  
  /// A list of window patterns.
  using WindowPatterns_t = std::vector<WindowPattern>;
  
  //----------------------------------------------------------------------------
  std::string to_string(WindowPattern const& pattern);

  //----------------------------------------------------------------------------

} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Specification of the requirement of sliding window firing pattern.
 *
 * This structure contains the requirements of a trigger on a ("main") window.
 * The main window is supposed to have neighbors upstream and downstream of it
 * (the "stream" being the neutrino beam's), and a single window opposite to it.
 * The requirements are a minimum activity (as a count of LVDS signals open at
 * the same time) for each of the windows. A minimum requirement of `0` is
 * considered to be always satisfied (even if effectively the window it refers
 * to does not exist).
 * 
 * Special requirements are whether it is mandatory for the window to have a
 * upstream and/or downstream window (the presence of an opposite window is
 * assumed to always be optional).
 */
struct icarus::trigger::WindowPattern {
  
  /// @name Minimum required number of open trigger primitives per window.
  /// @{
  unsigned int minInMainWindow = 0U;
  unsigned int minInUpstreamWindow = 0U;
  unsigned int minInDownstreamWindow = 0U;
  unsigned int minInOppositeWindow = 0U;
  unsigned int minSumInOppositeWindows = 0U;
  /// @}
  
  /// Whether a window location with no upstream window should be discarded.
  bool requireUpstreamWindow = false;
  
  /// Whether a window location with no downstream window should be discarded.
  bool requireDownstreamWindow = false;
  
  /// Returns whether the main requirement (`M`) contributes to specification.
  /// This is not the case when `S` is specified that is twice `M` or when its
  /// value is `0`.
  bool isMainRequirementRelevant() const;
  
  /// Returns whether the sum requirement (`S`) contributes to specification.
  /// This is not the case when `S` is not larger than the sum of `M` and `O`.
  bool isSumRequirementRelevant() const;
  
  
  /**
   * @brief Returns a tag summarizing the pattern.
   * 
   * The tag encodes the requirements on the main window (_R(M)_), its
   * opposite window (_R(O)_), their minimum sum (_R(S)_) and downstream
   * (_R(D)_) and upstream (_R(U)_) windows.
   * A requirement _R(X)_ is in the format `X##[req]`, where `X` is the tag
   * letter of the requirement, `##` is the requirement level for that window,
   * and the optional `req` tag means that if for a main window this window
   * does not exist, that main window is not considered (e.g. the downstream
   * window of a main window which is the most downstream in the detector).
   * 
   * For example, `M5O2D2reqU1` requires 5 openings in the main window (`M5`),
   * 2 in the window opposite to the main one (`O2`) and also 2 on the window
   * downstream of the main one (`D2req`) and also 1 on the window
   * upstream of the main one (`U1`); in addition, if the main window
   * has no downstream window (i.e. it's at the "far end" of the detector),
   * the downstream requirement is never satisfied and the trigger is
   * considered to never fire. Instead, if there is no upstream window (i.e.
   * the main window is in the "near end" of the detector) the upstream window
   * requirement is considered to be satisfied (or ignored).
   * 
   * Redundant requirements are omitted (except for `M0` if not superseded by a
   * sum requirement).
   */
  std::string tag() const;
  
  /// Returns a description of the pattern.
  std::string description() const;
  
  
}; // icarus::triggerWindowPattern


//------------------------------------------------------------------------------
//---  Inline definitions
//------------------------------------------------------------------------------
inline std::string icarus::trigger::to_string(WindowPattern const& pattern)
  { return pattern.tag(); }


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWPATTERN_H
