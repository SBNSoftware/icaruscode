/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.cxx
 * @brief  Applies sliding window trigger patterns.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.h
 */


// library header
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowPatternAlg.h"


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <algorithm> // std::binary_search()
#include <utility> // std::pair<>, std::move()
#include <cassert>


//------------------------------------------------------------------------------
icarus::trigger::SlidingWindowPatternAlg::SlidingWindowPatternAlg(
  WindowTopology_t windowTopology,
  WindowPattern_t windowPattern,
  icarus::trigger::ApplyBeamGateClass beamGate,
  std::string const& logCategory /* = "SlidingWindowPatternAlg" */
  )
  : icarus::ns::util::mfLoggingClass(logCategory)
  , fWindowTopology(std::move(windowTopology))
  , fWindowPattern(std::move(windowPattern))
  , fBeamGate(std::move(beamGate))
  {}


//------------------------------------------------------------------------------
icarus::trigger::SlidingWindowPatternAlg::SlidingWindowPatternAlg(
  WindowTopology_t windowTopology,
  WindowPattern_t windowPattern,
  std::string const& logCategory /* = "SlidingWindowPatternAlg" */
  )
  : icarus::ns::util::mfLoggingClass(logCategory)
  , fWindowTopology(std::move(windowTopology))
  , fWindowPattern(std::move(windowPattern))
  {}


//------------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowPatternAlg::simulateResponse
  (TriggerGates_t const& gates) const -> AllTriggerInfo_t
{
  
  // ensures input gates are in the same order as the configured windows
  verifyInputTopology(gates);
  
  auto const& inBeamGates = fBeamGate? fBeamGate->applyToAll(gates): gates;
  
  //
  // 2.   apply pattern:
  //
  std::size_t const nWindows = fWindowTopology.nWindows();
    
  //
  // 2.1.   for each main window, apply the pattern
  //
  WindowTriggerInfo_t triggerInfo; // start empty
  for (std::size_t const iWindow: util::counter(nWindows)) {
    
    TriggerInfo_t const windowResponse
      = applyWindowPattern(fWindowPattern, iWindow, inBeamGates);
    
    if (!windowResponse) continue;
    
    mfLogTrace()
      << "Pattern fired on window #" << iWindow
      << " at tick " << windowResponse.atTick()
      ;
    
    //
    // 2.2. pick the main window with the earliest successful response, if any;
    //      that defines location and time of the trigger
    //
    if (!triggerInfo || triggerInfo.info.atTick() > windowResponse.atTick())
      triggerInfo.emplace(iWindow, windowResponse);
    
  } // main window choice
  
  return { std::move(triggerInfo.info), MoreInfo_t{ triggerInfo.windowIndex } };
} // icarus::trigger::SlidingWindowPatternAlg::simulateResponse()


//------------------------------------------------------------------------------
bool icarus::trigger::SlidingWindowPatternAlg::hasBeamGate() const
  { return fBeamGate.has_value(); }


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowPatternAlg::setBeamGate
  (icarus::trigger::ApplyBeamGateClass beamGate)
  { fBeamGate.emplace(std::move(beamGate)); }


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowPatternAlg::clearBeamGate()
  { fBeamGate.reset(); }


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowPatternAlg::verifyInputTopology(TriggerGates_t const& gates) const {
  /*
   * Verifies that the `gates` are in the expected order and have
   * the expected channel content.
   */

  std::string errorMsg; // if this stays `empty()` there is no error
  for (auto const& [ iWindow, gate ]: util::enumerate(gates)) {
    
    std::string windowError; // if this stays `empty()` there is no error
    
    // more input gates than windows?
    if (iWindow >= fWindowTopology.nWindows()) {
      windowError = "unexpected input gate #" + std::to_string(iWindow) + " (";
      for (raw::Channel_t const channel: gate.channels())
        windowError += " " + std::to_string(channel);
      windowError += " )";
      errorMsg += windowError + '\n';
      continue; // we collect info of all spurious gates
    }
    
    WindowTopology_t::WindowInfo_t const& windowInfo
      = fWindowTopology.info(iWindow);
    
    auto const channelInWindow
      = [begin=windowInfo.channels.cbegin(),end=windowInfo.channels.cend()]
      (raw::Channel_t channel)
      { return std::binary_search(begin, end, channel); }
      ;
    
    for (raw::Channel_t const channel: gate.channels()) {
      if (channelInWindow(channel)) continue;
      if (windowError.empty()) {
        windowError =
          "channels not in window #" + std::to_string(windowInfo.index)
          + ":";
      } // if first error
      windowError += " " + std::to_string(channel);
    } // for all channels in gate
    
    if (!windowError.empty()) errorMsg += windowError + '\n';
  } // for gates

  // more input gates than windows?
  if (fWindowTopology.nWindows() > size(gates)) {
    errorMsg +=
      "Not enough input gates: " + std::to_string(size(gates)) + " gates for "
      + std::to_string(fWindowTopology.nWindows()) + " windows\n";
  }
  
  if (errorMsg.empty()) return;
  
  // put together the exception message and throw it.
  throw cet::exception("SlidingWindowPatternAlg")
    << "Some channels from trigger gates do not match"
    " the configured window allocation:\n"
    << errorMsg
    << "\n" // empty line
    << "Window configuration: "
    << fWindowTopology << "\n";
  
} // icarus::trigger::SlidingWindowPatternAlg::verifyInputTopology()


//------------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowPatternAlg::applyWindowPattern(
  WindowTopology_t::WindowInfo_t const& windowInfo,
  WindowPattern_t const& pattern,
  TriggerGates_t const& gates
  ) -> TriggerInfo_t
{
  
  /*
   * 1. check that the pattern can be applied; if not, return no trigger
   * 2. discriminate all the relevant gates against their required minimum count
   * 3. combine them in AND
   * 4. find the trigger time, fill the trigger information accordingly
   */
  TriggerInfo_t res; // no trigger by default
  assert(!res);

  //
  // 1. check that the pattern can be applied; if not, return no trigger
  //
  
  // check that the pattern centered into iWindow has all it needs:
  if (pattern.requireUpstreamWindow && !windowInfo.hasUpstreamWindow())
    return res;
  if (pattern.requireDownstreamWindow && !windowInfo.hasDownstreamWindow())
    return res;
  
  
  //
  // 2. discriminate all the relevant gates against their required minimum count
  // 3. combine them in AND
  //
  
  // main window
  TriggerGateData_t trigPrimitive
    = discriminate(gates[windowInfo.index], pattern.minInMainWindow);
  
  // add opposite window requirement (if any)
  if ((pattern.minInOppositeWindow > 0U) && windowInfo.hasOppositeWindow()) {
    trigPrimitive.Mul
      (discriminate(gates[windowInfo.opposite], pattern.minInOppositeWindow));
  } // if
  
  // add upstream window requirement (if any)
  if ((pattern.minInUpstreamWindow > 0U) && windowInfo.hasUpstreamWindow()) {
    trigPrimitive.Mul
      (discriminate(gates[windowInfo.upstream], pattern.minInUpstreamWindow));
  } // if
  
  // add downstream window requirement (if any)
  if ((pattern.minInDownstreamWindow > 0U) && windowInfo.hasDownstreamWindow())
  {
    trigPrimitive.Mul(
      discriminate(gates[windowInfo.downstream], pattern.minInDownstreamWindow)
      );
  } // if
  
  //
  // 4. find the trigger time, fill the trigger information accordingly
  //
  auto const trigTick = trigPrimitive.findOpen(); // first trigger
  if (trigTick != trigPrimitive.MaxTick) {
    res.emplace(detinfo::timescales::optical_tick{ trigTick });
    assert(res);
  }
  
  return res;
  
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::applyWindowPattern()


//------------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowPatternAlg::applyWindowPattern(
  WindowPattern_t const& pattern,
  std::size_t iWindow,
  TriggerGates_t const& gates
  ) const -> TriggerInfo_t
{
  WindowTopology_t::WindowInfo_t const& windowInfo
    = fWindowTopology.info(iWindow);
  assert(windowInfo.index == iWindow);
  
  return applyWindowPattern(windowInfo, pattern, gates);
} // icarus::trigger::SlidingWindowTriggerEfficiencyPlots::applyWindowPattern()


//------------------------------------------------------------------------------
