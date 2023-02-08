/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WindowTopologyAlg.cxx
 * @brief  Assembles the topology of trigger windows (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 25, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/WindowTopologyAlg.h
 */


// library header
#include "icaruscode/PMT/Trigger/Algorithms/WindowTopologyAlg.h"

// ICARUS libraries
#include "icarusalg/Utilities/sortBy.h" // also icarus::util::sortCollBy()

// LArSoft libraries
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // MiddlePointAccumulator
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"

// C/C++ standard libraries
#include <algorithm> // std::move(), std::sort()
#include <array> // std::move()
#include <utility> // std::move()
#include <cassert>


//------------------------------------------------------------------------------
namespace {
  
  template <typename T>
  std::vector<T>& append(std::vector<T>& destColl, std::vector<T>&& srcColl) {
    
    if (destColl.empty()) destColl = std::move(srcColl);
    else {
      destColl.reserve(destColl.size() + srcColl.size());
      std::move(srcColl.begin(), srcColl.end(), back_inserter(destColl));
      srcColl.clear();
    }
    
    return destColl;
  } // append()
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::trigger::WindowTopologyAlg
//------------------------------------------------------------------------------
icarus::trigger::WindowTopologyAlg::WindowTopologyAlg(
  geo::GeometryCore const& geom,
  std::string const& logCategory /* = "WindowTopologyAlg" */
)
  : icarus::ns::util::mfLoggingClass(logCategory)
  , fGeom(&geom)
{}


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::createFromGates
  (TriggerGatesPerCryostat_t const& gates) const -> WindowChannelMap
{
  
  // store the window topology information here:
  std::vector<WindowChannelMap::WindowInfo_t> windows;
  
  for (auto const& [ cryoGates, cryo ]
         : util::zip(gates, fGeom->Iterate<geo::CryostatGeo>()))
  {
    
    append(
      windows,
      createWindowsFromCryostat
        (extractGateChannels(cryoGates), cryo, *fGeom, windows.size())
      );
    
  } // for cryostats
  
  return emplaceAndDumpMap(std::move(windows));
  
} // icarus::trigger::WindowTopologyAlg::createFromGates()


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::createFromGates
  (TriggerGates_t const& gates) const -> WindowChannelMap
{
  
  // split gates by cryostat
  TriggerGatesPerCryostat_t gatesByCryostat { fGeom->Ncryostats() };
  for (InputTriggerGate_t const& gate: gates) {
    
    geo::CryostatID cid; // invalid
    assert(!cid.isValid);
    
    for (raw::Channel_t const channel: gate.channels()) {
      geo::OpDetGeo const& opDet = fGeom->OpDetGeoFromOpChannel(channel);
      geo::OpDetID const oid = opDet.ID();
      if (!cid) cid = oid;
      else if (cid != oid) { // just in case
        throw cet::exception("WindowTopologyAlg")
          << "Input gate includes gates from different cryostats!!\n";
      }
    } // for channels in gate
    if (!cid) continue; // gate with no channels does not contribute
    
    gatesByCryostat.at(cid.Cryostat).push_back(gate); // (copy)
  } // for gates
  
  return createFromGates(gatesByCryostat);
  
} // icarus::trigger::WindowTopologyAlg::createFromGates()


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::createFromCryostatGates
  (TriggerGates_t const& windowChannels, geo::CryostatGeo const& cryo) const
  -> WindowChannelMap
{
  return emplaceAndDumpMap(
    createWindowsFromCryostat(extractGateChannels(windowChannels), cryo, *fGeom)
    );
} // icarus::trigger::WindowTopologyAlg::createFromCryostatGates()


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::createFromCryostatGates
  (TriggerGates_t const& windowChannels, geo::CryostatID cryoID) const
  -> WindowChannelMap
{
  return createFromCryostatGates(windowChannels, fGeom->Cryostat(cryoID));
} // icarus::trigger::WindowTopologyAlg::createFromCryostatGates()


//------------------------------------------------------------------------------
template <typename... Args>
auto icarus::trigger::WindowTopologyAlg::emplaceAndDumpMap(Args&&... args) const
  -> WindowChannelMap
{
  WindowChannelMap const map { std::forward<Args>(args)... };
  mfLogTrace() << "Window map: << " << map;
  return map;
} // icarus::trigger::WindowTopologyAlg::emplaceAndDumpMap()


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::createWindowsFromCryostat(
  WindowChannelColl_t const& windowChannels,
  geo::CryostatGeo const& cryo,
  geo::GeometryCore const& geom,
  std::size_t firstWindowIndex /* = 0U */
  ) -> std::vector<WindowChannelMap::WindowInfo_t>
{
  /*
   * 1.     fill the window information with local information
   * 2.     sort the windows in drift plane (first cryostat TPC as reference)
   * 3.     split the windows per plane
   * 4.     for each window plane:
   * 4.1.     sort windows on beam direction (TPC width direction)
   * 4.2.     fill the neighbour information on each window
   */
  std::vector<WindowChannelMap::WindowInfo_t> windows;
  windows.reserve(windowChannels.size());
  
  // use the first TPC of the cryostat as reference for directions
  assert(cryo.NTPC() > 0U);
  geo::TPCGeo const& refTPC = cryo.TPC(0U);
  
  //
  // 1. fill the window information with local information
  //
  using WindowInfoPtrs_t = std::vector<WindowChannelMap::WindowInfo_t*>;
  
  WindowInfoPtrs_t cryoWindowInfo;
  cryoWindowInfo.reserve(windowChannels.size());
  
  std::size_t iWindow = firstWindowIndex;
  for (auto const& channels: windowChannels) {
    
    WindowChannelMap::WindowInfo_t wInfo;
    WindowChannelMap::WindowComposition_t& wComp = wInfo.composition;
    
    wInfo.topology.index = iWindow++;
    wComp.channels = channels;
    std::sort(wComp.channels.begin(), wComp.channels.end());
    wComp.cryoid = channels.empty()
      ? geo::CryostatID{}
      : geom.OpDetGeoFromOpChannel(channels.front()).ID().asCryostatID()
      ;
    
    geo::vect::MiddlePointAccumulator middlePoint;
    for (raw::Channel_t const channel: channels) {
      // documentation of OpDetGeoFromOpChannel() does not say what on error...
      geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpChannel(channel);
      middlePoint.add(opDet.GetCenter());
      if (opDet.ID() != wComp.cryoid) wComp.cryoid = geo::CryostatID{};
    } // for channel
    wComp.center = middlePoint.middlePoint();
    
    windows.push_back(std::move(wInfo));
    cryoWindowInfo.push_back(&windows.back());
    
  } // for windows
  
  //
  // 2. sort the windows in drift plane (first cryostat TPC as reference)
  //
  auto const normalProjection = [&refTPC](auto const* info)
    { return refTPC.DistanceFromReferencePlane(info->composition.center); };
  WindowInfoPtrs_t const windowsByNormal
    = util::sortCollBy(cryoWindowInfo, normalProjection);
  
  //
  // 3. split the windows per plane
  //
  // split the list in two; there is a good deal of faith here
  auto const beamCoordinate
    = [&refPlane=refTPC.ReferencePlane()](auto const* info)
      {
        return refPlane.PointWidthDepthProjection(info->composition.center).X();
      }
    ;

  //
  // 4. for each plane:
  //
  // 4.1. sort windows on beam direction (TPC width direction)
  //
  auto const iMiddleWindow
    = std::next(windowsByNormal.cbegin(), windowsByNormal.size() / 2U);
  std::array<WindowInfoPtrs_t, 2U> const windowsByPlane = {
    util::sortBy(windowsByNormal.cbegin(), iMiddleWindow, beamCoordinate),
    util::sortBy(iMiddleWindow, windowsByNormal.cend(), beamCoordinate)
  };
  
  for (auto const& [ iPlane, planeWindows ]: util::enumerate(windowsByPlane)) {
    //
    // 4.2.     fill the neighbour information on each window
    //
    auto const& otherPlaneWindows
      = windowsByPlane.at(windowsByPlane.size() - 1U - iPlane);
    std::size_t const iLastPlaneWindow = planeWindows.size() - 1U;
    for (auto [ iPlaneWindow, windowInfo ]: util::enumerate(planeWindows)) {
      
      // assumes all topology information is InvalidWindowIndex by default
      
      if (WindowChannelMap::WindowInfo_t const* closestOppositeWindow
        = findClosestWindow(otherPlaneWindows, windowInfo)
      ) {
        windowInfo->topology.opposite = closestOppositeWindow->topology.index;
      }

      if (iPlaneWindow > 0U) {
        windowInfo->topology.upstream
          = planeWindows[iPlaneWindow - 1U]->topology.index;
      }
      
      if (iPlaneWindow < iLastPlaneWindow) {
        windowInfo->topology.downstream
          = planeWindows[iPlaneWindow + 1U]->topology.index;
      }
      
    } // for window in plane
  } // for planes
  
  return windows;
  
} // icarus::trigger::WindowTopologyAlg::createWindowsFromCryostat()


// -----------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::extractGateChannels
  (TriggerGates_t const& gates) -> WindowChannelColl_t
{
  
  WindowChannelColl_t windowChannels;
  windowChannels.reserve(gates.size());
  
  for (InputTriggerGate_t const& gate: gates) {
    auto const& gateChannels = gate.channels();
    windowChannels.emplace_back(gateChannels.begin(), gateChannels.end());
  } // for gates
  
  return windowChannels;
  
} // icarus::trigger::WindowTopologyAlg::extractGateChannels()


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::findClosestWindow(
  std::vector<WindowChannelMap::WindowInfo_t*> const& windowList,
  WindowChannelMap::WindowInfo_t const* targetWindow
  ) -> WindowChannelMap::WindowInfo_t const*
{
  
  if (!targetWindow || windowList.empty()) return nullptr;
  
  WindowChannelMap::WindowInfo_t const* closest = nullptr;
  double minDistance2 = std::numeric_limits<double>::max();
  for (auto const* window: windowList) {
    if (!window) continue;
    
    double const d2
      = (window->composition.center - targetWindow->composition.center).mag2();
    if (minDistance2 <= d2) continue;
    minDistance2 = d2;
    closest = window;
  } // for
  
  return closest;
} // icarus::trigger::WindowTopologyAlg::findClosestWindow()


// -----------------------------------------------------------------------------
//--- icarus::trigger::WindowTopologyVerification
//------------------------------------------------------------------------------
std::string icarus::trigger::WindowTopologyVerification::verify
  (TriggerGates_t const& gates) const
{
  /*
   * Verifies that the `gates` are in the expected order and have
   * the expected channel content.
   */

  if (!hasTopology()) {
    throw cet::exception("WindowTopologyVerification")
      << "verify() called without any window topology set to be verified.\n";
  }
  
  std::size_t iWindow = 0U;
  std::string errorMsg; // if this stays `empty()` there is no error
  for (auto const& gate: gates) {
    
    // if error message is `empty()` there is no error
    std::string const windowError = verifyGate(iWindow++, gate);
    
    if (!windowError.empty()) errorMsg += windowError + '\n';
    
  } // for gates
  
  return errorMsg;
  
} // icarus::trigger::WindowTopologyVerification::verify()


//------------------------------------------------------------------------------
std::string icarus::trigger::WindowTopologyVerification::verify
  (TriggerGatesPerCryostat_t const& gates) const
{
  /*
   * Verifies that the `gates` are in the expected order and have
   * the expected channel content.
   */

  if (!hasTopology()) {
    throw cet::exception("WindowTopologyVerification")
      << "verify() called without any window topology set to be verified.\n";
  }
  
  std::size_t iWindow = 0U;
  std::string errorMsg; // if this stays `empty()` there is no error
  for (auto const& cryoGates: gates) {
    for (auto const& gate: cryoGates) {
      
      // if error message is `empty()` there is no error
      std::string const windowError = verifyGate(iWindow++, gate);
      
      if (!windowError.empty()) errorMsg += windowError + '\n';
      
    } // for gates in cryostat
  } // for cryostats
  
  return errorMsg;
  
} // icarus::trigger::WindowTopologyVerification::verifyTopologicalMap()


//------------------------------------------------------------------------------
std::string icarus::trigger::WindowTopologyVerification::verifyGate
  (std::size_t iWindow, InputTriggerGate_t const& gate) const
{
  std::string errors; // if this stays `empty()` there is no error
  
  WindowChannelMap::WindowInfo_t const& windowInfo = fWindowMap->info(iWindow);
  
  auto const channelInWindow = [
      begin=windowInfo.composition.channels.cbegin(),
      end=windowInfo.composition.channels.cend()
    ](raw::Channel_t channel)
    { return std::binary_search(begin, end, channel); }
    ;
  
  for (raw::Channel_t const channel: gate.channels()) {
    if (channelInWindow(channel)) continue;
    if (errors.empty()) {
      errors =
        "channels not in window #" + std::to_string(windowInfo.topology.index)
        + ":";
    } // if first error
    errors += " " + std::to_string(channel);
  } // for all channels in gate
  
  return errors;
} // icarus::trigger::WindowTopologyVerification::verifyGate()


//------------------------------------------------------------------------------
