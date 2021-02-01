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
{
  
} // icarus::trigger::WindowTopologyAlg::WindowTopologyAlg()


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::createFromGates
  (TriggerGatesPerCryostat_t const& gates) const -> WindowChannelMap
{
  
  // store the window topology information here:
  std::vector<WindowChannelMap::WindowInfo_t> windows;
  
  for (auto const& [ cryoGates, cryo ]: util::zip(gates, fGeom->IterateCryostats())) {
    
    append(
      windows,
      createWindowsFromCryostat
        (extractGateChannels(cryoGates), cryo, *fGeom, windows.size())
      );
    
  } // for cryostats
  
  return WindowChannelMap{ std::move(windows) };
  
} // icarus::trigger::WindowTopologyAlg::createFromGates()


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::createFromCryostatGates
  (TriggerGates_t const& windowChannels, geo::CryostatGeo const& cryo) const
  -> WindowChannelMap
{
  return WindowChannelMap{
    createWindowsFromCryostat(extractGateChannels(windowChannels), cryo, *fGeom)
    };
} // icarus::trigger::WindowTopologyAlg::createFromVolume()


//------------------------------------------------------------------------------
auto icarus::trigger::WindowTopologyAlg::createFromCryostatGates
  (TriggerGates_t const& windowChannels, geo::CryostatID cryoID) const
  -> WindowChannelMap
{
  return createFromCryostatGates(windowChannels, fGeom->Cryostat(cryoID));
} // icarus::trigger::WindowTopologyAlg::createFromVolume()


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
    
    wInfo.index = iWindow++;
    wInfo.channels = channels;
    std::sort(wInfo.channels.begin(), wInfo.channels.end());
    wInfo.cryoid = channels.empty()
      ? geo::CryostatID{}
      : geom.OpDetGeoFromOpChannel(channels.front()).ID().asCryostatID()
      ;
    
    geo::vect::MiddlePointAccumulator middlePoint;
    for (raw::Channel_t const channel: channels) {
      // documentation of OpDetGeoFromOpChannel() does not say what on error...
      geo::OpDetGeo const& opDet = geom.OpDetGeoFromOpChannel(channel);
      middlePoint.add(opDet.GetCenter());
      if (opDet.ID() != wInfo.cryoid) wInfo.cryoid = geo::CryostatID{};
    } // for channel
    wInfo.center = middlePoint.middlePoint();
    
    windows.push_back(std::move(wInfo));
    cryoWindowInfo.push_back(&windows.back());
    
  } // for windows
  
  //
  // 2. sort the windows in drift plane (first cryostat TPC as reference)
  //
  auto const normalProjection = [&refTPC](auto const* info)
    { return refTPC.DistanceFromReferencePlane(info->center); };
  WindowInfoPtrs_t const windowsByNormal
    = util::sortCollBy(cryoWindowInfo, normalProjection);
  
  //
  // 3. split the windows per plane
  //
  // split the list in two; there is a good deal of faith here
  auto const beamCoordinate
    = [&refPlane=refTPC.ReferencePlane()](auto const* info)
      { return refPlane.PointWidthDepthProjection(info->center).X(); }
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
      windowInfo->opposite = otherPlaneWindows[iPlaneWindow]->index;
      
      if (iPlaneWindow > 0U)
        windowInfo->upstream = planeWindows[iPlaneWindow - 1U]->index;
      
      if (iPlaneWindow < iLastPlaneWindow)
        windowInfo->downstream = planeWindows[iPlaneWindow + 1U]->index;
      
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
//--- icarus::trigger::WindowTopologyVerification
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
      
      std::string windowError; // if this stays `empty()` there is no error
      
      WindowChannelMap::WindowInfo_t const& windowInfo
        = fWindowMap->info(iWindow++);
      
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
    } // for gates in cryostat
  } // for cryostats
  
  return errorMsg;
  
} // icarus::trigger::WindowTopologyVerification::verifyTopologicalMap()


// -----------------------------------------------------------------------------
void icarus::trigger::WindowTopologyVerification::operator()
  (TriggerGatesPerCryostat_t const& gates) const
{
  
  std::string const errorMsg = verify(gates);
  if (errorMsg.empty()) return;
  
  // put together the exception message and throw it.
  throw cet::exception("WindowTopologyVerification")
    << "Some channels from trigger gates do not match the previous window allocation:\n"
    << errorMsg
    << "\n"
    << "Window allocation: "
    << fWindowMap.value()
    ;
  
} // icarus::trigger::WindowTopologyVerification::operator()


//------------------------------------------------------------------------------
//--- icarus::trigger::WindowTopologyManager
//------------------------------------------------------------------------------
std::string icarus::trigger::WindowTopologyManager::setOrVerify
  (TriggerGatesPerCryostat_t const& gates)
{
  if (fVerify.hasTopology()) {
    return fVerify.verify(gates);
  }
  else {
    extractTopology(gates);
    return {};
  }
} // icarus::trigger::WindowTopologyManager::setOrVerify()


//------------------------------------------------------------------------------
bool icarus::trigger::WindowTopologyManager::operator()
  (TriggerGatesPerCryostat_t const& gates)
{
  if (fVerify.hasTopology()) {
    fVerify(gates);
    return false;
  }
  else {
    extractTopology(gates);
    return true;
  }
} // icarus::trigger::WindowTopologyManager::operator()


//------------------------------------------------------------------------------
void icarus::trigger::WindowTopologyManager::extractTopology
  (TriggerGatesPerCryostat_t const& gates)
{
  icarus::trigger::WindowTopologyAlg const topoMaker
    { *fGeom, logCategory() + ":Extractor" };
  fVerify.setTopology(topoMaker.createFromGates(gates));
} // icarus::trigger::WindowTopologyManager::extractTopology()


//------------------------------------------------------------------------------
