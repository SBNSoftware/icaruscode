/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.cxx
 * @brief  Algorithm composing PMT sliding windows from geometry information.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 26, 2021
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.h
 * 
 */


// library header
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.h"

// ICARUS libraries
#include "icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.h"

// LArSoft libraries
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <ostream>
#include <optional>
#include <cassert>


// -----------------------------------------------------------------------------
// ---  icarus::trigger::SlidingWindowDefinitionAlg implementation
// -----------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowDefinitionAlg::makeWindows
  (unsigned int windowSize, unsigned int windowStride) const
  -> WindowDefs_t
{
  /*
   * 1. compute the vertical PMT towers in each separate optical detector plane
   * 2. fill the windows by counting channels (i.e. op. det.)
   */
  using icarus::trigger::PMTverticalSlicingAlg;

  //
  // 1. compute the vertical PMT towers in each separate optical detector plane
  //
  PMTverticalSlicingAlg slicerAlg(fLogCategory);
  PMTverticalSlicingAlg::Slices_t slices;
  for (geo::CryostatGeo const& cryo: fGeom.Iterate<geo::CryostatGeo>())
    slicerAlg.appendCryoSlices(slices, cryo);

  //
  // 2. fill the windows by counting channels (i.e. optical detectors)
  //
  WindowDefs_t windows;

  for (PMTverticalSlicingAlg::PMTtowerOnPlane_t const& planeSlices: slices) {

    auto itSlice = planeSlices.begin();
    auto const send = planeSlices.end();
    while (itSlice != send) {

      mf::LogTrace(fLogCategory) << "Assembling window #" << windows.size();

      WindowDefs_t::value_type window;
      window.reserve(windowSize);

      std::optional<decltype(itSlice)> nextStart;
      unsigned int nChannels = 0U;
      while (nChannels < windowSize) {
        if (itSlice == send) break;

        // aside: check if this is the right place to start the next window
        if (nChannels == windowStride) {
          mf::LogTrace(fLogCategory)
            << "  (next window will start from this slice)";
          nextStart = itSlice;
        }
        else if ((nChannels > windowStride) && !nextStart) {
          throw cet::exception("SlidingWindowDefinitionAlg")
            << "Unable to start a new window " << windowStride
            << " channels after window #" << windows.size()
            << " (next slice starts " << nChannels << " channels after)\n";
        }

        mf::LogTrace(fLogCategory)
          << "  adding " << itSlice->size() << " channels to existing "
          << nChannels;
        for (geo::OpDetGeo const* opDet: *itSlice) {
          geo::OpDetID const& id = opDet->ID();
          raw::Channel_t const channel
            = fGeom.OpDetFromCryo(id.OpDet, id.Cryostat);
          mf::LogTrace(fLogCategory)
            << "   * " << id << " (channel " << channel << ")";
          window.push_back(channel);
        } // for channels in slice
        nChannels += (itSlice++)->size();
      } // while
      if (nChannels == windowStride) nextStart = itSlice;
      assert(nextStart);

      if (nChannels < windowSize) {
        if (!windows.empty()) {
          // if the last window can't be completed, it's still fine for us
          mf::LogTrace(fLogCategory)
            << "  ... couldn't complete the last " << windowSize
            << "-channel window: only " << nChannels << " channels collected";
          break;
        }
      } // if missing windows
      if (nChannels != windowSize) {
        throw cet::exception("SlidingWindowDefinitionAlg")
          << "Definition of one window yielded " << nChannels
          << " elements (window should be of size " << windowSize
          << " and with stride " << windowStride << ").\n";
      }

      windows.push_back(std::move(window));

      itSlice = nextStart.value();
    } // for all slices
  } // for all windows
  mf::LogTrace(fLogCategory)
    << "SlidingWindowTrigger defined " << windows.size() << " windows.";

  return windows;
} // icarus::trigger::SlidingWindowDefinitionAlg::makeWindows()


// -----------------------------------------------------------------------------
