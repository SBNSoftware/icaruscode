/** ****************************************************************************
 * @file ChannelROI.cxx
 * @brief Definition of basic channel signal object.
 * @author brebel@fnal.gov
 * @see  ChannelROI.h
 *
 * ****************************************************************************/

#include "icaruscode/IcarusObj/ChannelROI.h"

// C/C++ standard libraries
#include <utility> // std::move()

namespace recob{

  //----------------------------------------------------------------------
  ChannelROI::ChannelROI()
    : fChannel(raw::InvalidChannelID)
    , fSignalROI()
    {}

  //----------------------------------------------------------------------
  ChannelROI::ChannelROI(
    RegionsOfInterest_t const& sigROIlist,
    raw::ChannelID_t channel
    )
    : fChannel(channel)
    , fSignalROI(sigROIlist)
    {}

  //----------------------------------------------------------------------
  ChannelROI::ChannelROI(
    RegionsOfInterest_t&& sigROIlist,
    raw::ChannelID_t channel
    )
    : fChannel(channel)
    , fSignalROI(std::move(sigROIlist))
    {}


  //----------------------------------------------------------------------
  std::vector<short int> ChannelROI::Signal() const {
    return { fSignalROI.begin(), fSignalROI.end() };
  } // ChannelROI::Signal()


}
////////////////////////////////////////////////////////////////////////

