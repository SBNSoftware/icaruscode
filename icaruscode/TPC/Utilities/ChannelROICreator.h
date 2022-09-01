/** ****************************************************************************
 * @file   ChannelROICreator.h
 * @brief  Helper functions to create a wire
 * @date   December 11, 2014
 * @author petrillo@fnal.gov
 * @see    Wire.h ChannelROICreator.cxx
 *
 * ****************************************************************************/

#ifndef ChannelROICreator_H
#define ChannelROICreator_H

// C/C++ standard library
#include <utility> // std::move()

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "icaruscode/IcarusObj/ChannelROI.h"

namespace raw { class RawDigit; }

/// Reconstruction base classes
namespace recob {

  /**
   * @brief Class managing the creation of a new recob::Wire object
   *
   * In order to be as simple as possible (Plain Old Data), data products like
   * recob::Wire need to be stripped of most of their functions, including the
   * ability to communicate whether a value we try to store is invalid
   * (that would require a art::Exception -- art -- or at least a message on the
   * screen -- MessageFacility) and the ability to read things from event,
   * services (e.g. geometry) etc.
   *
   * A Creator is a class that creates a temporary data product, and at the
   * end it yields it to the caller for storage.
   * This last step should be by move construction, although a copy method is
   * also provided.
   *
   * An example of creating a Wire object:
   *
   *     // let RoIsignal be a recob::Wire::RegionsOfInterest_t already filled
   *     // with the signal regions, and rawdigit the raw::RawDigit of the
   *     // channel; RoIsignal will become empty
   *     recob::ChannelROICreator wire(std::move(RoIsignal), rawdigit);
   *     wires.push_back(wire.move()); // wire content is not valid any more
   *
   * This is a one-step creation object: the wire is constructed at the same
   * time the ChannelROICreator is, and no facility is offered to modify the
   * constructed wire, or to create another one.
   */
  class ChannelROICreator {
    public:
      /// Alias for the type of regions of interest
      using RegionsOfInterest_t = ChannelROI::RegionsOfInterest_t;

      // destructor, copy and move constructor and assignment as default

      /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param rawdigit the raw digit this channel is associated to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       */
      ChannelROICreator
        (const RegionsOfInterest_t& sigROIlist, const raw::RawDigit& rawdigit);


      /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param rawdigit the raw digit this channel is associated to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       *
       * Signal information is moved from sigROIlist, that becomes empty.
       */
      ChannelROICreator
        (RegionsOfInterest_t&& sigROIlist, const raw::RawDigit& rawdigit);


      /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param channel the ID of the channel
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       */
      ChannelROICreator(
        RegionsOfInterest_t const& sigROIlist,
        raw::ChannelID_t channel
        );


      /**
       * @brief Constructor: uses specified signal in regions of interest
       * @param sigROIlist signal organized in regions of interest
       * @param channel the ID of the channel
       * @param view the view the channel belongs to
       *
       * The information used from the raw digit are the channel ID and the
       * length in samples (TDC ticks) of the original readout window.
       *
       * Signal information is moved from sigROIlist, that becomes empty.
       */
      ChannelROICreator(
        RegionsOfInterest_t&& sigROIlist,
        raw::ChannelID_t channel
        );

      /**
       * @brief Prepares the constructed wire to be moved away
       * @return a right-value reference to the constructed wire
       *
       * Despite the name, no move happens in this function.
       * Move takes place in the caller code as proper; for example:
       *
       *     // be wire a ChannelROICreator instance:
       *     std::vector<recob::Wire> Wires;
       *     wire.move();                          // nothing happens
       *     Wires.push_back(wire.move());         // here the copy happens
       *     recob::Wire single_wire(wire.move()); // wrong! wire is empty now
       *
       */
      ChannelROI&& move() { return std::move(channelROI); }


      /**
       * @brief Returns the constructed wire
       * @return a constant reference to the constructed wire
       *
       * Despite the name, no copy happens in this function.
       * Copy takes place in the caller code as proper; for example:
       *
       *     // be wire a ChannelROICreator instance:
       *     std::vector<recob::Wire> Wires;
       *     wire.copy();                          // nothing happens
       *     Wires.push_back(wire.copy());         // here a copy happens
       *     recob::Wire single_wire(wire.copy()); // wire is copied again
       *
       */
      const ChannelROI& copy() const { return channelROI; }

    protected:

      ChannelROI channelROI; ///< local instance of the wire being constructed

  }; // class ChannelROICreator

} // namespace recob

#endif // ChannelROICreator_H
