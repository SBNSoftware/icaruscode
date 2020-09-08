/**
 * @file   icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h
 * @brief  Utilities for the conversion of trigger gate data formats.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERDATAUTILS_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERDATAUTILS_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/Utilities/DataProductPointerMap.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/Exception.h"

// C/C++ standard libraries
#include <map>
#include <vector>
#include <tuple>
#include <string> // std::to_string()
#include <stdexcept> // std::out_of_range
#include <utility> // std::move()


namespace icarus::trigger {

  // ---------------------------------------------------------------------------

  /// Map `util::DataProductPointerMap_t` for `raw::OpDetWaveform` objects.
  using OpDetWaveformDataProductMap_t
    = util::DataProductPointerMap_t<raw::OpDetWaveform>;


  // ---------------------------------------------------------------------------
  /**
   * @brief Returns the trigger gates in serializable format.
   * @tparam Gates type of the source of trigger gate data
   * @param gates the data to be reformatted (*data will be stolen!*)
   * @return gate data product
   *
   * This function transfers the data from the original structured `gates` into
   * a data collection suitable for serialization with _art_, but
   * *not including the association of each gate with all its contributing
   * waveforms*.
   * It must be stressed that this causes information loss, because the trigger
   * gates will not be associated any more not only to the contributing
   * waveforms, but also to the number of the optical detector channel(s).
   * This happens because the trigger gates do not store that information,
   * which is instead conveyed by the associated waveforms.
   *
   * The return value is a collection of trigger gate data
   * (`icarus::trigger::OpticalTriggerGateData_t`), with data _stolen_
   * from `gates`.
   *
   * The trigger gates are processed in the same order as they are in `gates`.
   *
   * After the function returns, `gates` will have been depleted of all the gate
   * data; the waveform information will be still associated to each gate, whose
   * gate data will be in an invalid state anyway, only good for destruction.
   *
   * The object `gates` is a iterable collection with each element required to
   * respond to:
   * * `waveforms()` returning a collection of pointers to the waveforms
   *   associated to the gate;
   * * `gateLevels()` returning a moveable container with the data of the
   *   trigger gate.
   * Collections of objects of type `icarus::trigger::OpticalTriggerGate` or
   * derived from it (`icarus::trigger::SingleChannelOpticalTriggerGate`,
   * `icarus::trigger::MultiChannelOpticalTriggerGate`) satisfy this
   * requirement.
   *
   * In the following example, we start with trigger gates already serialized
   * in a _art_ event. The serialization splits a trigger gate object in two
   * components: the gate levels, and the associated waveforms.
   * In the first part of the example we recover the information from the event
   * and we assemble it into the standard trigger gate objects (of type
   * `icarus::trigger::SingleChannelOpticalTriggerGate`).
   * After some unspecified and optional processing, `gates` are disassembled
   * to be saved into the event: this is achieved by a call to
   * `transformIntoOpticalTriggerGate()` which produces the trigger gate data
   * and their associations to the waveforms.
   * In the last part, these components are both stored into the event.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using icarus::trigger::OpticalTriggerGateData_t; // for convenience
   * 
   * //
   * // somehow get/build a collection of trigger gates;
   * // here we use `FillTriggerGates()` to read existing data from the event
   * //
   * auto gates = icarus::trigger::FillTriggerGates<icarus::trigger::SingleChannelOpticalTriggerGate>
   *   (
   *   *(event.getValidHandle<std::vector<OpticalTriggerGateData_t>>("orig")),
   *   *(event.getValidHandle<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>("orig"))
   *   );
   * 
   * 
   * // ...
   * 
   * //
   * // use the created vector and associations (e.g. put them into art event)
   * //
   * event.put(std::make_unique<std::vector<TriggerGateData_t>>
   *   (icarus::trigger::transformIntoOpticalTriggerGate(std::move(gates))));
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Also note that in the omitted processing part the trigger gates might have
   * been combined into a different type of gates instead, e.g. into a
   * collection of `icarus::trigger::MultiChannelOpticalTriggerGate`.
   * To serialize that collection, the exact same procedure would be used.
   * Once more, remember that the modules reading those newly put trigger gates
   * will now know the channel or channels associated with them.
   */
  template <typename Gates>
  std::vector<icarus::trigger::OpticalTriggerGateData_t>
  transformIntoOpticalTriggerGate(Gates&& gates);


  // ---------------------------------------------------------------------------
  /**
   * @brief Returns the trigger gates in serializable format.
   * @tparam Gates type of the source of trigger gate data
   * @param gates the data to be reformatted (*data will be stolen!*)
   * @param makeGatePtr _art_ pointer maker for the gate data
   * @param opDetWavePtrs map of art pointers to optical waveforms
   * @return a pair: gate data product and associations to optical waveforms
   *
   * This function transfers the data from the original structured `gates` into
   * a data collection suitable for serialization with _art_, including the
   * association of each gate with all its contributing waveforms.
   *
   * The return value is a tuple of two elements:
   * * `0`: collection of trigger gate data
   *     (`icarus::trigger::OpticalTriggerGateData_t`), with data _stolen_
   *     from `gates`;
   * * `1`: association between trigger gate data and their optical waveforms.
   *
   * The trigger gates are processed in the same order as they are in `gates`,
   * and the associations to the waveforms are set gate by gate, in the same
   * order as they are reported by the `waveforms()` method of the input gate.
   *
   * After the function returns, `gates` will have been depleted of all the data
   * and left in an undefined state only good for destruction.
   *
   * The object `gates` is a iterable collection with each element required to
   * respond to:
   * * `waveforms()` returning a collection of pointers to the waveforms
   *   associated to the gate;
   * * `gateLevels()` returning a moveable container with the data of the
   *   trigger gate.
   * Collections of objects of type `icarus::trigger::OpticalTriggerGate` or
   * derived from it (`icarus::trigger::SingleChannelOpticalTriggerGate`,
   * `icarus::trigger::MultiChannelOpticalTriggerGate`) satisfy this
   * requirement.
   *
   *
   * In the following example, we start with trigger gates already serialized
   * in a _art_ event. The serialization splits a trigger gate object in two
   * components: the gate levels, and the associated waveforms.
   * In the first part of the example we recover the information from the event
   * and we assemble it into the standard trigger gate objects (of type
   * `icarus::trigger::SingleChannelOpticalTriggerGate`).
   * After some unspecified and optional processing, `gates` are disassembled
   * to be saved into the event: this is achieved by a call to
   * `transformIntoOpticalTriggerGate()` which produces the trigger gate data
   * and their associations to the waveforms.
   * In the last part, these components are both stored into the event.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using icarus::trigger::OpticalTriggerGateData_t; // for convenience
   * 
   * //
   * // somehow get/build a collection of trigger gates;
   * // here we use `FillTriggerGates()` to read existing data from the event
   * //
   * auto gates = icarus::trigger::FillTriggerGates<icarus::trigger::SingleChannelOpticalTriggerGate>
   *   (
   *   *(event.getValidHandle<std::vector<OpticalTriggerGateData_t>>("orig")),
   *   *(event.getValidHandle<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>("orig"))
   *   );
   * 
   * 
   * // ...
   * 
   * //
   * // optical waveform to pointer map is required to create associations
   * // between the trigger gates and their waveforms
   * //
   * auto const& opDetWavePtrs = util::mapDataProductPointers
   *   (event, event.getValidHandle<std::vector<raw::OpDetWaveform>>("opdaq"));
   * // transform the data; after this line, `gates` is not usable any more
   * auto thresholdData = icarus::trigger::transformIntoOpticalTriggerGate
   *   (std::move(gates), makeGatePtr, opDetWavePtrs);
   * 
   * //
   * // use the created vector and associations (e.g. put them into art event)
   * //
   * event.put(
   *   std::make_unique<std::vector<TriggerGateData_t>>
   *     (std::move(std::get<0U>(thresholdData)))
   *   );
   * event.put(
   *   std::make_unique<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
   *     (std::move(std::get<1U>(thresholdData)))
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Also note that in the omitted processing part the trigger gates might have
   * been combined into a different type of gates instead, e.g. into a
   * collection of `icarus::trigger::MultiChannelOpticalTriggerGate`.
   * To serialize that collection, the exact same procedure would be used,
   * with the different outcome that now each trigger gate may be associated
   * to waveforms from different optical detector channels.
   *
   */
  template <typename Gates>
  std::tuple<
    std::vector<icarus::trigger::OpticalTriggerGateData_t>,
    art::Assns<icarus::trigger::OpticalTriggerGateData_t, raw::OpDetWaveform>
    >
  transformIntoOpticalTriggerGate(
    Gates&& gates,
    art::PtrMaker<icarus::trigger::OpticalTriggerGateData_t> const& makeGatePtr,
    OpDetWaveformDataProductMap_t const& opDetWavePtrs
    );


  // ---------------------------------------------------------------------------
  /**
   * @brief Creates a gate object out of trigger gate data products.
   * @tparam GateObject type of gate object to create
   * @tparam TriggerGateData the type of data stored
   * @param gates collection of data from all the gates
   * @param gateToWaveforms association between each gate and its waveforms
   * @return a STL vector of `GateObject` with copy of data from `gates`
   *
   * Objects like `icarus::trigger::OpticalTriggerGate` are complex enough that
   * they are not saved directly into an _art_ event. Rather, they are diced
   * into pieces and the pieces are stored.
   * This function stitches the pieces and returns back an object like
   * `icarus::trigger::OpticalTriggerGate`.
   *
   * The supported types for `GateObject` need to expose a
   * `icarus::trigger::OpticalTriggerGate`-like interface (including also e.g.
   * `icarus::trigger::SingleChannelOpticalTriggerGate` and
   * `icarus::trigger::MultiChannelOpticalTriggerGate`).
   * They must accept:
   * * to be constructed with a `raw::OpDetWaveform` reference;
   * * to be added a `raw::OpDetWaveform` reference via `add()` call;
   * * to be assigned a `TriggerGateData` type.
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * icarus::trigger::TriggerGateBuilder::TriggerGates readTriggerGates(
   *   art::Event const& event,
   *   icarus::trigger::ADCCounts_t const threshold,
   *   art::InputTag const& dataTag
   * ) {
   *
   *   auto const& gates
   *     = *(event.getValidHandle<std::vector<TriggerGateData_t>>(dataTag));
   *   auto const& gateToWaveforms = *(
   *     event.getValidHandle<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
   *       (dataTag)
   *     );
   *   return {
   *     threshold,
   *     icarus::trigger::FillTriggerGates
   *       <icarus::trigger::SingleChannelOpticalTriggerGate>
   *       (gates, gateToWaveforms)
   *     };
   *
   * } // readTriggerGates()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will return a `icarus::trigger::TriggerGateBuilder::TriggerGates` object
   * for the specified threshold, reading the information from `dateTag` data
   * products in an _art_ `event`.
   *
   *
   * Return value
   * -------------
   *
   * The returned collection contains one element for each `TriggerGateData`
   * object in `gates`, in the same order. Each of these elements is of type
   * `GateObject`, contains a copy of the data of the corresponding gate,
   * and a list of optical waveforms (`raw::OpDetWaveform` objects) it is
   * associated to.
   *
   *
   * Requirements
   * -------------
   *
   * The requirements bind the gates to their association to waveforms:
   * * each gate must be associated to at least one waveform
   * * the associations must be grouped so that all the association pairs
   *   pertaining a gate are contiguous
   *     * within each of these groups, which is made of at least one waveform,
   *       the waveforms must be ordered by increasing timestamp
   *     * the groups must be in the same order as their associated gates
   *   This constitutes the requirement of
   *   @ref LArSoftProxyDefinitionOneToManySeqAssn "one-to-many sequential association"
   *   with the addition that each element in `gates` must have at least one
   *   associated waveform.
   *
   */
  template <typename GateObject, typename TriggerGateData>
  std::vector<GateObject> FillTriggerGates(
    std::vector<TriggerGateData> const& gates,
    art::Assns<TriggerGateData, raw::OpDetWaveform> const& gateToWaveforms
    );


  // ---------------------------------------------------------------------------
  /// Associates each optical detector channel to a gate.
  template <typename GateObject>
  class TriggerGateIndex {

      public:

    using TriggerGate_t = GateObject; ///< Type of gate the index lists.

    /// Constructs an empty index expecting to be filled by `expectedChannels`
    /// registered channels.
    TriggerGateIndex(std::size_t expectedChannels = 0U);

    /// Initializes the index from all the `gates` in the specified collection.
    TriggerGateIndex(std::vector<TriggerGate_t> const& gates);

    /// Adds a new gate to the index.
    /// @throws cet::exception (`"TriggerGateIndex"`) if it conflicts with a
    ///         gate already added
    void add(TriggerGate_t const& gate);

    /// Returns the total number of registered channels.
    unsigned int nChannels() const;

    /// Returns the gate corresponding to the specified `channel`.
    TriggerGate_t const* find(raw::Channel_t const channel) const;

    /// Returns the gate corresponding to the specified `channel`.
    /// @throws std::out_of_range if the channel is not registered
    TriggerGate_t const& operator[](raw::Channel_t const channel) const;

      private:

    /// Index of gates by channel number (the same gate may appear many times).
    std::vector<TriggerGate_t const*> fGates;

    /// Extends internal map with `nullptr` to hold at least `chIndex`.
    /// @return whether an extension was performed
    bool expandToHold(std::size_t chIndex);

    /// Converts a channel number into an index.
    static std::size_t channelIndex(raw::Channel_t channel);

  }; // class TriggerGateIndex


  // ---------------------------------------------------------------------------

} // namespace icarus::trigger


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename Gates>
std::vector<icarus::trigger::OpticalTriggerGateData_t>
icarus::trigger::transformIntoOpticalTriggerGate(Gates&& gates) {
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGateData_t;

  //
  // create objects for the data products
  //
  std::vector<TriggerGateData_t> gateData;

  for (auto& gate: gates) {
    assert(gate.hasChannels());
    
    // we steal the data from the gate
    gateData.push_back(std::move(gate.gateLevels()));
    
    // now we add all channels;
    // the cast is a cheat for objects like `SingleChannelOpticalTriggerGate`
    // which hide the `channels()` method
    auto& newGate = gateData.back();
    for (auto channel: static_cast<TriggerGateData_t const&>(gate).channels())
      newGate.addChannel(channel);

  } // for all channels

  return gateData;

} // icarus::trigger::transformIntoOpticalTriggerGate()


// -----------------------------------------------------------------------------
template <typename Gates>
std::tuple<
  std::vector<icarus::trigger::OpticalTriggerGateData_t>,
  art::Assns<icarus::trigger::OpticalTriggerGateData_t, raw::OpDetWaveform>
  >
icarus::trigger::transformIntoOpticalTriggerGate(
  Gates&& gates,
  art::PtrMaker<icarus::trigger::OpticalTriggerGateData_t> const& makeGatePtr,
  OpDetWaveformDataProductMap_t const& opDetWavePtrs
  )
{
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGateData_t;

  //
  // create objects for the data products
  //
  std::vector<TriggerGateData_t> gateData
    = transformIntoOpticalTriggerGate(std::move(gates));
  
  // note that we rely on the fact that `gates` still contains trigger gates,
  // and those trigger gates still know about the associated waveforms;
  // attempting to access their levels, though, would be an error
  // (if really needed, we can still find them in `gateData`)
  art::Assns<TriggerGateData_t, raw::OpDetWaveform> gateToWaveforms;

  for (auto const& [ iGate, gate ]: util::enumerate(gates)) {

    if (gate.waveforms().empty()) continue;

    //
    // produce the associations
    //

    // pointer to the gate data we have just added:
    art::Ptr<TriggerGateData_t> const& gatePtr = makeGatePtr(iGate);
    for (raw::OpDetWaveform const* waveform: gate.waveforms()) {
      
      gateToWaveforms.addSingle(gatePtr, opDetWavePtrs.at(waveform));

    } // for waveforms

  } // for all channels

  return { std::move(gateData), std::move(gateToWaveforms) };

} // icarus::trigger::transformIntoOpticalTriggerGate()


// -----------------------------------------------------------------------------
template <typename GateObject, typename TriggerGateData>
std::vector<GateObject> icarus::trigger::FillTriggerGates(
  std::vector<TriggerGateData> const& gates,
  art::Assns<TriggerGateData, raw::OpDetWaveform> const& gateToWaveforms
) {
  using GateData_t = GateObject; ///< Type of object being returned,

  std::vector<GateData_t> allGates;

  auto iGateToWaveform = gateToWaveforms.begin();
  auto const gwend = gateToWaveforms.end();

  for (auto const& [ iGate, gate ]: util::enumerate(gates)) {

    //
    // find the relevant waveforms for this gate;
    // match by index of the element in the data product collection
    //
    while (iGateToWaveform != gwend) {
      if (iGateToWaveform->first.key() == iGate) break;
      ++iGateToWaveform;
    } // while
    auto const iFirstWaveform = iGateToWaveform;
    while (iGateToWaveform != gwend) {
      if (iGateToWaveform->first.key() != iGate) break;
      ++iGateToWaveform;
    } // while
    if (iFirstWaveform == iGateToWaveform) {
      throw cet::exception("FillTriggerGates")
        << " Could not find any waveform associated to trigger gate #" << iGate
        << "\n";
    }

    // NOTE we do not control that all waveforms come from the same channel
    GateData_t gateData(*(iFirstWaveform->second)); // add the first waveform
    auto iEndWaveform = iFirstWaveform;
    while (++iEndWaveform != iGateToWaveform) // add the other waveforms
      gateData.add(*(iEndWaveform->second));

    gateData = gate; // copy the gate data from the data product

    allGates.push_back(std::move(gateData));

  } // for gates

  return allGates;
} // icarus::trigger::FillTriggerGates()


// -----------------------------------------------------------------------------
// --- icarus::trigger::TriggerGateIndex
// -----------------------------------------------------------------------------
template <typename GateObject>
icarus::trigger::TriggerGateIndex<GateObject>::TriggerGateIndex
  (std::size_t expectedChannels /* = 0U */)
{
  if (expectedChannels > 0) fGates.reserve(expectedChannels);
} // icarus::trigger::TriggerGateIndex::TriggerGateIndex(size_t)

// -----------------------------------------------------------------------------
template <typename GateObject>
icarus::trigger::TriggerGateIndex<GateObject>::TriggerGateIndex
  (std::vector<TriggerGate_t> const& gates)
{
  for (TriggerGate_t const& gate: gates) add(gate);
} // icarus::trigger::TriggerGateIndex::TriggerGateIndex(gates)


// -----------------------------------------------------------------------------
template <typename GateObject>
void icarus::trigger::TriggerGateIndex<GateObject>::add
  (TriggerGate_t const& gate)
{

  for (raw::Channel_t const channel: gate.channels()) {
    auto const chIndex = static_cast<std::size_t>(channel);

    if (expandToHold(chIndex) || !fGates[chIndex]) {
      fGates[chIndex] = &gate;
    }
    else if (fGates[chIndex] != &gate) {
      throw cet::exception("TriggerGateIndex") << "TriggerGateIndex::add(): "
        << "A gate was already present for channel " << channel << ".\n";
    }

  } // for

} // icarus::trigger::TriggerGateIndex::add()


// -----------------------------------------------------------------------------
template <typename GateObject>
unsigned int icarus::trigger::TriggerGateIndex<GateObject>::nChannels() const
  { return static_cast<unsigned int>(fGates.size()); }


// -----------------------------------------------------------------------------
template <typename GateObject>
auto icarus::trigger::TriggerGateIndex<GateObject>::find
  (raw::Channel_t const channel) const -> TriggerGate_t const*
{
  auto const chIndex = channelIndex(channel);
  return (chIndex < nChannels())? fGates[chIndex]: nullptr;
} // icarus::trigger::TriggerGateIndex::find()


// -----------------------------------------------------------------------------
template <typename GateObject>
auto icarus::trigger::TriggerGateIndex<GateObject>::operator[]
  (raw::Channel_t const channel) const -> TriggerGate_t const&
{
  auto const gate = find(channel);
  if (gate) return *gate;
  throw std::out_of_range(std::to_string(channel));
} // icarus::trigger::TriggerGateIndex::operator[]()


// -----------------------------------------------------------------------------
template <typename GateObject>
bool icarus::trigger::TriggerGateIndex<GateObject>::expandToHold
  (std::size_t chIndex)
{

  if (chIndex < fGates.size()) return false;

  fGates.resize(chIndex + 1U, nullptr);
  return true;

} // icarus::trigger::TriggerGateIndex::expandToHold()


// -----------------------------------------------------------------------------
template <typename GateObject>
std::size_t icarus::trigger::TriggerGateIndex<GateObject>::channelIndex
  (raw::Channel_t channel)
  { return static_cast<std::size_t>(channel); }


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERDATAUTILS_H
