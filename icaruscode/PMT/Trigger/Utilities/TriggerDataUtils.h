/**
 * @file   icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h
 * @brief  Utilities for the conversion of trigger gate data formats.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERDATAUTILS_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERDATAUTILS_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h"
#include "icaruscode/Utilities/DataProductPointerMap.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <map>
#include <vector>
#include <tuple>
#include <string> // std::to_string()
#include <stdexcept> // std::out_of_range
#include <utility> // std::move()
#include <type_traits> // std::enable_if_t()


//------------------------------------------------------------------------------
namespace icarus::trigger {

  // ---------------------------------------------------------------------------

  /// Map `util::DataProductPointerMap_t` for `raw::OpDetWaveform` objects.
  using OpDetWaveformDataProductMap_t
    = util::DataProductPointerMap_t<raw::OpDetWaveform>;

  /// Map `util::DataProductPointerMap_t` for `sbn::OpDetWaveformMeta` objects.
  using OpDetWaveformMetaDataProductMap_t
    = util::DataProductPointerMap_t<sbn::OpDetWaveformMeta>;

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
   * gates will not be associated any more to the contributing waveforms.
   *
   * The return value is a collection of trigger gate data
   * (`icarus::trigger::OpticalTriggerGateData_t`), with data _stolen_
   * from `gates`.
   *
   * The trigger gates are processed in the same order as they are in `gates`.
   *
   * After the function returns, `gates` will have been depleted of all the gate
   * data; any tracking information will be still associated to each gate, whose
   * gate data will be in an invalid state anyway, only good for destruction.
   *
   * The object `gates` is a iterable collection of gate objects: types derived
   * from `icarus::trigger::ReadoutTriggerGate` and from
   * `TrackedOpticalTriggerGate` are supported.
   *
   * In the following example, we start with trigger gates already serialized
   * in a _art_ event. The serialization splits a trigger gate object in two
   * components: the gate levels, and the associated waveforms.
   * In the first part of the example we recover the information from the event
   * and we assemble it into the standard trigger gate objects (of type
   * `icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>`).
   * After some unspecified and optional processing, `gates` are disassembled
   * to be saved into the event: this is achieved by a call to
   * `transformIntoOpticalTriggerGate()` which produces the trigger gate data
   * and their associations to the waveforms.
   * In the last part, these components are both stored into the event.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * //
   * // somehow get/build a collection of trigger gates;
   * // here we use `FillTriggerGates()` to read existing data from the event
   * //
   * std::vector<icarus::trigger::TrackedOpticalTriggerGate<>> gates
   *   = icarus::trigger::ReadTriggerGates(event, "orig");
   * 
   * // ...
   * 
   * //
   * // use the created vector and associations (e.g. put them into art event)
   * //
   * event.put(std::make_unique<std::vector<icarus::trigger::OpticalTriggerGateData_t>>
   *   (icarus::trigger::transformIntoOpticalTriggerGate(std::move(gates))));
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Also note that in the omitted processing part the trigger gates might have
   * been combined into a different type of gates instead.
   * To serialize that collection, the exact same procedure would be used.
   * Note that in this example the association to the source waveform metadata
   * is not serialised (and presumably lost). Another version of 
   * `transformIntoOpticalTriggerGate()` also provides support for the
   * association to waveform sources.
   */
  template <typename Gates>
  std::vector<icarus::trigger::OpticalTriggerGateData_t>
  transformIntoOpticalTriggerGate(Gates&& gates);
  

  // ---------------------------------------------------------------------------
  /**
   * @brief Returns the trigger gates in serializable format.
   * @tparam OpDetInfo type of the associated waveforms
   * @param gates the data to be reformatted (*data will be stolen!*)
   * @param makeGatePtr _art_ pointer maker for the gate data
   * @param opDetInfoPtrs map of art pointers to optical waveform information
   * @return a pair: gate data product and associations to `OpDetInfo`
   *
   * This function transfers the data from the original structured `gates` into
   * a data collection suitable for serialization with _art_, including the
   * association of each gate with all its contributing waveforms.
   *
   * The return value is a tuple of two elements:
   * * `0`: collection of trigger gate data
   *     (`icarus::trigger::OpticalTriggerGateData_t`), with data _stolen_
   *     from `gates`;
   * * `1`: association between trigger gate data and their optical waveform
   *     information, as reported by `opDetInfoPtrs`.
   *
   * The trigger gates are processed in the same order as they are in `gates`,
   * and the associations to the waveform information are set gate by gate,
   * in the same order as they are reported by the tracking information of the
   * input gate.
   *
   * After the function returns, `gates` will have been depleted of all the data
   * and left in an unspecified state only good for destruction.
   *
   * The object `gates` is a sequence of
   * `icarus::trigger::TrackedOpticalTriggerGate` objects, each of them
   * containing a trigger gate object
   * (`icarus::trigger::OpticalTriggerGateData_t`) and a list of tracked sources
   * in the form of `OpDetInfo` objects.
   *
   *
   * In the following example, we start with trigger gates already serialized
   * in a _art_ event. The serialization splits a trigger gate object in two
   * components: the gate levels, and the associated waveform information.
   * In the first part of the example we recover the information from the event
   * and we assemble it into the standard trigger gate objects (of type
   * `icarus::trigger::OpticalTriggerGateData_t`).
   * After some unspecified and optional processing, `gates` are disassembled
   * to be saved into the event: this is achieved by a call to
   * `transformIntoOpticalTriggerGate()` which produces the trigger gate data
   * and their associations to the waveform information.
   * In the last part, these components are both stored into the event.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using icarus::trigger::OpticalTriggerGateData_t; // for convenience
   * 
   * //
   * // somehow get/build a collection of trigger gates;
   * // here we use `FillTriggerGates()` to read existing data from the event
   * //
   * std::vector<icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>> gates
   *   = icarus::trigger::FillTriggerGates<sbn::OpDetWaveformMeta>
   *   (
   *   event.getProduct<std::vector<OpticalTriggerGateData_t>>("orig"),
   *   event.getProduct<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>("orig")
   *   );
   * 
   * 
   * // ...
   * 
   * //
   * // optical waveform information to pointer map is required to create
   * // associations between the trigger gates and their waveforms
   * //
   * auto const& opDetInfoPtrs = util::mapDataProductPointers
   *   (event, event.getValidHandle<std::vector<raw::OpDetWaveformMeta>>("opdaq"));
   * // transform the data; after this line, `gates` is not usable any more
   * auto thresholdData = icarus::trigger::transformIntoOpticalTriggerGate
   *   (std::move(gates), makeGatePtr, opDetInfoPtrs);
   * 
   * //
   * // use the created vector and associations (e.g. put them into art event)
   * //
   * event.put(
   *   std::make_unique<std::vector<OpticalTriggerGateData_t>>
   *     (std::move(std::get<0U>(thresholdData)))
   *   );
   * event.put(
   *   std::make_unique<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformInfo>>
   *     (std::move(std::get<1U>(thresholdData)))
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this example the most convenient waveform information is used,
   * `sbn::OpDetWaveformMeta`; the same syntax would hold for using the more
   * traditional, but encumbering, `raw::OpDetWaveform`.
   */
  template <typename OpDetInfo = sbn::OpDetWaveformMeta>
  std::tuple<
    std::vector<icarus::trigger::OpticalTriggerGateData_t>,
    art::Assns<icarus::trigger::OpticalTriggerGateData_t, OpDetInfo>
    >
  transformIntoOpticalTriggerGate(
    std::vector<icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>>&& gates,
    art::PtrMaker<icarus::trigger::OpticalTriggerGateData_t> const& makeGatePtr,
    util::DataProductPointerMap_t<OpDetInfo> const& opDetInfoPtrs
    );
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Creates gate objects out of trigger gate data products.
   * @tparam OpDetInfo type of the associated waveforms
   * @param gates collection of data from all the gates
   * @param gateToWaveformInfo association between each gate and its waveforms
   * @return a STL vector of `TrackedOpticalTriggerGate` with copy of data from
   *         `gates`
   *
   * Objects like `icarus::trigger::TrackedOpticalTriggerGate` include the gate
   * level information plus the connections to the source objects (waveforms).
   * These objects are complicate enough that they are not saved directly into
   * an _art_ event. Rather, they are diced into pieces and the pieces are
   * stored.
   * 
   * This function stitches the pieces and returns back an object like
   * `icarus::trigger::TrackedOpticalTriggerGate`.
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>>
   * readTriggerGates(art::Event const& event, art::InputTag const& dataTag)
   * {
   *
   *   auto const& gates
   *     = event.getProduct<std::vector<TriggerGateData_t>>(dataTag);
   *   auto const& gateToWaveforms
   *     = event.getProduct<art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta>>
   *       (dataTag);
   *   return icarus::trigger::FillTriggerGates<sbn::OpDetWaveformMeta>
   *     (gates, gateToWaveforms);
   *
   * } // readTriggerGates()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will read the information from `dateTag` data products in an _art_ `event`
   * and return a collection of `icarus::trigger::TrackedOpticalTriggerGate`
   * objects, each one with a _copy_ of the data of the original data product
   * and a pointer to each of the associated waveform metadata.
   * Note that this specific functionality is now wrapped into
   * `icarus::trigger::ReadTriggerGates()`.
   *
   *
   * Return value
   * -------------
   *
   * The returned collection contains one element for each
   * `OpticalTriggerGateData_t` object in `gates`, in the same order.
   * Each of these elements is of type
   * `icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>`, contains a copy
   * of the data of the corresponding gate, and a list of optical waveform
   * information pointers (`OpDetInfo` objects) it is associated to.
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
  template <typename OpDetInfo>
  std::vector<icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>>
  FillTriggerGates(
    std::vector<icarus::trigger::OpticalTriggerGateData_t> const& gates,
    art::Assns<icarus::trigger::OpticalTriggerGateData_t, OpDetInfo> const& gateToWaveformInfo
    );

  
  // ---------------------------------------------------------------------------
  template <typename OpDetInfo = sbn::OpDetWaveformMeta>
  class TriggerGateReader; // documentation is at definition
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Assembles and returns trigger gates from serialized data.
   * @tparam OpDetInfo type of object associated to the gates
   * @tparam Event type of object to read the needed information from
   * @param event object to read all the needed information from
   * @param dataTag tag of the data to read from the `event`
   * @return a collection of "tracking" trigger gates
   * @throw cet::exception wrapping any `cet::exception` thrown internally
   *                       (typically a data product not found)
   * @see `icarus::trigger::TriggerGateReader`, `icarus::trigger::FillTriggerGates()`
   * 
   * This function returns "tracking" trigger gates from data read from the data
   * read from `event`.
   * It is a one-shot replacement for `TriggerGateReader` to be used when
   * data to be consumed needn't or can't be declared, or when that declaration
   * is performed separately.
   */
  template <typename OpDetInfo = sbn::OpDetWaveformMeta, typename Event>
  std::vector<icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>>
  ReadTriggerGates(Event const& event, art::InputTag const& dataTag);
  
  
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
/**
 * @brief Assembles and returns trigger gates from serialized data.
 * @tparam OpDetInfo type of object associated to the gates
 * @tparam Event type of object to read the needed information from
 * @param event object to read all the needed information from
 * @param dataTag tag of the data to read from the `event`
 * @return a collection of "tracking" trigger gates
 * @throw cet::exception wrapping any `cet::exception` thrown internally
 *                       (typically a data product not found)
 * @see `icarus::trigger::FillTriggerGates()`
 * 
 * This class manages reading from a data source trigger gates
 * (vector of `icarus::trigger::OpticalTriggerGateData_t` elements) and their
 * associations to a tracking information (`OpDetInfo`: `sbn::OpDetWaveformMeta`
 * by default, but it may also be e.g. `raw::OpDetWaveform`).
 * 
 * All this information is merged into the "tracking" trigger gates, which
 * wrap _ a copy_ of the trigger gates and a reference to the tracked
 * information (by simple C pointer).
 * 
 * This class is a framework interface to `FillTriggerGates()`;
 * interaction  with the data source and the framework is driven by the
 * `Event` object in the `read()` method, which is expected to expose an
 * `interface similar to _art_'s `art::Event`, and in particular:
 * * `Event::getProduct<T>(art::InputTag)` to read and return directly
 *   a reference to the persistent, framework-owned data products of type `T`.
 * 
 * In addition, this class provides a shortcut to declare the consumed data
 * products, again in a _art_-like interface.
 * A full usage example in a module includes:
 * 
 * 1. a reader data member in the module class reading the data:
 *     
 *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 *     icarus::trigger::TriggerGateReader<> const fGateReader;
 *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *     
 * 2. construction in the module constructor initializer list from a configured
 *     input tag representing the data to be read:
 *     
 *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 *       , fGateReader{ config().TriggerGateTag() }
 *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 3. also in the module constructor body, declaration of what is going to be
 *    read:
 *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 *       fGateReader.declareConsumes(consumesCollector());
 *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 4. actual reading in a context where the current `event` is known:
 *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 *     std::vector<icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>>
 *       gates = fGateReader(event);
 *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * In a context where the consume declaration is not relevant/desired, the first
 * three steps can be omitted and the reading can happen with a temporary
 * `TriggerGateReader` object; or the shortcut function
 * `icarus::trigger::ReadTriggerGates()` can be used instead:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::vector<icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>>
 *   gates = icarus::trigger::ReadTriggerGates<sbn::OpDetWaveformMeta>(event, dataTag);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * is equivalent to the last step above. Note that it is still possible to use
 * temporaries also for the consume declaration: for example, a program using
 * `icarus::trigger::ReadTriggerGates()` may in the proper place (e.g. a module
 * constructor) call explicitly the consume declaration:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 *   icarus::trigger::TriggerGateReader<>{ config().TriggerGateTag() }
 *     .declareConsumes(consumesCollector());
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * is equivalent to the three steps above. The only advantage of the approach
 * in multiple steps of the first example is that consistency of declaration
 * and use is guaranteed by construction.
 * 
 */
template <typename OpDetInfo /* = sbn::OpDetWaveformMeta */>
class icarus::trigger::TriggerGateReader {
  
  art::InputTag fDataTag; ///< Data tag of the data products being read.
  
    public:
  
  using OpDetWaveformMeta_t = OpDetInfo; ///< Type of associated information.
  
  /// Type of data objects returned.
  using TriggerGates_t
    = icarus::trigger::TrackedOpticalTriggerGate<OpDetWaveformMeta_t>;
  
  
  /// Constructor: configure to read the data with the specified tag.
  TriggerGateReader(art::InputTag dataTag): fDataTag(std::move(dataTag)) {}
  
  /// Declares to the `collector` the data products that are going to be read.
  template <typename ConsumesCollector>
  void declareConsumes(ConsumesCollector& collector) const;
  
  /// Returns the input tag the data is read from.
  art::InputTag const& inputTag() const { return fDataTag; }
  
  // @{
  /// Reads the configured data product from the specified `event`.
  template <typename Event>
  std::vector<TriggerGates_t> read(Event const& event) const;
  template <typename Event>
  std::vector<TriggerGates_t> operator() (Event const& event) const
    { return read(event); }
  // @}
  
  
}; // class icarus::trigger::TriggerGateReader



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename Gates>
std::vector<icarus::trigger::OpticalTriggerGateData_t>
icarus::trigger::transformIntoOpticalTriggerGate(Gates&& gates) {
  
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGateData_t;
  
  std::vector<TriggerGateData_t> gateData;

  // this function also captures a `Gate const& gates`, which would yield to
  // const `gate` in this loop: `auto&&` matches it as exactly as it can
  for (auto&& gate: icarus::trigger::gatesIn(gates))
    gateData.push_back(std::move(gate));
  
  return gateData;
  
} // icarus::trigger::transformIntoOpticalTriggerGate(Gates)


// -----------------------------------------------------------------------------
template <typename OpDetInfo /* = sbn::OpDetWaveformMeta */>
std::tuple<
  std::vector<icarus::trigger::OpticalTriggerGateData_t>,
  art::Assns<icarus::trigger::OpticalTriggerGateData_t, OpDetInfo>
  >
icarus::trigger::transformIntoOpticalTriggerGate(
  std::vector<icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>>&& gates,
  art::PtrMaker<icarus::trigger::OpticalTriggerGateData_t> const& makeGatePtr,
  util::DataProductPointerMap_t<OpDetInfo> const& opDetInfoPtrs
) {
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGateData_t;

  // note that we rely on the fact that `gates` still contains trigger gates,
  // and those trigger gates still know about the associated waveforms;
  // attempting to access their levels, though, would be an error
  // (if really needed, we can still find them in `gateData`)
  art::Assns<TriggerGateData_t, OpDetInfo> gateToWaveforms;

  for (auto const& [ iGate, gate ]: util::enumerate(gates)) {

    if (!gate.tracking().hasTracked()) continue;

    //
    // produce the associations
    //

    // pointer to the gate data we have just added:
    art::Ptr<TriggerGateData_t> const& gatePtr = makeGatePtr(iGate);
    for (OpDetInfo const* waveform: gate.tracking().getTracked()) {
      
      gateToWaveforms.addSingle(gatePtr, opDetInfoPtrs.at(waveform));

    } // for all waveform info

  } // for all gates

  //
  // create objects for the data products
  //
  std::vector<TriggerGateData_t> gateData
    = transformIntoOpticalTriggerGate(icarus::trigger::gatesIn(gates));
  
  return { std::move(gateData), std::move(gateToWaveforms) };

} // icarus::trigger::transformIntoOpticalTriggerGate()


// -----------------------------------------------------------------------------
template <typename OpDetInfo>
std::vector<icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>>
icarus::trigger::FillTriggerGates(
  std::vector<icarus::trigger::OpticalTriggerGateData_t> const& gates,
  art::Assns<icarus::trigger::OpticalTriggerGateData_t, OpDetInfo> const& gateToWaveformInfo
) {

  std::vector<icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>>
    allGates; // the output

  auto iGateToWaveform = gateToWaveformInfo.begin();
  auto const gwend = gateToWaveformInfo.end();

  for (auto const& [ iGate, gate ]: util::enumerate(gates)) {

    icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo> trackedGate { gate };
    
    //
    // find the relevant waveforms for this gate;
    // match by index of the element in the data product collection
    //
    while (iGateToWaveform != gwend) {
      if (iGateToWaveform->first.key() == iGate) break;
      ++iGateToWaveform;
    } // while

    while (iGateToWaveform != gwend) {
      if (iGateToWaveform->first.key() != iGate) break;
      
      OpDetInfo const& wfInfo = *(iGateToWaveform->second);
      
      // add the associated information to the gate:
      trackedGate.gate().addChannel(wfInfo.ChannelNumber());
      trackedGate.tracking().add(&wfInfo);
      
      ++iGateToWaveform;
    } // while
    
    allGates.push_back(std::move(trackedGate));
    
  } // for gates

  return allGates;
} // icarus::trigger::FillTriggerGates()


// -----------------------------------------------------------------------------
// --- icarus::trigger::TriggerGateReader
// -----------------------------------------------------------------------------
template <typename OpDetInfo /* = sbd::OpDetWaveformMeta */>
template <typename ConsumesCollector>
void icarus::trigger::TriggerGateReader<OpDetInfo>::declareConsumes
  (ConsumesCollector& collector) const
{

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  collector.template consumes<std::vector<OpticalTriggerGateData_t>>(fDataTag);
  collector.template consumes<art::Assns<OpticalTriggerGateData_t, OpDetInfo>>
    (fDataTag);
  
} // icarus::trigger::TriggerGateReader::declareConsumes()


// -----------------------------------------------------------------------------
template <typename OpDetInfo /* = sbd::OpDetWaveformMeta */>
template <typename Event>
std::vector<icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>>
icarus::trigger::TriggerGateReader<OpDetInfo>::read(Event const& event) const {

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // currently the associations are a waste of time memory...
  auto const& gates
    = event.template getProduct<std::vector<OpticalTriggerGateData_t>>(fDataTag);
  auto const& gateToWaveforms = event.template getProduct
    <art::Assns<OpticalTriggerGateData_t, OpDetInfo>>(fDataTag);
  
  try {
    return icarus::trigger::FillTriggerGates(gates, gateToWaveforms);
  }
  catch (cet::exception const& e) {
    throw cet::exception("readTriggerGates", "", e)
      << "readTriggerGates() encountered an error while reading data products from '"
      << fDataTag.encode() << "'\n";
  }

} // icarus::trigger::TriggerGateReader::read()


// -----------------------------------------------------------------------------
template <typename OpDetInfo /* = sbd::OpDetWaveformMeta */, typename Event>
std::vector<icarus::trigger::TrackedOpticalTriggerGate<OpDetInfo>>
icarus::trigger::ReadTriggerGates
  (Event const& event, art::InputTag const& dataTag)
{
  return TriggerGateReader<OpDetInfo>{ dataTag }.read(event);
} // icarus::trigger::TriggerSimulationOnGates::ReadTriggerGates()


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
