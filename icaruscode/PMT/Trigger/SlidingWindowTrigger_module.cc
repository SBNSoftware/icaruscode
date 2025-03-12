/**
 * @file   SlidingWindowTrigger_module.cc
 * @brief  Combines trigger channels into sliding windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 6, 2020
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowCombinerAlg.h"
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefinitionAlg.h"
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/OpDetWaveformMetaMatcher.h"
#include "icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateDataFormatting.h"
#include "icaruscode/PMT/Algorithms/PMTverticalSlicingAlg.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h" // sbn::OpDetWaveformMeta

// LArSoft libraries
#include "lardata/Utilities/NestedIterator.h" // lar::double_fwd_const_iterator
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/get_elements.h" // util::get_elements()
// #include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::::static_assert_on<>
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"


// C/C++ standard libraries
#include <numeric> // std::iota()
#include <map>
#include <vector>
#include <string>
#include <optional>
#include <memory> // std::unique_ptr
#include <utility> // std::move()
#include <cassert>
#if 0
#include <set>
#include <cstddef> // std::size_t

#endif // 0

//------------------------------------------------------------------------------
namespace icarus::trigger { class SlidingWindowTrigger; }

/**
 * @brief Combines trigger channels (V1730 LVDS gates) into sliding windows.
 * 
 * This module combines the input discriminated waveforms into a set of new
 * multi-level discriminated gates according to the geometric criterion of a
 * sliding window.
 *
 * The combination is a sum of all the discriminated waveforms within the
 * window.
 * Different optical detector planes are treated independently and separately.
 * It is expected that the input is from V1730 LVDS gates
 * (see `icarus::trigger::LVDSgates` module) which combine (usually) two optical
 * detector channels each.
 *
 *
 * Sliding window configuration
 * -----------------------------
 *
 * The size of the window is specified in terms on number of optical detector
 * channels per window and of stride. The windows are allowed to overlap, but
 * they have all the same number of channels.
 * If the parameters yield an invalid configuration, the module will throw an
 * exception.
 *
 * The definition of a valid window configuration is quite restrictive.
 * Optical detector channels are grouped by the TPC width coordinate (_z_) into
 * "towers". In ICARUS, a tower is made of either two or three channels,
 * depending on the position in _z_. A window can't split a tower, but rather
 * has to include whole towers.
 * In addition, the input discriminated waveforms may include information from
 * multiple channels (that is the common case if the input is LVDS gates).
 * A window can't cover only part of the channels of any input, but rather must
 * extend wide enough so that all the channels in the included gates are
 * covered. Channels may be exempted by marking them as "missing".
 * Finally, all the channels which are not "missing" must be covered by at least
 * one window.
 *
 * For the sliding window configuration to be valid, all the windows starting
 * with an offset multiple of the stride must be valid.
 *
 * For the standard pairing of LVDS gates in ICARUS, this translates into the
 * only allowed window sizes are by 10 and 15 channels or their multiples,
 * with a stride also in multiples of 10 or 15 respectively.
 * For example, the following configurations are all valid:
 * * 10 channels per window, with stride of 10 channels (9 windows per plane);
 * * 20 channels per window, with stride of 10 channels (8 windows per plane);
 * * 30 channels per window, with stride of 10 channels (7 windows per plane);
 * * 15 channels per window, with stride of 15 channels (6 windows per plane);
 * * 30 channels per window, with stride of 30 channels (3 windows per plane);
 * * 30 channels per window, with stride of 15 channels (5 windows per plane);
 * * 45 channels per window, with stride of 45 channels (2 windows per plane).
 *
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (labels out of
 *   `TriggerGatesTag` and `Thresholds`): full sets of discriminated waveforms,
 *   each waveform possibly covering multiple optical channels,
 *   and their associations to waveform metadata (`sbn::OpDetWaveformMeta`).
 *   One set per threshold.
 *
 * Requirements
 * -------------
 *
 * The algorithms require that in each set all optical channels are covered.
 *
 *
 * Output data products
 * =====================
 *
 * * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (instance name:
 *   same as the input gates): sets of gates combined according to the window
 *   configuration; one set per input threshold.
 * * `art::Assns<icarus::trigger::OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>`
 *   (instance name: same as the input gates): associations between each
 *   produced gate and the metadata of optical waveforms providing the original
 *   data.
 * * `art::Assns<icarus::trigger::OpticalTriggerGateData_t, raw::OpDetWaveform>`
 *   (instance name: same as the input gates): associations between each
 *   produced gate and the optical waveforms providing the original data.
 *   It is produced only if `ProduceWaveformAssns` configuration parameter is
 *   `true`, and it relies on the assumption that there is an association
 *   available between each `sbn::OpDetWaveformMeta` and its
 *   `raw::OpDetWaveform`, produced by the same module (i.e. with the same input
 *   tag) as the one of the original `sbn::OpDetWaveformMeta` data product
 *   itself.
 *
 * If window selection is requested (with `EnableOnlyWindows` or
 * `DisableWindows` configuration parameters), only the selected windows will
 * have an output entry. While each trigger data object comes with all the
 * channels it covers, there is no explicit information of the _index_ of the
 * surviving windows.
 *
 * @note At the moment, associations between input and output gates is not
 *       produced. `art::Assns` would not support it, requiring distinct
 *       types for the associated objects.
 *
 *
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description SlidingWindowTrigger`.
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module
 *     instance which produced the discriminated waveforms; it must not include
 *     any instance name, as the instance names will be automatically added from
 *     `Thresholds` parameter.
 * * `Thresholds` (list of strings, mandatory): list of the discrimination
 *     thresholds to consider. A data product containing a digital signal is
 *     read for each one of the thresholds, and the tag of the data product is
 *     expected to be the module label `TriggerGatesTag` with as instance name
 *     the value of the entry in this list. While it is common for this entry to
 *     be the actual threshold value in ADC (e.g. `"400"`), this is not
 *     required.
 * * `WindowSize` (integral, mandatory): the number of _optical detector
 *     channels_ that every window must include. Note that this is not
 *     equivalent to the number of input objects from the `TriggerGatesTag` data
 *     products: for example, in the standard ICARUS configuration, a window
 *     spanning 15 channels is generated out of only 8 LVDS input gates.
 * * `Stride` (integral, default: as `WindowSize`): the base offset, in number
 *     of optical detector channels, where to start a new window. Windows are
 *     started at all multiples of `Stride`. For example, windows of 30 channels
 *     with `Stride` 15 will start after 0, 15, 30, 45 and 60 channels on each
 *     optical detector plane.
 * * `EnableOnlyWindows` (list of integers, default: omitted): if specified,
 *     only the windows with index in this list will be processed and output;
 *     mutually exclusive with `DisableWindows`.
 * * `DisableWindows` (list of integers, default: omitted): if specified,
 *     the windows with index in this list will be excluded from processing and
 *     from output; mutually exclusive with `EnableOnlyWindows`.
 * * `MissingChannels` (list of integers, default: empty): the channels whose ID
 *     is included in this list are expected and required not to be present in
 *     the input (i.e. no input gate should include them).
 * * `ProduceWaveformAssns` (flag, default: `true`): produce also associations
 *     between each gate and the `raw::OpDetWaveform` which contributed to it.
 * * `IgnoreErrors` (flag, default: `false`): catch processing exceptions and,
 *     in case, just stop processing the current event. Some of the data
 *     products may be written, while others not; _art_ won't be happy with
 *     this, and module's `errorOnFailureToPut` option will need to be set to
 *     `false`.
 * * `LogCategory` (string): name of the output stream category for console
 *     messages (managed by MessageFacility library).
 *
 * *Note*: when using `EnableOnlyWindows` or `DisableWindows` it may be useful
 * to check the output of the module (at construction) to verify that the
 * selected indices match the intended windows.
 * 
 */
class icarus::trigger::SlidingWindowTrigger: public art::EDProducer {
  
    public:
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;


    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    fhicl::Sequence<std::string> Thresholds {
      Name("Thresholds"),
      Comment
        ("thresholds to consider (instance names of `TriggerGatesTag` input)")
      };

    fhicl::Atom<unsigned int> WindowSize {
      Name("WindowSize"),
      Comment("number of optical channels to be included in each window")
      };
    fhicl::OptionalAtom<unsigned int> Stride {
      Name("Stride"),
      Comment(
        "number of optical channel used as offset for sliding window [as WindowSize]"
        )
      };
    
    fhicl::OptionalSequence<std::size_t> DisableWindows {
      Name("DisableWindows"),
      Comment("ignores the windows with the specified index"),
      [this](){ return !EnableOnlyWindows.hasValue(); }
      };

    fhicl::OptionalSequence<std::size_t> EnableOnlyWindows {
      Name("EnableOnlyWindows"),
      Comment("only enables the windows with the specified index"),
      [this](){ return !DisableWindows.hasValue(); }
      };
    
    fhicl::Sequence<raw::Channel_t> MissingChannels {
      Name("MissingChannels"),
      Comment("list of ID of channels missing from the input"),
      std::vector<raw::Channel_t>{}
      };

    fhicl::Atom<bool> ProduceWaveformAssns {
      Name("ProduceWaveformAssns"),
      Comment
        ("also produce gate/waveform associations together with gate/metadata"),
      true
      };
    
    fhicl::Atom<bool> IgnoreErrors {
      Name("IgnoreErrors"),
      Comment("in case of fatal error just stop processing the event"),
      false
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "SlidingWindowTrigger" // default
      };

  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit SlidingWindowTrigger(Parameters const& config);
  
  // Plugins should not be copied or assigned.
  SlidingWindowTrigger(SlidingWindowTrigger const&) = delete;
  SlidingWindowTrigger(SlidingWindowTrigger&&) = delete;
  SlidingWindowTrigger& operator=(SlidingWindowTrigger const&) = delete;
  SlidingWindowTrigger& operator=(SlidingWindowTrigger&&) = delete;
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Fills the plots. Also extracts the information to fill them with.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  /// Definition of all windows.
  using WindowDefs_t = icarus::trigger::TriggerWindowDefs_t;
  
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGate::GateData_t;
  
  /// Reconstituted trigger gate type internally used.
  using TrackedTriggerGate_t
    = icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>;
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// ADC thresholds to read, and the input tag connected to their data.
  std::map<std::string, art::InputTag> fADCthresholds;

  unsigned int const fWindowSize; ///< Sliding window size in number of channels.
  unsigned int const fWindowStride; ///< Sliding window base offset.
  
  /// Channel content of each window.
  WindowDefs_t const fWindowChannels;

  /// List of windows to be included.
  std::vector<std::size_t> const fEnabledWindows;

  /// Whether to produce gate/waveform associations.
  bool const fProduceWaveformAssns;
  
  bool const fIgnoreErrors; /// Whether to survive errors in the processing.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  /// Combining algorithm.
  icarus::trigger::SlidingWindowCombinerAlg const fCombiner;


  /// Returns the number of disabled windows.
  unsigned int nDisabledWindows() const noexcept
    { return fWindowChannels.size() - fEnabledWindows.size(); }
  
  /// Performs the combination for data with a specified threshold.
  void produceThreshold(
    art::Event& event,
    icarus::trigger::OpDetWaveformMetaDataProductMap_t& waveformMap,
    std::string const& threshold,
    art::InputTag const& dataTag
    ) const;
  
  /// Reads a set of input gates from the `event` and updates `waveformMap`.
  std::vector<icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>> ReadTriggerGates(
    art::Event const& event,
    art::InputTag const& dataTag,
    icarus::trigger::OpDetWaveformMetaDataProductMap_t& waveformMap
    ) const;
  
  /// Defines the channels falling in each window.
  WindowDefs_t defineWindows() const;
  
  /// Returns a list of enabled window indices.
  static std::vector<std::size_t> makeEnabledWindowIndices(
      std::size_t nWindows
    , fhicl::OptionalSequence<std::size_t> const& enabled
    , fhicl::OptionalSequence<std::size_t> const& disabled
    );

  /// Adds the waveforms in the specified association to the waveform `map`.
  static void UpdateWaveformMap(
    icarus::trigger::OpDetWaveformMetaDataProductMap_t& map,
    art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta> const& assns
    );

}; // class icarus::trigger::SlidingWindowTrigger



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Returns a collection of type CONT with copies of only the elements with
  /// the specified `indices`.
  template <typename Cont, typename Indices>
  Cont filter(Cont const& coll, Indices const& indices) {
    
    Cont selected;
    for (auto const index: indices) selected.push_back(coll.at(index));
    return selected;
    
  } // filter()
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::trigger::SlidingWindowTrigger
//------------------------------------------------------------------------------
icarus::trigger::SlidingWindowTrigger::SlidingWindowTrigger
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fWindowSize(config().WindowSize())
  , fWindowStride(config().Stride().value_or(fWindowSize))
  , fWindowChannels(defineWindows())
  , fEnabledWindows(makeEnabledWindowIndices(
     fWindowChannels.size(),
     config().EnableOnlyWindows, config().DisableWindows
     ))
  , fProduceWaveformAssns(config().ProduceWaveformAssns())
  , fIgnoreErrors(config().IgnoreErrors())
  , fLogCategory(config().LogCategory())
    // demand full PMT coverage only if no window was disabled:
  , fCombiner(
      filter(fWindowChannels, fEnabledWindows),
      config().MissingChannels(),
      nDisabledWindows() == 0U
    )
{
  //
  // more complex parameter parsing
  //
  std::string const discrModuleLabel = config().TriggerGatesTag();
  for (std::string const& threshold: config().Thresholds())
    fADCthresholds[threshold] = art::InputTag{ discrModuleLabel, threshold };

  //
  // configuration report (short)
  //
  {
    mf::LogInfo log(fLogCategory);
    log <<   "Trigger configuration: " << fWindowChannels.size() << " windows";
    if (fEnabledWindows.size() < fWindowChannels.size())
      log << " (only " << fEnabledWindows.size() << " enabled)";
    log << ":";
    for (auto const& [ iWindow, window ]: util::enumerate(fWindowChannels)) {
      log << "\n #" << iWindow << ": "
        << icarus::trigger::dumpTriggerWindowChannels(window);
      if (std::find(fEnabledWindows.begin(), fEnabledWindows.end(), iWindow)
        == fEnabledWindows.end())
      {
        log << " (disabled)";
      }
    } // for
    
    log << "\nConfigured " << fADCthresholds.size() << " thresholds (ADC):";
    for (auto const& [ threshold, dataTag ]: fADCthresholds)
      log << "\n * " << threshold << " (from '" << dataTag.encode() << "')";
  } // local block


  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  //
  // input data declaration
  //
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    consumes<std::vector<OpticalTriggerGateData_t>>(inputDataTag);
    consumes<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>
      (inputDataTag);
  } // for


  //
  // output data declaration
  //
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    produces<std::vector<OpticalTriggerGateData_t>>(inputDataTag.instance());
    produces<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>
      (inputDataTag.instance());
    if (fProduceWaveformAssns) {
      produces<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
        (inputDataTag.instance());
    }
  } // for

  
} // icarus::trigger::SlidingWindowTrigger::SlidingWindowTrigger()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTrigger::produce(art::Event& event) {
  
  icarus::trigger::OpDetWaveformMetaDataProductMap_t waveformMap;

  for (auto const& [ thresholdStr, dataTag ]: fADCthresholds) {
    try {
      produceThreshold(event, waveformMap, thresholdStr, dataTag);
    }
    catch (cet::exception const& e) {
      if (!fIgnoreErrors) throw;
      
      mf::LogError{ fLogCategory }
        << moduleDescription().moduleName()
        << "[" << moduleDescription().moduleLabel()
        << "]: Exception thrown while processing threshold " << thresholdStr
        << " of input '" << dataTag.encode() << "':\n"
        << e.what()
        << "\nException is being ignored, but no data product will be added.\n";
      // wise?
    }
    
  } // for all thresholds
  
} // icarus::trigger::SlidingWindowTrigger::produce()


//------------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowTrigger::defineWindows() const
  -> WindowDefs_t
{
  // delegated to `SlidingWindowDefinitionAlg`
  
  icarus::trigger::SlidingWindowDefinitionAlg const algo {
    *(lar::providerFrom<geo::Geometry>()),
    "SlidingWindowTrigger/SlidingWindowDefinitionAlg"
    };
  
  return algo.makeWindows(fWindowSize, fWindowStride);
  
} // icarus::trigger::SlidingWindowTrigger::defineWindows()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTrigger::produceThreshold(
  art::Event& event,
  icarus::trigger::OpDetWaveformMetaDataProductMap_t& waveformMap,
  std::string const& threshold,
  art::InputTag const& dataTag
) const {
  
  mf::LogDebug(fLogCategory)
    << "Processing threshold " << threshold
    << " from '" << dataTag.encode() << "'"
    ;

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience
  using TrackedTriggerGate_t
    = icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>;
  
  std::vector<TrackedTriggerGate_t> const& gates
    = ReadTriggerGates(event, dataTag, waveformMap);
  
  std::vector<TrackedTriggerGate_t> combinedGates = fCombiner.combine(gates);

  // transform the data; after this line, `gates` is not usable any more
  art::PtrMaker<OpticalTriggerGateData_t> const makeGatePtr
    (event, dataTag.instance());
  auto [ outputGates, outputAssns ]
    = icarus::trigger::transformIntoOpticalTriggerGate
    (std::move(combinedGates), makeGatePtr, waveformMap);

  {
    mf::LogTrace log { fLogCategory };
    log
      << "Threshold " << threshold << " ('" << dataTag.encode() << "'): "
      << gates.size() << " input channels, "
      << outputGates.size() << " output channels (+"
      << outputAssns.size() << " associations to waveforms) into '"
      << moduleDescription().moduleLabel() << ":" << dataTag.instance() << "'"
      ;
    for (auto const& [ iGate, gate ]: util::enumerate(outputGates)) {
      log << "\n [#" << iGate << "] " << compactdump(gate);
    } // for
  } // local scope

  if (fProduceWaveformAssns) {
    
    // produce one gate-waveform association for each gate-metadata one;
    // do it now while the gate/metadata association is still locally available
    icarus::trigger::OpDetWaveformMetaMatcher waveformMetaMatcher{ event };
    art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform> outputWaveAssns;
    for (auto const &[ gatePtr, metaPtr ]: outputAssns)
      outputWaveAssns.addSingle(gatePtr, waveformMetaMatcher(metaPtr));
    
    event.put(
      std::make_unique<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
        (std::move(outputWaveAssns)),
      dataTag.instance()
      );
  } // if fProduceWaveformAssns
  
  event.put(
    std::make_unique<std::vector<OpticalTriggerGateData_t>>
      (std::move(outputGates)),
    dataTag.instance()
    );
  event.put(
    std::make_unique<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>
      (std::move(outputAssns)),
    dataTag.instance()
    );

} // icarus::trigger::SlidingWindowTrigger::produceThreshold()


//------------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowTrigger::makeEnabledWindowIndices(
    std::size_t nWindows
  , fhicl::OptionalSequence<std::size_t> const& enabled
  , fhicl::OptionalSequence<std::size_t> const& disabled
) -> std::vector<std::size_t> {
  
  std::vector<std::size_t> selection; // indices of the enabled windows
  
  // add all enabled windows (or all if no enabled list is specified)
  if (!enabled(selection)) {
    selection.resize(nWindows);
    std::iota(selection.begin(), selection.end(), 0U);
  }
  
  // remove all disabled windows (if any)
  std::vector<std::size_t> removeWindows;
  disabled(removeWindows);
  auto const find = [&selection](std::size_t index)
    { return std::find(selection.begin(), selection.end(), index); };
  for (std::size_t const index: removeWindows) {
    if (auto const it = find(index); it != selection.end()) selection.erase(it);
  } // for
  
  return selection;
} // icarus::trigger::SlidingWindowTrigger::makeEnabledWindowIndices()


//------------------------------------------------------------------------------
void icarus::trigger::SlidingWindowTrigger::UpdateWaveformMap(
  icarus::trigger::OpDetWaveformMetaDataProductMap_t& map,
  art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta> const& assns
) {
  
  for (art::Ptr<sbn::OpDetWaveformMeta> const& wave: util::get_elements<1U>(assns))
    map.emplace(wave.get(), wave);
  
} // icarus::trigger::SlidingWindowTrigger::UpdateWaveformMap()


//------------------------------------------------------------------------------
auto icarus::trigger::SlidingWindowTrigger::ReadTriggerGates(
  art::Event const& event,
  art::InputTag const& dataTag,
  icarus::trigger::OpDetWaveformMetaDataProductMap_t& waveformMap
) const -> std::vector<TrackedTriggerGate_t> {
  using icarus::trigger::OpticalTriggerGateData_t;
  
  auto const& assns =
    event.getProduct<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>
    (dataTag);
  
  UpdateWaveformMap(waveformMap, assns);
  
  return icarus::trigger::FillTriggerGates
    (event.getProduct<std::vector<OpticalTriggerGateData_t>>(dataTag), assns);
  
} // icarus::trigger::SlidingWindowTrigger::ReadTriggerGates()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::SlidingWindowTrigger)


//------------------------------------------------------------------------------
