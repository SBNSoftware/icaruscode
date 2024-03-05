/**
 * @file   LVDSgates_module.cc
 * @brief  Combines discriminated PMT outputs into V1730 LVDS gates.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 13, 2019
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Utilities/TriggerGateOperations.h" // OpGates
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h"
#include "icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/OpDetWaveformMetaMatcher.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h" // sbn::OpDetWaveformMeta
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/Utilities/MultipleChoiceSelection.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/get_elements.h" // util::get_elements()
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
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"


// C/C++ standard libraries
#include <map>
#include <set>
#include <vector>
#include <regex>
#include <string>
#include <memory> // std::unique_ptr
#include <utility> // std::move()
#include <cassert>
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
namespace icarus::trigger { class LVDSgates; }

/**
 * @brief Combines discriminated PMT outputs into V1730 LVDS gates.
 * 
 * This module simulates the combination of discriminated PMT outputs
 * by the V1730 board LVDS output.
 * 
 * 
 * Input data products
 * ====================
 *
 * * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (labels out of
 *   `TriggerGatesTag` and `Thresholds`): full sets of discriminated waveforms,
 *   each waveform possibly covering multiple optical channels,
 *   and their associations to optical waveform metadata. One set per threshold.
 *
 *
 * Requirements
 * -------------
 * 
 * The algorithms require that the list of discriminated waveforms has in each
 * position the waveform from the optical detector channel number matching the
 * position (e.g. entry `0` will be the discriminated waveform from channel 0,
 * and so forth). No gaps are allowed: if no signal is present above threshold
 * for a certain channel, a discriminated waveform should appear for that
 * channel, with a gate always closed.
 * 
 *
 * Output data products
 * =====================
 *
 * * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (instance name:
 *   concatenation of input module label and instance name of the input gates;
 *   the former is omitted if it is `TriggerGatesTag`): sets of gates combined
 *   according to the configuration; one set per input threshold.
 * * `art::Assns<icarus::trigger::OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>`
 *   (instance name: same as above): associations between each
 *   produced gate and the metadata of the optical waveforms providing the
 *   original data.
 * * `art::Assns<icarus::trigger::OpticalTriggerGateData_t, raw::OpDetWaveform>`
 *   (instance name: same as above; optional): associations between each
 *   produced gate and the optical waveforms providing the original data.
 *   It is produced only if `ProduceWaveformAssns` configuration parameter is
 *   `true`, and it relies on the assumption that there is an association
 *   available between each `sbn::OpDetWaveformMeta` and its
 *   `raw::OpDetWaveform`, produced by the same module (i.e. with the same input
 *   tag) as the one of the original `sbn::OpDetWaveformMeta` data product
 *   itself.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description LVDSgates`.
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module
 *     instance which produced the discriminated waveforms; it must not include
 *     any instance name, as the instance names will be automatically added from
 *     `Thresholds` parameter. This value is used as module name for all
 *     thresholds which do not specify one.
 * * `Thresholds` (list of tags, mandatory): list of the discrimination
 *     thresholds to consider. A data product containing a digital signal is
 *     read for each one of the thresholds; each threshold entry is expected
 *     to be a input tag in the form `"ModuleLabel:InstanceName"`; both
 *     `ModuleName` and `InstanceName` may be omitted, in which case the colon
 *     `:` may also be omitted. In case the colon is omitted, the parameter will
 *     be considered an instance name if it represents an integer, and a module
 *     name otherwise. When the parameter is considered an instance name, the
 *     module is taken from the `TriggerGatesTag` parameter; when the instance
 *     name is omitted, it is considered empty. In this way, thresholds can be
 *     specified simply by a number, e.g. `[ 200, 400, 600 ]`, which will be
 *     translated into e.g. `"pmtfixedthr:200"`, `"pmtfixedthr:400"` and
 *     `"pmtfixedthr:600"`, while a full specification can be also used:
 *     `[ 200, "pmtthr:config" ]`, which will turn into `"pmtfixedthr:200"` and
 *     `"pmtthr:config"`; or a label-only one: `[ 200, "pmtthr" ]`, which will
 *     turn into `"pmtfixedthr:200"` and `"pmtthr"`.
 * * `ChannelPairing` (list of integral number lists): channels to combine;
 *     each element of this list is itself a list of optical detector channel
 *     numbers. The channels within each group are combined according to
 *     `CombinationMode`, resulting in a single combined gate associated to
 *     all the channels in the group.
 *     All optical detector channels *must* appear in the configuration, and
 *     each channel can be combined only in one group. Groups must not be empty.
 * * `IgnoreChannels` (list of integral numbers, optional): ID of the optical
 *     detector channels to skip.
 * * `CombinationMode` (either: `"disable"`, `"input1"`, `"input2"`, `"AND"`
 *     or `"OR"`): for each group of channels defined in `ChannelPairing`,
 *     the selected operation is used to combine all the channels in the group:
 *      * `disable`: all channels are discarded, and the resulting gates are
 *          all closed all the way;
 *      * `input1`: the first channel in the pairing is kept, and the others
 *          are discarded;
 *      * `input2`: the second channel in the pairing is kept, and the others
 *          are discarded;
 *      * `AND`: all channels are combined in a AND, resulting in a gate open
 *          only at the time when all the channel gates in the group are open;
 *          if there is a single channel in the group, the combined gate is a
 *          copy of it;
 *      * `OR`: all channels are combined in a OR, resulting in a gate open
 *          at every time when any of the channel gates in the group is open;
 *          if there is a single channel in the group, the combined gate is a
 *          copy of it.
 *     All pairings undergo the same operation (it is not possible, for example,
 *     disabling only one group of channels and combining the other groups with
 *     a `AND`).
 * * `ProduceWaveformAssns` (flag, default: `true`): produce also associations
 *     between each gate and the `raw::OpDetWaveform` which contributed to it.
 * * `LogCategory` (string): name of the output stream category for console
 *     messages (managed by MessageFacility library).
 *
 *
 * 
 */
class icarus::trigger::LVDSgates: public art::EDProducer {
  
    public:
  
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGate::GateData_t;
  
  /// Available operations for the combination of channels.
  enum class ComboMode { disable, Input1, Input2, AND, OR };
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    /// Selector for `CombinationMode` parameter.
    util::MultipleChoiceSelection<ComboMode> const CombinationModeSelector
      {
        { ComboMode::disable, "disable", "off" },
        { ComboMode::Input1,  "input1" },
        { ComboMode::Input2,  "input2" },
        { ComboMode::AND,     "AND" },
        { ComboMode::OR,      "OR" }
      };
    
    fhicl::OptionalAtom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of trigger gate extraction module (no instance name)")
      };

    fhicl::Sequence<std::string> Thresholds {
      Name("Thresholds"),
      Comment
        ("thresholds to consider (full tag or instance name of TriggerGatesTag)")
      };
    
    fhicl::Sequence<fhicl::Sequence<raw::Channel_t>> ChannelPairing {
      Name("ChannelPairing"),
      Comment
        ("grouping of optical detector channels (e.g.: [ [ 0, 2 ], [ 4, 6 ], [ 8 ], ... ])")
      };

    fhicl::Sequence<raw::Channel_t> IgnoreChannels {
      Name("IgnoreChannels"),
      Comment
        ("optical detector channels to ignore (e.g.: [ 54, 58, 67, 76, ... ])")
      };    

    fhicl::Atom<std::string> CombinationMode {
      Name("CombinationMode"),
      Comment("channel combination mode: " + CombinationModeSelector.optionListString())
      };
    
    fhicl::Atom<bool> ProduceWaveformAssns {
      Name("ProduceWaveformAssns"),
      Comment
        ("also produce gate/waveform associations together with gate/metadata"),
      true
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "LVDSgates" // default
      };
    
    
    ComboMode getCombinationMode() const
      {
        try {
          return CombinationModeSelector.parse(CombinationMode());
        }
        catch (util::MultipleChoiceSelectionBase::UnknownOptionError const& e) {
          throw art::Exception(art::errors::Configuration)
            << "Invalid value for 'CombinationMode' parameter: '" << e.label()
            << "'; valid options: "
            << CombinationModeSelector.optionListString() << ".\n";
        }
      } // getCombinationMode()
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit LVDSgates(Parameters const& config);
  
  // Plugins should not be copied or assigned.
  LVDSgates(LVDSgates const&) = delete;
  LVDSgates(LVDSgates&&) = delete;
  LVDSgates& operator=(LVDSgates const&) = delete;
  LVDSgates& operator=(LVDSgates&&) = delete;
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Fills the plots. Also extracts the information to fill them with.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  /// Type of member function to combine the current gate with another one.
  using BinaryCombinationFunc_t
    = OpticalTriggerGate& (OpticalTriggerGate::*)(OpticalTriggerGate const&);
  
  /// Reconstituted trigger gate type internally used.
  using TrackedTriggerGate_t
    = icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta>;
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Information for each source.
  struct SourceInfo_t {
    art::InputTag inputTag;
    std::string outputInstanceName;
  };
  
  /// ADC thresholds to read, and the input tag connected to their data.
  std::map<std::string, SourceInfo_t> fADCthresholds;
  
  /// Pairing of optical detector channels.
  std::vector<std::vector<raw::Channel_t>> fChannelPairing;
  ComboMode fComboMode; ///< The operation used to combinate channels.

  std::vector<raw::Channel_t> fIgnoreChannels;
  
  /// Whether to produce gate/waveform associations.
  bool fProduceWaveformAssns;
  
  std::string fLogCategory; ///< Message facility stream category for output.
  
  // --- END Configuration variables -------------------------------------------
  
  
  /// Completely processes the data for the specified threshold.
  void produceThreshold(
    art::Event& event,
    icarus::trigger::OpDetWaveformMetaDataProductMap_t& waveformMap,
    std::string const& thresholdStr,
    SourceInfo_t const& srcInfo
    ) const;

  /**
   * @brief Removes all the channel in `ignoreChannels` from `channelPairing`.
   * @return a copy of `channelPairing` with no element from `ignoreChannels`
   * 
   * Channel lists in the returned list may become empty.
   */
  std::vector<std::vector<raw::Channel_t>> removeChannels(
    std::vector<std::vector<raw::Channel_t>> channelPairing,
    std::vector<raw::Channel_t> const& ignoreChannels
    ) const;
  
  
  // --- BEGIN -- Combination modes --------------------------------------------
  /// Combination function for the `disable` mode.
  TrackedTriggerGate_t discardChannels(
    std::vector<TrackedTriggerGate_t> const&,
    std::vector<raw::Channel_t> const& pairing
    ) const;
  
  /// Combination function for the `input1` and `input2` modes.
  TrackedTriggerGate_t selectChannel(
    std::vector<TrackedTriggerGate_t> const&,
    std::vector<raw::Channel_t> const& pairing, std::size_t chosenIndex
    ) const;
  
  /// Combination function for the `AND` and `OR` modes; `op` is from `GateOps`.
  template <typename Op>
  TrackedTriggerGate_t binaryCombineChannel(
    std::vector<TrackedTriggerGate_t> const& gates,
    std::vector<raw::Channel_t> const& pairing,
    Op combine
    ) const;
  
  /// Performs the combination of a group of channels.
  TrackedTriggerGate_t combineChannels(
    std::vector<TrackedTriggerGate_t> const& gates,
    std::vector<raw::Channel_t> const& pairing,
    ComboMode comboMode
    ) const;
  // --- END -- Combination modes ----------------------------------------------
  
  
  
  // --- BEGIN -- Checks -------------------------------------------------------
  /// Checks the pairing configuration, throws on error.
  /// Returns the number of configured channels.
  unsigned int checkPairings() const;

  /// Input requirement check. Throws if requirements fail.
  void checkInput(std::vector<TrackedTriggerGate_t> const&  gates) const;
  
  // --- END -- Checks ---------------------------------------------------------
  
  
  /// Adds the associated waveforms into the map.
  static void UpdateWaveformMap(
    icarus::trigger::OpDetWaveformMetaDataProductMap_t& map,
    art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta> const& assns
    );
  
  /// Assembles trigger gates from `dataTag` data products in `event`.
  static std::vector<TrackedTriggerGate_t> ReadTriggerGates(
    art::Event const& event, art::InputTag const& dataTag,
    icarus::trigger::OpDetWaveformMetaDataProductMap_t& waveformMap
    );

  /// Converts a threshold string into an input tag.
  static art::InputTag makeTag
    (std::string const& thresholdStr, std::string const& defModule);

  /// Converts an input tag into an instance name for the corresponding output.
  static std::string makeOutputInstanceName
    (art::InputTag const& inputTag, std::string const& defModule);
  
}; // class icarus::trigger::LVDSgates



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
icarus::trigger::LVDSgates::LVDSgates
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fChannelPairing(config().ChannelPairing())
  , fComboMode(config().getCombinationMode())
  , fIgnoreChannels(config().IgnoreChannels())
  , fProduceWaveformAssns(config().ProduceWaveformAssns())
  , fLogCategory(config().LogCategory())
{
  //
  // more complex parameter parsing
  //
  std::string const discrModuleLabel = config().TriggerGatesTag().value_or("");
  for (std::string const& thresholdStr: config().Thresholds()) {
    art::InputTag const inputTag = makeTag(thresholdStr, discrModuleLabel);
    fADCthresholds[thresholdStr]
      = { inputTag, makeOutputInstanceName(inputTag, discrModuleLabel) };
  } // for all thresholds
  
  //
  // configuration validation
  //
  unsigned int nConfiguredChannels = checkPairings();

  //
  // configuration report (short)
  //
  {
    mf::LogInfo log(fLogCategory);
    log << nConfiguredChannels
      << " optical channels configured to be combined into "
      << fChannelPairing.size() << " LVDS outputs.";
    log << "\nIgnored " << fIgnoreChannels.size() << " channels.";
    log << "\nConfigured " << fADCthresholds.size() << " thresholds (ADC):";
    for (auto const& [ threshold, srcInfo ]: fADCthresholds) {
      log << "\n * " << threshold
        << " (from '" << srcInfo.inputTag.encode() << "')";
    }
  } // local block

  
  using icarus::trigger::OpticalTriggerGateData_t; // for convenience
  
  //
  // input data declaration
  //
  for (auto const& [ inputTag, outName ]: util::const_values(fADCthresholds)) {
    consumes<std::vector<OpticalTriggerGateData_t>>(inputTag);
    consumes<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>
      (inputTag);
  } // for
  
  
  //
  // output data declaration
  //
  for (auto const& [ inputTag, outName ]: util::const_values(fADCthresholds)) {
    produces<std::vector<OpticalTriggerGateData_t>>(outName);
    produces<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>(outName);
    if (fProduceWaveformAssns) {
      produces<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>
        (outName);
    }
  } // for
  
  
} // icarus::trigger::LVDSgates::LVDSgates()


//------------------------------------------------------------------------------
void icarus::trigger::LVDSgates::produce(art::Event& event) {
  
  icarus::trigger::OpDetWaveformMetaDataProductMap_t waveformMap;
  
  for (auto const& [ thresholdStr, srcInfo ]: fADCthresholds) {
    
    produceThreshold(event, waveformMap, thresholdStr, srcInfo);
    
  } // for all thresholds
  
} // icarus::trigger::LVDSgates::produce()


//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void icarus::trigger::LVDSgates::produceThreshold(
  art::Event& event,
  icarus::trigger::OpDetWaveformMetaDataProductMap_t& waveformMap,
  std::string const& thresholdStr,
  SourceInfo_t const& srcInfo
) const {
  
  auto const& [ dataTag, outputInstanceName ] = srcInfo;
  
  mf::LogDebug(fLogCategory)
    << "Processing threshold " << thresholdStr
    << " from '" << dataTag.encode() << "'";
  
  using icarus::trigger::OpticalTriggerGateData_t; // for convenience
  
  std::vector<TrackedTriggerGate_t> const& gates
    = ReadTriggerGates(event, dataTag, waveformMap);

  checkInput(gates);
  
  std::vector<TrackedTriggerGate_t> combinedGates;
  
  std::vector<std::vector<raw::Channel_t>> const cleanedChannels
    = removeChannels(fChannelPairing, fIgnoreChannels);

  for (std::vector<raw::Channel_t> const& pairing: cleanedChannels) {
    if (pairing.empty()) continue; // ???
    
    combinedGates.push_back(combineChannels(gates, pairing, fComboMode));
    
  } // for
  
  // transform the data; after this line, `gates` is not usable any more
  art::PtrMaker<OpticalTriggerGateData_t> const makeGatePtr
    (event, outputInstanceName);
  auto [ outputGates, outputAssns ]
    = icarus::trigger::transformIntoOpticalTriggerGate
    (std::move(combinedGates), makeGatePtr, waveformMap);
  
  mf::LogTrace(fLogCategory)
    << "Threshold " << thresholdStr << " ('" << dataTag.encode() << "'): "
    << gates.size() << " input channels, "
    << outputGates.size() << " output channels (+"
    << outputAssns.size() << " associations to waveforms) into '"
    << moduleDescription().moduleLabel() << ":" << outputInstanceName << "'"
    ;

  
  if (fProduceWaveformAssns) {
    
    // produce one gate-waveform association for each gate-metadata one;
    // do it now while the gate/metadata association is still locally available
    icarus::trigger::OpDetWaveformMetaMatcher waveformMetaMatcher{ event };
    art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform> outputWaveAssns;
    for (auto const & [ gatePtr, metaPtr ]: outputAssns)
      outputWaveAssns.addSingle(gatePtr, waveformMetaMatcher(metaPtr));
    
    event.put(
      std::make_unique<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
        (std::move(outputWaveAssns)),
      outputInstanceName
      );
    
  } // if fProduceWaveformAssns
  
  event.put(
    std::make_unique<std::vector<OpticalTriggerGateData_t>>
      (std::move(outputGates)),
    outputInstanceName
    );
  event.put(
    std::make_unique<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>
      (std::move(outputAssns)),
    outputInstanceName
    );
  
} // icarus::trigger::LVDSgates::produceThreshold()


//------------------------------------------------------------------------------
auto icarus::trigger::LVDSgates::combineChannels(
  std::vector<TrackedTriggerGate_t> const& gates,
  std::vector<raw::Channel_t> const& pairing,
  ComboMode comboMode
) const -> TrackedTriggerGate_t {
  
  namespace GateOps = icarus::trigger::GateOps;
  
  switch (comboMode) {
    case ComboMode::disable:
      return discardChannels(gates, pairing);
    case ComboMode::AND:
      return binaryCombineChannel(gates, pairing, GateOps::Min);
    case ComboMode::OR:
      return binaryCombineChannel(gates, pairing, GateOps::Max);
    case ComboMode::Input1:
      return selectChannel(gates, pairing, 0U);
    case ComboMode::Input2:
      return selectChannel(gates, pairing, 1U);
    default:
      throw art::Exception(art::errors::LogicError)
        << "Unexpected combination mode (#" << static_cast<int>(comboMode)
        << ").\n";
  } // switch(comboMode)
} // icarus::trigger::LVDSgates::combineChannels()


//------------------------------------------------------------------------------
std::vector<std::vector<raw::Channel_t>> icarus::trigger::LVDSgates::removeChannels(
  std::vector<std::vector<raw::Channel_t>> channelPairing,
  std::vector<raw::Channel_t> const& ignoreChannels
  ) const {

  for (std::vector<raw::Channel_t>& channels: channelPairing) {
    for (unsigned int j = channels.size(); j-- > 0; ) {
      std::vector<raw::Channel_t>::const_iterator const it
        = find(ignoreChannels.begin(), ignoreChannels.end(), channels[j]);
      if (it != fIgnoreChannels.end()) {
        channels.erase(channels.begin() + j);
      }
    }
  }

  return channelPairing;

} // icarus::trigger::LVDSgates::removeChannels()


//------------------------------------------------------------------------------
auto icarus::trigger::LVDSgates::discardChannels(
  std::vector<TrackedTriggerGate_t> const& gates,
  std::vector<raw::Channel_t> const& pairing
) const -> TrackedTriggerGate_t {
  TrackedTriggerGate_t gate;
  for (raw::Channel_t channel: pairing) gate.gate().addChannel(channel);
  for (auto const& srcGate: gates) gate.tracking().add(srcGate.tracking());
  return gate;
} // icarus::trigger::LVDSgates::combineChannels()


//------------------------------------------------------------------------------
auto icarus::trigger::LVDSgates::selectChannel(
  std::vector<TrackedTriggerGate_t> const& gates,
  std::vector<raw::Channel_t> const& pairing,
  std::size_t chosenIndex
) const -> TrackedTriggerGate_t {
  // requiring that gates are at an index matching their channel number
  TrackedTriggerGate_t gate { gates[pairing[chosenIndex]] };
  for (raw::Channel_t channel: pairing) {
    mf::LogTrace(fLogCategory) << "Input:  " << gates[channel].gate();
    gate.gate().addChannel(channel);
  }
  mf::LogTrace(fLogCategory) << "Output: " << gate.gate();
  for (auto const& srcGate: gates) gate.tracking().add(srcGate.tracking());
  return gate;
} // icarus::trigger::LVDSgates::selectChannel()


//------------------------------------------------------------------------------
template <typename Op>
auto icarus::trigger::LVDSgates::binaryCombineChannel(
  std::vector<TrackedTriggerGate_t> const& gates,
  std::vector<raw::Channel_t> const& pairing,
  Op combine
) const -> TrackedTriggerGate_t {
  if (pairing.empty()) return discardChannels(gates, pairing);

#if 1
  
  auto byIndex = [&gates](std::size_t index) -> TrackedTriggerGate_t const&
    { return gates[index]; };
  
# if 0
  // C++20: the following loses the debug output:
  
  return icarus::trigger::OpGateColl
    (combine, pairing | std::ranges::views::transform(byIndex));
# else
  
  return icarus::trigger::OpGateColl
    (combine, util::make_transformed_span(pairing, byIndex));
  
# endif // 0
  
#else
  auto iChannel = pairing.begin();
  auto cend = pairing.end();

  // requiring that gates are at an index matching their channel number
  TrackedTriggerGate_t gate { gates[*iChannel] };

  mf::LogTrace(fLogCategory) << "Input:  " << gates[*iChannel].gate();
  while (++iChannel != cend) {

    mf::LogTrace(fLogCategory) << "Input:  " << gates[*iChannel].gate();
    gate = icarus::trigger::OpGates(combine, gate, gates[*iChannel]);

  } // while
  mf::LogTrace(fLogCategory) << "Output: " << gate.gate();

  return gate;
#endif // 0
} // icarus::trigger::LVDSgates::binaryCombineChannel()


//------------------------------------------------------------------------------
unsigned int icarus::trigger::LVDSgates::checkPairings() const {

  //
  // collect all the errors first
  //
  auto const& geom = *(lar::providerFrom<geo::Geometry>());


  // collect configured and duplicate channel numbers
  unsigned int nEmptyGroups = 0U;
  std::set<raw::Channel_t> configuredChannels, duplicateChannels;
  for (auto const& pairing: fChannelPairing) {
    if (pairing.empty()) {
      ++nEmptyGroups;
      continue;
    }
    for (raw::Channel_t channel: pairing) {
      if (!configuredChannels.insert(channel).second)
        duplicateChannels.insert(channel);
    } // for channel
  } // for pairing

  // collect invalid channel numbers
  std::vector<raw::Channel_t> invalidChannels;
  for (raw::Channel_t channel: configuredChannels) {
    if (!geom.IsValidOpChannel(channel)) invalidChannels.push_back(channel);
  } // for configured channels

  // collect missing channels
  auto const endChannel = geom.MaxOpChannel() + 1U;
  std::vector<raw::Channel_t> missingChannels;
  for (auto channel: util::counter<raw::Channel_t>(endChannel)) {
    if (!geom.IsValidOpChannel(channel)) continue;
    if (configuredChannels.count(channel) > 0U) continue;
    missingChannels.push_back(channel);
  } // for channel

  // no error? good we are
  if (invalidChannels.empty()
    && missingChannels.empty()
    && duplicateChannels.empty()
    && (nEmptyGroups == 0U)
  ) {
    return configuredChannels.size();
  }

  // error? good we are not
  art::Exception e(art::errors::Configuration);
  e << "The pairing configuration of " << configuredChannels.size()
    << " optical channels in `ChannelPairing` presents the following errors.\n";
  if (!invalidChannels.empty()) {
    auto iChannel = invalidChannels.cbegin();
    auto const cend = invalidChannels.cend();
    e
      << invalidChannels.size()
      << " channel numbers do not represent valid channels: " << *iChannel;
    while (++iChannel != cend) e << ", " << *iChannel;
    e << ".\n";
  } // if invalid channels

  if (!missingChannels.empty()) {
    auto iChannel = missingChannels.cbegin();
    auto const cend = missingChannels.cend();
    e
      << missingChannels.size()
      << " channels do not appear in the configuration: " << *iChannel;
    while (++iChannel != cend) e << ", " << *iChannel;
    e << ".\n";
  } // if missing channels

  if (!duplicateChannels.empty()) {
    auto iChannel = duplicateChannels.cbegin();
    auto const cend = duplicateChannels.cend();
    e
      << duplicateChannels.size()
      << " channels are specified more than once: " << *iChannel;
    while (++iChannel != cend) e << ", " << *iChannel;
    e << ".\n";
  } // if duplicate channels
  
  if (nEmptyGroups > 0U) {
    e << nEmptyGroups << " channel groups are empty.";
  }
  
  throw e;
} // icarus::trigger::LVDSgates::checkPairings()


//------------------------------------------------------------------------------
void icarus::trigger::LVDSgates::checkInput
  (std::vector<TrackedTriggerGate_t> const& gates) const
{
  /*
   * check that 
   *  * each gate has a single channel
   *  * each gate is in the position matching its channel number
   */
  
  using icarus::trigger::gatesIn;
  
  for (auto const& [ iGate, gate ]: util::enumerate(gatesIn(gates))) {
    // gate is a OpticalTriggerGateData_t
    
    auto const expectedChannel = static_cast<raw::Channel_t>(iGate);
    
    if (gate.nChannels() != 1U) {
      cet::exception e("LVDSgates");
      e << "Requirement failure: gate #" << iGate << " includes ";
      if (!gate.hasChannels())
        e << "no channels";
      else {
        auto const& channels = gate.channels();
        auto iChannel = channels.begin();
        auto const cend = channels.end();
        e << gate.nChannels() << " channels (" << *iChannel;
        while (++iChannel != cend) e << ", " << *iChannel;
        e << ")";
      }
      e << " should have exactly one channel (" << expectedChannel << ").\n";
      throw e;
    }
    
    if (gate.channel() != expectedChannel) {
      throw cet::exception("LVDSgates")
        << "Requirement failure: gate #" << iGate << " gate covers channel "
        << gate.channel() << ", should cover " << expectedChannel
        << " instead.\n";
    }
    
  } // for
  
  MF_LOG_TRACE(fLogCategory)
    << "LVDSgates[" << moduleDescription().moduleLabel() << "]: input checked.";
  
} // icarus::trigger::LVDSgates::checkInput()


//------------------------------------------------------------------------------
void icarus::trigger::LVDSgates::UpdateWaveformMap(
  icarus::trigger::OpDetWaveformMetaDataProductMap_t& map,
  art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta> const& assns
) {
  
  for (art::Ptr<sbn::OpDetWaveformMeta> const& wave: util::get_elements<1U>(assns))
    map.emplace(wave.get(), wave);
  
} // icarus::trigger::LVDSgates::UpdateWaveformMap()


//------------------------------------------------------------------------------
auto icarus::trigger::LVDSgates::ReadTriggerGates(
  art::Event const& event,
  art::InputTag const& dataTag,
  icarus::trigger::OpDetWaveformMetaDataProductMap_t& waveformMap
) -> std::vector<TrackedTriggerGate_t> {
  using icarus::trigger::OpticalTriggerGateData_t;
  
  auto const& assns =
    event.getProduct<art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>
    (dataTag);
  
  UpdateWaveformMap(waveformMap, assns);
  
  return icarus::trigger::FillTriggerGates
    (event.getProduct<std::vector<OpticalTriggerGateData_t>>(dataTag), assns);
  
} // icarus::trigger::LVDSgates::ReadTriggerGates()


//------------------------------------------------------------------------------
art::InputTag icarus::trigger::LVDSgates::makeTag
  (std::string const& thresholdStr, std::string const& defModule)
{
  auto const isNumber = [pattern=std::regex{ "[+-]?[0-9]+" }]
    (std::string const& s) -> bool
    { return std::regex_match(s, pattern); };
  return 
    ((thresholdStr.find(':') != std::string::npos) || !isNumber(thresholdStr))
    ? art::InputTag{ thresholdStr }
    : defModule.empty()
      ? throw (art::Exception(art::errors::Configuration)
        << "No default module label (`TriggerGatesTag`) specified"
           ", and it's needed for threshold '"
        << thresholdStr << "'.\n")
      : art::InputTag{ defModule, thresholdStr }
    ;
} // icarus::trigger::LVDSgates::makeTag()


//------------------------------------------------------------------------------
std::string icarus::trigger::LVDSgates::makeOutputInstanceName
  (art::InputTag const& inputTag, std::string const& defModule)
{
  return (inputTag.label() == defModule)
    ? inputTag.instance()
    : inputTag.instance().empty()
      ? inputTag.label(): inputTag.label() + inputTag.instance()
    ;
} // icarus::trigger::LVDSgates::makeTag()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::LVDSgates)


//------------------------------------------------------------------------------
