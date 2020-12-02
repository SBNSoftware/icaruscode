/**
 * @file   LVDSgates_module.cc
 * @brief  Combines discriminated PMT outputs into V1730 LVDS gates.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 13, 2019
 */

// ICARUS libraries
#include "sbnobj/ICARUS/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/Utilities/MultipleChoiceSelection.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/get_elements.h" // util::get_elements()
#include "lardataobj/RawData/OpDetWaveform.h"
// #include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::::static_assert_on<>

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"


// C/C++ standard libraries
#include <map>
#include <set>
#include <vector>
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
 *   and their associations to optical waveforms. One set per threshold.
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
 *   same as the input gates): sets of gates combined according to the
 *   configuration; one set per input threshold.
 * * `art::Assns<icarus::trigger::OpticalTriggerGateData_t, raw::OpDetWaveform>`
 *   (instance name: same as the input gates): associations between each
 *   produced gate and the optical waveforms providing the original data.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description LVDSgates`.
 * 
 * * `TriggerGatesTag` (string, default: `discrimopdaq`): name of the module
 *     instance which produced the discriminated waveforms; it must not include
 *     any instance name, as the instance names will be automatically added from
 *     `Thresholds` parameter.
 * * `Thresholds` (list of integers, mandatory): list of the discrimination
 *     thresholds to consider, in ADC counts. A data product containing a
 *     digital signal is read for each one of the thresholds, and the tag of the
 *     data product is expected to be the module label `TriggerGatesTag` with as
 *     instance name the value of the threshold (e.g. for a threshold of 6 ADC
 *     counts the data product tag might be `discrimopdaq:6`).
 * * `ChannelPairing` (list of integral number lists): channels to combine;
 *     each element of this list is itself a list of optical detector channel
 *     numbers. The channels within each group are combined according to
 *     `CombinationMode`, resulting in a single combined gate associated to
 *     all the channels in the group.
 *     All optical detector channels *must* appear in the configuration, and
 *     each channel can be combined only in one group. Groups must not be empty.
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
    
    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of trigger gate extraction module (no instance name)"),
      "discrimopdaq"
      };

    fhicl::Sequence<raw::ADC_Count_t> Thresholds {
      Name("Thresholds"),
      Comment("thresholds to consider [ADC counts]")
      };
    
    fhicl::Sequence<fhicl::Sequence<raw::Channel_t>> ChannelPairing {
      Name("ChannelPairing"),
      Comment
        ("grouping of optical detector channels (e.g.: [ [ 0, 2 ], [ 4, 6 ], [ 8 ], ... ])")
      };

    fhicl::Sequence<raw::Channel_t> IgnoreChannels {
      Name("IgnoreChannels"),
      Comment
        ("grouping of optical detector channels to ignore (e.g.: [ 54, 58, 67, 76, ... ])")
      };    

    fhicl::Atom<std::string> CombinationMode {
      Name("CombinationMode"),
      Comment("channel combination mode: " + CombinationModeSelector.optionListString())
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
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// ADC thresholds to read, and the input tag connected to their data.
  std::map<icarus::trigger::ADCCounts_t, art::InputTag> fADCthresholds;
  
  /// Pairing of optical detector channels.
  std::vector<std::vector<raw::Channel_t>> fChannelPairing;
  ComboMode fComboMode; ///< The operation used to combinate channels.

  std::vector<raw::Channel_t> fIgnoreChannels;
  
  std::string fLogCategory; ///< Message facility stream category for output.
  
  // --- END Configuration variables -------------------------------------------
  
  
  /// Completely processes the data for the specified threshold.
  void produceThreshold(
    art::Event& event,
    icarus::trigger::OpDetWaveformDataProductMap_t& waveformMap,
    icarus::trigger::ADCCounts_t const threshold,
    art::InputTag const& dataTag
    ) const;

  std::vector<std::vector<raw::Channel_t>> removeChannels(
    std::vector<std::vector<raw::Channel_t>> channelPairing,
    std::vector<raw::Channel_t> fIgnoreChannels
    ) const;
  
  
  // --- BEGIN -- Combination modes --------------------------------------------
  /// Combination function for the `disable` mode.
  icarus::trigger::MultiChannelOpticalTriggerGate discardChannels(
    std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const&,
    std::vector<raw::Channel_t> const& pairing
    ) const;
  
  /// Combination function for the `input1` and `input2` modes.
  icarus::trigger::MultiChannelOpticalTriggerGate selectChannel(
    std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const&,
    std::vector<raw::Channel_t> const& pairing, std::size_t chosenIndex
    ) const;
  
  /// Combination function for the `AND` and `OR` modes.
  icarus::trigger::MultiChannelOpticalTriggerGate binaryCombineChannel(
    std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const& gates,
    std::vector<raw::Channel_t> const& pairing,
    BinaryCombinationFunc_t combine
    ) const;
  
  /// Performs the combination of a group of channels.
  icarus::trigger::MultiChannelOpticalTriggerGate combineChannels(
    std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const& gates,
    std::vector<raw::Channel_t> const& pairing,
    ComboMode comboMode
    ) const;
  // --- END -- Combination modes ----------------------------------------------
  
  
  
  // --- BEGIN -- Checks -------------------------------------------------------
  /// Checks the pairing configuration, throws on error.
  /// Returns the number of configured channels.
  unsigned int checkPairings() const;

  /// Input requirement check. Throws if requirements fail.
  void checkInput(
    std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const&  gates
    ) const;
  
  // --- END -- Checks ---------------------------------------------------------
  
  
  /// Adds the associated waveforms into the map.
  static void UpdateWaveformMap(
    icarus::trigger::OpDetWaveformDataProductMap_t& map,
    art::Assns<TriggerGateData_t, raw::OpDetWaveform> const& assns
    );
  
  /// Assembles trigger gates from `dataTag` data products in `event`.
  static std::vector<icarus::trigger::SingleChannelOpticalTriggerGate>
  ReadTriggerGates(
    art::Event const& event, art::InputTag const& dataTag,
    icarus::trigger::OpDetWaveformDataProductMap_t& waveformMap
    );

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
  , fLogCategory(config().LogCategory())
{
  //
  // more complex parameter parsing
  //
  std::string const discrModuleLabel = config().TriggerGatesTag();
  for (raw::ADC_Count_t threshold: config().Thresholds()) {
    fADCthresholds[icarus::trigger::ADCCounts_t{threshold}]
      = art::InputTag{ discrModuleLabel, util::to_string(threshold) };
  }
  
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
    log << "\nignored " << fIgnoreChannels.size() << " channels.";
    log << "\nConfigured " << fADCthresholds.size() << " thresholds:";
    for (auto const& [ threshold, dataTag ]: fADCthresholds)
      log << "\n * " << threshold << " ADC (from '" << dataTag.encode() << "')";
  } // local block

  
  using icarus::trigger::OpticalTriggerGateData_t; // for convenience
  
  //
  // input data declaration
  //
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    consumes<std::vector<OpticalTriggerGateData_t>>(inputDataTag);
    consumes<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (inputDataTag);
  } // for
  
  
  //
  // output data declaration
  //
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    produces<std::vector<OpticalTriggerGateData_t>>(inputDataTag.instance());
    produces<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (inputDataTag.instance());
  } // for
  
  
} // icarus::trigger::LVDSgates::LVDSgates()


//------------------------------------------------------------------------------
void icarus::trigger::LVDSgates::produce(art::Event& event) {
  
  icarus::trigger::OpDetWaveformDataProductMap_t waveformMap;
  
  for (auto const& [ threshold, dataTag ]: fADCthresholds) {
    
    produceThreshold(event, waveformMap, threshold, dataTag);
    
  } // for all thresholds
  
} // icarus::trigger::LVDSgates::produce()

std::vector<std::vector<raw::Channel_t>> icarus::trigger::LVDSgates::removeChannels(
  std::vector<std::vector<raw::Channel_t>> channelPairing,
  std::vector<raw::Channel_t> fIgnoreChannels
  ) const {

  for (unsigned int i = 0; i < channelPairing.size(); i++) {
    for (unsigned int j = channelPairing[i].size(); j-- > 0; ) {
      std::vector<raw::Channel_t>::iterator it = find(fIgnoreChannels.begin(), fIgnoreChannels.end(), channelPairing[i][j]);
      if (it != fIgnoreChannels.end()) {
        channelPairing[i].erase(channelPairing[i].begin() + j);
      }
    }
  }

  return channelPairing;

}

//------------------------------------------------------------------------------
void icarus::trigger::LVDSgates::produceThreshold(
  art::Event& event,
  icarus::trigger::OpDetWaveformDataProductMap_t& waveformMap,
  icarus::trigger::ADCCounts_t const threshold,
  art::InputTag const& dataTag
) const {
  
  mf::LogDebug(fLogCategory)
    << "Processing threshold " << threshold
    << " from '" << dataTag.encode() << "'";

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience
  
  std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const& gates
    = ReadTriggerGates(event, dataTag, waveformMap);
  
  checkInput(gates);
  
  std::vector<icarus::trigger::MultiChannelOpticalTriggerGate> combinedGates;
  std::vector<std::vector<raw::Channel_t>> cleanedChannels = removeChannels(fChannelPairing, fIgnoreChannels);

  for (std::vector<raw::Channel_t> const& pairing: cleanedChannels) {
    if (pairing.empty()) continue; // ???
    
    combinedGates.push_back(combineChannels(gates, pairing, fComboMode));
    
  } // for
  
  // transform the data; after this line, `gates` is not usable any more
  art::PtrMaker<OpticalTriggerGateData_t> const makeGatePtr
    (event, dataTag.instance());
  auto&& [ outputGates, outputAssns ]
    = icarus::trigger::transformIntoOpticalTriggerGate
    (combinedGates, makeGatePtr, waveformMap);
  
  mf::LogTrace(fLogCategory)
    << "Threshold " << threshold << " ('" << dataTag.encode() << "'): "
    << gates.size() << " input channels, "
    << outputGates.size() << " output channels (+"
    << outputAssns.size() << " associations to waveforms) into '"
    << moduleDescription().moduleLabel() << ":" << dataTag.instance() << "'"
    ;

  event.put(
    std::make_unique<std::vector<OpticalTriggerGateData_t>>
      (std::move(outputGates)),
    dataTag.instance()
    );
  event.put(
    std::make_unique<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (std::move(outputAssns)),
    dataTag.instance()
    );
  
} // icarus::trigger::LVDSgates::produceThreshold()


//------------------------------------------------------------------------------
icarus::trigger::MultiChannelOpticalTriggerGate
icarus::trigger::LVDSgates::combineChannels(
  std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const& gates,
  std::vector<raw::Channel_t> const& pairing,
  ComboMode comboMode
) const {
  switch (comboMode) {
    case ComboMode::disable:
      return discardChannels(gates, pairing);
    case ComboMode::AND:
      return binaryCombineChannel(gates, pairing, &OpticalTriggerGate::Min);
    case ComboMode::OR:
      return binaryCombineChannel(gates, pairing, &OpticalTriggerGate::Max);
    case ComboMode::Input1:
      return selectChannel(gates, pairing, 0U);
    case ComboMode::Input2:
      return selectChannel(gates, pairing, 1U);
    default:
      throw art::Exception(art::errors::LogicError)
        << "Unexpected combination mode (#" << static_cast<int>(comboMode)
        << ").\n";
  } // switch(comboMode)
} // icarus::trigger::LVDSgates::discardChannels()


//------------------------------------------------------------------------------
icarus::trigger::MultiChannelOpticalTriggerGate
icarus::trigger::LVDSgates::discardChannels(
  std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const&,
  std::vector<raw::Channel_t> const& pairing
) const {
  icarus::trigger::MultiChannelOpticalTriggerGate gate;
  for (raw::Channel_t channel: pairing) gate.addChannel(channel);
  return gate;
} // icarus::trigger::LVDSgates::combineChannels()


//------------------------------------------------------------------------------
icarus::trigger::MultiChannelOpticalTriggerGate
icarus::trigger::LVDSgates::selectChannel(
  std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const& gates,
  std::vector<raw::Channel_t> const& pairing,
  std::size_t chosenIndex
) const {
  // requiring that gates are at an index matching their channel number
  icarus::trigger::MultiChannelOpticalTriggerGate gate {
    static_cast<icarus::trigger::OpticalTriggerGate const&>
      (gates[pairing[chosenIndex]])
    };
  for (raw::Channel_t channel: pairing) {
    mf::LogTrace(fLogCategory) << "Input:  " << gates[channel];
    gate.addChannel(channel);
  }
  mf::LogTrace(fLogCategory) << "Output: " << gate;
  return gate;
} // icarus::trigger::LVDSgates::selectChannel()


//------------------------------------------------------------------------------
icarus::trigger::MultiChannelOpticalTriggerGate
icarus::trigger::LVDSgates::binaryCombineChannel(
  std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const& gates,
  std::vector<raw::Channel_t> const& pairing,
  BinaryCombinationFunc_t combine
  ) const
{
  if (pairing.empty()) return discardChannels(gates, pairing);

  auto iChannel = pairing.begin();
  auto cend = pairing.end();

  // requiring that gates are at an index matching their channel number
  icarus::trigger::MultiChannelOpticalTriggerGate gate
    {static_cast<icarus::trigger::OpticalTriggerGate const&>(gates[*iChannel])};

  mf::LogTrace(fLogCategory) << "Input:  " << gates[*iChannel];
  while (++iChannel != cend) {

    mf::LogTrace(fLogCategory) << "Input:  " << gates[*iChannel];
    (gate.*combine)(gates[*iChannel]);

  } // while
  mf::LogTrace(fLogCategory) << "Output: " << gate;

  return gate;
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
  (std::vector<icarus::trigger::SingleChannelOpticalTriggerGate> const&  gates)
  const
{
  /*
   * check that 
   *  * each gate has a single channel
   *  * each gate is in the position matching its channel number
   */
  for (auto const& [ iGate, gate ]: util::enumerate(gates)) {
    
    auto const expectedChannel = static_cast<raw::Channel_t>(iGate);
    
    if (gate.nChannels() != 1U) {
      cet::exception e("LVDSgates");
      e << "Requirement failure: gate #" << iGate << " includes ";
      if (!gate.hasChannels())
        e << "no channels";
      else {
        // we are cheating here, because `SingleChannelOpticalTriggerGate`
        // actually does not want more than one channel;
        // having them is sort of a logic error
        auto const& channels = gate.OpticalTriggerGate::channels();
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
  icarus::trigger::OpDetWaveformDataProductMap_t& map,
  art::Assns<TriggerGateData_t, raw::OpDetWaveform> const& assns
) {
  
  for (art::Ptr<raw::OpDetWaveform> const& wave: util::get_elements<1U>(assns))
    map.emplace(wave.get(), wave);
  
} // icarus::trigger::LVDSgates::UpdateWaveformMap()


//------------------------------------------------------------------------------
std::vector<icarus::trigger::SingleChannelOpticalTriggerGate>
icarus::trigger::LVDSgates::ReadTriggerGates(
  art::Event const& event,
  art::InputTag const& dataTag,
  icarus::trigger::OpDetWaveformDataProductMap_t& waveformMap
) {
  using icarus::trigger::OpticalTriggerGateData_t;
  
  auto const& assns =
    *(event.getValidHandle<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>(dataTag));
  
  UpdateWaveformMap(waveformMap, assns);
  
  return icarus::trigger::FillTriggerGates<icarus::trigger::SingleChannelOpticalTriggerGate>
    (*(event.getValidHandle<std::vector<OpticalTriggerGateData_t>>(dataTag)), assns);
  
} // icarus::trigger::LVDSgates::ReadTriggerGates()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::LVDSgates)


//------------------------------------------------------------------------------
