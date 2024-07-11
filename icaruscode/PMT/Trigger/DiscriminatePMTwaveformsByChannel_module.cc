/**
 * @file   DiscriminatePMTwaveformsByChannel_module.cc
 * @brief  Module producing discriminated waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 25, 2021
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h"
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/TriggerGateData.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h"
#include "icaruscode/Utilities/DataProductPointerMap.h"
#include "icarusalg/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/Utilities/quantities_fhicl.h" // for ADCCounts_t parameters
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <sstream>
#include <map>
#include <optional>
#include <vector>
#include <string>


//------------------------------------------------------------------------------
namespace icarus::trigger { class DiscriminatePMTwaveformsByChannel; }

/**
 * @brief Produces discriminated optical waveforms.
 * 
 * This module produces a "discriminated waveform" for each optical detector
 * channel and for a single set of discrimination thresholds.
 * The thresholds are specified one per channel.
 * 
 * 
 * Output data products
 * =====================
 * 
 * The output data products can be assigned an instance name
 * (`OutputInstanceName` configuration parameter) if desired. They are:
 * 
 * * `std::vector<icarus::trigger::OpticalTriggerGate::GateData_t>`,
 *   with one entry per optical channel, including _all_ optical channels
 *   (even when no signal above threshold is ever present); each gate object
 *   has as index in the vector the number of the optical channel;
 * * an _art_ association, one trigger gate data to many optical waveforms,
 *   of each trigger gate data (as above) with all the optical waveforms
 *   which contributed to it; the order of the association pairs is the same
 *   as the one of the trigger gate data, i.e. by channel, and some trigger
 *   gate data entries may be not associated to any waveform and therefore
 *   they can not be listed in the association; association pairs within
 *   the same optical channel are sorted by optical waveform timestamp;
 *   the type of the association is
 *   `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, raw::OpDetWaveform>`;
 *   similarly, if `SavePMTcoverage` parameter was set, also an association
 *   `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, sbn::OpDetWaveformMeta>`
 *   is produced.
 * * `std::vector<sbn::OpDetWaveformMeta>` parallel to the input optical
 *   detector waveforms, if `SavePMTcoverage` parameter is set; each entry in
 *   the collection is matched with the corresponding one in the input waveform
 *   collection; an association
 *   `art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>` is also produced.
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<raw::OpDetWaveform>`: a single waveform for each recorded
 *      optical detector activity; the activity belongs to a single channel, but
 *      there may be multiple waveforms on the same channel. The time stamp is
 *      expected to be from the
 *      @ref DetectorClocksElectronicsTime "electronics time scale"
 *      and therefore expressed in microseconds.
 * 
 * 
 * Service requirements
 * ---------------------
 * 
 * The following services are _required_:
 * 
 * * an implementation of `detinfo::DetectorClocksService`
 * * _art_ message facility
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * There are several (too many?) ways to specify the discrimination thresholds.
 * In the end, all PMT channels must be assigned a threshold: if at any time,
 * on any event, a channel is found whose threshold was not set, an exception
 * is thrown.
 * All threshold specifications are relative to the baseline. There are
 * currently four levels of threshold setting:
 * 
 * 1. if `NChannels` (positive integer) and `DefaultThreshold` (integral number)
 *     are both specified, the first `NChannels` PMT channels are assigned the
 *     `DefaultThreshold`; note that they both need to be specified at the same
 *     time, and when they are, channels with ID equal or larger than
 *     `NChannels` are still not configured, i.e. they will not be assigned the
 *     `DefaultThreshold` and if encountered an exception will be thrown;
 * 2. `ThresholdList` (list of integers) will define a threshold for each
 *     channel starting with ID `0` on; thresholds already set with the previous
 *     step will be overridden;
 * 3. `ThresholdsFromPMTconfig` (input tag) points to a data product with
 *     a PMT configuration (`sbn::PMTconfiguration`) which contains
 *     information on each channel, including the thresholds; thresholds are set
 *     for each channel from the value returned by
 *     `sbn::V1730channelConfiguration::relativeThreshold()` for that channel;
 *     only channel configuration entries that have a offline channel ID
 *     (`sbn::V1730channelConfiguration::hasChannelID()`) are configured in this
 *     way; thresholds already set with the previous steps will be overridden;
 *     when using this parameter, `NChannels` is also mandatory and channels
 *     with ID greater or equal than `NChannels` will trigger an exception.
 * 4. `Thresholds` (list of tables): in each entry, a threshold is specified for
 *     a single channel; the format of the entry is:
 *     * `Channel` (integral, mandatory): channel number;
 *     * `Threshold` (integral, mandatory): the threshold for that channel.
 *     .
 *     Thresholds already set with the previous steps will be overridden.
 * 
 * All these are optional (but if the configuration has none, not much can be
 * done since there will be no channel configured).
 * 
 * A terse description of the remaining parameters follows.
 * 
 * * `OpticalWaveforms` (input tag): the data product containing all optical
 *   detector waveforms.
 * * `Baselines` (input tag): the data product containing one baseline per
 *   waveform in data product (from `OpticalWaveforms`).
 * * `NChannels` (integral): if specified, only channels with ID from `0` to
 *   `NChannels - 1` are processed; waveforms from channels with higher ID will
 *   be ignored. If `DefaultThreshold` or `ThresholdsFromConfiguration`
 *   parameters are specified (see above), then `NChannels` is mandatory;
 *   otherwise, its value is deduced from the highest channel configured with
 *   one of the other threshold parameters.
 * * `TriggerGateBuilder` (tool configuration): template configuration of the
 *   _art_ tool used to discriminate the optional waveforms; the tool interface
 *   is `icarus::trigger::TriggerGateBuilder`. Note that the threshold settings
 *   of this configuration is going to be ignored.
 * * `SavePMTcoverage` (flag, default: `true`): also produces a collection of
 *   `sbn::OpDetWaveformMeta` representing each of the input waveforms; trigger
 *   tools can use this information in place for the more space-hungry waveforms
 *   for further processing.
 * * `OutputCategory` (string, default: `"DiscriminatePMTwaveformsByChannel"`):
 *   label for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service.
 *   configuration.
 * 
 * 
 */
class icarus::trigger::DiscriminatePMTwaveformsByChannel: public art::EDProducer
{
  
    public:
  
  /// The type of data produced for a single channel.
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGate::GateData_t;
  
  /// Collection of settings for one channel.
  struct ChannelInfo_t {
    raw::Channel_t channel;
    std::optional<raw::ADC_Count_t> baseline;
    std::optional<raw::ADC_Count_t> threshold;
    
    bool operator== (ChannelInfo_t const& other) const;
    bool operator< (ChannelInfo_t const& other) const;
  }; // ChannelInfo_t
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  
  struct ChannelConfig {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<raw::Channel_t> Channel {
      Name("Channel"),
      Comment("off-line channel ID")
      };
    
    fhicl::OptionalAtom<raw::ADC_Count_t> Baseline {
      Name("Baseline"),
      Comment("the baseline for this channel [ADC]")
      };
    
    fhicl::OptionalAtom<raw::ADC_Count_t> Threshold {
      Name("Threshold"),
      Comment("threshold relative to the baseline for this channel [ADC]")
      };
    
  }; // ChannelConfig
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> OpticalWaveforms {
      Name("OpticalWaveforms"),
      Comment("label of input digitized optical waveform data product"),
      "opdaq" // tradition demands
      };
    
    fhicl::OptionalAtom<art::InputTag> Baselines {
      Name("Baselines"),
      Comment("label of input waveform baselines (parallel to the waveforms)")
      };
    
    fhicl::OptionalAtom<float> Baseline {
      Name("Baseline"),
      Comment("constant baseline for all waveforms, in ADC counts")
      };
    
    fhicl::DelegatedParameter TriggerGateBuilder_ {
      Name("TriggerGateBuilder"),
      Comment
        ("parameters for generating trigger gates from optical channel output")
      };
    
    fhicl::OptionalAtom<raw::ADC_Count_t> DefaultThreshold {
      Name("DefaultThreshold"),
      Comment("threshold for all channels up to NChannel [ADC]")
      };
    
    fhicl::Sequence<raw::ADC_Count_t> ThresholdList {
      Name("ThresholdList"),
      Comment("relative thresholds for each channel from #0 on [ADC]"),
      std::vector<raw::ADC_Count_t>{}
      };
    
    fhicl::OptionalAtom<art::InputTag> ThresholdsFromPMTconfig {
      Name("ThresholdsFromPMTconfig"),
      Comment("read thresholds from this PMT configuration run by run")
      };
    
    fhicl::Sequence<fhicl::TableAs<ChannelInfo_t, ChannelConfig>> Thresholds {
      Name("Thresholds"),
      Comment
        ("thresholds for the specified channels; overrides `ThresholdList`"),
      std::vector<ChannelInfo_t>{}
      };
    
    
    fhicl::OptionalAtom<unsigned int> NChannels {
      Name("NChannels"),
      Comment("minimum number of channels to provide (default: all)")
      };
    
    fhicl::Atom<bool> SavePMTcoverage {
      Name("SavePMTcoverage"),
      Comment("write also a sbn::OpDetWaveformMeta collection"),
      true
      };
    
    fhicl::Atom<std::string> OutputInstanceName {
      Name("OutputInstanceName"),
      Comment("instance name for the output products (none by default)"),
      "" // default
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("tag of the module output to console via message facility"),
      "DiscriminatePMTwaveformsByChannel"
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit DiscriminatePMTwaveformsByChannel(Parameters const& config);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  /// Prepares the plots to be filled.
  /// Prepares the plots to be filled.
  virtual void beginRun(art::Run& run) override;
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  using ADCCounts_t = icarus::trigger::ADCCounts_t; // alias
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fOpDetWaveformTag; ///< Input optical waveform tag.
  
  ///< Input waveform baseline tag.
  std::optional<art::InputTag> const fBaselineTag;
  
  std::optional<float> const fBaseline; ///< A constant baseline level.

  unsigned int const fNOpDetChannels; ///< Number of optical detector channels.
  
  /// Default threshold.
  std::optional<raw::ADC_Count_t> const fDefaultThreshold;
  
  /// The thresholds configured by a sequence parameter.
  std::vector<raw::ADC_Count_t> const fThresholdList;
  
  /// Tag of `PMTconfiguration` data product to read thresholds from.
  std::optional<art::InputTag> const fThresholdsFromPMTconfig;
  
  /// A map of thresholds: channel ID -> relative threshold value [ADC]
  std::vector<ChannelInfo_t> const fChannelInfos;
  
  std::string const fOutputInstanceName; ///< Instance name for output.
  
  bool const fSavePMTcoverage; ///< Whether to save also `sbn::OpDetWaveformMeta`.
  
  /// Category name for the console output stream.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  /// PMT readout configuration for the current run.
  std::optional<sbn::PMTconfiguration> fPMTconfig;
  
  /// Current baseline list from run (one per channel).
  std::vector<icarus::WaveformBaseline> fBaselinesRun;
  
  enum class Source_t
    { Default, List, Unset, Run, Event, PMTconfig, ChannelSpec };
  
  /// Information on the current setting of a threshold.
  struct ADCvalueSetting_t {
    Source_t source = Source_t::Unset; ///< Where the value comes from.
    icarus::trigger::ADCCounts_t value {}; ///< Threshold value.
    bool isUnset() const noexcept { return source == Source_t::Unset; }
  }; // ADCvalueSetting_t
  
  /// A map of current thresholds: channel ID -> relative threshold value [ADC]
  std::vector<ADCvalueSetting_t> fCurrentThresholds;
  
  /// A map of current baselines: channel ID -> baseline [ADC]
  std::vector<ADCvalueSetting_t> fCurrentBaselines;
  
  /// Algorithms to simulate trigger gates out of optical channel output.
  std::unique_ptr<icarus::trigger::TriggerGateBuilder> fTriggerGateBuilder;
  
  /// Updates `fCurrentThresholds` with the current status of the object.
  void refreshCurrentThresholds();
  
  /// Updates `fCurrentBaselines` with the current status of the object.
  void refreshCurrentBaselines();
  
  /// Prints the current thresholds to the message logger.
  template <typename Logger = mf::LogInfo>
  void printCurrentThresholdsAndBaselines(Logger& log) const;
  
  /// Prints the current thresholds to the message logger.
  void printCurrentThresholdsAndBaselines() const
    { mf::LogInfo l{ fLogCategory }; printCurrentThresholdsAndBaselines(l); }
  
  // --- END Algorithms --------------------------------------------------------
  
  
  /**
   * @brief Creates a collection with one gate per channel and no gaps.
   * @param gates the trigger gates to be put in the collection
   * @return a collection with with one gate per channel and no gaps
   * 
   * The returned collection includes one gate per channel, and at least
   * `fNOpDetChannels` channels (more if there are channels with higher number
   * than that).
   * For each channel, a gate is set in the item of the collection with index
   * the channel ID (i.e. the first gate in the collection will be the one for
   * channel `0`, the next the one for channel `1` and so on). It is assumed
   * that there is only at most one gate per channel.
   * The `gates` specified in the argument are _moved_ to the proper item in the
   * returned collection. The other gates are assigned the proper channel number
   * but are closed gates and associated with no waveform.
   */
  /// Fetches the optical detector count from geometry unless in `param`.
  static unsigned int getNOpDetChannels
    (fhicl::OptionalAtom<unsigned int> const& param);
  
  static std::vector<ChannelInfo_t> readChannelInfoSpecs
    (std::vector<ChannelInfo_t> const& config);
  
}; // icarus::trigger::DiscriminatePMTwaveformsByChannel


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::DiscriminatePMTwaveformsByChannel::ChannelInfo_t
//------------------------------------------------------------------------------
bool icarus::trigger::DiscriminatePMTwaveformsByChannel::ChannelInfo_t::operator==
  (ChannelInfo_t const& other) const
{
  return (channel == other.channel)
    && (baseline == other.baseline)
    && (threshold == other.threshold)
    ;
} // icarus::...::DiscriminatePMTwaveformsByChannel::ChannelInfo_t::operator==()


//------------------------------------------------------------------------------
namespace icarus::trigger {
  // conversion function: ChannelConfig -> ChannelInfo_t
  DiscriminatePMTwaveformsByChannel::ChannelInfo_t convert
    (DiscriminatePMTwaveformsByChannel::ChannelConfig const& config)
  {
    return DiscriminatePMTwaveformsByChannel::ChannelInfo_t{
      config.Channel(),                   // channel
      util::fhicl::getOptionalValue(config.Baseline),  // baseline
      util::fhicl::getOptionalValue(config.Threshold)  // threshold
      };
  } // icarus::trigger::convert()
} // namespace icarus::trigger


//------------------------------------------------------------------------------
bool icarus::trigger::DiscriminatePMTwaveformsByChannel::ChannelInfo_t::operator<
  (ChannelInfo_t const& other) const
  { return channel < other.channel; }


//------------------------------------------------------------------------------
//--- icarus::trigger::DiscriminatePMTwaveformsByChannel
//------------------------------------------------------------------------------
icarus::trigger::DiscriminatePMTwaveformsByChannel::DiscriminatePMTwaveformsByChannel
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fOpDetWaveformTag(config().OpticalWaveforms())
  , fBaselineTag(util::fhicl::getOptionalValue(config().Baselines))
  , fBaseline(util::fhicl::getOptionalValue(config().Baseline))
  , fNOpDetChannels(getNOpDetChannels(config().NChannels))
  , fDefaultThreshold(util::fhicl::getOptionalValue(config().DefaultThreshold))
  , fThresholdList(config().ThresholdList())
  , fThresholdsFromPMTconfig
    (util::fhicl::getOptionalValue(config().ThresholdsFromPMTconfig))
  , fChannelInfos(readChannelInfoSpecs(config().Thresholds()))
  , fOutputInstanceName(config().OutputInstanceName())
  , fSavePMTcoverage(config().SavePMTcoverage())
  , fLogCategory(config().OutputCategory())
  , fTriggerGateBuilder
    (
      art::make_tool<icarus::trigger::TriggerGateBuilder>
        (config().TriggerGateBuilder_.get<fhicl::ParameterSet>())
    )
{
  
  if (config().DefaultThreshold.hasValue()
    || config().ThresholdsFromPMTconfig.hasValue()
  ) {
    if (!config().NChannels.hasValue()) {
      throw art::Exception(art::errors::Configuration)
        << "Configuration parameter `" << config().NChannels.name()
        << "` not specified: it is mandatory when either `"
        << config().DefaultThreshold.name() << "` or `"
        << config().ThresholdsFromPMTconfig.name() << "` are specified.\n";
    }
  }
  
  //
  // optional configuration parameters
  //
  if (bool(fBaseline) + bool(fBaselineTag) != 1)
  {
    throw art::Exception(art::errors::Configuration)
      << "Exactly one among options `" << config().Baselines.name() << "` ("
      << (fBaselineTag? ("'" + fBaselineTag->encode() + "'"): "not set")
      << ") and `" << config().Baseline.name() << "` ("
      << (fBaseline? ("'" + std::to_string(*fBaseline) + "'"): "not set")
      << ") must be specified.\n";
  }
  
  
  refreshCurrentThresholds();
  refreshCurrentBaselines();
  
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  
  if (fThresholdsFromPMTconfig)
    consumes<sbn::PMTconfiguration, art::InRun>(*fThresholdsFromPMTconfig);
  
  
  //
  // declaration of output
  //
  
  produces<std::vector<TriggerGateData_t>>(fOutputInstanceName);
  produces<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
    (fOutputInstanceName);
  if (fSavePMTcoverage) {
    produces<std::vector<sbn::OpDetWaveformMeta>>(fOutputInstanceName);
    produces<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
      (fOutputInstanceName);
    produces<art::Assns<sbn::OpDetWaveformMeta, TriggerGateData_t>>
      (fOutputInstanceName);
  }
  
  //
  // configuration dump
  //
  
  {
    mf::LogInfo log { fLogCategory };
    if (!fThresholdsFromPMTconfig) printCurrentThresholdsAndBaselines(log);
    log << "Algorithm configuration:\n";
    std::ostringstream sstr;
    fTriggerGateBuilder->dumpConfiguration(sstr, /* indent */ "  ");
    log << sstr.str().c_str() // C strings are treated special (and faster)
      << "\nDiscrimination thresholds will be overwritten with DAQ values.";
  }
  
} // icarus::trigger::DiscriminatePMTwaveformsByChannel::DiscriminatePMTwaveformsByChannel()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatePMTwaveformsByChannel::beginRun
  (art::Run& run)
{
  if (fThresholdsFromPMTconfig) {
    auto const& PMTconfig
      = run.getProduct<sbn::PMTconfiguration>(fThresholdsFromPMTconfig.value());
    
    if (!fPMTconfig || (fPMTconfig != PMTconfig)) {
      fPMTconfig = PMTconfig;
      refreshCurrentThresholds();
      printCurrentThresholdsAndBaselines();
    }
  } // if thresholds from config
  
} // icarus::trigger::DiscriminatePMTwaveformsByChannel::beginRun()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatePMTwaveformsByChannel::produce(art::Event& event) {
  
  //
  // set up the algorithm to create the trigger gates
  //
  detinfo::DetectorTimings const detTimings {
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };

  //
  // fetch input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  auto const& waveforms = *waveformHandle;
  
  // map address of waveform to art pointer to that waveform
  auto const& opDetWavePtrs
    = util::mapDataProductPointers(event, waveformHandle);
  
  //
  // retrieve the baseline information if provided event by event
  //
  std::vector<icarus::WaveformBaseline> baselines;
  if (fBaselineTag) {
    // read and assign from the data product; configured baselines are ignored
    
    baselines = event.getProduct<std::vector<icarus::WaveformBaseline>>
      (fBaselineTag.value());
    
  } // if baselines from event
  else {
    // assign baselines from configuration
    baselines.reserve(waveforms.size());
    for (auto const& waveform: waveforms) {
      auto const channelSlot
        = static_cast<std::size_t>(waveform.ChannelNumber());
      if (channelSlot >= fCurrentBaselines.size()) {
        throw cet::exception("DiscriminatePMTwaveformsByChannel")
          << "Waveforms on channel " << waveform.ChannelNumber()
          << " have no baseline configured.\n";
      }
      baselines.emplace_back(fCurrentBaselines[channelSlot].value.value());
    } // for all waveforms
  } // if ... else
  
  {
    mf::LogDebug log { fLogCategory };
    log << "Discrimination of " << waveforms.size() << " PMT waveforms from '"
      << fOpDetWaveformTag.encode() << "'";
    if (fBaselineTag)
      log << " with baselines from '" << fBaselineTag->encode() << "'";
  }
  
  
  // We use an algorithm designed for applying to all channels the same
  // thresholds, and many of them. But here we have a different threshold per
  // channel, and a single threshold in that.
  // So we'll reconfigure the algorithm for each channel, and process only
  // the waveforms (and all of them) on that channel
  // (a possible alternative is to set all the thresholds that appear in any
  // of the channels, and after processing pick the right one for each channel;
  // this normally works well because the threshold values are actually very few
  // (a nominal one, and a few exceptions) but it may not scale.
  
  // group the waveforms by channel; ignore the channels we'll not process:
  // channel ID -> all input waveforms and baselines on that channel
  std::vector<std::vector<icarus::trigger::WaveformWithBaseline>>
    waveformInfoPerChannel(fNOpDetChannels);

  for (auto const& [ waveform, baseline ]: util::zip(waveforms, baselines)) {
    
    auto const channelSlot = static_cast<std::size_t>(waveform.ChannelNumber());
    if (channelSlot >= waveformInfoPerChannel.size()) continue; // ignore
    
    waveformInfoPerChannel[channelSlot].emplace_back
      (&waveform, &baseline);
    
  } // for all channels
  
  assert(fCurrentThresholds.size() >= fNOpDetChannels);
  icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t gates;
  for (auto const& [ channelSlot, waveInfo ]:
    util::enumerate(waveformInfoPerChannel)
  ) {
    auto const channel = static_cast<raw::Channel_t>(channelSlot);
    if (waveInfo.empty()) {
      // gate associated to channel, always closed
      gates.emplace_back(icarus::trigger::OpticalTriggerGateData_t{ channel });
      continue;
    }
    
    // set the threshold
    auto const& threshold = fCurrentThresholds[channelSlot];
    if (threshold.isUnset()) {
      throw cet::exception("DiscriminatePMTwaveformsByChannel")
        << "No threshold set up for PMT channel #" << channel << ".\n";
    }
    fTriggerGateBuilder->resetup
      (detTimings, { ADCCounts_t::castFrom(threshold.value) });
    assert(fTriggerGateBuilder->nChannelThresholds() == 1U);
    mf::LogTrace(fLogCategory)
      << "Processing PMT channel #" << channel
      << " (first waveform baseline: " << waveInfo.front().baseline()
      << ") and threshold " << fTriggerGateBuilder->channelThresholds().front()
      ;
    
    // if there are waveforms, run the algorithm;
    // save the first collection, i.e. the first (and only) threshold;
    // extract the trigger gates from TriggerGates object immediately,
    // GateData_t is a collection of TrackedTriggerGate objects,
    // one per channel; that is, only one entry in this case
    icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t channelGates
      = std::move(fTriggerGateBuilder->build(waveInfo).front()).gates();
    assert(channelGates.size() == 1);
    gates.push_back(std::move(channelGates.front()));
    
  } // for all channels

  { // nameless block
    mf::LogTrace log(fLogCategory);
    log << "Trigger gates:\n";
    unsigned int nOpenGates = gates.size();
    for (auto const& gate: gates) if (gate.gate().alwaysClosed()) --nOpenGates;
    log << nOpenGates << "/" << gates.size() << " trigger gates";
    if (nOpenGates) {
      log << ":";
      for (auto const& gate: gates) {
        if (gate.gate().alwaysClosed()) continue;
        log << "\n  " << gate;
      }
    }
  } // nameless block
  
  //
  // prepare output
  //
  art::PtrMaker<TriggerGateData_t> const makeGatePtr
    { event, fOutputInstanceName };
  
  //
  // reformat the results for the threshold
  //
  auto [ data, assns ] = icarus::trigger::transformIntoOpticalTriggerGate
    (std::move(gates), makeGatePtr, opDetWavePtrs);

  if (fSavePMTcoverage) {
    
    art::PtrMaker<sbn::OpDetWaveformMeta> const makePMTinfoPtr
      { event, fOutputInstanceName };
    
    // PMT info and associations to input waveforms
    sbn::OpDetWaveformMetaMaker const makePMTinfo { detTimings };
    
    std::vector<sbn::OpDetWaveformMeta> PMTinfo;
    art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> WavePMTinfoAssns;
    art::PtrMaker<raw::OpDetWaveform> const makeWavePtr
      { event, waveformHandle.id() };
    
    for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
      
      PMTinfo.push_back(makePMTinfo(waveform));
      
      WavePMTinfoAssns.addSingle
        (makeWavePtr(iWaveform), makePMTinfoPtr(iWaveform));
      
    } // for
    
    // replica of discriminated gate-waveform association replacing the latter
    // with the PMT coverage with the same index as the waveform
    art::Assns<sbn::OpDetWaveformMeta, TriggerGateData_t> PMTinfoGateAssns;
    for (auto [ gatePtr, wavePtr ]: assns)
      PMTinfoGateAssns.addSingle(makePMTinfoPtr(wavePtr.key()), gatePtr);
    
    event.put(
      std::make_unique<std::vector<sbn::OpDetWaveformMeta>>(std::move(PMTinfo)),
      fOutputInstanceName
      );
    event.put(
      std::make_unique<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
        (std::move(WavePMTinfoAssns)),
      fOutputInstanceName
      );
    event.put(
      std::make_unique<art::Assns<sbn::OpDetWaveformMeta, TriggerGateData_t>>
        (std::move(PMTinfoGateAssns)),
      fOutputInstanceName
      );
    
  } // if save PMT coverage
  
  //
  // move the reformatted data into the event
  //
  event.put(
    std::make_unique<std::vector<TriggerGateData_t>>(std::move(data)),
    fOutputInstanceName
    );
  event.put(
    std::make_unique<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
      (std::move(assns)),
    fOutputInstanceName
    );
  
} // icarus::trigger::DiscriminatePMTwaveformsByChannel::produce()


//------------------------------------------------------------------------------
void
icarus::trigger::DiscriminatePMTwaveformsByChannel::refreshCurrentThresholds()
{
  //
  // start from scratch
  //
  fCurrentThresholds.clear();
  
  //
  // first set the default values for the first time
  //
  if (fNOpDetChannels > 0U) {
    fCurrentThresholds.resize(
      fNOpDetChannels,
      fDefaultThreshold.has_value()
        ? ADCvalueSetting_t
         { Source_t::Default, ADCCounts_t::castFrom(fDefaultThreshold.value()) }
        : ADCvalueSetting_t{ Source_t::Unset }
      );
  } // if
  
  //
  // override from the list
  //
  if (fCurrentThresholds.size() < fThresholdList.size())
    fCurrentThresholds.resize(fThresholdList.size());
  for (auto [ iThr, thr ]: util::enumerate(fThresholdList))
    fCurrentThresholds[iThr] = { Source_t::List, ADCCounts_t::castFrom(thr) };
  
  //
  // override from PMT configuration
  //
  if (fPMTconfig) {
    for (sbn::V1730Configuration const& boardConfig: fPMTconfig->boards) {
      for (sbn::V1730channelConfiguration const& channelConfig
        : boardConfig.channels
      ) {
        if (!channelConfig.hasChannelID()) continue;
        
        auto const channelSlot
          = static_cast<std::size_t>(channelConfig.channelID);
        if (fCurrentThresholds.size() <= channelSlot)
          fCurrentThresholds.resize(channelSlot + 1U);
        
        fCurrentThresholds[channelSlot] = {
          Source_t::PMTconfig,
          ADCCounts_t::castFrom(channelConfig.relativeThreshold())
          };
      
      } // for channel
    } // for board
  } // if we have PMT configuration
  
  //
  // override from channel list
  //
  for (auto const& chInfo: fChannelInfos) {
    if (!chInfo.threshold.has_value()) continue;
    auto const channelSlot = static_cast<std::size_t>(chInfo.channel);
    if (fCurrentThresholds.size() <= channelSlot) // make room
      fCurrentThresholds.resize(channelSlot + 1U);
    fCurrentThresholds[channelSlot] = {
      Source_t::ChannelSpec, ADCCounts_t::castFrom(chInfo.threshold.value())
      };
    
  } // for
  
} // icarus::trigger::DiscriminatePMTwaveformsByChannel::refreshCurrentThresholds()

//------------------------------------------------------------------------------
void
icarus::trigger::DiscriminatePMTwaveformsByChannel::refreshCurrentBaselines()
{
  //
  // start from scratch
  //
  fCurrentBaselines.clear();
  
  //
  // first set the default values for the first time
  //
  if ((fNOpDetChannels > 0U) && fBaseline.has_value()) {
    fCurrentBaselines.resize(
      fNOpDetChannels,
      ADCvalueSetting_t
        { Source_t::Default, ADCCounts_t::castFrom(fBaseline.value()) }
      );
  } // if
  
  //
  // override from channel list
  //
  for (auto const& chInfo: fChannelInfos) {
    if (!chInfo.baseline.has_value()) continue;
    auto const channelSlot = static_cast<std::size_t>(chInfo.channel);
    if (fCurrentBaselines.size() <= channelSlot) // make room
      fCurrentBaselines.resize(channelSlot + 1U);
    fCurrentBaselines[channelSlot] = {
      Source_t::ChannelSpec, ADCCounts_t::castFrom(chInfo.baseline.value())
      };
    
  } // for
  
} // icarus::trigger::DiscriminatePMTwaveformsByChannel::refreshCurrentBaselines()


//------------------------------------------------------------------------------
template <typename Logger /* mf::LogInfo */>
void
icarus::trigger::DiscriminatePMTwaveformsByChannel::printCurrentThresholdsAndBaselines
  (Logger& log) const
{
  
  auto const printSourceTag = [&log](Source_t source)
    {
      switch (source) {
        case Source_t::Default:     log << "(default value)";
          break;
        case Source_t::List:        log << "(from threshold list)";
          break;
        case Source_t::Unset:       log << "not set";
          break;
        case Source_t::PMTconfig:   log << "(from PMT readout configuration)";
          break;
        case Source_t::Run:         log << "(from run information)";
          break;
        case Source_t::Event:       log << "(from event information)";
          break;
        case Source_t::ChannelSpec: log << "(specified for this channel)";
          break;
      } // switch
    };
  
  std::size_t const nChannels
    = std::max(fCurrentThresholds.size(), fCurrentBaselines.size());
  log << "discrimination settings for " << nChannels << " channels";
  if ((fNOpDetChannels > 0U) && (fNOpDetChannels < nChannels))
    log << " (but only the first " << fNOpDetChannels << " will be processed)";
  log << ':';
  for (std::size_t const channel: util::counter(nChannels)) {
    
    log << "\n  channel " << channel
      << ":  threshold ";
    if (channel < fCurrentThresholds.size()) {
      auto const& thrInfo = fCurrentThresholds[channel];
      if (thrInfo.source != Source_t::Unset) log << thrInfo.value << ' ';
      printSourceTag(thrInfo.source);
    }
    else log << "n/a";
    
    log << ";  baseline ";
    if (channel < fCurrentBaselines.size()) {
      auto const& blineInfo = fCurrentBaselines[channel];
      if (blineInfo.source != Source_t::Unset) log << blineInfo.value << ' ';
      printSourceTag(blineInfo.source);
    }
    else log << "n/a";
    
  } // for channel
  
} // icarus::trigger::DiscriminatePMTwaveformsByChannel::printCurrentThresholdsAndBaselines()


//------------------------------------------------------------------------------
unsigned int icarus::trigger::DiscriminatePMTwaveformsByChannel::getNOpDetChannels
  (fhicl::OptionalAtom<unsigned int> const& param)
{
  unsigned int nChannels;
  return
    param(nChannels)? nChannels: lar::providerFrom<geo::Geometry>()->NOpDets();
} // icarus::trigger::DiscriminatePMTwaveformsByChannel::getNOpDetChannels()


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatePMTwaveformsByChannel::readChannelInfoSpecs
  (std::vector<ChannelInfo_t> const& config) -> std::vector<ChannelInfo_t>
{
  std::vector<ChannelInfo_t> chMap = config;
  
  auto const compareByChannel
    = [](ChannelInfo_t const& a, ChannelInfo_t const& b)
      { return a.channel == b.channel; }
    ;
  
  // sort and check for duplicates
  std::sort(chMap.begin(), chMap.end());
  auto iDup = std::unique(chMap.begin(), chMap.end(), compareByChannel);
  if (iDup != chMap.end()) {
    art::Exception e { art::errors::Configuration };
    e << "Found " << std::distance(iDup, chMap.end())
      << " duplicate channels in the configuration:";
    while (iDup != chMap.end()) e << " " << (iDup++)->channel;
    throw e;
  } // if duplicates
  
  return chMap;
} // icarus::trigger::DiscriminatePMTwaveformsByChannel::readChannelInfoSpecs()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::DiscriminatePMTwaveformsByChannel)


//------------------------------------------------------------------------------
