/**
 * @file   DiscriminatedAdderSignal_module.cc
 * @brief  Module to discriminate adder signals.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 6, 2023; rewritten on December 5, 2025
 */

#undef NDEBUG // FIXME

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Algorithms/AdderSignalSimulationTypes.h" // ADCsettings_t
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Algorithms/WaveformWithBaseline.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // transformIntoOpticalTriggerGate()
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icaruscode/IcarusObj/ChannelToChannelMap.h"
#include "icaruscode/Utilities/DataProductPointerMap.h" // mapDataProductPointers()
#include "icarusalg/Utilities/TimeIntervalConfig.h" // TimeInterval in FHiCL
#include "icarusalg/Utilities/TimeInterval.h"
#include "icarusalg/Utilities/mfLoggingClass.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib/maybe_ref.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <algorithm> // std::find(), std::remove_if()
#include <cassert>
#include <memory> // std::unique_ptr
#include <map>
#include <optional>
#include <stdexcept> // std::logic_error
#include <string>
#include <utility> // std::move()
#include <vector>


//------------------------------------------------------------------------------
namespace icarus::trigger { class DiscriminatedAdderSignal; }
/**
 * @brief Produces discriminated optical waveforms for PMT adders.
 * 
 * This module produces a "discriminated waveform" for each adder board and
 * for each configured threshold, via an algorithm of the class of
 * `icarus::trigger::TriggerGateBuilder` (`TriggerGateBuilder` configuration
 * parameter).
 * 
 * The input is adder waveforms (for example from simulation via
 * `SimulateAdderSignal` module / `AdderSignalSimulation` algorithm, or directly
 * from data with special detector setup).
 * The signals are discriminated against thresholds and saved.
 * 
 * Thresholds are ultimately chosen by the tool in charge of actually running
 * the discrimination algorithm. Out of the thresholds that this algorithm
 * produces, it is possible to choose only a subset of them
 * (`SelectThresholds`).
 * 
 * 
 * Output data products
 * =====================
 * 
 * * for each threshold _thr_ (in ADC counts), data products are created with
 *   an instance name equivalent to the threshold (e.g. for a threshold of 60 mV
 *   the data product instance name will be `"491"`, that is the 60 millivolt
 *   translated in ADC counts). The data products are:
 *     * a data product of type
 *       `std::vector<icarus::trigger::OpticalTriggerGate::GateData_t>`,
 *       with one entry per adder module; each gate object is associated to
 *       the adder channel ID and to the IDs of all contributing PMT channels;
 *     * _art_ associations: one trigger gate data to the adder waveforms
 *       it originates from. The order of the association pairs is the same
 *       as the one of the trigger gate data, i.e. by channel, and some trigger
 *       gate data entries may be not associated to any waveform and therefore
 *       they may not appear in the association; association pairs within
 *       the same optical channel are sorted by optical waveform timestamp;
 *       the type of the association is
 *       `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, raw::OpDetWaveform>`.
 *       If `SaveMetaAssns` is set, also an association
 *       art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, sbn::OpDetWaveformMeta>`
 *       is produced, just like the previous one but associating to the adder
 *       waveform metadata instead of the waveform itself.
 * 
 * 
 * Input data products
 * ====================
 * 
 * From run data:
 * * `icarus::ChannelToChannelMap<raw::Channel_t>` (`AdderChannelMapTag`):
 *   mapping of the adder channels to the PMT ones. In this module this mapping
 *   is not fundamental, but it allows to let the module know which adder
 *   channels are expected to be present and to associate them to the original
 *   PMT channels. This data product should be reflecting the input adder
 *   waveform sources, and ideally is produced by the same module that also
 *   produced those (e.g. `SimulateAdderSignal`).
 * 
 * From event data:
 * * `std::vector<raw::OpDetWaveform>` (`AdderWaveformTag`): a single waveform
 *   for each recorded adder activity; the activity belongs to a single
 *   channel, but there may be multiple waveforms on the same channel. The time
 *   stamp is expected to be from the
 *   @ref DetectorClocksElectronicsTime "electronics time scale" and therefore
 *   expressed in microseconds. The waveforms are expected to be in order of
 *   increasing adder channel ID, and then of time stamp.
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>` (`BaselineTag`,
 *   optional): the baseline associated to each of the input adder waveforms.
 *   If missing, baselines are assumed to be 0.
 * * `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>` (`WaveformTag`;
 *   only if `SaveMetaAssns` is set): associations between the input waveforms
 *   and their metadata (will be used to create associations to the
 *   discriminated waveforms).
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
 * A terse description of the parameters is printed by running
 * `lar --print-description DiscriminatedAdderSignal`.
 * 
 * * `AdderWaveformTag` (input tag): the data product containing all adder
 *   waveforms. The time windows not covered by any waveforms on one channel are
 *   considered to have no signal.
 * * `BaselineTag` (input tag, optional): the data product containing the
 *   association between each input waveform (`AdderWaveformTag`) and its
 *   estimated baseline. If omitted, input waveforms are considered already
 *   baseline-subtracted.
 * * `AdderWaveformMetaTag` (input tag, optional): the data product containing
 *   the association between each input waveform (`AdderWaveformTag`) and its
 *   metadata; if omitted, the tag from `AdderWaveformTag` is used.
 * * `AdderChannelMapTag` (input tag, default: same as `AdderWaveformTag`): the
 *   mapping of adder channels to PMT channels used for the adder waveforms;
 *   possible sources include `SimulateAdderSignal` and `WriteAdderChannelMap`
 *   modules.
 * * `TimeInterval` (table, optional): limits the time range of the adders to
 *   the specified time interval; if omitted, no limit is applied. The time is
 *   specified in @ref DetectorClocksElectronicsTime "electronics time scale".
 *   No more than two of the following must be specified:
 *     * `Start` (time, default: ages ago): start from this time point.
 *     * `End` (time, default: ages in the future): end at this time point.
 *     * `Duration` (time, default: ages): duration of the time interval.
 * * `TriggerGateBuilder` (tool configuration): configuration of the _art_ tool
 *   used to discriminate the optional waveforms; the tool interface is
 *   `icarus::trigger::TriggerGateBuilder`.
 * * `SelectThresholds` (sequence of ADC thresholds): if specified, only the
 *   waveforms discriminated with the thresholds in this list will be saved
 *   by the module; if not specified, all the thresholds produced by the
 *   algorithm will be saved. @note If a requested threshold is not eventually
 *   produced by the algorithm, an exception will be thrown by the module.
 * * `SaveMetaAssns` (flag, default: `true`): if set, an association between the
 *   metadata of the input waveforms (from `AdderWaveformMetaTag`) and the
 *   discriminated waveforms is also saved.
 * * `LogCategory` (string, default: `"DiscriminatedAdderSignal"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * 
 * 
 */
class icarus::trigger::DiscriminatedAdderSignal
  : public art::SharedProducer, private icarus::ns::util::mfLoggingClass
{
  
    public:
  
  using electronics_time = detinfo::timescales::electronics_time;
  
  /// The type of data produced for a single channel.
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGate::GateData_t;
  
  /// Type of channel map read from the run data.
  using ChannelMap_t = icarus::ChannelToChannelMap<raw::Channel_t>;
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> AdderWaveformTag{
      Name{ "AdderWaveformTag" },
      Comment{ "input tag of adder board waveform data product" }
      };
    
    fhicl::OptionalAtom<art::InputTag> BaselineTag{
      Name{ "BaselineTag" },
      Comment{ "input tag for associations between waveforms and baselines" }
      };
    
    fhicl::OptionalAtom<art::InputTag> AdderWaveformMetaTag{
      Name{ "AdderWaveformMetaTag" },
      Comment{
        "input tag for associations between waveforms and metadata"
        " [default: as AdderWaveformTag]"
        }
      };
    
    fhicl::OptionalAtom<art::InputTag> AdderChannelMapTag{
      Name{ "AdderChannelMapTag" },
      Comment
        { "input tag for the adder channel map [default: as AdderWaveformTag]" }
      };
    
    fhicl::DelegatedParameter TriggerGateBuilder_{
      Name{ "TriggerGateBuilder" },
      Comment{
        "parameters for generating trigger gates from optical channel output"
        }
      };
    
    icarus::ns::fhicl::TimeIntervalOptionalTable<electronics_time> TimeInterval{
      Name{ "TimeInterval" },
      Comment{ "limit the adder time to this interval, relative to beam gate" }
      };
    
    fhicl::OptionalSequence<std::string> SelectThresholds{
      Name{ "SelectThresholds" },
      Comment{ "thresholds to save (default: all those produced by algorithm)" }
      };
    
    fhicl::Atom<bool> SaveMetaAssns{
      Name{ "SaveMetaAssns" },
      Comment{
        "write associations between discriminated waveforms"
        " and their input metadata"
      },
      true
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "tag of the module output to console via message facility" },
      "DiscriminatedAdderSignal"
      };
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  explicit DiscriminatedAdderSignal
    (Parameters const& config, art::ProcessingFrame const& frame);
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Loads run-specific adder channel mapping.
  virtual void beginRun(art::Run& run, art::ProcessingFrame const&) override;
  
  /// Creates the data products.
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  /// Time interval in electronics time.
  using TimeInterval_t = icarus::ns::util::TimeInterval<electronics_time>;
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fAdderWaveformTag; ///< Input adder waveform tag.
  
  art::InputTag const fBaselineTag; ///< Input adder waveform baseline tag.
  
  art::InputTag const fAdderWaveformMetaTag; ///< Input waveform metadata tag.
  
  art::InputTag const fAdderChannelMapTag; ///< Input adder channel mapping tag.
  
  TimeInterval_t const fTimeInterval; ///< Range of time to create adders on.
  
  /// Thresholds selected for saving, and their instance name.
  std::map<icarus::trigger::ADCCounts_t, std::string> fSelectedThresholds;
  
  bool const fSaveMetaAssns; ///< Whether to save input waveform associations.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  /// Algorithms to simulate trigger gates out of optical channel output.
  std::unique_ptr<icarus::trigger::TriggerGateBuilder> fTriggerGateBuilder;
  
  // --- END Algorithms --------------------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  /// Duration of an optical tick.
  util::quantities::intervals::microseconds const fOpticalTick;
  
  adder::types::ADCsettings_t const fADCsettings; ///< ADC conversion utility.
  
  /// Adder/PMT channel mapping (changes every run).
  ChannelMap_t fAdderChannelMap;
  
  // --- END Service variables -------------------------------------------------
  
  
  /**
   * @brief Fetches all the adder waveforms for input.
   * @param waveforms the adder waveforms to be matched
   * @param toBaselines association finder fro waveforms to baselines
   * @param timeInterval the time interval we are interested in
   * @return the adder waveforms and their baselines
   * 
   * The data products (adder waveforms and optionally associations to baseline)
   * are read from the configured input tags.
   * If no baseline tag is configured, the baseline pointers are left null and
   * it is up to the user to react appropriately.
   * 
   * The returned list includes all and only the waveforms which have at least
   * one sample in the requested `timeInterval`.
   */
  std::vector<WaveformWithBaseline> matchWaveformsAndBaselines(
    std::vector<raw::OpDetWaveform> const& waveforms,
    art::FindOne<icarus::WaveformBaseline> const* toBaselines,
    TimeInterval_t const& timeInterval
    ) const;
  
  /**
   * @brief Fills missing channel gates and adds PMT channel numbers.
   * @param gates the result of discrimination (for a single threshold)
   * @return the same `gates`, enriched
   * 
   * The output collection of tracking gates includes one gate for each channel
   * in the current adder map (`fAdderChannelMap`).
   * If a channel is featured in the input `gates`, the tracking trigger gate
   * for that channel is _moved_ to the returned collection. Otherwise, a new
   * empty gate is created for that channel (always closed and with no tracked
   * waveforms). In both cases, the PMT channels for the adder channel mapping
   * are also added to the channel list of the gate, which will have first the
   * adder channel number and then all the PMT channels contributing to it
   * (whether there is any actual contribution in this gate from them).
   */
  TriggerGateBuilder::TriggerGates::GateData_t fillChannelGaps
    (TriggerGateBuilder::TriggerGates::GateData_t gates) const;
  
  
  //@{
  /// Returns the time interval covered by the `waveform`.
  TimeInterval_t waveformInterval(raw::OpDetWaveform const& waveform) const;
  TimeInterval_t waveformInterval(WaveformWithBaseline const& waveform) const;
  //@}
  
  /// Returns an `ADCsettings_t` object with the appropriate sampling time.
  adder::types::ADCsettings_t makeADCsettings() const;
  
  
}; // icarus::trigger::DiscriminatedAdderSignal



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  template <typename Coll, typename T>
  bool contains(Coll const& coll, T&& value) {
    using std::begin, std::end;
    auto const b = begin(coll), e = end(coll);
    return std::find(b, e, value) != e;
  }
  
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& obj)
    { return std::make_unique<T>(std::move(obj)); }
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::trigger::DiscriminatedAdderSignal
//------------------------------------------------------------------------------
icarus::trigger::DiscriminatedAdderSignal::DiscriminatedAdderSignal
  (Parameters const& config, art::ProcessingFrame const& frame)
  : art::SharedProducer{ config }
  , icarus::ns::util::mfLoggingClass{ config().LogCategory() }
  // configuration
  , fAdderWaveformTag  { config().AdderWaveformTag() }
  , fBaselineTag       { config().BaselineTag().value_or(art::InputTag{}) }
  , fAdderWaveformMetaTag
    { config().AdderWaveformMetaTag().value_or(fAdderWaveformTag) }
  , fAdderChannelMapTag
    { config().AdderChannelMapTag().value_or(fAdderWaveformTag) }
  , fTimeInterval
    {
      icarus::ns::fhicl::makeTimeInterval(config().TimeInterval())
        .value_or(TimeInterval_t{})
    }
  , fSaveMetaAssns     { config().SaveMetaAssns()    }
  // algorithms
  , fTriggerGateBuilder{
      art::make_tool<icarus::trigger::TriggerGateBuilder>
        (config().TriggerGateBuilder_.get<fhicl::ParameterSet>())
    }
  // service cache
  , fOpticalTick  {
    detinfo::makeDetectorClocksWithUnits
      (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob())
      .OpticalClockPeriod()
    }
  , fADCsettings{ makeADCsettings() }
{
  async<art::InEvent>();
  
  //
  // configuration post-processing
  //
  
  // fill or check the selected thresholds
  std::vector<ADCCounts_t> const& allThresholds
    = fTriggerGateBuilder->channelThresholds(); // guaranteed sorted
  if (allThresholds.empty()) {
    throw art::Exception{ art::errors::Configuration }
      << "Trigger building algorithm reports no threshold!\n"
      << "Check the configuration of `TriggerGateBuilder` tool.\n";
  }
  
  std::vector<ADCCounts_t> selectedThresholds;
  if (!config().SelectThresholds()) selectedThresholds = allThresholds;
  else {
    selectedThresholds = TriggerGateBuilder::parseThresholds
      (config().SelectThresholds().value(), fADCsettings);
    std::sort(begin(selectedThresholds), end(selectedThresholds));
    std::vector<ADCCounts_t> spurious;
    for (ADCCounts_t threshold: selectedThresholds) {
      if (!contains(allThresholds, threshold)) spurious.push_back(threshold);
    }
    if (!spurious.empty()) {
      art::Exception e{ art::errors::Configuration };
      e << "Selected " << spurious.size() << "/" << selectedThresholds.size()
        << " thresholds that are not configured in `TriggerGateBuilder` tool:";
      for (ADCCounts_t threshold: spurious) e << " " << threshold;
      e << "\nThe configured tool thresholds are:";
      for (ADCCounts_t threshold: allThresholds) e << " " << threshold;
      throw e << ".";
    }
  }
  if (selectedThresholds.empty()) {
    throw art::Exception{ art::errors::Configuration }
      << "Configured to save no thresholds!\n"
      << "Add values to the `SelectThresholds` configuration parameter.\n";
  }
  
  for (auto const& threshold: selectedThresholds)
    fSelectedThresholds[threshold] = util::to_string(threshold.value());
  
  //
  // configuration dump
  //
  {
    auto log = mfLogInfo();
    log << "Configuration:";
    if (fTimeInterval.empty()) log << "\n * no time restriction";
    else {
      log << "\n * time interval: [ "
        << fTimeInterval.start << " - " << fTimeInterval.stop << " ]";
    }
    
    bool const allSelected = selectedThresholds.size() == allThresholds.size();
    log << "\n * selected";
    if (!allSelected)
      log << " " << selectedThresholds.size() << "/" << allThresholds.size();
    log << " thresholds";
    if (!allSelected) log << " (and [ignored]):";
    for (auto const& thr: allThresholds) {
      bool const selected = (fSelectedThresholds.count(thr) > 0);
      log << " " << (selected? "": "[")
        << thr << " ADC = " << fADCsettings.to_mV(thr);
      if (selected) log << " (\"" << fSelectedThresholds[thr] << "\")";
      else log << "]";
    }
  } // nameless block
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fAdderWaveformTag);
  if (!fBaselineTag.empty()) {
    consumes<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>
      (fBaselineTag);
  }
  if (fSaveMetaAssns && !fAdderWaveformMetaTag.empty()) {
    consumes<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
      (fAdderWaveformMetaTag);
  }
  
  //
  // declaration of output
  //
  for (std::string const& instanceName: util::const_values(fSelectedThresholds))
  {
    produces<std::vector<TriggerGateData_t>>(instanceName);
    produces<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>(instanceName);
    if (fSaveMetaAssns) {
      produces<art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta>>
        (instanceName);
    }
  } // for thresholds
  
} // icarus::trigger::DiscriminatedAdderSignal::DiscriminatedAdderSignal()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatedAdderSignal::beginRun
  (art::Run& run, art::ProcessingFrame const&)
{
  
  fAdderChannelMap = run.getProduct<ChannelMap_t>(fAdderChannelMapTag);
  
} // icarus::trigger::DiscriminatedAdderSignal::beginRun()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatedAdderSignal::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  //
  // set up the discrimination algorithm
  //
  detinfo::DetectorTimings const detTimings {
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  fTriggerGateBuilder->resetup(detTimings);
  mfLogDebug()
    << "Trigger at " << detTimings.TriggerTime() << ", "
    << (detTimings.TriggerTime() - detTimings.BeamGateTime())
    << " after beam gate opened at " << detTimings.BeamGateTime();
  
  //
  // prepare input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fAdderWaveformTag);
  auto const& waveforms = *waveformHandle;
  
  std::optional const toBaselines = fBaselineTag.empty()? std::nullopt
    : std::make_optional(art::FindOne<icarus::WaveformBaseline>
        (waveformHandle, event, fBaselineTag)
      );
  
  // attach to each waveform additional information: baseline
  std::vector<WaveformWithBaseline> const waveformInfo
    = matchWaveformsAndBaselines
      (waveforms, toBaselines? &*toBaselines: nullptr, fTimeInterval)
    ;

  mfLogTrace()
    << "Processing " << waveformInfo.size() << "/" << waveforms.size()
    << " adder waveforms from '" << fAdderWaveformTag.encode() << "'";
  
  //
  // discrimination (all thresholds at once)
  //
  
  // this is a collection where each entry (of type `TriggerGates`) contains
  // the complete set of trigger gates for an event for a given threshold
  // (the `TriggerGates` type knows which one).
  // Each trigger gate is tracking the waveform(s) it comes from, that here is
  // all the adder waveforms from a single channel in the requested interval
  using TriggerGates_t = icarus::trigger::TriggerGateBuilder::TriggerGates;
  std::vector<TriggerGates_t> triggerGatesByThreshold
    = fTriggerGateBuilder->build(waveformInfo);
  
  { // nameless block
    auto log = mfLogTrace();
    log << "Trigger gates from " << triggerGatesByThreshold.size()
      << " thresholds:\n";
    for (TriggerGates_t const& triggerGates: triggerGatesByThreshold)
      triggerGates.dump(log);
    
  } // nameless block
  
  //
  // prepare output
  //

  // map address of waveform to art pointer to that waveform
  auto const& opDetWavePtrs
    = util::mapDataProductPointers(event, waveformHandle);
  
  std::optional<art::FindOneP<sbn::OpDetWaveformMeta>> waveformToMeta; // cache
  
  for (TriggerGates_t const& gates: triggerGatesByThreshold) {
    
    ADCCounts_t const thr = gates.threshold();
    
    //
    // find the threshold and its label, if we care of it
    //
    auto iThr = fSelectedThresholds.find(thr);
    if (iThr == fSelectedThresholds.end()) continue; // not interested
    
    std::string const& instanceName = iThr->second;
    art::PtrMaker<TriggerGateData_t> const makeGatePtr(event, instanceName);
    
    
    //
    // reformat the results for the threshold
    //
    auto [ discrGates, discrToAdder ]
      = icarus::trigger::transformIntoOpticalTriggerGate
        (fillChannelGaps(std::move(gates).gates()), makeGatePtr, opDetWavePtrs)
      ;
    
    if (fSaveMetaAssns) {
      if (!waveformToMeta) { // fill cache
        // more input: needed to retrieve input waveform metadata pointers
        waveformToMeta = std::make_optional<art::FindOneP<sbn::OpDetWaveformMeta>>
          (waveformHandle, event, fAdderWaveformMetaTag);
      }
      
      // replica of discriminated gate-waveform association replacing the latter
      // with the metadata with the same index as the waveform
      art::Assns<TriggerGateData_t, sbn::OpDetWaveformMeta> discrToAdderMeta;
      for (auto [ gatePtr, wavePtr ]: discrToAdder)
        discrToAdderMeta.addSingle(gatePtr, waveformToMeta->at(wavePtr.key()));
      event.put(moveToUniquePtr(discrToAdderMeta), instanceName);
    } // if save metadata associations

    //
    // move the reformatted data into the event
    //
    event.put(moveToUniquePtr(discrGates), instanceName);
    event.put(moveToUniquePtr(discrToAdder), instanceName);
    
  } // for all extracted thresholds
  
} // icarus::trigger::DiscriminatedAdderSignal::produce()


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatedAdderSignal::matchWaveformsAndBaselines(
  std::vector<raw::OpDetWaveform> const& waveforms,
  art::FindOne<icarus::WaveformBaseline> const* toBaselines,
  TimeInterval_t const& timeInterval
) const -> std::vector<WaveformWithBaseline> {

  std::vector<WaveformWithBaseline> waveformInfo;
  
  for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
    
    icarus::WaveformBaseline const* baseline = nullptr;
    try {
      baseline = toBaselines? &(toBaselines->at(iWaveform).ref()): nullptr;
    }
    catch (std::logic_error const&) { // thrown by `maybe_ref::ref()`
      throw art::Exception{ art::errors::NotFound }
        << "No baseline ('" << fBaselineTag.encode()
        << "') found for waveform #" << iWaveform << " (channel="
        << waveform.ChannelNumber() << " time=" << waveform.TimeStamp()
        << ").\n";
    }
    
    waveformInfo.emplace_back(&waveform, baseline);
  } // for waveforms
  
  // remove the waveforms that do not overlap at all with the requested interval
  if (!timeInterval.empty()) {
    
    auto notInTime = [this, &timeInterval](WaveformWithBaseline const& wi)
      { return !waveformInterval(wi).overlaps(timeInterval); };
    
    waveformInfo.erase(
      std::remove_if(begin(waveformInfo), end(waveformInfo), notInTime),
      end(waveformInfo)
      );
    
    mfLogTrace()
      << "Preserved " << waveformInfo.size() << "/" << waveforms.size()
      << " adder waveforms overlapping " << timeInterval;
  }
  
  return waveformInfo;
  
} // icarus::trigger::DiscriminatedAdderSignal::matchWaveformsAndBaselines()


//------------------------------------------------------------------------------
icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t
icarus::trigger::DiscriminatedAdderSignal::fillChannelGaps
  (TriggerGateBuilder::TriggerGates::GateData_t gates) const
{
  
  using GateDataColl_t = TriggerGateBuilder::TriggerGates::GateData_t;
  using Gate_t = GateDataColl_t::value_type; // TrackedOpticalTriggerGate
  
  //
  // fill a map channel -> gate
  //
  std::map<raw::Channel_t, Gate_t*> gateMap;
  for (Gate_t& gate: gates) gateMap[gate.gate().channel()] = &gate;
  
  //
  // fill
  //
  GateDataColl_t allGates;
  allGates.reserve(fAdderChannelMap.size());
  for (auto const& [ adderChannel, PMTchannels ]: fAdderChannelMap) {
    
    std::optional<Gate_t> gate;
    if (auto it = gateMap.find(adderChannel); it != end(gateMap)) {
      // we have the gate
      gate = std::move(*(it->second));
      assert(gate->gate().channel() == adderChannel);
    }
    else {
      // we have to make up a gate: empty, with no tracking and this channel:
      gate = Gate_t{ icarus::trigger::OpticalTriggerGateData_t{ adderChannel } };
    }
    assert(gate);
    for (raw::Channel_t const PMTchannel: PMTchannels)
      gate->gate().addChannel(PMTchannel);
    
    allGates.push_back(std::move(*gate));
    
  } // for all channels
  
  return allGates;
} // icarus::trigger::DiscriminatePMTwaveforms::fillChannelGaps()


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatedAdderSignal::waveformInterval
  (raw::OpDetWaveform const& waveform) const -> TimeInterval_t
{
  electronics_time const start{ waveform.TimeStamp() };
  return { start, start + fOpticalTick * waveform.size() };
} // icarus::trigger::DiscriminatedAdderSignal::waveformInterval(WaveformWithBaseline)


auto icarus::trigger::DiscriminatedAdderSignal::waveformInterval
  (WaveformWithBaseline const& waveform) const -> TimeInterval_t
  { return waveformInterval(waveform.waveform()); }


//------------------------------------------------------------------------------
auto icarus::trigger::DiscriminatedAdderSignal::makeADCsettings() const
  -> adder::types::ADCsettings_t
{
  adder::types::ADCsettings_t ADCsettings;
  ADCsettings.samplingTime = fOpticalTick.quantity();
  return ADCsettings;
}


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::DiscriminatedAdderSignal)


//------------------------------------------------------------------------------
