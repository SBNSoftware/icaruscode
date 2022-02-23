/**
 * @file   DiscriminatePMTwaveforms_module.cc
 * @brief  Module producing discriminated waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 4, 2019
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h"
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"
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
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
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
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <map>
#include <vector>
#include <string>
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::trigger { class DiscriminatePMTwaveforms; }

/**
 * @brief Produces discriminated optical waveforms.
 * 
 * This module produces a "discriminated waveform" for each optical detector
 * channel and for each configured threshold, via an algorithm of the class of
 * `icarus::trigger::TriggerGateBuilder` (`TriggerGateBuilder` configuration
 * parameter).
 * 
 * Thresholds are ultimately chosen by the tool in charge of actually run the
 * discrimination algorithm. Out of the thresholds that this algorithm produces,
 * it is possible to choose only a subset of them (`SelectThresholds`).
 * 
 * 
 * Output data products
 * =====================
 * 
 * * `std::vector<sbn::OpDetWaveformMeta>` parallel to the input optical
 *   detector waveforms, if `SavePMTcoverage` parameter is set; each entry in
 *   the collection is matched with the corresponding one in the input waveform
 *   collection; an association
 *   `art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>` is also produced.
 * * for each threshold _thr_ (in ADC counts), data products are created with
 *   an instance name equivalent to the threshold (e.g. for a threshold of 70
 *   ADC counts, the data product instance name will be `"70"`):
 *     * a data product of type
 *       `std::vector<icarus::trigger::OpticalTriggerGate::GateData_t>`,
 *       with one entry per optical channel, including _all_ optical channels
 *       (even when no signal above threshold is ever present); each gate object
 *       has as index in the vector the number of the optical channel;
 *     * an _art_ association, one trigger gate data to many optical waveforms,
 *       of each trigger gate data (as above) with all the optical waveforms
 *       which contributed to it; the order of the association pairs is the same
 *       as the one of the trigger gate data, i.e. by channel, and some trigger
 *       gate data entries may be not associated to any waveform and therefore
 *       they can not be listed in the association; association pairs within
 *       the same optical channel are sorted by optical waveform timestamp;
 *       the type of the association is
 *       `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, raw::OpDetWaveform>`;
 *       similarly, if `SavePMTcoverage` parameter was set, also an association
 *       `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, sbn::OpDetWaveformMeta>`
 *       is produced.
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
 * A terse description of the parameters is printed by running
 * `lar --print-description DiscriminatePMTwaveforms`.
 * 
 * * `OpticalWaveforms` (input tag): the data product containing all optical
 *   detector waveforms
 * * `Baselines` (input tag): the data product containing one baseline per
 *   waveform in data product (from `OpticalWaveforms`)
 * * `TriggerGateBuilder` (tool configuration): configuration of the _art_ tool
 *   used to discriminate the optional waveforms; the tool interface is
 *   `icarus::trigger::TriggerGateBuilder`.
 * * `SavePMTcoverage` (flag, default: `true`): also produces a collection of
 *   `sbn::OpDetWaveformMeta` representing each of the input waveforms; trigger
 *   tools can use this information in place for the more space-hungry waveforms
 *   for further processing.
 * * `OutputCategory` (string, default: `"DiscriminatePMTwaveforms"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * * `SelectThresholds` (sequence of ADC thresholds): if specified, only the
 *   waveforms discriminated with the thresholds in this list will be saved
 *   by the module; if not specified, all the thresholds produced by the
 *   algorithm will be saved. @note If a requested threshold is not eventually
 *   produced by the algorithm, an exception will be thrown by the module.
 * 
 * 
 */
class icarus::trigger::DiscriminatePMTwaveforms: public art::EDProducer {
  
    public:
  
  /// The type of data produced for a single channel.
  using TriggerGateData_t = icarus::trigger::OpticalTriggerGate::GateData_t;
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> OpticalWaveforms{
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
    
    fhicl::OptionalAtom<unsigned int> NChannels {
      Name("NChannels"),
      Comment("minimum number of channels to provide (default: all)")
      };
    
    fhicl::Atom<bool> SavePMTcoverage {
      Name("SavePMTcoverage"),
      Comment("write also a sbn::OpDetWaveformMeta collection"),
      true
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("tag of the module output to console via message facility"),
      "DiscriminatePMTwaveforms"
      };
    
    fhicl::OptionalSequence<raw::ADC_Count_t> SelectThresholds {
      Name("SelectThresholds"),
      Comment("thresholds to save (default: all those produced by algorithm)")
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit DiscriminatePMTwaveforms(Parameters const& config);
  
  // Plugins should not be copied or assigned.
  DiscriminatePMTwaveforms(DiscriminatePMTwaveforms const&) = delete;
  DiscriminatePMTwaveforms(DiscriminatePMTwaveforms&&) = delete;
  DiscriminatePMTwaveforms& operator=(DiscriminatePMTwaveforms const&) = delete;
  DiscriminatePMTwaveforms& operator=(DiscriminatePMTwaveforms&&) = delete;
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fOpDetWaveformTag; ///< Input optical waveform tag.
  
  ///< Input waveform baseline tag.
  std::optional<art::InputTag> const fBaselineTag;
  
  std::optional<float> const fBaseline; ///< A constant baseline level.
  
  unsigned int const fNOpDetChannels; ///< Number of optical detector channels.
  
  /// Thresholds selected for saving, and their instance name.
  std::map<icarus::trigger::ADCCounts_t, std::string> fSelectedThresholds;
  
  bool const fSavePMTcoverage; ///< Whether to save also `sbn::OpDetWaveformMeta`.
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  /// Algorithms to simulate trigger gates out of optical channel output.
  std::unique_ptr<icarus::trigger::TriggerGateBuilder> fTriggerGateBuilder;
  
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
  icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t fillChannelGaps
    (icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t gates) const;
  
  /// Fetches the optical detector count from geometry unless in `param`.
  static unsigned int getNOpDetChannels
    (fhicl::OptionalAtom<unsigned int> const& param);
  
}; // icarus::trigger::DiscriminatePMTwaveforms



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::DiscriminatePMTwaveforms
//------------------------------------------------------------------------------
icarus::trigger::DiscriminatePMTwaveforms::DiscriminatePMTwaveforms
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fOpDetWaveformTag(config().OpticalWaveforms())
  , fBaselineTag(util::fhicl::getOptionalValue(config().Baselines))
  , fBaseline(util::fhicl::getOptionalValue(config().Baseline))
  , fNOpDetChannels(getNOpDetChannels(config().NChannels))
  , fSavePMTcoverage(config().SavePMTcoverage())
  , fLogCategory(config().OutputCategory())
  , fTriggerGateBuilder
    (
      art::make_tool<icarus::trigger::TriggerGateBuilder>
        (config().TriggerGateBuilder_.get<fhicl::ParameterSet>())
    )
{
  //
  // optional configuration parameters
  //
  if (fBaseline && fBaselineTag) {
    throw art::Exception(art::errors::Configuration)
      << "Both `Baselines` ('" << fBaselineTag->encode()
      << "') and `Baseline` (" << *fBaseline
      << ") parameters specified, but they are exclusive!\n";
  }
  if (!fBaseline && !fBaselineTag) {
    throw art::Exception(art::errors::Configuration)
      << "Either `Baselines` or `Baseline` parameters is required.\n";
  }
  
  std::vector<raw::ADC_Count_t> selectedThresholds;
  if (!config().SelectThresholds(selectedThresholds)) {
    std::vector<ADCCounts_t> const& allThresholds
      = fTriggerGateBuilder->channelThresholds();
    
    if (allThresholds.empty()) {
      throw art::Exception(art::errors::Configuration)
        << "Trigger building algorithm reports no threshold!\n"
        << "Check the configuration of `TriggerGateBuilder` tool.\n";
    }
    std::transform(
      allThresholds.begin(), allThresholds.end(),
      std::back_inserter(selectedThresholds),
      std::mem_fn(&icarus::trigger::ADCCounts_t::value)
      );
    
  }
  else if (selectedThresholds.empty()) {
    throw art::Exception(art::errors::Configuration)
      << "Configured to save no thresholds!\n"
      << "Add values to the `SelectThresholds` configuration parameter.\n";
  }
  
  for (auto const& threshold: selectedThresholds) {
    fSelectedThresholds[icarus::trigger::ADCCounts_t{threshold}]
      = util::to_string(threshold);
  }
  
  {
    mf::LogInfo log(fLogCategory);
    log << "Selected thresholds:";
    for (auto const& [ threshold, instance ]: fSelectedThresholds)
      log << " " << threshold << " (\"" << instance << "\")";
  } // nameless block
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  
  //
  // declaration of output
  //
  
  if (fSavePMTcoverage) {
    produces<std::vector<sbn::OpDetWaveformMeta>>();
    produces<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>();
  }
  for (std::string const& instanceName: util::const_values(fSelectedThresholds))
  {
    produces<std::vector<TriggerGateData_t>>(instanceName);
    produces<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>(instanceName);
    if (fSavePMTcoverage) {
      produces<art::Assns<sbn::OpDetWaveformMeta, TriggerGateData_t>>
        (instanceName);
    }
  } // for
  
} // icarus::trigger::DiscriminatePMTwaveforms::DiscriminatePMTwaveforms()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatePMTwaveforms::produce(art::Event& event) {
  
  //
  // set up the algorithm to create the trigger gates
  //
  detinfo::DetectorTimings const detTimings {
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  fTriggerGateBuilder->resetup(detTimings);
  
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
  // retrieve the baseline information
  //
  std::vector<icarus::WaveformBaseline> fixedBaselines;
  std::vector<icarus::WaveformBaseline> const* baselines = nullptr;
  if (fBaselineTag) {
    baselines =
      &(event.getProduct<std::vector<icarus::WaveformBaseline>>(*fBaselineTag));
  }
  else {
    fixedBaselines.resize
      (waveforms.size(), icarus::WaveformBaseline{ *fBaseline });
    baselines = &fixedBaselines;
  }
  
  //
  // provide each waveform with additional information: baseline
  //
  if (baselines->size() != waveforms.size()) {
    assert(fBaselineTag);
    throw cet::exception("DiscriminatePMTwaveforms")
      << "Incompatible baseline information for the waveforms: "
      << waveforms.size() << " waveforms (" << fOpDetWaveformTag.encode()
      << ") for " << baselines->size() << " baselines ("
      << fBaselineTag->encode() << ")!\n";
  }
  std::vector<icarus::trigger::WaveformWithBaseline> waveformInfo;
  waveformInfo.reserve(waveforms.size());
  for (auto const& [ waveform, baseline ]: util::zip(waveforms, *baselines))
    waveformInfo.emplace_back(&waveform, &baseline);
  
  {
    mf::LogDebug log { fLogCategory };
    log << "Discrimination of " << waveforms.size() << " PMT waveforms from '"
      << fOpDetWaveformTag.encode() << "'";
    if (fBaselineTag)
      log << " with baselines from '" << fBaselineTag->encode() << "'";
  }
  
  
  //
  // define channel-level trigger gate openings as function on threshold
  //
  
  // this is a collection where each entry (of type `TriggerGates`) contains
  // the complete set of trigger gates for an event.
  std::vector<icarus::trigger::TriggerGateBuilder::TriggerGates> const&
    triggerGatesByThreshold = fTriggerGateBuilder->build(waveformInfo);
  
  { // nameless block
    mf::LogTrace log(fLogCategory);
    log << "Trigger gates from " << triggerGatesByThreshold.size()
      << " thresholds:\n";
    for (auto const& triggerGates: triggerGatesByThreshold)
      triggerGates.dump(log);
    
  } // nameless block
  
  
  //
  // prepare output
  //
  std::optional<art::PtrMaker<sbn::OpDetWaveformMeta>> const makePMTinfoPtr
    = fSavePMTcoverage
    ? std::optional<art::PtrMaker<sbn::OpDetWaveformMeta>>(event)
    : std::nullopt
    ;
  
  for (icarus::trigger::TriggerGateBuilder::TriggerGates const& gates
    : triggerGatesByThreshold
  ) {
    
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
    auto [ discrGates, waveGateAssns ]
      = icarus::trigger::transformIntoOpticalTriggerGate
        (fillChannelGaps(std::move(gates).gates()), makeGatePtr, opDetWavePtrs)
      ;
    
    if (fSavePMTcoverage) {
      // replica of discriminated gate-waveform association replacing the latter
      // with the PMT coverage with the same index as the waveform
      assert(makePMTinfoPtr);
      art::Assns<sbn::OpDetWaveformMeta, TriggerGateData_t> assns;
      for (auto [ gatePtr, wavePtr ]: waveGateAssns)
        assns.addSingle((*makePMTinfoPtr)(wavePtr.key()), gatePtr);
      event.put(
        std::make_unique<art::Assns<sbn::OpDetWaveformMeta, TriggerGateData_t>>
          (std::move(assns)),
        instanceName
        );
    } // if save PMT coverage

    //
    // move the reformatted data into the event
    //
    event.put(
      std::make_unique<std::vector<TriggerGateData_t>>(std::move(discrGates)),
      instanceName
      );
    event.put(
      std::make_unique<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
        (std::move(waveGateAssns)),
      instanceName
      );
    
  } // for all extracted thresholds
  
  
  // add a simple one-to-one PMT coverage - waveform association just in case
  if (fSavePMTcoverage) {
    assert(makePMTinfoPtr);
    
    sbn::OpDetWaveformMetaMaker const makePMTinfo { detTimings };
    
    std::vector<sbn::OpDetWaveformMeta> PMTinfo;
    art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> assns;
    art::PtrMaker<raw::OpDetWaveform> const makeWavePtr
      { event, waveformHandle.id() };
    
    for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
      
      PMTinfo.push_back(makePMTinfo(waveform));
      
      assns.addSingle(makeWavePtr(iWaveform), (*makePMTinfoPtr)(iWaveform));
      
    } // for
    
    event.put
      (std::make_unique<std::vector<sbn::OpDetWaveformMeta>>(std::move(PMTinfo)));
    event.put(
      std::make_unique<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
        (std::move(assns))
      );
    
  } // if save PMT coverage
  
  
} // icarus::trigger::DiscriminatePMTwaveforms::produce()


//------------------------------------------------------------------------------
icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t
icarus::trigger::DiscriminatePMTwaveforms::fillChannelGaps
  (icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t gates) const
{
  
  using GateDataColl_t
    = icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t;
  using Gate_t = GateDataColl_t::value_type; // TrackedOpticalTriggerGate
  
  //
  // fill a map channel -> gate (missing channels have a nullptr gate)
  //
  std::vector<Gate_t const*> gateMap(fNOpDetChannels, nullptr);
  for (Gate_t const& gate: gates) {
    assert(!gate.channels().empty());
    
    auto const channel = gate.channels().front();
    if (static_cast<std::size_t>(channel) >= gateMap.size())
      gateMap.resize(channel + 1U, nullptr);
    assert(gateMap[channel] == nullptr);
    gateMap[channel] = &gate;
    
  } // for all gates
  
  //
  // fill
  //
  GateDataColl_t allGates;
  allGates.reserve(gateMap.size());
  for (auto const& [ channelNo, gate ]: util::enumerate(gateMap)) {
    
    if (gate) {
      assert(gate->channels().front() == channelNo);
      allGates.push_back(std::move(*gate));
    }
    else {
      allGates.emplace_back(
        icarus::trigger::OpticalTriggerGateData_t{ raw::Channel_t(channelNo) }
        );
    }
    
  } // for
  
  return allGates;
} // icarus::trigger::DiscriminatePMTwaveforms::fillChannelGaps()


//------------------------------------------------------------------------------
unsigned int icarus::trigger::DiscriminatePMTwaveforms::getNOpDetChannels
  (fhicl::OptionalAtom<unsigned int> const& param)
{
  unsigned int nChannels;
  return
    param(nChannels)? nChannels: lar::providerFrom<geo::Geometry>()->NOpDets();
} // icarus::trigger::DiscriminatePMTwaveforms::getNOpDetChannels()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::DiscriminatePMTwaveforms)


//------------------------------------------------------------------------------
