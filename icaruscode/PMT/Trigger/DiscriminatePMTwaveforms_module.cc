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
#include "icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Data/TriggerGateData.h"
#include "icaruscode/Utilities/DataProductPointerMap.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/Utilities/quantities_fhicl.h" // for ADCCounts_t parameters
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/enumerate.h"
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
 *       `art::Assns<icarus::trigger::OpticalTriggerGate::GateData_t, raw::OpDetWaveform>`.
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<raw::OpDetWaveform>`: a single waveform for each recorded
 *      optical detector activity; the activity belongs to a single channel, but
 *      there may be multiple waveforms on the same channel. The time stamp is
 *      expected to be from the
 *      @anchor DetectorClocksElectronicsTime "electronics time scale"
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
 * * `TriggerGateBuilder` (tool configuration): configuration of the _art_ tool
 *   used to discriminate the optional waveforms; the tool interface is
 *   `icarus::trigger::TriggerGateBuilder`.
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
    
    fhicl::DelegatedParameter TriggerGateBuilder_ {
      Name("TriggerGateBuilder"),
      Comment
        ("parameters for generating trigger gates from optical channel output")
      };
    
    fhicl::OptionalAtom<unsigned int> NChannels {
      Name("NChannels"),
      Comment("minimum number of channels to provide (default: all)")
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
  /// Prepares the plots to be filled.
  virtual void beginJob() override;
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fOpDetWaveformTag; ///< Input optical waveform tag.
  
  unsigned int const fNOpDetChannels; ///< Number of optical detector channels.
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  /// Thresholds selected for saving, and their instance name.
  std::map<icarus::trigger::ADCCounts_t, std::string> fSelectedThresholds;
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  detinfo::DetectorClocksData fClockData;   // FIXME: this assumes that the cached clock data is valid for the whole job.
  detinfo::DetectorTimings fDetTimings;
  
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
   * channel `0`, the next the one for channel `1` and so on).
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
  , fNOpDetChannels(getNOpDetChannels(config().NChannels))
  , fLogCategory(config().OutputCategory())
  , fClockData{art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob()}
  , fDetTimings{fClockData}
  , fTriggerGateBuilder
    (
      art::make_tool<icarus::trigger::TriggerGateBuilder>
        (config().TriggerGateBuilder_.get<fhicl::ParameterSet>())
    )
{
  //
  // optional configuration parameters
  //
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
  
  for (std::string const& instanceName: util::const_values(fSelectedThresholds))
  {
    produces<std::vector<TriggerGateData_t>>(instanceName);
    produces<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>(instanceName);
  } // for
  
} // icarus::trigger::DiscriminatePMTwaveforms::DiscriminatePMTwaveforms()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatePMTwaveforms::beginJob() {
  
  //
  // set up the algorithm to create the trigger gates
  //
  fTriggerGateBuilder->setup(fDetTimings);
  
} // icarus::trigger::DiscriminatePMTwaveforms::beginJob()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatePMTwaveforms::produce(art::Event& event) {
  
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
  // define channel-level trigger gate openings as function on threshold
  //
  
  // this is a collection where each entry (of type `TriggerGates`) contains
  // the complete set of trigger gates for an event.
  std::vector<icarus::trigger::TriggerGateBuilder::TriggerGates> const&
    triggerGatesByThreshold = fTriggerGateBuilder->build(waveforms);
  
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
    auto thresholdData = icarus::trigger::transformIntoOpticalTriggerGate
      (fillChannelGaps(std::move(gates).gates()), makeGatePtr, opDetWavePtrs);

    //
    // move the reformatted data into the event
    //
    event.put(
      std::make_unique<std::vector<TriggerGateData_t>>
        (std::move(std::get<0U>(thresholdData))),
      instanceName
      );
    event.put(
      std::make_unique<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>
        (std::move(std::get<1U>(thresholdData))),
      instanceName
      );
    
  } // for all extracted thresholds
  
} // icarus::trigger::DiscriminatePMTwaveforms::produce()


//------------------------------------------------------------------------------
icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t
icarus::trigger::DiscriminatePMTwaveforms::fillChannelGaps
  (icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t gates) const
{
  
  using GateDataColl_t
    = icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t;
  using Gate_t = GateDataColl_t::value_type; // SingleChannelOpticalTriggerGate
  
  //
  // fill a map channel -> gate (missing channels have a nullptr gate)
  //
  std::vector<Gate_t const*> gateMap(fNOpDetChannels, nullptr);
  for (Gate_t const& gate: gates) {
    assert(gate.hasChannels());
    
    auto const channel = gate.channel();
    if (static_cast<std::size_t>(channel) >= gateMap.size())
      gateMap.resize(channel + 1U, nullptr);
    gateMap[channel] = &gate;
    
  } // for all gates
  
  //
  // fill
  //
  GateDataColl_t allGates;
  allGates.reserve(gateMap.size());
  for (auto const& [ channelNo, gate ]: util::enumerate(gateMap)) {
    
    if (gate) {
      assert(gate->channel() == channelNo);
      allGates.push_back(std::move(*gate));
    }
    else allGates.emplace_back(Gate_t::ChannelID_t(channelNo));
    
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
