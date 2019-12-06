/**
 * @file   DiscriminatePMTwaveforms_module.cc
 * @brief  Module producing discriminated waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 4, 2019
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Data/TriggerGateData.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
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


namespace icarus::trigger {
  
  class DiscriminatePMTwaveforms;
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
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
 * * `TriggerGateBuilder` (tool configuration): configuration of the _art_ tool
 *   used to discriminate the optional waveforms; the tool interface is
 *   `icarus::trigger::TriggerGateBuilder`.
 * * `OutputCategory` (string, default: `"DiscriminatePMTwaveforms"`): label
 *     for the category of messages in the console output; this is the label
 *     that can be used for filtering messages via MessageFacility service
 *     configuration.
 * * `SelectThresholds` (sequence of ADC thresholds): if specified, only the
 *     waveforms discriminated with the thresholds in this list will be saved
 *     by the module; if not specified, all the thresholds produced by the
 *     algorithm will be saved. @note If a specified threshold is not eventually
 *     produced by the algorithm, an exception will be thrown by the module.
 *     @todo To be implemented
 * 
 * @todo Complete the documentation
 * 
 * * `PlotSets` (list of selection criteria, see below):
 *     for each element in the list, a set of plots is created. The plots are
 *     contributed only by the events fulfilling the specified category.
 *     If no plot set is specified at all, a single plot set is created in the
 *     main ROOT directory, contributed by all the events in the input.
 * 
 * 
 * Selection criteria for plot sets
 * ---------------------------------
 * 
 * * `CategoryKey` (string; default: from the activated criteria):
 *     name used for the ROOT subdirectory where the plots are stored.
 *     It may be used in ROOT object names as needed.
 * * `CategoryName` (string; default: from the activated criteria):
 *     string used in plot titles and other ROOT objects.
 * * `CurrentType` (string; default: `"any"`): plot only the events from a weak
 *     current interaction of this type; valid types are `"charged"` and
 *     `"neutral"`. All events are plotted when `"any"` is instead specified.
 * 
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
    
#if 0
    fhicl::Atom<art::InputTag> PropagatedParticles{
      Name("PropagatedParticles"),
      Comment("label of input simulated particles in the detector"),
      "largeant" // tradition demands
      };
#endif // 0
    
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
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("tag of the module output to console via message facility"),
      "DiscriminatePMTwaveforms"
      };
    
    fhicl::OptionalSequence<raw::ADC_Count_t> SelectThresholds {
      Name("SelectThresholds"),
      Comment("thresholds to save (default: all those produced by algorithm)")
      };
    
#if 0
    fhicl::OptionalTable<icarus::trigger::RegionFinder::Config> RegionFinder {
      Name("RegionFinder"),
      Comment
        ("parameters for determining regions of interest on optical waveforms")
      };
    
    fhicl::Atom<microsecond::value_t> TriggerGateDuration {
      Name("TriggerGateDuration"),
      Comment("length of time interval when optical triggers are accepted [us]")
      };
    
    fhicl::Atom<unsigned int> ADCbits {
      Name("ADCbits"),
      Comment("number of bits in the ADC; count is assumed signed"),
      14U
      };
    
    fhicl::Atom<bool> ApplyStyle {
      Name("ApplyStyle"),
      Comment("apply a ROOT style to all plots (this style still needs tuning!)"), // TODO
      true
      };
    
    fhicl::OptionalSequence
      <icarus::trigger::PlotDispatcherByGenerator::Parameters>
    PlotSets {
      Name("PlotSets"),
      Comment("list of plot categories (default: one with all events)")
      };
    
    fhicl::DelegatedParameter Patterns {
      Name("Patterns"),
      Comment("configuration of the trigger patterns to be simulated (only one pattern at a time supported so far)")
      };
    
#endif // 0
    
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
  
#if 0
  /// Sets the run-dependent parameters up.
  virtual void beginRun(art::Run const&) override;
#endif // 0
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
#if 0
  /// Type of data histogram (non-ROOT, just the binned data).
  using BinnedData_t = util::DataHistogram<>;
  
  /// Type of list of output ROOT directories for plots.
  using PlotDestinations_t
    = icarus::trigger::TriggerPattern::PlotDestinations_t;
  using ConstPlotDestinations_t
    = icarus::trigger::TriggerPattern::ConstPlotDestinations_t;
  
  
  /// Name of the plot box (and ROOT directory) for region of interest plots.
  static std::string const RegionsPlotboxName;
  
#endif // 0
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag fOpDetWaveformTag; ///< Input optical waveform tag.
  std::string fLogCategory; ///< Category name for the console output stream.
  
  /// Thresholds selected for saving, and their instance name.
  std::map<icarus::trigger::ADCCounts_t, std::string> fSelectedThresholds;
  
#if 0
  art::InputTag fPropagatedParticlesTag; ///< Input simulated particles.
  
  /// Duration of the gate during with global optical triggers are accepted.
  microsecond fTriggerGateDuration;
  
  // more obscure entries here
  unsigned int fMaxADC; ///< Maximum ADC value.
  
  bool fApplyStyle; ///< Whether to apply ROOT style to plots.
  
#endif // 0
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  detinfo::DetectorTimings fDetTimings;
  
#if 0
  /// ROOT directory where all the plots are written.
  art::TFileDirectory fOutputDir;
  
  geo::GeometryCore const& fGeom;
  detinfo::DetectorClocks const& fDetClocks;
  
  /// Total number of optical channels (PMTs).
  unsigned int fNOpDetChannels = 0;
  microsecond const fOpDetTickDuration; ///< Optical detector sampling period.
  
#endif // 0
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
#if 0
  /// Algorithm to find the interesting regions in a optical waveform.
  std::optional<RegionFinder> fRegionFinder;
#endif // 0
  
  /// Algorithms to simulate trigger gates out of optical channel output.
  std::unique_ptr<icarus::trigger::TriggerGateBuilder> fTriggerGateBuilder;
  
  
#if 0
  // --- BEGIN -- Trigger patterns ---------------------------------------------
  std::vector<std::unique_ptr<icarus::trigger::TriggerPattern>>
    fTriggerPatterns;
  // --- END -- Trigger patterns -----------------------------------------------
#endif // 0
  
  // --- END Algorithms --------------------------------------------------------
  
  
  // --- BEGIN Internal variables ----------------------------------------------
#if 0
  std::vector<std::unique_ptr<PlotDispatcherBase>> fPlotSets; ///< Set of plots.
  
  TStyle* fStyle = nullptr; ///< ROOT style used to create the plots
  
  // --- END Internal variables ------------------------------------------------
  
  /// Returns whether a region finder algorithm is configured.
  bool hasRegionFinder() const { return fRegionFinder.has_value(); }
  
  //@{
  /// Returns the current region finding algorithm.
  RegionFinder& regionFinder() { return fRegionFinder.value(); }
  RegionFinder const& regionFinder() const { return fRegionFinder.value(); }
  //@}
  
  /**
   * @brief Runs the algorithm to find the regions of interest of the waveforms.
   * @param waveforms all PMT waveforms
   * @return number of interestig regions found in each channel
   * 
   * The result is expressed as function of the optical detector channel number,
   * containing the number of detected interesting regions in that channel.
   */
  std::vector<unsigned int> findRegions(
    std::vector<raw::OpDetWaveform> const& waveforms,
    PlotDestinations_t& plotSets
    );
  
  /// Dispatches the beam gate setup to all configured patterns.
  void updateBeamGates();
  
  /// Constructs all the plot sets.
  std::vector<std::unique_ptr<PlotDispatcherBase>> initPlotSets(
    art::TFileDirectory baseOutputDir,
    std::vector<icarus::trigger::PlotDispatcherByGenerator::Config> const&
      setConfig
    );
  
  /// Initializes all the plots.
  void preparePlots();
  
  /// Initializes all the region of interest plots in the specified plot box.
  void initRegionPlots(PlotSandbox& plotSet);
  
  /// Creates a ROOT style for our plots.
  TStyle* preparePlotStyle(art::TFileDirectory outputDir);
  
  
  /**
   * @brief Fills an histogram with the energy loss evolution in time.
   * @tparam BIter type of iterator to the first of the particles
   * @tparam EIter type of iterator after the last of the particles
   * @param begin iterator to the first of the particles
   * @param end iterator after the last of the particles
   * @param maxTime duration of the interval where energy loss is profiled [us]
   * @param nTimeBins number of bins in the profile
   * @return a `BinnedData_t` object containing the energy loss profile
   * 
   * A list of particles forming a cascade is provided between `begin` and
   * `end`.
   * The returned object contains the amount of energy [GeV] lost during each
   * of the time bins by all the particles in the list.
   * The time is plot relative to the start of the interaction
   * (first particle of the list).
   * If the particles survive with energy, that energy is not deposited but it
   * is considered to have escaped.
   * 
   * No constraint is imposed on the volume where energy deposition happens.
   */
  template <typename BIter, typename EIter>
  BinnedData_t plotEnergyLossOverTime
    (BIter begin, EIter end, double maxTime, unsigned int nTimeBins) const;
  
  /**
   * @brief Fills an histogram with the energy loss evolution in time.
   * @tparam Coll type of collection of particles
   * @param coll collection of particles
   * @param maxTime duration of the interval where energy loss is profiled [us]
   * @param nTimeBins number of bins in the profile
   * @return a `BinnedData_t` object containing the energy loss profile
   * 
   * This is the equivalent of the version with two iterators, but accepting
   * a whole collection.
   */
  template <typename Coll>
  BinnedData_t plotEnergyLossOverTime
    (Coll const& coll, double maxTime, unsigned int nTimeBins) const;
  
  
  // --- BEGIN -- Event filtering ----------------------------------------------
  /// @name Event classification for category plots
  /// @{
  
  /// Common implementation for `collectPlotSets()`.
  template <typename PlotSetType>
  static auto collectPlotSetsImpl
    (PlotSetType& plotSets, art::Event const* event = nullptr);
  
  //@{
  /// Returns a list of all plot sets.
  ConstPlotDestinations_t collectPlotSets() const;
  PlotDestinations_t collectPlotSets();
  //@}
  
  //@{
  /// Returns a list of all plot sets the `event` is selected for.
  ConstPlotDestinations_t collectPlotSets
    (art::Event const& event) const;
  PlotDestinations_t collectPlotSets(art::Event const& event);
  //@}
  
  auto plotDestinations(art::Event const& event) const;
  
  /// @}
  // --- END -- Event filtering ------------------------------------------------
  
#endif // 0
  
  /// Creates a map from address of data product element to _art_ pointer to it.
  template <typename Handle>
  static auto mapDataProductPointers
    (art::Event const& event, Handle const& handle);
  
  
}; // icarus::trigger::DiscriminatePMTwaveforms


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::DiscriminatePMTwaveforms
//------------------------------------------------------------------------------
#if 0
std::string const icarus::trigger::TriggerDesignPlots::RegionsPlotboxName
  = "regions";

#endif // 0
//------------------------------------------------------------------------------
icarus::trigger::DiscriminatePMTwaveforms::DiscriminatePMTwaveforms
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fOpDetWaveformTag(config().OpticalWaveforms())
  , fLogCategory(config().OutputCategory())
  , fDetTimings(*lar::providerFrom<detinfo::DetectorClocksService>())
  , fTriggerGateBuilder
    (
      art::make_tool<icarus::trigger::TriggerGateBuilder>
        (config().TriggerGateBuilder_.get<fhicl::ParameterSet>())
    )
#if 0
  , fPropagatedParticlesTag(config().PropagatedParticles())
  , fTriggerGateDuration   (config().TriggerGateDuration())
  , fMaxADC                ((1UL << (config().ADCbits() - 1)) - 1)
  , fApplyStyle            (config().ApplyStyle())
  // services
  , fOutputDir             (*art::ServiceHandle<art::TFileService>())
  , fGeom                  (*lar::providerFrom<geo::Geometry>())
  , fDetClocks             (*lar::providerFrom<detinfo::DetectorClocksService>())
  , fNOpDetChannels        (fGeom.NOpDets())
  , fOpDetTickDuration     (fDetClocks.OpticalClock().TickPeriod())
  // others
  , fPlotSets
    (
      initPlotSets(
        fOutputDir,
        util::fhicl::getOptionalValue(
          config().PlotSets,
          std::vector<icarus::trigger::PlotDispatcherByGenerator::Config>{}
        )
      )
    )
#endif // 0
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
  
  // TODO
  for (std::string const& instanceName: util::values(fSelectedThresholds))
  {
    produces<std::vector<TriggerGateData_t>>(instanceName);
    produces<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>(instanceName);
  } // for
  
#if 0
  consumes<std::vector<simb::MCParticle>>(fPropagatedParticlesTag);
  
  for (auto const& plotSet: fPlotSets)
    plotSet->declareConsume(consumesCollector());
  
  //
  // set up the algorithm for regions of interest, if requested
  //
  icarus::trigger::RegionFinder::Config regionFinderConfig;
  if (config().RegionFinder(regionFinderConfig))
    fRegionFinder.emplace(regionFinderConfig);
  
  //
  // trigger pattern algorithm creation
  //
  
  /*
   * waiting for art experts to explain me how to have a
   * `Sequence<DelegatedParameter>` or similar to manage a list of tools...
   */
  for (auto const& patternCfg: { config().Patterns }) {
    
    fTriggerPatterns.push_back(
      art::make_tool<icarus::trigger::TriggerPattern>
        (patternCfg.get<fhicl::ParameterSet>())
      );
    
  } // for
  
  mf::LogInfo log(LogCategory);
  log
    <<   "The detector has " << fNOpDetChannels << " optical channels."
    << "\nOptical detector digitization clock period is: " << fOpDetTickDuration
    ;
  
#endif // 0
  
} // icarus::trigger::DiscriminatePMTwaveforms::DiscriminatePMTwaveforms()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatePMTwaveforms::beginJob() {
  
  //
  // set up the algorithm to create the trigger gates
  //
  fTriggerGateBuilder->setup(fDetTimings);
  
#if 0
  //
  // set up of trigger pattern algorithms
  // (this includes also their plots)
  //
  
  // create a list of all possible plot set directory destinations
  // and inform the trigger pattern algorithms about it
  PlotDestinations_t allDestinations = collectPlotSets();
  
  std::optional<util::ROOT::TStyleChanger> styleChanger;
  if (fApplyStyle) styleChanger.emplace(preparePlotStyle(fOutputDir));
  
  icarus::trigger::BeamGateMaker const makeBeamGate { fDetClocks };
  auto beamGate = makeBeamGate(fTriggerGateDuration);
  
  for (auto& triggerPattern:
    util::make_transformed_span(fTriggerPatterns, util::dereference())
    )
  {
    triggerPattern.setup(allDestinations);
  } // for
  
  preparePlots();
  
#endif // 0
} // icarus::trigger::DiscriminatePMTwaveforms::beginJob()


//------------------------------------------------------------------------------
#if 0
void icarus::trigger::TriggerDesignPlots::beginRun(art::Run const&) {
  
  // clocks settings might change on every run
  updateBeamGates();
  
} // icarus::trigger::TriggerDesignPlots::beginRun()


#endif // 0

//------------------------------------------------------------------------------
template <typename Handle>
auto icarus::trigger::DiscriminatePMTwaveforms::mapDataProductPointers
  (art::Event const& event, Handle const& handle)
{
  using Data_t = typename Handle::element_type::value_type;
  using Map_t = std::map<Data_t const*, art::Ptr<Data_t>>;
  
  static_assert(
    std::is_same_v<std::vector<Data_t>, typename Handle::element_type>,
    "mapDataProductPointers() requires handles to STL vectors of data"
    );
  
  Map_t map;
  art::PtrMaker<Data_t> makePtr { event, handle.id() };
  for (auto const& [ i, item ]: util::enumerate(*handle))
    map[&item] = makePtr(i);
  return map;
} // icarus::trigger::DiscriminatePMTwaveforms::mapDataProductPointers()


//------------------------------------------------------------------------------
void icarus::trigger::DiscriminatePMTwaveforms::produce(art::Event& event) {
  
  //
  // fetch input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  auto const& waveforms = *waveformHandle;
  
  // map address of waveform to art pointer to that waveform
  auto const& opDetWavePtrs = mapDataProductPointers(event, waveformHandle);
  
  
#if 0
  auto const& particles
    = *(event.getValidHandle<std::vector<simb::MCParticle>>(fPropagatedParticlesTag));
  
  // find out which plot sets are interested to this event
  PlotDestinations_t destinations = collectPlotSets(event);
  
  //
  // 1. Truth information analysis
  // 1.1 group truth into interactions
  // 1.2 extract relevant information for each interaction (time, size; generated energy? lost energy?)
  // 1.2.1 tag each particle as entering, exiting or contained
  // 1.3 compute the deposited energy for each interaction
  // 
  // 2. Optical waveform analysis
  // 2.1 find regions of interest with minimal threshold
  // 2.2 determine the photoelectron threshold equivalent
  //     (this is likely non-linear, more refined parametrisation may be needed)
  // 2.3 define channel-level trigger gate openings as function on threshold
  // 2.4 define combination of trigger gates
  //
  // 3. Plots
  // 3.1 sanity check plots
  // 3.1.1 on simulated particles
  // 3.1.2 on waveforms
  // 3.1.2.1 on peak amplitude distribution
  
  //
  // 1. Truth information analysis
  //
  
  // 1.1 group truth into interactions
  
  sim::ParticlePartitioner particleSorter;
  sim::ParticlePartitioner::Interactions_t interactions
    = particleSorter.makeSortedPartition(particles);
  
  MF_LOG_TRACE("TriggerDesignPlots")
    << "Event contains " << particles.size() << " particles grouped in "
    << interactions.interactions.size() << " interactions.";
  if (!interactions.others.empty()) {
    MF_LOG_TRACE("TriggerDesignPlots")
      << "Also, event contains " << interactions.others.size()
      << " particles not grouped in any interaction.";
  } // if unsorted
  
  for (auto const& interaction: interactions.interactions) {
    assert(!interaction.empty());
    if (!sim::ParticlePartitioner::isPrimary(*(interaction.front()))) {
      MF_LOG_TRACE("TriggerDesignPlots") << "Particle tree starting with "
        << *(interaction.front()) << " (not a primary particle)";
    }
    
    BinnedData_t const energyLossOverTime = plotEnergyLossOverTime(
      util::make_transformed_span(interaction, util::dereference()),
      100.0, 100000
      );
    
    // TODO need to extract statistics
    
  } // for
  
  
  // 1.2 extract relevant information for each interaction (time, size; generated energy? lost energy?)
  
  // TODO
  
  // 1.2.1 tag each particle as entering, exiting or contained
  
  // TODO
  
  // 1.3 compute the deposited energy for each interaction
  
  // TODO
  
  // 
  
  
  //
  // 2.1 find regions of interest with minimal threshold
  //
  // Thresholds come from region finder algorithm configuration.
  //
  
  if (hasRegionFinder()) {
    
    findRegions(waveforms, destinations); // ignoring the result so far
    
  } // if finding regions of interest
  
  //
  // 2.2 determine the photoelectron threshold equivalent
  //
  
  // TODO
  
#endif // 0
  //
  // 2.3 define channel-level trigger gate openings as function on threshold
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
    
    //
    // collect the input
    //
    icarus::trigger::TriggerGateBuilder::TriggerGates::GateData_t channelGates
      = std::move(gates).gates(); // stolen
    
    //
    // create objects for the data products
    //
    auto gateData
      = std::make_unique<std::vector<TriggerGateData_t>>();
    auto gateToWaveforms
      = std::make_unique<art::Assns<TriggerGateData_t, raw::OpDetWaveform>>();
    
    art::PtrMaker<TriggerGateData_t> const makeGatePtr { event, instanceName };
    
    for (icarus::trigger::SingleChannelOpticalTriggerGate& channelGate
      : channelGates
    ) {
      
      if (channelGate.waveforms().empty()) {
        // special case: no waveforms, no gate... we actually don't even know
        // which channel, but we assume it's the next one
        
        gateData->emplace_back(); // empty vector, no associations
        
        continue;
      }
      
      // check sorting; if this consistently fails,
      // we might be forced to change this algorithm or the gate builder
      raw::Channel_t const channel = channelGate.channel();
      if (channel != static_cast<raw::Channel_t>(gateData->size())) {
        throw art::Exception(art::errors::LogicError)
          << "Discriminated channels contain gaps: got " << channel
          << " instead of " << gateData->size() << "!\n";
      }
      
      // we steal the data from the gate
      gateData->push_back(std::move(channelGate.gateLevels()));
      
      // pointer to the gate data we have just added:
      art::Ptr<TriggerGateData_t> gatePtr = makeGatePtr(gateData->size() - 1U);
      
      //
      // produce the associations
      //
      for (raw::OpDetWaveform const* waveform: channelGate.waveforms()) {
        
        gateToWaveforms->addSingle(gatePtr, opDetWavePtrs.at(waveform));
        
      } // for waveforms
      
      
    } // for all channels
    
    event.put(std::move(gateData), instanceName);
    event.put(std::move(gateToWaveforms), instanceName);
    
  } // for all extracted thresholds
  
  
  
  

#if 0
  //
  // 3.1 "Channels with signal over XXX p.e. threshold"
  //     (thresholds picked from configuration file)
  //
  
  // 3.1.1 group trigger gates and sum them according to the pattern
  
  using detinfo::timescales::optical_tick;
  
  icarus::trigger::BeamGateMaker const makeBeamGate { fDetClocks };
  icarus::trigger::OpticalTriggerGate const openTriggerGate
    = makeBeamGate(fTriggerGateDuration);
  MF_LOG_TRACE(LogCategory)
    << "Open gate for optical triggers: " << openTriggerGate;
  
  // simulate all patterns and tell them which plot sets to fill
  for (auto& pattern: fTriggerPatterns)
    pattern->simulateResponseByThreshold(triggerGatesByThreshold, destinations);
  
#endif // 0
  
} // icarus::trigger::DiscriminatePMTwaveforms::produce()


//------------------------------------------------------------------------------
#if 0
std::vector<unsigned int> icarus::trigger::TriggerDesignPlots::findRegions(
  std::vector<raw::OpDetWaveform> const& waveforms, PlotDestinations_t& plotSets
) {
  if (!hasRegionFinder()) return {};
  
  auto regionPlotbox
    = [](auto* plotSet){ return plotSet->findSandbox(RegionsPlotboxName); };
  
  std::vector<unsigned int> NRegionsPerChannel(fNOpDetChannels, 0U);
  
  for (raw::OpDetWaveform const& waveform: waveforms) {
    
    // find the regions of interest
    // find the peak of each
    icarus::trigger::RegionFinder::Regions_t const& regions
      = regionFinder().extractRegions(waveform);
    
    NRegionsPerChannel.at(waveform.ChannelNumber()) += regions.size();
    
    for (icarus::trigger::PlotSandbox* plotSet:
      util::make_transformed_span(plotSets, regionPlotbox)
      )
    {
      
      auto* pHRegionDuration = plotSet->use<TH1F>("RegionDuration");
      auto* pHWaveformPeaks = plotSet->use<TH1F>("PeakSignal");
      auto* pHPeakDelay = plotSet->use<TH1F>("PeakTimeDelay");
      auto* pHPeakDelayVsAmplitude
        = plotSet->use<TH2F>("PeakTimeDelayVsAmplitude");
      
      for (auto const& region: regions) {
        // add it to the histogram
        pHWaveformPeaks->Fill(region.peakAmplitude().value());
        pHRegionDuration->Fill(region.length());
        pHPeakDelay->Fill(region.peakDelay());
        pHPeakDelayVsAmplitude->Fill
          (region.peakAmplitude().value(), region.peakDelay());
      } // for regions
    } // for plot sets
    
  } // for waveforms
  
  unsigned int const nRegions = 
    std::accumulate(NRegionsPerChannel.begin(), NRegionsPerChannel.end(), 0U);
  for (icarus::trigger::PlotSandbox* plotSet:
    util::make_transformed_span(plotSets, regionPlotbox)
    )
  {
    
    auto* pHNRegionsPerChannel = plotSet->use<TH1F>("NRegionsPerChannel"); 
    for (unsigned int count: NRegionsPerChannel)
      pHNRegionsPerChannel->Fill(count);
    
    plotSet->use<TH1F>("NRegions")->Fill(nRegions);
    
  } // for plot sets
  
  return NRegionsPerChannel;
  
} // icarus::trigger::TriggerDesignPlots::findRegions()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerDesignPlots::updateBeamGates() {
  
  icarus::trigger::BeamGateMaker const makeBeamGate { fDetClocks };
  auto const beamGate = makeBeamGate(fTriggerGateDuration);
  
  for (auto& triggerPattern:
    util::make_transformed_span(fTriggerPatterns, util::dereference())
    )
  {
    triggerPattern.setBeamGate(beamGate);
  } // for
  
} // icarus::trigger::TriggerDesignPlots::updateBeamGates()


//------------------------------------------------------------------------------
std::vector<std::unique_ptr<icarus::trigger::PlotDispatcherBase>>
icarus::trigger::TriggerDesignPlots::initPlotSets(
  art::TFileDirectory baseOutputDir,
  std::vector<icarus::trigger::PlotDispatcherByGenerator::Config> const&
    setConfig
  )
{
  std::vector<std::unique_ptr<PlotDispatcherBase>> plotSets;
  std::set<std::string> keys;
  plotSets.reserve(setConfig.size());
  for (auto const& config: setConfig) {
    plotSets.push_back(
      std::make_unique<icarus::trigger::PlotDispatcherByGenerator>
        (baseOutputDir, config)
      );
    std::string const key = plotSets.back()->categoryKey("");
    if (keys.count(key) > 1) {
      throw art::Exception(art::errors::Configuration)
        << "Multiple plot sets are configured with key '" << key << "'!\n";
    }
    keys.insert(key);
  } // for
  
  if (plotSets.empty()) {
    plotSets.push_back(
      std::make_unique<icarus::trigger::ConsentingPlotDispatcher>(baseOutputDir)
      );
  }
  
  return plotSets;
} // icarus::trigger::TriggerDesignPlots::PlotSets::initPlotSets()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerDesignPlots::preparePlots() {
  
  // set our style as current if available, for all and only this method
  std::optional<util::ROOT::TStyleChanger> styleChanger;
  if (fApplyStyle && fStyle) styleChanger.emplace(fStyle);
  
  //
  // create empty histograms
  //
  PlotDestinations_t destinations = collectPlotSets();
  for (PlotSandbox* plotSet: destinations) {
    initRegionPlots
      (plotSet->addSubSandbox(RegionsPlotboxName, "Regions of interest"));
  } // for
  
} // icarus::trigger::TriggerDesignPlots::preparePlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerDesignPlots::initRegionPlots
  (PlotSandbox& plotSet)
{
  unsigned int const maxADC = fMaxADC;
  
  // plots are left in the box and retrieved via ROOT when needed.
  
  // * regions of interest per channel
  plotSet.make<TH1F>(
    "NRegionsPerChannel",
    "Regions of interest per channel;regions of interest;channels",
    200, 0.0, 200.0
    );
  
  // * regions of interest per event
  plotSet.make<TH1F>(
    "NRegions",
    "Regions of interest per event;regions of interest;event",
    1'000 - 1, 1.0, 1'000.0
    );
  
  //   * region of interest length
  plotSet.make<TH1F>(
    "RegionDuration",
    "Duration of each region of interest;duration  [ ticks ];regions of interest",
    1'000, 0.0, 10'000.0
    );
  
  
  // * peak amplitude
  plotSet.make<TH1F>(
    "PeakSignal",
    "Signal peak;peak level  [ ADC ];waveform regions",
    maxADC, 0.0, double(maxADC)
    );
  
  
  // * peak delay from above threshold
  plotSet.make<TH1F>(
    "PeakTimeDelay",
    "Signal peak time from raise above threshold;#Delta t  [ ticks ];waveform regions",
    100, 0.0, 100.0
    );
  plotSet.make<TH2F>(
    "PeakTimeDelayVsAmplitude",
    "Signal peak time from raise above threshold vs. peak amplitude;peak amplitude  [ / 5 ADC counts ];#Delta t  [ ticks ]",
    maxADC / 5, 0.0, double(maxADC),
    100, 0.0, 100.0
    );
  
  //   * peak time (or above threshold time) TODO
  //   * peak amplitude vs. peak time TODO
  
  
} // icarus::trigger::TriggerDesignPlots::initRegionPlots()


//------------------------------------------------------------------------------
TStyle* icarus::trigger::TriggerDesignPlots::preparePlotStyle
  (art::TFileDirectory outputDir)
{
  if (fStyle) return fStyle; // already there!
  
  auto const& modDesc = moduleDescription(); // from `art::ModuleBase`
  std::string const title
    = "Style for " + modDesc.moduleName() + "[" + modDesc.moduleLabel() + "]";
  std::string const name
    = "Style_" + modDesc.moduleName() + "_" + modDesc.moduleLabel();
  
  fStyle = outputDir.makeAndRegister<TStyle>(name.c_str(), title.c_str(),
    name.c_str(), title.c_str()
    );
  
  fStyle->SetFillColor(10); // the white that does not make it B&W
  fStyle->SetFillStyle(kSolid);
  fStyle->SetLineWidth(2);
  fStyle->SetMarkerStyle(kFullDotLarge);
//  fStyle->SetMarkerSize();
//  fStyle->SetTextFont(42); // sans serif, medium thickness, normal
  fStyle->SetTextSize(0.07);
  
  // canvases and pads
  fStyle->SetCanvasBorderMode(0);
  fStyle->SetCanvasColor(kWhite);
  fStyle->SetPadBorderMode(0);
  fStyle->SetPadColor(kWhite);
  fStyle->SetFrameBorderMode(0);
  fStyle->SetFrameFillColor(kWhite);
  fStyle->SetFrameFillStyle(kSolid);
  fStyle->SetFrameLineWidth(1);
  fStyle->SetPadGridX(1);
  fStyle->SetPadGridY(1);
  fStyle->SetPadLeftMargin(0.12);
  fStyle->SetPadRightMargin(0.05);
  fStyle->SetPadBottomMargin(0.15);
//  fStyle->SetPadTopMargin(0.1);
  
  // grid and ticks
  fStyle->SetGridColor(kGray);
  fStyle->SetGridStyle(kDotted);
  fStyle->SetGridWidth(1);
  fStyle->SetPadTickX(1);
  fStyle->SetPadTickY(1);
  
  // plots
  fStyle->SetOptTitle(1);
  fStyle->SetTitleAlign(22); // centered ("2") both horizontally and vertically
  fStyle->SetTitleBorderSize(0);
  fStyle->SetTitleFont(42, "XYZ");
  fStyle->SetTitleOffset(0.85, "XYZ");
  fStyle->SetTitleSize(0.08, "XYZ");
  fStyle->SetHistLineWidth(3);
  fStyle->SetHistMinimumZero();
  fStyle->SetStripDecimals(kFALSE);
  fStyle->SetLabelFont(42, "XYZ");
//  fStyle->SetLabelOffset(0.005, "XYZ");
  fStyle->SetLabelSize(0.06, "XYZ"); // large
  
  // legend and statistics
  fStyle->SetLegendBorderSize(0); // no border at all
  fStyle->SetLegendFillColor(kWhite);
  fStyle->SetLegendFont(42);
  fStyle->SetLegendTextSize(0.8);
  fStyle->SetOptStat(1110);
  
  return fStyle;
} // icarus::trigger::TriggerDesignPlots::preparePlotStyle()


//------------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto icarus::trigger::TriggerDesignPlots::plotEnergyLossOverTime
  (BIter begin, EIter end, double maxTime, unsigned int nTimeBins) const
  -> BinnedData_t
{
  BinnedData_t hist { 0.0, maxTime, nTimeBins };
  
  auto pos
    = [](simb::MCTrajectory::const_iterator const& it){ return it->first; };
  auto mom
    = [](simb::MCTrajectory::const_iterator const& it){ return it->second; };
  
  double const startTime = begin->T(); // time of first point of first particle
  for (auto iPart = begin; iPart != end; ++iPart) {
    
    simb::MCTrajectory const& traj = iPart->Trajectory();
    if (traj.size() < 2) continue;
    
    auto beforeStep = traj.begin();
    auto afterStep = beforeStep;
    auto const endPoint = traj.end();
    while (++afterStep != endPoint) {
      double const lostEnergy = mom(afterStep).E() - mom(beforeStep).E();
      double const time = pos(afterStep).T() - startTime;
      hist.fill(time, lostEnergy);
      beforeStep = afterStep;
    } // while
  } // for iPart
  
  return hist;
} // icarus::trigger::TriggerDesignPlots::plotEnergyLossOverTime(BIter, EIter)


//------------------------------------------------------------------------------
template <typename Coll>
auto icarus::trigger::TriggerDesignPlots::plotEnergyLossOverTime
  (Coll const& coll, double maxTime, unsigned int nTimeBins) const
  -> BinnedData_t
{
  return plotEnergyLossOverTime
    (util::begin(coll), util::end(coll), maxTime, nTimeBins);
} // icarus::trigger::TriggerDesignPlots::plotEnergyLossOverTime(Coll)


//------------------------------------------------------------------------------
template <typename PlotSetType>
auto icarus::trigger::TriggerDesignPlots::collectPlotSetsImpl
  (PlotSetType& plotSets, art::Event const* event /* = nullptr */)
{
  // PlotSetType is some std::vector<std::unique_ptr<PlotDispatcherBase>>;
  // we store pointers to PlotSandbox with the same constantness as it.
  using PlotSet_t
    = util::with_const_as_t<icarus::trigger::PlotSandbox, PlotSetType>;
  
  std::vector<PlotSet_t*> destinations;
  destinations.reserve(plotSets.size()); // may be too much
  for (auto&& plotSet: plotSets) { // unique_ptr<PlotDispatcherBase> [const]
    if (event && !plotSet->passes(*event)) continue;
    destinations.push_back(&(plotSet->sandbox()));
  }
  return destinations;
} // icarus::trigger::TriggerDesignPlots::collectPlotSetsImpl()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerDesignPlots::collectPlotSets() const
  -> ConstPlotDestinations_t
  { return collectPlotSetsImpl(fPlotSets); }

//------------------------------------------------------------------------------
auto icarus::trigger::TriggerDesignPlots::collectPlotSets()
  -> PlotDestinations_t
  { return collectPlotSetsImpl(fPlotSets); }


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerDesignPlots::collectPlotSets
  (art::Event const& event) const
  -> ConstPlotDestinations_t
  { return collectPlotSetsImpl(fPlotSets, &event); }

//------------------------------------------------------------------------------
auto icarus::trigger::TriggerDesignPlots::collectPlotSets
  (art::Event const& event)
  -> PlotDestinations_t
  { return collectPlotSetsImpl(fPlotSets, &event); }

#endif // 0

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::DiscriminatePMTwaveforms)


//------------------------------------------------------------------------------
