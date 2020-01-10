/**
 * @file   TriggerEfficiencyPlots_module.cc
 * @brief  Plots of trigger efficiency.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 9, 2020
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // FillTriggerGates()
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icaruscode/PMT/Trigger/Utilities/ROOTutils.h" // util::ROOT
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"
#if 0
#include "icaruscode/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()
#endif // 0

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds, ...
#include "lardataalg/Utilities/quantities/energy.h" // gigaelectronvolt // FIXME
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/CoreUtils/get_elements.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#if 0
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
// #include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::::static_assert_on<>
#include "lardataobj/RawData/OpDetWaveform.h"
#endif // 0

// nutools libraries
#include "nusimdata/SimulationBase/MCTruth.h"
// #include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"
#if 0
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"
#endif // 0

// ROOT libraries
#include "TEfficiency.h"
#include "TH1F.h"
#include "TH2F.h"

// C/C++ standard libraries
#include <ostream>
#include <map>
#include <vector>
#include <string>
#include <functional> // std::function<>, std::reference_wrapper<>
#include <memory> // std::unique_ptr
#include <utility> // std::pair<>, std::move()
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
// using GeV = double;
// FIXME: replace by
using GeV = util::quantities::gigaelectronvolt;
// when `lardataalg/Utilities/quantities/energy.h` is released

using namespace util::quantities::time_literals; // ""_ns ...


//------------------------------------------------------------------------------
/// Information about the event.
struct EventInfo_t {

  // --- BEGIN Query interface -----------------------------------------------
  /// @name Query interface
  /// @{

  /// Returns whether the event is generated as a neutrino interaction.
  bool isWeakChargedCurrent() const { return fInteractionTypeMask & itWCC; }

  /// Returns whether the event is generated as a neutrino interaction.
  bool isWeakNeutralCurrent() const { return fInteractionTypeMask & itWNC; }

  /// Returns whether the event is generated as a neutrino interaction.
  bool isNeutrino() const { return fInteractionTypeMask & itWC; }
  
  /// Returns the total energy deposited in the detector during the event [GeV]
  GeV DepositedEnergy() const { return fEnergyDepTotal; }
  
  /// Returns the total energy deposited in the detector during beam [GeV]
  GeV DepositedEnergyInSpill() const { return fEnergyDepSpill; }

  /// @}
  // --- END Query interface -------------------------------------------------


  // --- BEGIN Set interface -------------------------------------------------
  /// @name Set interface
  /// @{

  /// Marks this event as including a weak charged current interaction.
  void SetWeakChargedCurrent() { fInteractionTypeMask |= itWCC; }

  /// Marks this event as including a weak neutral current interaction.
  void SetWeakNeutralCurrent() { fInteractionTypeMask |= itWNC; }

  /// Sets the total deposited energy of the event [GeV]
  void SetDepositedEnergy(GeV e) { fEnergyDepTotal = e; }

  /// Sets the energy of the event deposited during beam gate [GeV]
  void SetDepositedEnergyInSpill(GeV e) { fEnergyDepSpill = e; }

  /// @}
  // --- END Set interface ---------------------------------------------------

  /// Prints the content of the object into a stream.
  void dump(std::ostream& out) const
    {
      out << "Event contains:";
      if (isNeutrino()) {
        if (isWeakChargedCurrent()) out << " CC";
        if (isWeakNeutralCurrent()) out << " NC";
      }
      else {
        out << " no neutrino interaction";
      }
      out << "\nTotal deposited energy: " << DepositedEnergy()
        << ", of which in spill " << DepositedEnergyInSpill();
      out << "\n";
      out.flush();
    } // dump()

    private:

  // --- BEGIN interaction type constants ------------------------------------
  using InteractionType_t = std::size_t;
  static constexpr InteractionType_t itWCC = 0x01; ///< Charged weak current.
  static constexpr InteractionType_t itWNC = 0x02; ///< Neutral weak current.

  static constexpr InteractionType_t itWC = itWNC | itWCC; ///< Neutral current.

  // --- END interaction type constants --------------------------------------

  InteractionType_t fInteractionTypeMask {}; ///< Type of interaction.

  GeV fEnergyDepTotal { 0.0 }; ///< Total deposited energy [GeV]
  GeV fEnergyDepSpill { 0.0 }; ///< Total deposited energy [GeV]


}; // struct EventInfo_t

std::ostream& operator<< (std::ostream& out, EventInfo_t const& info)
  { info.dump(out); return out; }



//------------------------------------------------------------------------------
class PlotCategory {
  
    public:
  
  /// Type of test function.
  using QualifyFunc_t = std::function<bool(EventInfo_t const&)>;
  
  
  /// Constructor from category name and test function.
  PlotCategory(
    std::string name, std::string descr = {},
    QualifyFunc_t&& test = [](EventInfo_t const&){ return true; }
    )
    : fName(std::move(name)), fDescr(std::move(descr)), fTest(std::move(test))
    {}
  
  /// Returns the name of the category.
  std::string const& name() const { return fName; }
  
  /// Returns the description of the category.
  std::string const& description() const { return fDescr; }
  
  /// Returns whether the event belong to this category.
  bool test(EventInfo_t const& info) const { return fTest(info); }
  
  operator std::string() const { return name(); }
  bool operator() (EventInfo_t const& info) const { return test(info); }
  
    private:
  
  std::string fName;
  std::string fDescr;
  QualifyFunc_t fTest;
  
}; // class PlotCategory
using PlotCategories_t = std::vector<PlotCategory>;


PlotCategories_t const PlotCategories {

  PlotCategory{
    "All"
    },

  PlotCategory{
    "NuCC", "CC",
    [](EventInfo_t const& info){ return info.isWeakChargedCurrent(); }
    },

  PlotCategory{
    "NuNC", "NC",
    [](EventInfo_t const& info){ return info.isWeakNeutralCurrent(); }
    }

}; // PlotCategories[]


//------------------------------------------------------------------------------
namespace icarus::trigger { class TriggerEfficiencyPlots; }
/**
 * @brief Produces plots to inform trigger design.
 * 
 * This module produces sets of plots based on the configured trigger settings.
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (labels out of
 *   `TriggerGatesTag` and `Thresholds`): full sets of discriminated waveforms,
 *   each waveform possibly covering multiple optical channels,
 *   and their associations to optical waveforms. One set per threshold.
 * * `std::vector<simb::MCTruth>`: generator information, used for categorising
 *      the events for plot sets
 * * `std::vector<simb::MCParticle>`: simulation information, used for
 *      categorizing the events for plotting
 * 
 * 
 * Output plots
 * =============
 * 
 * For each event category, a set of plots is left into a ROOT subdirectory.
 * 
 * @todo Document which plots these are!
 * 
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description TriggerEfficiencyPlots`.
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module
 *     instance which produced the trigger primitives to be used as input;
 *     it must not include any instance name, as the instance names will be
 *     automatically added from `Thresholds` parameter.
 *     The typical trigger primitives used as input may be LVDS discriminated
 *     output (e.g. from `icarus::trigger::LVDSgates` module) or combinations
 *     of them (e.g. from `icarus::trigger::SlidingWindowTrigger` module).
 * * `Thresholds` (list of integers, mandatory): list of the discrimination
 *     thresholds to consider, in ADC counts. A data product containing a
 *     digital signal is read for each one of the thresholds, and the tag of the
 *     data product is expected to be the module label `TriggerGatesTag` with as
 *     instance name the value of the threshold (e.g. for a threshold of 60 ADC
 *     counts the data product tag might be `LVDSgates:60`).
 * * `LogCategory` (string, default `TriggerEfficiencyPlots`): name of category
 *     used to stream messages from this module into message facility.
 * * `PropagatedParticles` (input tag, default: `largeant`): data product
 *     containing the list of particles going through the argon in the detector.
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
class icarus::trigger::TriggerEfficiencyPlots: public art::EDAnalyzer {

  using microseconds = util::quantities::intervals::microseconds;
  using nanoseconds = util::quantities::intervals::nanoseconds;

    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> GeneratorTags {
      Name("GeneratorTags"),
      Comment("labels of the event generators"),
      std::vector<art::InputTag>{ "generator" }
      };

    fhicl::Atom<art::InputTag> DetectorParticleTag {
      Name("DetectorParticleTag"),
      Comment("label of simulated particles through the detector"),
      "largeant" // tradition demands
      };

    fhicl::Sequence<art::InputTag> EnergyDepositTags {
      Name("EnergyDeposits"),
      Comment("label of energy deposition data product(s) in the detector"),
      std::vector<art::InputTag>{ "largeant:TPCActive" }
      };

    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    fhicl::Sequence<raw::ADC_Count_t> Thresholds {
      Name("Thresholds"),
      Comment("thresholds to consider [ADC counts]")
      };

    fhicl::Atom<microseconds> BeamGateDuration {
      Name("BeamGateDuration"),
      Comment("length of time interval when optical triggers are accepted [us]")
      };

    fhicl::Sequence<unsigned int> MinimumPrimitives {
      Name("MinimumPrimitives"),
      Comment("minimum required number of trigger primitives for a trigger")
      };
    
    fhicl::Atom<nanoseconds> TriggerTimeResolution {
      Name("TriggerTimeResolution"),
      Comment("resolution of trigger in time [ns]"),
      8_ns
      };

    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "SlidingWindowTrigger" // default
      };

  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  explicit TriggerEfficiencyPlots(Parameters const& config);

  // Plugins should not be copied or assigned.
  TriggerEfficiencyPlots(TriggerEfficiencyPlots const&) = delete;
  TriggerEfficiencyPlots(TriggerEfficiencyPlots&&) = delete;
  TriggerEfficiencyPlots& operator=(TriggerEfficiencyPlots const&) = delete;
  TriggerEfficiencyPlots& operator=(TriggerEfficiencyPlots&&) = delete;

  // --- END Constructors ------------------------------------------------------


  // --- BEGIN Framework hooks -------------------------------------------------

  /// Fills the plots. Also extracts the information to fill them with.
  virtual void analyze(art::Event const& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  /// Type of gate data without channel information.
  using TriggerGateData_t
    = icarus::trigger::OpticalTriggerGateData_t::GateData_t;
  
  using PlotSandbox = icarus::trigger::PlotSandbox; ///< Import type.
  
  /// List of references to plot sandboxes.
  using PlotSandboxRefs_t
    = std::vector<std::reference_wrapper<PlotSandbox const>>;
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  std::vector<art::InputTag> fGeneratorTags; ///< Generator data product tags.
  art::InputTag fDetectorParticleTag; ///< Input simulated particles.
  /// Energy deposition data product tags.
  std::vector<art::InputTag> fEnergyDepositTags;
  
  /// ADC thresholds to read, and the input tag connected to their data.
  std::map<icarus::trigger::ADCCounts_t, art::InputTag> fADCthresholds;
  
  /// Duration of the gate during with global optical triggers are accepted.
  microseconds fBeamGateDuration;
  
  /// Minimum number of trigger primitives for a trigger to happen.
  std::vector<unsigned int> fMinimumPrimitives;
  
  nanoseconds fTriggerTimeResolution; ///< Trigger resolution in time.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;

  // --- END Configuration variables -------------------------------------------


  // --- BEGIN Service variables -----------------------------------------------

  detinfo::DetectorClocks const& fDetClocks;
  detinfo::DetectorTimings fDetTimings;

  /// ROOT directory where all the plots are written.
  art::TFileDirectory fOutputDir;

  // --- END Service variables -------------------------------------------------

  
  // --- BEGIN Internal variables ----------------------------------------------
  
  /// Gate representing the time we expect light from beam interactions.
  icarus::trigger::OpticalTriggerGate const fBeamGate;
  
  /// Beam gate start and stop tick in optical detector scale.
  std::pair
    <detinfo::timescales::optical_tick, detinfo::timescales::optical_tick>
    const fBeamGateOpt;

  /// Beam gate start and stop time in simulation scale.
  std::pair
    <detinfo::timescales::simulation_time, detinfo::timescales::simulation_time>
    const fBeamGateSim;

  /// All plots, one set per ADC threshold.
  std::vector<PlotSandbox> fThresholdPlots;

  // --- END Internal variables ------------------------------------------------


  /// Initializes all the plot sets, one per threshold.
  void initializePlots(PlotCategories_t const& categories);

  /// Initializes a single plot set (ADC threshold + category) into `plots`.
  void initializePlotSet(PlotSandbox& plots) const;

  /// Returns the names of the plot categories event qualifies for.
  std::vector<std::string> selectPlotCategories
    (EventInfo_t const& info, PlotCategories_t const& categories) const;

  /// Extracts from `event` the relevant information on the physics event.
  EventInfo_t extractEventInfo(art::Event const& event) const;
  
  /// Computes the trigger response from primitives with the given `threshold`.
  TriggerGateData_t combineTriggerPrimitives(
    art::Event const& event,
    icarus::trigger::ADCCounts_t const threshold,
    art::InputTag const& dataTag
    ) const;

  /// Applies the beam gate coincidence to the specified trigger primitive.
  template <typename GateObject>
  GateObject applyBeamGate(GateObject gate) const;

  /// Adds all the `responses` (one per threshold) to the plots.
  void plotResponses(
    icarus::trigger::ADCCounts_t const threshold,
    PlotSandboxRefs_t const& plots, EventInfo_t const& info,
    TriggerGateData_t const& primitiveCount
    ) const;

  /// Reads a set of input gates from the `event`.
  std::vector<icarus::trigger::MultiChannelOpticalTriggerGate> ReadTriggerGates
    (art::Event const& event, art::InputTag const& dataTag) const;
  
  static std::string thrAndCatName
    (std::string const& boxName, std::string const& category)
    { return boxName + "_" + category; }
  static std::string thrAndCatName
    (PlotSandbox const& box, std::string const& category)
    { return thrAndCatName(box.name(), category); }
  
}; // icarus::trigger::TriggerEfficiencyPlots


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::TriggerEfficiencyPlots
//------------------------------------------------------------------------------
icarus::trigger::TriggerEfficiencyPlots::TriggerEfficiencyPlots
  (Parameters const& config)
  : art::EDAnalyzer(config)
  // configuration
  , fGeneratorTags        (config().GeneratorTags())
  , fDetectorParticleTag  (config().DetectorParticleTag())
  , fEnergyDepositTags    (config().EnergyDepositTags())
  , fBeamGateDuration     (config().BeamGateDuration())
  , fMinimumPrimitives    (config().MinimumPrimitives())
  , fTriggerTimeResolution(config().TriggerTimeResolution())
  , fLogCategory          (config().LogCategory())
  // services
  , fDetClocks (*lar::providerFrom<detinfo::DetectorClocksService>())
  , fDetTimings(fDetClocks)
  , fOutputDir (*art::ServiceHandle<art::TFileService>())
  // cached
  , fBeamGate(icarus::trigger::BeamGateMaker{ fDetClocks }(fBeamGateDuration))
  , fBeamGateOpt(
      fDetTimings.toOpticalTick(fDetTimings.BeamGateTime()),
      fDetTimings.toOpticalTick(fDetTimings.BeamGateTime() + fBeamGateDuration)
    )
  , fBeamGateSim(
      fDetTimings.toSimulationTime(fBeamGateOpt.first),
      fDetTimings.toSimulationTime(fDetTimings.BeamGateTime())
        + fBeamGateDuration
    )
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
  // input data declaration
  //
  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // event information
  for (art::InputTag const& inputTag: fGeneratorTags)
    consumes<std::vector<simb::MCTruth>>(inputTag);
  for (art::InputTag const& inputTag: fEnergyDepositTags)
    consumes<std::vector<sim::SimEnergyDeposit>>(inputTag);
//   consumes<std::vector<simb::MCParticle>>(fDetectorParticleTag);
  
  // trigger primitives
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    consumes<std::vector<OpticalTriggerGateData_t>>(inputDataTag);
    consumes<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (inputDataTag);
  } // for

  //
  // initialization of plots and plot categories
  //
  initializePlots(::PlotCategories);

  {
    mf::LogInfo log(fLogCategory);
    log << "\nConfigured " << fADCthresholds.size() << " thresholds:";
    for (auto const& [ threshold, dataTag ]: fADCthresholds)
      log << "\n * " << threshold << " ADC (from '" << dataTag.encode() << "')";
    log << "\nBeam gate is " << fBeamGate << " (" << fBeamGateSim.first
      << " -- " << fBeamGateSim.second << ")";
  } // local block

} // icarus::trigger::TriggerEfficiencyPlots::TriggerEfficiencyPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::analyze(art::Event const& event) {

  /*
   * 1. find out the features of the event and the categories it belongs to
   * 2. for each threshold:
   *   1. read the trigger primitives
   *   2. apply the beam gate on the primitives
   *   3. (optional) combine the trigger primitives
   *   4. generate the trigger response
   *   5. pick the plots to be filled
   *   6. add the response to all the plots
   *
   */

  //
  // 1. find out the features of the event and the categories it belongs to
  //
  EventInfo_t const eventInfo = extractEventInfo(event);
  
  std::vector<std::string> selectedPlotCategories
    = selectPlotCategories(eventInfo, ::PlotCategories);
  {
    mf::LogTrace log(fLogCategory);
    log
      << "Event " << event.id() << " falls in " << selectedPlotCategories.size()
      << " categories:"
      ;
    for (std::string const& name: selectedPlotCategories)
      log << " \"" << name << "\"";
    
  } // local block

  //
  // 2. for each threshold:
  //
  for (auto&& [ thrPair, thrPlots ]: util::zip(fADCthresholds, fThresholdPlots))
  {

    auto const& [ threshold, dataTag ] = thrPair;

    //
    // 1.1. read the trigger primitives
    // 1.2. (optional) combine the trigger primitives
    // 1.3. apply the beam gate on the primitives
    // 1.4. generate the trigger response
    //
    TriggerGateData_t const primitiveCount
      = applyBeamGate(combineTriggerPrimitives(event, threshold, dataTag));

    //
    // 1.5. pick the plots to be filled
    //
    PlotSandboxRefs_t selectedPlots;
    
    for (std::string const& name: selectedPlotCategories)
      selectedPlots.emplace_back(*(thrPlots.findSandbox(name)));
    
    // 1.6. add the response to the appropriate plots
    plotResponses(threshold, selectedPlots, eventInfo, primitiveCount);
    
  } // for


} // icarus::trigger::TriggerEfficiencyPlots::analyze()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::initializePlots
  (PlotCategories_t const& categories)
{
  using namespace std::string_literals;
  
  for (icarus::trigger::ADCCounts_t const threshold
    : util::get_elements<0U>(fADCthresholds))
  {
    // create a plot sandbox inside `fOutputDir` with a name/prefix `Thr###`
    auto const thr = threshold.value();
    icarus::trigger::PlotSandbox thrPlots { fOutputDir,
      "Thr"s + util::to_string(thr), "(thr: "s + util::to_string(thr) + ")"s };
    
    // create a subbox for each plot category
    for (PlotCategory const& category: categories) {
      PlotSandbox& plots = thrPlots.addSubSandbox(
        category.name(),
        category.description()
        );
      
      initializePlotSet(plots);
    }
    fThresholdPlots.push_back(std::move(thrPlots));
  } // for thresholds
  
  mf::LogTrace log(fLogCategory);
  log << "Created " << fThresholdPlots.size() << " plot boxes:\n";
  for (auto const& box: fThresholdPlots) {
    box.dump(log, "  ");
  } // for
  
} // icarus::trigger::TriggerEfficiencyPlots::initializePlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::initializePlotSet
  (PlotSandbox& plots) const
{
  
  // a variable binning for the required number of trigger primitives
  auto [ minimumPrimBinning, minimumPrimBinningLabels ]
    = util::ROOT::makeVariableBinningAndLabels(fMinimumPrimitives);
  assert(minimumPrimBinning.size() == minimumPrimBinningLabels.size() + 1U);

  {
    mf::LogTrace log(fLogCategory);
    log << "TriggerEfficiencyPlots (plots '"
      << plots.name() << "') variable binning including the "
      << fMinimumPrimitives.size() << " points {";
    for (auto value: fMinimumPrimitives) log << " " << value;
    log << " } => " << minimumPrimBinningLabels.size() << " bins: ";
    for (auto const& [ value, label ]
      : util::zip<1U>(minimumPrimBinning, minimumPrimBinningLabels))
    {
      log << " " << value << " (\"" << label << "\") =>";
    } // for
    log << " " << minimumPrimBinning.back();
  } // debug output block

  //
  // Triggering efficiency vs. threshold.
  //
  detinfo::timescales::optical_time_ticks const triggerResolutionTicks
    { fDetTimings.toOpticalTicks(fTriggerTimeResolution) };
  auto* TrigTime = plots.make<TH2F>(
    "TriggerTick",
    "Trigger time tick"
      ";minimum requested number of trigger primitives"
      ";optical time tick [ /" + util::to_string(triggerResolutionTicks) + " ]",
    minimumPrimBinning.size() - 1U, minimumPrimBinning.data(),
//    fMinimumPrimitives.back(), 0, fMinimumPrimitives.back() + 1
    (fBeamGateOpt.second - fBeamGateOpt.first) / triggerResolutionTicks,
    fBeamGateOpt.first.value(), fBeamGateOpt.second.value()
    );
  
  util::ROOT::applyAxisLabels(TrigTime->GetXaxis(), minimumPrimBinningLabels);
  
  auto* Eff = plots.make<TEfficiency>(
    "Eff",
    "Efficiency of triggering"
      ";minimum requested number of trigger primitives"
      ";trigger efficiency",
    minimumPrimBinning.size() - 1U, minimumPrimBinning.data()
//    fMinimumPrimitives.back(), 0, fMinimumPrimitives.back() + 1
    );
  
  // people are said to have earned hell for things like this;
  // but TEfficiency really does not expose the interface to assign labels to
  // its axes, which supposedly could be done had we chosen to create it by
  // histograms instead of directly as recommended.
  util::ROOT::applyAxisLabels(
    const_cast<TH1*>(Eff->GetTotalHistogram())->GetXaxis(),
    minimumPrimBinningLabels
    );
  
  //
  // plots independent of the trigger primitive requirements
  //
  plots.make<TH1F>(
    "NPrimitives",
    "Number of trigger primitives (\"channels firing at once\")"
    ";maximum trigger primitives at the same time"
    ";events",
    192, 0.0, 192.0 // large number, zoom in presentations!
    );
  
  //
  // Selection-related plots
  //
  plots.make<TH1F>(
    "EnergyInSpill",
    "Energy eposited during the beam gate opening"
      ";deposited energy  [ GeV ]"
      ";events  [ / 50 MeV ]",
    120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
    );
  
} // icarus::trigger::TriggerEfficiencyPlots::initializePlotSet()


//------------------------------------------------------------------------------
std::vector<std::string>
icarus::trigger::TriggerEfficiencyPlots::selectPlotCategories
  (EventInfo_t const& info, PlotCategories_t const& categories) const
{
  std::vector<std::string> selected;
  
  for (auto const& category: categories)
    if (category(info)) selected.push_back(category);
  
  return selected;
  
} // icarus::trigger::TriggerEfficiencyPlots::selectPlotCategories()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlots::extractEventInfo
  (art::Event const& event) const -> EventInfo_t
{
  
  EventInfo_t info;
  
  //
  // generator information
  //
  for (art::InputTag const& inputTag: fGeneratorTags) {
  
    auto const& truthRecords
      = *(event.getValidHandle<std::vector<simb::MCTruth>>(inputTag));
    
    for (simb::MCTruth const& truth: truthRecords) {
      
      if (truth.NeutrinoSet()) {
        //
        // interaction type (CC, NC)
        //
        switch (truth.GetNeutrino().CCNC()) {
          case simb::kCC: info.SetWeakChargedCurrent(); break;
          case simb::kNC: info.SetWeakNeutralCurrent(); break;
          default:
            mf::LogWarning(fLogCategory)
              << "Event " << event.id() << " has unexpected NC/CC flag ("
              << truth.GetNeutrino().CCNC() << ")";
        } // switch
        
      } // if neutrino event
      
    } // for truth records
    
  } // for generators
  
  //
  // propagation in the detector
  //
  GeV totalEnergy { 0.0 }, inSpillEnergy { 0.0 };
  
  for (art::InputTag const& edepTag: fEnergyDepositTags) {
    
    auto const& energyDeposits
      = *(event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(edepTag));
    mf::LogTrace(fLogCategory)
      << "Event " << event.id() << " has " << energyDeposits.size()
      << " energy deposits recorded in '" << edepTag.encode() << "'";
    
    for (sim::SimEnergyDeposit const& edep: energyDeposits) {
      
      GeV const e { edep.Energy() }; // assuming it's stored in GeV
      
      detinfo::timescales::simulation_time const t { edep.Time() };
      bool const inSpill
        = (t >= fBeamGateSim.first) && (t <= fBeamGateSim.second);
      
      totalEnergy += e;
      if (inSpill) inSpillEnergy += e;
      
    } // for all energy deposits in the data product
    
  } // for all energy deposit tags
  
  info.SetDepositedEnergy(totalEnergy);
  info.SetDepositedEnergyInSpill(inSpillEnergy);
  
  mf::LogTrace(fLogCategory) << "Event " << event.id() << ": " << info;
  
  return info;
} // icarus::trigger::TriggerEfficiencyPlots::extractEventInfo()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlots::combineTriggerPrimitives(
  art::Event const& event,
  icarus::trigger::ADCCounts_t const threshold,
  art::InputTag const& dataTag
) const -> TriggerGateData_t {

  //
  // simple count
  // TODO replace by a art tool
  //

  std::vector<icarus::trigger::MultiChannelOpticalTriggerGate> const& gates
    = ReadTriggerGates(event, dataTag);

  mf::LogTrace(fLogCategory)
    << "Simulating trigger response with ADC threshold " << threshold
    << " from '" << dataTag.encode() << "' (" << gates.size() << " primitives)";

  TriggerGateData_t combinedGate;

  auto itGate = gates.begin();
  auto const gend = gates.end();
  if (itGate == gend) { // this is unexpected...
    mf::LogWarning(fLogCategory)
      << "  No trigger primitive found for threshold " << threshold
      << " ('" << dataTag.encode() << "')";
    return {};
  } // if no gates

  combinedGate = *itGate;
  while (++itGate != gend) combinedGate.Sum(*itGate);

  return combinedGate;
} // icarus::trigger::TriggerEfficiencyPlots::combineTriggerPrimitives()


//------------------------------------------------------------------------------
template <typename GateObject>
GateObject icarus::trigger::TriggerEfficiencyPlots::applyBeamGate
  (GateObject gate) const
{
  return std::move(gate.Mul(fBeamGate));
} // icarus::trigger::TriggerEfficiencyPlots::applyBeamGate()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::plotResponses(
  icarus::trigger::ADCCounts_t const threshold,
  PlotSandboxRefs_t const& plotSets,
  EventInfo_t const& eventInfo,
  TriggerGateData_t const& primitiveCount
) const {
  
  /*
   * This function plots according to the configured minimum number of trigger
   * primitives: for each requirement of minimum number of primitives, the
   * earliest time where that requirement is met is found, and that is
   * considered as the trigger time.
   * 
   * The following quantities are drawn per ADC threshold and per plot category:
   * 
   * * per minimum number of primitives:
   *    * trigger efficiency (i.e. whether there was _any time_ a number of
   *      primitives fulfilling the requirement)
   *    * trigger time in ticks (distribution as 2D histogram)
   * * maximum number of trigger primitives present at any time
   * * deposited energy during beam spill
   * 
   */
  using namespace std::string_literals;
  
  using ClockTick_t = TriggerGateData_t::ClockTick_t;
  using OpeningCount_t = TriggerGateData_t::OpeningCount_t;
  
  using PrimitiveCount_t = std::pair<ClockTick_t, OpeningCount_t>;
  
  // largest number of trigger primitives at any time
  detinfo::DetectorTimings::optical_tick const maxPrimitiveTime
    { primitiveCount.findMaxOpen() };
  PrimitiveCount_t const maxPrimitives { 
    maxPrimitiveTime.value(),
    primitiveCount.openingCount(maxPrimitiveTime.value())
    };

  mf::LogTrace(fLogCategory)
    << "Max primitive count in " << threshold << ": "
    << maxPrimitives.second << " at tick " << maxPrimitives.first
    << " (" << fDetTimings.toElectronicsTime(maxPrimitiveTime) << ")"
    ;
  
  /*
   * Fill all the histograms for all the minimum primitive requirements
   * (filling the information whether or not the trigger fired),
   * for all the qualifying categories.
   * Note that in this type of plots each event appears in all bins
   * (may be with "fired" or "not fired" on each bin)
   */
  PrimitiveCount_t lastMinCount { TriggerGateData_t::MinTick, 0 };
  bool fired = true;
  for (OpeningCount_t minCount: fMinimumPrimitives) {
    
    if (fired && (lastMinCount.second < minCount)) {
      // if we haven't passed this minimum yet
      ClockTick_t const time = primitiveCount.findOpen(minCount);
      if (time == TriggerGateData_t::MaxTick) {
        mf::LogTrace(fLogCategory)
          << "Never got at " << minCount << " primitives or above.";
        fired = false;
      }
      lastMinCount = { time, primitiveCount.openingCount(time) };
      mf::LogTrace(fLogCategory)
        << "Reached " << minCount << " primitives or above ("
        << lastMinCount.second << ") at " << lastMinCount.first << ".";
    } // if
    
    // at this point we know we have minCount or more trigger primitives,
    // and the time of this one is in lastMinCount.first (just in case)
    
    // go through all the plot categories this event qualifies for
    for (icarus::trigger::PlotSandbox const& plotSet: plotSets) {
      
      auto getEff = [&plotSet](std::string const& name)
        { return plotSet.use<TEfficiency>(name); };
      auto getHist2D = [&plotSet](std::string const& name)
        { return plotSet.use<TH2>(name); };
      
      // simple efficiency
      getEff("Eff"s)->Fill(fired, minCount);

      // trigger time (if any)
      if (fired)
        getHist2D("TriggerTick"s)->Fill(minCount, lastMinCount.first);

    } // for all qualifying plot categories
    
  } // for all thresholds
  
  /*
   * Now fill the plots independent of the minimum trigger primitives:
   * the same value is plotted in all plot sets.
   * 
   */
  for (icarus::trigger::PlotSandbox const& plotSet: plotSets) {
    auto getHist = [&plotSet](std::string const& name)
      { return plotSet.use<TH1>(name); };
    
    // selection-related plots:
    getHist("EnergyInSpill"s)
      ->Fill(double(eventInfo.DepositedEnergyInSpill()));
    
    // number of primitives
    getHist("NPrimitives"s)->Fill(maxPrimitives.second);
    
  } // for 

} // icarus::trigger::TriggerEfficiencyPlots::plotResponses()


//------------------------------------------------------------------------------
std::vector<icarus::trigger::MultiChannelOpticalTriggerGate>
icarus::trigger::TriggerEfficiencyPlots::ReadTriggerGates
  (art::Event const& event, art::InputTag const& dataTag) const
{

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // currently the associations are a waste (we include it not to have
  // potentially broken `MultiChannelOpticalTriggerGate` objects, but at the
  // time this module was written that was not important)
  auto const& gates
    = *(event.getValidHandle<std::vector<OpticalTriggerGateData_t>>(dataTag));
  auto const& gateToWaveforms = *(
    event.getValidHandle
      <art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>(dataTag)
    );
  try {
    return icarus::trigger::FillTriggerGates
      <icarus::trigger::MultiChannelOpticalTriggerGate>(gates, gateToWaveforms);
  }
  catch (cet::exception const& e) {
    throw cet::exception("TriggerEfficiencyPlots", "", e)
      << "Error encountered while reading data products from '"
      << dataTag.encode() << "'\n";
  }

} // icarus::trigger::TriggerEfficiencyPlots::ReadTriggerGates()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::TriggerEfficiencyPlots)


//------------------------------------------------------------------------------
