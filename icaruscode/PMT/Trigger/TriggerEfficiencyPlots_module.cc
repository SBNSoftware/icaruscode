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

/**
 * @brief List of event categories.
 * 
 * category name  | condition
 * -------------- | ------------------------------------------------------------
 * `All`          | any event
 * `NuCC`         | at least one generated charged current neutrino interaction
 * `NuNC`         | at least one generated neutral current neutrino interaction
 * 
 * 
 */
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
 * @brief Produces plots about trigger simulation and trigger efficiency.
 * 
 * This module produces sets of plots based on trigger primitives given in
 * input.
 * 
 * The following documentation mainly deals with the most standard configuration
 * and operations in ICARUS. This module and the ones upstream of it have quite
 * some knobs that can be manipulated to test unorthodox configurations.
 * 
 * 
 * Overview of trigger simulation steps
 * =====================================
 * 
 * This simulation only trigger primitives derived from the optical detector
 * via hardware (V1730B boards).
 * The trigger simulation branches from the standard simulation of optical
 * detector waveforms (`icarus::SimPMTIcarus` module).
 * From there, multiple steps are needed.
 * 
 * 1. Produce single-PMT-channel discriminated waveforms: the discrimination
 *    threshold can be chosen ("ADC threshold" or just "threshold" in the
 *    following), and the result is one binary discriminated waveform that
 *    has value `0` when the waveform is under threshold and `1` when it's
 *    above threshold, with the same time discretization as the PMT waveform.
 *    Each discriminated waveform spans the whole event readout time and
 *    therefore it may merge information from multiple readout waveforms
 *    (`raw::OpDetWaveform`), but all from the same optical detector channel.
 *    This step can be performed with the module
 *    `icarus::trigger::DiscriminatePMTwaveforms`.
 * 2. Combine the discriminated waveforms in pairs, in the same fashion as the
 *    V1730B readout board does to produce LVDS output. The output of this step
 *    is roughly half as many discriminated outputs as there were from the
 *    previous step (some channels are not paired). This output will be called
 *    _trigger primitives_ because it is what the trigger logics is based on.
 *    We say that a trigger primitive is "on" when its level is `1`.
 *    This step can be performed with the module `icarus::trigger::LVDSgates`.
 * 3. Simulate the trigger logic based on the trigger primitives _(see below)_.
 *    This usually includes the requirement of coincidence with the beam gate.
 * 
 * Trigger logic may be complex, being implemented in a FPGA.
 * Many options are available, including:
 * 
 * * coincidence of a minimum number of trigger primitives: the event is trigger
 *   if at least _N_ trigger primitives are on at the same time;
 * * sliding window: primitives are grouped depending on their location, and
 *   there is a requirement on the minimum number of primitives in "on" level
 *   within one or more of these groups. For example, any sliding window of 30
 *   channels (corresponding to 16 trigger primitives in standard ICARUS PMT
 *   readout installation) has to contain at least 10 trigger primitives "on"
 *   at the same time; or there must be two sliding windows with that condition;
 *   etc. Sliding windows can be obtained from further processing of the LVDS
 *   trigger primitives, for example with the module
 *   `icarus::trigger::SlidingWindowTrigger`.
 * 
 * This module (`icarus::trigger::TriggerEfficiencyPlots`) is designed for
 * implementing a trigger logic based on coincidence of trigger primitives,
 * and to plot trigger efficiency and related information based on a requirement
 * of a minimum number of trigger primitives, in that enabling the simulation
 * in the fashion of the first option above. More details on the steps performed
 * by this module are documented @ref TriggerEfficiencyPlots_Algorithm "later".
 * 
 * For the implementation of sliding window, different logic needs to be
 * implemented, possibly keeping track of the location of each primitive, and
 * different plots are needed as well. This module may provide a template for
 * such implementation, but it hasn't been designed with enough flexibility to
 * accommodate that directly.
 * 
 * 
 * Data objects for discriminated waveforms
 * -----------------------------------------
 * 
 * A discriminated waveform is the information whether the level of a waveform
 * is beyond threshold, as function of time.
 * A discriminated waveform may be binary, i.e. with only levels `0` and `1`
 * based on a single threshold, or with multiple levels.
 * Also the numerical _addition_ of two binary discriminated waveforms
 * yields a multi-level waveform (in fact, three levels -- `0`, `1` and `2`).
 * 
 * We represent this data in the forms of "events": an event is a change of
 * level happening at a certain time. The class holding this information,
 * `icarus::trigger::TriggerGateData`, covers the whole time, starting with a
 * ground level `0`. The next event will be a "opening" that increases the
 * level, usually to `1`. Other changing events may follow, and typically the
 * last one will bring the level back to `0`.
 * 
 * This information is joined by a list of _channel numbers_ in order to
 * represent a discriminated waveform e.g. from the optical detector.
 * There may be one or more channels associated to a discriminated waveform,
 * but for us there should always be at least one.
 * The discriminated waveforms from PMT readout (`raw::OpDetWaveform`) are
 * associated to a single channel, while LVDS trigger primitives are most often
 * associated to two channels (some channels are not paired and will have only
 * one channel). A sliding window will have as many channels as the PMT it
 * covers. The global trigger should have _all_ channels, while in ICARUS each
 * of the two discriminated wabeforms from a cryostat should have half the
 * channels.
 * This information is represented in the class
 * `icarus::trigger::ReadoutTriggerGate`, which inherits from
 * `icarus::trigger::TriggerGateData`.
 * This class is generic and can hold any representation for the time of the
 * level changing events, for the levels, and for the identifiers of the
 * channels. ICARUS specifies a data type for each of these quantities, and
 * the resulting `icarus::trigger::ReadoutTriggerGate` class instance is called
 * `icarus::trigger::OpticalTriggerGateData_t`.
 * 
 * @note The class `icarus::trigger::OpticalTriggerGateData_t` is the one that
 *       gets written to disk in the _art_ ROOT files.
 *       That is not a class by itself, but rather an alias of
 *       `icarus::trigger::ReadoutTriggerGate`, and many coding tools will call
 *       it in the latter way.
 * 
 * The class `icarus::trigger::OpticalTriggerGate` is currently the most
 * commonly used in the code. It adds to the information of
 * `icarus::trigger::ReadoutTriggerGate`, from which it derives, a list of
 * optical waveforms (`raw::OpDetWaveform`) it originates from.
 * 
 * Finally, the classes `icarus::trigger::SingleChannelOpticalTriggerGate` and
 * `icarus::trigger::MultiChannelOpticalTriggerGate` do not add any information
 * to the `icarus::trigger::OpticalTriggerGate` they derive from, but they
 * have an interface explicitly tuned for discriminated waveforms from a single
 * channel or from multiple channels, respectively (for example, the former
 * provides a `channel()` method returning a single channel identifier, while
 * the latter provides a `channels()` method returning a list of channels).
 * 
 * These three higher level classes, `icarus::trigger::OpticalTriggerGate` and
 * derivatives, _can't be directly saved_ in _art_ ROOT files.
 * There are utilities available in
 * `icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h` that can convert them
 * into a collection of `icarus::trigger::OpticalTriggerGateData_t` objects
 * plus a collection of _art_ associations (for writing), and the other way
 * around (for reading). The module `icarus::trigger::LVDSgates` uses both sides
 * and can be used as an illustration of the functionality.
 * 
 * A module is provided, called `icarus::trigger::DumpTriggerGateData`, which
 * dumps on screen or on text file the information from a collection of
 * discriminated waveforms in a (sort-of) human readable format.
 * An example configuration for this module is provided in `icaruscode`, called
 * `dump_triggergatedata_icarus.fcl`.
 * 
 * The main functions to manipulate the trigger gates are defined in the very
 * base class, `icarus::trigger::TriggerGateData`: these allow e.g. to find
 * events and query the level at a given time.
 * Another important set of features is also present in
 * `icarus::trigger::TriggerGateData` and replicated in the higher levels:
 * the combination of trigger gates by sum (`OR`), multiplication (_AND_),
 * minimum and maximum value.
 * 
 * @note Combining a multi-level gate via `Min()` with a binary gate results
 *       into a binary gate which is logic AND of the two gates.
 * 
 * 
 * On the terminology
 * -------------------
 * 
 * In this documentation, the "discriminated waveforms" are sometimes called
 * "trigger gates" and sometimes "trigger primitives".
 * Although these concepts are, strictly, different, usually the difference
 * does not affect the meaning and they end up being exchanged carelessly.
 * 
 * 
 * Trigger logic algorithm
 * ========================
 * 
 * @anchor TriggerEfficiencyPlots_Algorithm
 * 
 * This section describes the trigger logic algorithm used in
 * `icarus::trigger::TriggerEfficiencyPlots` and its assumptions.
 * 
 * The algorithm treats all the trigger primitives equally, whether they
 * originate from one or from two channels (or 10 or 30), and wherever their
 * channels are in the detector.
 * The trigger primitives are combined in a multi-level gate by adding them,
 * so that the level of the resulting gate describes at any time matches how
 * many trigger primitives are on at that time.
 * 
 * This multi-level gate is set in coincidence with the beam gate by multiplying
 * the multi-level and the beam gates.
 * The beam gate opens at a time configured in `DetectorClocks` service provider
 * (`detinfo::DetectorClocks::BeamGateTime()`) and has a duration configured
 * in this module (`BeamGateDuration`).
 * 
 * At this point, the trigger gate is a multi-level gate suppressed everywhere
 * except than during the beam gate.
 * The algorithm handles multiple trigger definitions, or requirements.
 * Each requirement is simply how many trigger primitives must be open at the
 * same time for the trigger to fire. The values of these requirements are
 * set in the configuration (`MinimumPrimitives`).
 * To determine whether a trigger with a given requirement, i.e. with a required
 * minimum number of trigger primitives open at the same time, has fired, the
 * gate is scanned to find _the first tick_ where the level of the gate reaches
 * or passes this minimum required number. If such tick exists, the trigger is
 * considered to have fired, and at that time.
 * 
 * As a consequence, there are for each event as many different trigger
 * responses as how many different requirements are configured in
 * `MinimumPrimitives`, _times_ how many ADC thresholds are provided in input,
 * configured in `Thresholds`.
 * 
 * While there _is_ a parameter describing the time resolution of the trigger
 * (`TriggerTimeResolution`), this is currently only used for aesthetic purposes
 * to choose the binning of some plots: the resolution is _not_ superimposed
 * to the gates (yet).
 * 
 * This algorithm is currently implemented in
 * `detinfo::trigger::TriggerEfficiencyPlots::plotResponses()`.
 * 
 * 
 * Event categories
 * =================
 * 
 * Each event is assigned to a set of categories, and its information
 * contributes to the plots of all those categories.
 * The categories are defined in the `PlotCategories` list (hard-coded).
 * 
 * 
 * Module usage
 * =============
 * 
 * 
 * Input data products
 * --------------------
 * 
 * This module uses the following data products as input:
 * 
 * * trigger primitives:
 *     * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (labels out of
 *       `TriggerGatesTag` and `Thresholds`): full sets of discriminated
 *       waveforms, each waveform possibly covering multiple optical channels,
 *       and their associations to optical waveforms. One set per threshold;
 *     * associations with `raw::OpDetWaveform` (currently superfluous);
 * * event characterization:
 *     * `std::vector<simb::MCTruth>`: generator information (from
 *       `GeneratorTags`; multiple generators are allowed at the same time);
 *       it is used mostly to categorize the type of event (background,
 *       weak charged current, electron neutrino, etc.);
 *     * `std::vector<simb::MCParticle>`: particles propagating in the detector
 *       (from `PropagatedParticles`); currently not used;
 *     * `std::vector<sim::SimEnergyDeposit>`: energy deposited in the active
 *       liquid argon volume (from `EnergyDepositTags`); it is used to
 *       quantify the energy available to be detected in the event.
 * 
 * 
 * Output plots
 * -------------
 * 
 * The module produces a standard set of plots for each configured ADC threshold
 * and for each event category.
 * The plots are saved via _art_ service `TFileService` in a ROOT file and they
 * are organized in nested ROOT directories under the module label directory
 * which is assigned by _art_:
 * 
 * * `Thr###` (outer level) describes the ADC threshold on the discriminated
 *   waveforms: the threshold is `###` from the baseline;
 * * `\<Category\>` (inner level) describes the category of events included in
 *   the plots (e.g. `All`, `NuCC`; see `PlotCategories`).
 * 
 * Each of the inner ROOT directories contains a full set of plots, whose name
 * is the standard plot name followed by its event category and threshold
 * (e.g. `Eff_NuCC_Thr15` for the trigger efficiency plot on neutrino charged
 * current events with 15 ADC counts as PMT threshold).
 * 
 * Each set of plots, defined in
 * `icarus::trigger::TriggerEfficiencyPlots::initializePlotSet()`, is
 * contributed only by the events in the set category.
 * 
 * Some plots are specific to the trigger, and they are different according to
 * the ADC threshold:
 * 
 * * `Eff`: trigger efficiency defined as number of triggered events over the
 *   total number of events, as function of the minimum number of trigger
 *   primitives (`MinimumPrimitives`) to define a firing trigger; uncertainties
 *   are managed by `TEfficiency`.
 * * `TriggerTick`: distribution of the time of the earliest trigger for the
 *   event, as function of the minimum number of trigger primitives (as in
 *   `Eff`). It may happen that the event is such that there is e.g. a
 *   20-primitive flash, then subsiding, and then another 30-primitive flash.
 *   In such a case, in the trigger requirement "&geq; 15 primitives" such event
 *   will show at the time of the 20-primitive flash, while in the trigger
 *   requirement "&geq; 25 primitives" it will show at the time of the 
 *   30-primitive flash. Each event appears at most once for each trigger
 *   requirement, and it may not appear at all if does not fire a trigger.
 * * `NPrimitives`: the maximum number of primitives "on" at any time.
 * 
 * Some plots are independent of the ADC threshold and they have the same
 * content in each of the plot set for a given event category:
 * 
 * * `EnergyInSpill`: total energy deposited in the detector during the time the
 *   beam gate is open. It is proportional to the amount of scintillation light
 *   in the event.
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
 * * `GeneratorTags` (list of input tags, default: `[ generator ]`): a list of
 *     data products containing the particles generated by event generators;
 * * `DetectorParticleTag` (input tag, default: `largeant`): data product
 *     containing the list of particles going through the argon in the detector;
 * * `EnergyDepositTags`
 *     (list of input tags, default: `[ "largeant:TPCActive" ]`): a list of
 *     data products with energy depositions;
 * * `BeamGateDuration` (time, _mandatory_): the duration of the beam
 *     gate; _the time requires the unit to be explicitly specified_: use
 *     `"1.6 us"` for BNB, `9.5 us` for NuMI (also available as
 *     `BNB_settings.spill_duration` and `NuMI_settings.spill_duration` in
 *     `trigger_icarus.fcl`);
 * * `MinimumPrimitives` (list of integers, _mandatory_): a list of alternative
 *     requirements for the definition of a trigger; each value is the number
 *     of trigger primitives needed to be "on" at the same time for the trigger
 *     to fire;
 * * `TriggerTimeResolution` (time, default: `8 ns`): time resolution for the
 *     trigger primitives;
 * * `LogCategory` (string, default `TriggerEfficiencyPlots`): name of category
 *     used to stream messages from this module into message facility.
 * 
 * An example job configuration is provided as `maketriggerplots_icarus.fcl`.
 * 
 * 
 * Technical description of the module
 * ====================================
 * 
 * @todo Meh. Include:
 *       
 *       * how the plots are organised, and how that reflects the algorithm
 *       * how response and plotting are entangled
 *       * the plot sand box mechanism
 *       * how the selection of event categories happens
 *       * how the result is distributed to different event categories
 *       * how to add new plots
 *       * how to change the algorithm
 * 
 * 
 * Code style: quantities with units
 * ----------------------------------
 * 
 * To avoid issues with missing or wrong conversions, this code often uses
 * LArSoft quantities library. A variable with a `Quantity` type is represented
 * with a specific unit (e.g. microseconds) and can undergo only limited
 * manipulation. The allowed manipulation should guarantee that the unit of
 * the quantity is correctly preserved. For example, it is not possible to
 * add to a `microseconds` interval a pure number (e.g. `9.0`), but rather
 * it is only possible to add time quantities (e.g. another `microseconds`
 * variable, or a `nanoseconds` variable, or a literal value with unit, like
 * `9.0_ns`), and those quantities are properly converted, with the caveat that
 * rounding errors may happen that a more careful explicit handling could avoid.
 * Also, there are two types of variables that can feature the same unit,
 * intervals and points. A point can't be scaled (e.g. you can't "double" the
 * time the beam gate opens, while you can double the beam gate _duration_)
 * and it can't be added to another point (the difference between two points
 * is an interval).
 * 
 * To avoid mistakes in the conversion between different time scales, the
 * LArSoft utility library `detinfo::DetectorTimings` is used, which is a
 * wrapper of `DetectorClocks` service provider that makes the interface to
 * convert times easier and uniform. Here we use it to convert time points
 * from the simulation (which are expressed in nanoseconds and are starting
 * at trigger time) into optical detector readout ticks, and vice versa.
 * The values returned by this library have encoded in them which time scale
 * they belong to, and in which unit they are measured.
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
      Comment("length of time interval when optical triggers are accepted")
      };

    fhicl::Sequence<unsigned int> MinimumPrimitives {
      Name("MinimumPrimitives"),
      Comment("minimum required number of trigger primitives for a trigger")
      };
    
    fhicl::Atom<nanoseconds> TriggerTimeResolution {
      Name("TriggerTimeResolution"),
      Comment("resolution of trigger in time"),
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
  bool fired = true; // the final trigger response (changes with requirement)
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
   * Now fill the plots independent of the trigger response:
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
