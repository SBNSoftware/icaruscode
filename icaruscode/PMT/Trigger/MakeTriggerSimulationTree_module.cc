/**
 * @file   MakeTriggerSimulationTree_module.cc
 * @brief  Creates a ROOT tree with trigger information.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateStruct.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfoTree.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventIDTree.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/TreeHolder.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h"
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/Utilities/DetectorClocksHelpers.h" // makeDetTimings()...
#include "icarusalg/Utilities/ChangeMonitor.h" // ThreadSafeChangeMonitor
#include "icaruscode/PMT/Trigger/Algorithms/TriggerGateBuilder.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h"
#include "icaruscode/PMT/Trigger/Data/SingleChannelOpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Data/TriggerGateData.h"
#include "icaruscode/PMT/Data/WaveformBaseline.h"
#include "icaruscode/Utilities/DataProductPointerMap.h"
#include "icarusalg/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time, ...
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/Geometry/geo_vectors_utils.h" // MiddlePointAccumulator
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataalg/Utilities/quantities_fhicl.h" // for ADCCounts_t parameters
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TTree.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <string>
#include <optional>
#include <utility> // std::pair<>
#include <limits> // std::numeric_limits<>


//------------------------------------------------------------------------------
using detinfo::timescales::simulation_time;

// --- BEGIN -- Trigger information --------------------------------------------

/// All information from the trigger gates.
struct TriggerGatesInfo {
  
  /// Information from a single trigger gate.
  struct TriggerGateInfo {
    
    /// Middle point of the center of all contributing channels.
    geo::Point_t center;
    
    /// Number of times the gate was "open".
    unsigned int nOpenings = 0U;
    
    /// The time of the first opening on the channel, in simulation time [ns]
    simulation_time firstOpenTime = std::numeric_limits<simulation_time>::max();

    double Amplitude; 
    
  }; // struct TriggerGateInfo
  
  /// Collection of all available trigger gates.
  std::vector<TriggerGateInfo> TriggerGates;
  
}; // struct TriggerGatesInfo

std::ostream& operator<<
  (std::ostream& out, TriggerGatesInfo::TriggerGateInfo const& info);
std::ostream& operator<< (std::ostream& out, TriggerGatesInfo const& info);

// --- END -- Trigger information ----------------------------------------------


// --- BEGIN -- ROOT tree helpers ----------------------------------------------
/**
 * @brief Class managing the serialization of trigger gates in a simple ROOT
 *        tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 * Information may be loaded into the tree via `assignTriggerGatesInfo()`
 * method, which sets each branch according to a value of the data structure
 * `TriggerGatesInfo`, as described below.
 *
 * The branch structure is:
 * `NChannels/i:OpDetPos_[]/D:NOpenings[]/i:OpeningTime[]/D`
 *
 *  * `NChannels` (from the size of `TriggerGatesInfo::info`): number of
 *    channels in this event;
 *  * `OpDetPos` (one 3D point per channel, from `TriggerGateInfo::center`):
 *    the location of the centroid of the optical detectors contributing to the
 *    channel, in world coordinates [cm]; it's a vector of GenVector 3D points
 *    (can access coordinates as `OpDetPos.X()` or `OpDetPos.fCoordinates.fX`);
 * * `NOpenings` (one count per channel, from `TriggerGateInfo::nOpenings`):
 *    the number of times the gate opened;
 * * `OpeningTime` (one time per channel, based on
 *    `TriggerGateInfo::firstOpenTime`): first time that channel opened, or
 *    very large number (`std::numeric_limits<double>::max()`)
 *    if no opening happened at all (note that in this case the value in
 *    `TriggerGateInfo::firstOpenTime` is ignored).
 * 
 * All branches come from a variable size vector with size the number of
 * available trigger gates.
 */
struct TriggerGateTree: public icarus::trigger::details::TreeHolder {

  /// Constructor.
  TriggerGateTree(TTree& tree);
  
  /// Copies the information from the `info` record into the ROOT tree buffers.
  void assignTriggerGatesInfo(TriggerGatesInfo const& info);

    private:
  
  UInt_t fNChannels; ///< Number of channels.
  
  std::vector<geo::Point_t> fOpDetPos; ///< Coordinates of the optical detector.
  std::vector<UInt_t> fNOpenings; ///< Number of openings (`0` if never opens).
  std::vector<Double_t> fOpeningTime; ///< Time of first opening.
  std::vector<Double_t> fAmplitude; ///< PMT amplitude.
  
  /// Internal check: all branch buffers have the same size.
  void checkSizes() const;

}; // struct TriggerGateTree


// --- END -- ROOT tree helpers ------------------------------------------------


//------------------------------------------------------------------------------
namespace icarus::trigger { class MakeTriggerSimulationTree; }
/**
 * @brief 
 * 
 * This module produces a ROOT tree with information about the event and about
 * the trigger input.
 * 
 * 
 * Input
 * ======
 * 
 * The following data products are required in input:
 * 
 * * `std::vector<simb::MCTruth>` (tag: configured in `GeneratorTags`):
 *   any number of truth information records;
 * * `std::vector<sim::SimEnergyDeposits>` (tag: configured in
 *   `EnergyDepositTags`): energy depositions to be counted in the event
 *   information;
 * * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (tag: configured
 *   in `TriggerGatesTag`): list of trigger gate outputs to be put into the
 *   tree.
 * 
 * 
 * Output
 * =======
 * 
 * The output tree is written via `TFileService` in the ROOT directory assigned
 * to this module.
 * 
 * The tree contains three "sections":
 * 
 *  * event identification: see `icarus::trigger::details::EventIDTree`
 *    for the details;
 *  * event information: see `icarus::trigger::details::EventInfoTree`
 *    for the details;
 *  * information on the trigger gates (data product configured via
 *    `TriggerGatesTag` parameter): information on the state of each trigger
 *    gate _during the beam gate_ is extracted, including the number of times
 *    the gate is opened and when; see `TriggerGateTree` for details on the
 *    structure of the tree, and `extractTriggerInfo()` for the details of how
 *    each value is computed.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * * `GeneratorTags` (list of input tags, default: `[ generator ]`): a list of
 *     data products containing the particles generated by event generators;
 * * `ParticleTag` (input tag, default: `largeant`): tag of the data product
 *     with particles propagating in the detector, as produced by `LArG4`;
 * * `EnergyDepositTags`
 *     (list of input tags, default: `[ "largeant:TPCActive" ]`): a list of
 *     data products with energy depositions;
 * * `EnergyDepositSummaryTag` (input tag): alternative to `EnergyDepositTags`,
 *     uses energy deposition information from a summary data product instead
 *     of extracting the information;
 * * `TriggerGatesTag` (string, mandatory): name of the module
 *     instance which produced the trigger primitives to be used as input.
 *     The typical trigger primitives used as input may be LVDS discriminated
 *     output (e.g. from `icarus::trigger::LVDSgates` module) or combinations
 *     of them (e.g. from `icarus::trigger::SlidingWindowTrigger` module);
 * * `BeamGateDuration` (time, _mandatory_): the duration of the beam
 *     gate; _the time requires the unit to be explicitly specified_: use
 *     `"1.6 us"` for BNB, `9.5 us` for NuMI (also available as
 *     `BNB_settings.spill_duration` and `NuMI_settings.spill_duration` in
 *     `trigger_icarus.fcl`);
 * * `BeamGateStart` (time, default: 0 &micro;s): open the beam gate this long
 *      after the nominal beam gate time;
 * * `PreSpillWindow` (time, default: 10 &micro;s): duration of the pre-spill
 *      window;
 * * `PreSpillWindowGap` (time, default: 0 &micro;s): gap from the end of
 *      pre-spill window to the start of beam gate;
 * * `EventTreeName` (string, default: "Treegger"): the name of the ROOT tree
 *     being written on disk;
 * * `LogCategory` (string, default `MakeTriggerSimulationTree`): name of
 *     category used to stream messages from this module into message facility.
 * 
 * 
 * Technical description of the module
 * ====================================
 * 
 * 
 */
class icarus::trigger::MakeTriggerSimulationTree: public art::EDAnalyzer {

    public:
  
  using microseconds = util::quantities::intervals::microseconds;
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> GeneratorTags {
      Name("GeneratorTags"),
      Comment("labels of the event generators"),
      std::vector<art::InputTag>{ "generator" }
      };

    fhicl::Atom<art::InputTag> ParticleTag {
      Name("ParticleTag"),
      Comment("label of particles propagating through the detector"),
      "largeant" // default
      };

    fhicl::OptionalSequence<art::InputTag> EnergyDepositTags {
      Name("EnergyDepositTags"),
      Comment("label of energy deposition data product(s) in the detector")
      };

    fhicl::OptionalAtom<art::InputTag> EnergyDepositSummaryTag {
      Name("EnergyDepositSummaryTag"),
      Comment("label of energy deposition summary data product")
      };

    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    
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
    fhicl::Atom<microseconds> BeamGateDuration {
      Name("BeamGateDuration"),
      Comment("length of time interval when optical triggers are accepted")
      };

    fhicl::Atom<microseconds> BeamGateStart {
      Name("BeamGateStart"),
      Comment("open the beam gate this long after the nominal beam gate time"),
      microseconds{ 0.0 }
      };

    fhicl::Atom<microseconds> PreSpillWindow {
      Name("PreSpillWindow"),
      Comment("duration of the pre-spill window"),
      microseconds{ 10.0 }
      };

    fhicl::Atom<microseconds> PreSpillWindowGap {
      Name("PreSpillWindowGap"),
      Comment("gap from the end of pre-spill window to the start of beam gate"),
      microseconds{ 0.0 }
      };

    fhicl::Atom<std::string> EventTreeName {
      Name("EventTreeName"),
      Comment("name of the ROOT tree"),
      "Treegger" // default
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "MakeTriggerSimulationTree" // default
      };

  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  explicit MakeTriggerSimulationTree(Parameters const& config);

  // Plugins should not be copied or assigned.
  MakeTriggerSimulationTree(MakeTriggerSimulationTree const&) = delete;
  MakeTriggerSimulationTree(MakeTriggerSimulationTree&&) = delete;
  MakeTriggerSimulationTree& operator=(MakeTriggerSimulationTree const&) = delete;
  MakeTriggerSimulationTree& operator=(MakeTriggerSimulationTree&&) = delete;

  // --- END Constructors ------------------------------------------------------


  // --- BEGIN Framework hooks -------------------------------------------------

  /// Fills the plots. Also extracts the information to fill them with.
  virtual void analyze(art::Event const& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Duration of the gate during with global optical triggers are accepted.
  microseconds fBeamGateDuration;
  
  /// Start of the beam gate with respect to `BeamGate()`.
  microseconds fBeamGateStart;
  
  microseconds fPreSpillWindow; ///< Duration of the pre-spill gate.
  
  microseconds fPreSpillStart; ///< Start of the pre-spill gate.
  
  /// Tag for optical trigger gate data product.
  art::InputTag fTriggerGatesTag;
 
 
  // --- PMT AMPLITUDE Configuration variables -----------------------------------------
  art::InputTag const fOpDetWaveformTag; ///< Input optical waveform tag.
  
  ///< Input waveform baseline tag.
  std::optional<art::InputTag> const fBaselineTag;
  
  std::optional<float> const fBaseline; ///< A constant baseline level.
  
 // unsigned int const fNOpDetChannels; ///< Number of optical detector channels.
  
  /// Thresholds selected for saving, and their instance name.
  std::map<icarus::trigger::ADCCounts_t, std::string> fSelectedThresholds;
  
  std::string fLogCategory; ///< Name of output stream for message facility.
  
  // --- END Configuration variables -------------------------------------------


  // --- BEGIN Setup variables -------------------------------------------------
  
  geo::GeometryCore const& fGeom;
  
  // --- END Setup variables ---------------------------------------------------


  // --- BEGIN Internal variables ----------------------------------------------
  
  /// Main ROOT tree: event ID.
  details::EventIDTree fIDTree;
  
  // these trees are "optional" so that in the future they can be disabled;
  // but currently (v08_52_00) no option for disabling them is implemented
  
  /// ROOT tree: event information.
  std::optional<details::EventInfoTree> fEventTree;
  
  /// ROOT tree: trigger gates.
  std::optional<TriggerGateTree> fTriggerGateTree;
  
  /// Helper to fill a `EventInfo_t` from an _art_ event.
  details::EventInfoExtractorMaker const fEventInfoExtractorMaker;
  
  /// Functor returning whether a gate has changed.
  icarus::ns::util::ThreadSafeChangeMonitor<icarus::trigger::BeamGateStruct>
    fBeamGateChangeCheck;

  // --- END Internal variables ------------------------------------------------
  
  
  /**
   * @brief Returns a `TriggerGatesInfo` with trigger information for the tree.
   * @param event the _art_ event where to find the trigger gates
   * @return a `TriggerGatesInfo` with trigger information for the tree
   * 
   * Information from the trigger gates configured in the module is extracted
   * and stored into a `TriggerGatesInfo` object, which is returned.
   * The data fields are defined as follows:
   * 
   *  * a `TriggerGateInfo` element is stored in `TriggerGatesInfo::info` for
   *    each and every trigger gate in the input, in the same order as the
   *    input data product;
   *  * for each gate:
   *      * `TriggerGateInfo::center`: for each channel associated to the
   *        trigger gate, the geometric center of the optical detector connected
   *        to that channel is taken; `center` is set to the middle point
   *        (unweighted 3D centroid) of all those points; the center is in world
   *        coordinates and in centimeters;it is assumed that the trigger gate
   *        is associated to at least one channel, otherwise behavior is
   *        undefined;
   *      * `TriggerGateInfo::nOpenings`: the number of times in coincidence
   *        with the beam gate, when the gate changes state from closed (opening
   *        level `0`) to open (opening level `1` or larger);
   *      * `TriggerGateInfo::firstOpenTime`: the time the first of the openings
   *        (as defined in `nOpenings`) happens; the time is in
   *        @ref DetectorClocksSimulationTime "simulation time scale"; if there
   *        was no opening during the beam gate (`nOpenings` zero), the value
   *        is undefined;
   *      * 'TriggerGateInfo::Amplitude': the amplitude of the PMT
   * 
   */
  TriggerGatesInfo extractTriggerInfo(art::Event const& event) const;
  
  
  /// Returns the centroid (middle point) of all the centers of optical detectors
  /// contributing to the specified `gate`, in world coordinates [cm]
  geo::Point_t gateChannelCentroid
    (icarus::trigger::OpticalTriggerGateData_t const& gate) const;
  
  
  /// Creates a `EdepTags_t` from two optional parameters.
  static icarus::trigger::details::EventInfoExtractor::EdepTags_t makeEdepTag(
    fhicl::OptionalSequence<art::InputTag> const& EnergyDepositTags,
    fhicl::OptionalAtom<art::InputTag> const& EnergyDepositSummaryTag
    );
  
}; // icarus::trigger::MakeTriggerSimulationTree



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- TriggerGatesInfo
//------------------------------------------------------------------------------
std::ostream& operator<<
  (std::ostream& out, TriggerGatesInfo::TriggerGateInfo const& info)
{
  out << "center at " << info.center;
  if (info.nOpenings == 0) out << " never opened";
  else {
    out << " opened " << info.nOpenings << " times, first at "
      << info.firstOpenTime;
  }
  return out;
} // operator<< (TriggerGatesInfo::TriggerGateInfo const&)


//------------------------------------------------------------------------------
std::ostream& operator<< (std::ostream& out, TriggerGatesInfo const& info) {
  out << size(info.TriggerGates) << " gates:";
  for (auto const& [ iGate, gate ]: util::enumerate(info.TriggerGates))
    out << "\n [" << iGate << "] " << gate;
  out << "\n";
  return out;
} // operator<< (TriggerGatesInfo const&)



//------------------------------------------------------------------------------
//--- TriggerGateTree
//------------------------------------------------------------------------------
TriggerGateTree::TriggerGateTree(TTree& tree)
  : TreeHolder(tree)
{
  
  this->tree().Branch("NChannels",   &fNChannels);
  this->tree().Branch("OpDetPos",    &fOpDetPos);
  this->tree().Branch("NOpenings",   &fNOpenings);
  this->tree().Branch("OpeningTime", &fOpeningTime);
  this->tree().Branch("Amplitude",   &fAmplitude);
  
} // TriggerGateTree::TriggerGateTree()


//------------------------------------------------------------------------------
void TriggerGateTree::checkSizes() const {

  auto checkSize = [this](auto const& v){ return v.size() == fNChannels; };
  
  if (!checkSize(fOpDetPos)) {
    throw cet::exception("TriggerGateTree") << __func__
      << ": Internal error: unexpected buffer size (" << fOpDetPos.size()
      << ") : fOpDetPos\n";
  }
  if (!checkSize(fNOpenings)) {
    throw cet::exception("TriggerGateTree") << __func__
      << ": Internal error: unexpected buffer size (" << fNOpenings.size()
      << ") : fNOpenings\n";
  }
  if (!checkSize(fOpeningTime)) {
    throw cet::exception("TriggerGateTree") << __func__
      << ": Internal error: unexpected buffer size (" << fOpeningTime.size()
      << ") : fOpeningTime\n";
  }
  if (!checkSize(fAmplitude)) {
    throw cet::exception("TriggerGateTree") << __func__
      << ": Internal error: unexpected buffer size (" << fAmplitude.size()
      << ") : fOpeningTime\n";
  }
  
} // TriggerGateTree::checkSizes()


//------------------------------------------------------------------------------
void TriggerGateTree::assignTriggerGatesInfo(TriggerGatesInfo const& info) {
  
  fNChannels = info.TriggerGates.size();
  fOpDetPos.clear();
  fNOpenings.clear();
  fOpeningTime.clear();
  fAmplitude.clear();
  for
    (auto const& [ iChannel, channelInfo ]: util::enumerate(info.TriggerGates))
  {
    
    fOpDetPos.push_back(channelInfo.center);
    fNOpenings.push_back(channelInfo.nOpenings);
    
    // accepting the fallback value when there is no interaction
    // (that is `max()`)
    fOpeningTime.push_back(channelInfo.firstOpenTime.value());
    fAmplitude.push_back(channelInfo.Amplitude); 
  } // for
  
  checkSizes();
  
} // TriggerGateTree::assignTriggerGatesInfo()


//------------------------------------------------------------------------------
//--- icarus::trigger::MakeTriggerSimulationTree
//------------------------------------------------------------------------------
icarus::trigger::MakeTriggerSimulationTree::MakeTriggerSimulationTree
  (Parameters const& config)
  : art::EDAnalyzer(config)
  // persistent configuration
  , fBeamGateDuration     (config().BeamGateDuration())
  , fBeamGateStart        (config().BeamGateStart())
  , fPreSpillWindow       (config().PreSpillWindow())
  , fPreSpillStart
      (fBeamGateStart - config().PreSpillWindowGap() - fPreSpillWindow)
  , fTriggerGatesTag(config().TriggerGatesTag())
  , fOpDetWaveformTag(config().OpticalWaveforms())
  , fBaselineTag(util::fhicl::getOptionalValue(config().Baselines))
  , fBaseline(util::fhicl::getOptionalValue(config().Baseline))
  //, fNOpDetChannels(getNOpDetChannels(config().NChannels))
  , fLogCategory(config().LogCategory())
  // setup
  , fGeom(*lar::providerFrom<geo::Geometry>())
  // other data
  , fIDTree(*(art::ServiceHandle<art::TFileService>()
      ->make<TTree>(config().EventTreeName().c_str(), "Event information")
      ))
  , fEventTree(fIDTree.tree())
  , fTriggerGateTree(fIDTree.tree())
  , fEventInfoExtractorMaker(
    config().GeneratorTags(),              // truthTags
    config().ParticleTag(),                // particleTag
    makeEdepTag(config().EnergyDepositTags, config().EnergyDepositSummaryTag),
                                           // edepTags
    fGeom,                                 // geom
    fLogCategory,                          // logCategory
    consumesCollector()                    // consumesCollector
    )
{
  
  if (config().GeneratorTags().empty()) {
    throw art::Exception(art::errors::Configuration)
      << "No event generator data product specified (`"
        << config().GeneratorTags.name() << "`)."
      << "\nEvent generation information is mandatory."
      << "\n";
  }
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  
} // icarus::trigger::MakeTriggerSimulationTree::MakeTriggerSimulationTree()


//------------------------------------------------------------------------------
void icarus::trigger::MakeTriggerSimulationTree::analyze
  (art::Event const& event)
{
  
  // we need to convert the two relevant gates with the proper parameters
  auto const detTimings = icarus::ns::util::makeDetTimings(event);
  
  auto const beamGate = icarus::trigger::makeBeamGateStruct
    (detTimings, fBeamGateDuration, fBeamGateStart);
  
  if (auto oldGate = fBeamGateChangeCheck(beamGate); oldGate) {
    mf::LogWarning(fLogCategory)
      << "Beam gate has changed from " << oldGate->asOptTickRange()
      << " to " << beamGate.asOptTickRange() << " (optical tick)!";
  }
  
  details::EventInfo_t const eventInfo = fEventInfoExtractorMaker(
    beamGate.asSimulationRange(),
    icarus::trigger::makeBeamGateStruct
      (detTimings, fPreSpillStart, fPreSpillStart + fPreSpillWindow)
      .asSimulationRange()
    )(event);
  
  TriggerGatesInfo const triggerInfo = extractTriggerInfo(event);
  
  mf::LogDebug(fLogCategory) << event.id() << " trigger info: " << triggerInfo;
  
  fIDTree.assignID(event.id());
  if (fEventTree) fEventTree->assignEvent(eventInfo);
  if (fTriggerGateTree) fTriggerGateTree->assignTriggerGatesInfo(triggerInfo);
  
  fIDTree.tree().Fill();
  
} // icarus::trigger::MakeTriggerSimulationTree::analyze()


//------------------------------------------------------------------------------
TriggerGatesInfo icarus::trigger::MakeTriggerSimulationTree::extractTriggerInfo
  (art::Event const& event) const
{
  /*
   * 1. get the data product from the event
   * 2. fill one "channel" of `TriggerGateInfo` per entry in the data product
   *   1. get centroid of associated detectors
   *   2. find and convert to `simulation_time` the first opening
   *   3. find the number of openings
   *   4. find the amplitude
   *   4. fill the information
   * 
   */
  
  using detinfo::timescales::optical_tick;
  using icarus::trigger::OpticalTriggerGateData_t;
  
  //
  // 0. construct the beam gate for this event
  //
  auto const detTimings = icarus::ns::util::makeDetTimings(event);
  auto const beamGate = icarus::trigger::BeamGateMaker{ detTimings }
      .make(fBeamGateDuration, fBeamGateStart)
    ;
  
  //
  // 1. get the data product from the event
  //
  auto const& gates
    = event.getByLabel<std::vector<OpticalTriggerGateData_t>>(fTriggerGatesTag);
  
  //
  // 2. fill one "channel" of `TriggerGateInfo` per entry in the data product
  //
  TriggerGatesInfo info;
  
  unsigned int nOpenChannels = 0U; // count for debugging message
  for (OpticalTriggerGateData_t const& gate: gates) {
    
    // the gate, in coincidence with the beam gate:
    auto const beamAndGate = OpticalTriggerGateData_t::Mul(gate, beamGate);
    
    //
    // 2.2. find and convert to `simulation_time` the first opening
    //
    OpticalTriggerGateData_t::ClockTick_t const firstOpenTick
      = beamAndGate.findOpen();
    
    //
    // 2.3. find the number of openings
    //
    unsigned int nOpenings = 0U;
    OpticalTriggerGateData_t::ClockTick_t openTick = firstOpenTick;
    while (openTick != OpticalTriggerGateData_t::MaxTick) {
      ++nOpenings;
      // move out of this open region, then find the next one
      openTick = beamAndGate.findClose(1U, openTick);
      openTick = beamAndGate.findOpen(1U, openTick);
    } // while
    if (nOpenings > 0) ++nOpenChannels;
    //
    //2.35 if found openings>0 get the waveform amplitude
    //
    //

  auto const& waveforms
  = event.getByLabel<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  //
  // retrieve the baseline information
  //
  std::vector<icarus::WaveformBaseline> fixedBaselines;
  std::vector<icarus::WaveformBaseline> const* baselines = nullptr;
  if (fBaselineTag) {
    baselines =
      &(event.getByLabel<std::vector<icarus::WaveformBaseline>>(*fBaselineTag));
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
  std::vector<double> WaveformAndBaseline;

  for (auto const& [ waveform, baseline ]: util::zip(waveforms, *baselines)){
    //waveformInfo.emplace_back(&waveform, &baseline);
    for(auto const& w:waveform){
	WaveformAndBaseline.push_back(w - baseline.baseline());
}
 }

std::vector<double> pmtAmplitudes;

   if (nOpenings > 0){
   auto minAmplitudes = std::min_element(WaveformAndBaseline.begin(), WaveformAndBaseline.end());
   double a = static_cast<double>(minAmplitudes[0]);
   pmtAmplitudes.push_back(a);
 
//for(auto const& wf:waveformInfo){

//now wf is the type of const icarus::trigger::WaveformWithBaseline

//std::cout<<wf<<std::endl;

//std::vector<auto> waveformAndBaseline; 

/*for(auto const& v : wf.waveform())
WaveformAndBaseline.push_back(v - wf.baseline());
std::min(waveformAndBaseline.cbegin(), waveformAndBaseline.waveform.cend()); */
} else pmtAmplitudes.push_back(0.0);
    
   auto minAmplitude = std::min_element(pmtAmplitudes.begin(), pmtAmplitudes.end());
   double amplitude = static_cast<double>(minAmplitude[0]); 
    //
    // 2.4. fill the information
    //
    info.TriggerGates.push_back({
      gateChannelCentroid(gate), // 2.1. get centroid of associated detectors
      nOpenings,
        // users should ignore this if `nOpenings == 0`:
      (nOpenings == 0U)
        ? std::numeric_limits<simulation_time>::max() //std::numeric_limits<T>::max Returns the maximum finite value representable by the 
						      //numeric type T. Meaningful for all bounded types.
        : detTimings.toSimulationTime
          (detinfo::timescales::optical_tick{ firstOpenTick }),
       amplitude
      });
    
  } // for all gates
  
  mf::LogTrace(fLogCategory)
    << "Information from '" << fTriggerGatesTag.encode() << "' ("
    << gates.size() << " trigger gates, " << nOpenChannels
    << ") written to tree."
    ;
  
  return info;
} // icarus::trigger::MakeTriggerSimulationTree::extractTriggerInfo()


//------------------------------------------------------------------------------
geo::Point_t icarus::trigger::MakeTriggerSimulationTree::gateChannelCentroid
  (icarus::trigger::OpticalTriggerGateData_t const& gate) const
{
  
  geo::vect::MiddlePointAccumulator centroid;
  
  for (auto const channel: gate.channels()) {
    
    try {
      centroid.add(fGeom.OpDetGeoFromOpChannel(channel).GetCenter());
    }
    catch (cet::exception const& e) {
      throw cet::exception("MakeTriggerSimulationTree", "", e)
        << "Error while accessing position of optical detector with channel "
        << channel << "\n";
    }
    
  } // for
  
  return centroid.middlePoint();
  
} // icarus::trigger::MakeTriggerSimulationTree::gateChannelCentroid()


//------------------------------------------------------------------------------
icarus::trigger::details::EventInfoExtractor::EdepTags_t
icarus::trigger::MakeTriggerSimulationTree::makeEdepTag(
  fhicl::OptionalSequence<art::InputTag> const& EnergyDepositTags,
  fhicl::OptionalAtom<art::InputTag> const& EnergyDepositSummaryTag
) {
  
  if (EnergyDepositSummaryTag.hasValue()) {
    if (EnergyDepositTags.hasValue()) {
      throw art::Exception(art::errors::Configuration)
        << "MakeTriggerSimulationTree: "
        << "both EnergyDepositTags and EnergyDepositSummaryTag "
        << "have been specified, but they are exclusive: @erase one.\n";
    }
    
    art::InputTag tag;
    EnergyDepositSummaryTag(tag);
    return { 
      icarus::trigger::details::EventInfoExtractor::SimEnergyDepositSummaryInputTag
        { tag }
      };
  }
  else {
    
    std::vector<art::InputTag> tags;
    if (!EnergyDepositTags(tags)) tags.push_back("largeant:TPCActive");
    return { std::move(tags) };
    
  }
  
} // icarus::trigger::MakeTriggerSimulationTree::makeEdepTag()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::MakeTriggerSimulationTree)


//------------------------------------------------------------------------------
