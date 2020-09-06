/**
 * @file   MakeTriggerSimulationTree_module.cc
 * @brief  Creates a ROOT tree with trigger information.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfoTree.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventIDTree.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/TreeHolder.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h"
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time, ...
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/Geometry/geo_vectors_utils.h" // MiddlePointAccumulator
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
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
  std::vector<UInt_t> fNOpenings; ///< Number of openings (`0` if never open).
  std::vector<Double_t> fOpeningTime; ///< Time of first opening.
  
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
 * * `EnergyDepositTags`
 *     (list of input tags, default: `[ "largeant:TPCActive" ]`): a list of
 *     data products with energy depositions;
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
 * * `EventTreeName` (string, default: "Treegger"): the name of the ROOT tree
 *     being written on disk.
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

    fhicl::Sequence<art::InputTag> EnergyDepositTags {
      Name("EnergyDepositTags"),
      Comment("label of energy deposition data product(s) in the detector"),
      std::vector<art::InputTag>{ "largeant:TPCActive" }
      };

    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    fhicl::Atom<microseconds> BeamGateDuration {
      Name("BeamGateDuration"),
      Comment("length of time interval when optical triggers are accepted")
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
  
  /// Tag for optical trigger gate data product.
  art::InputTag fTriggerGatesTag;
  
  std::string fLogCategory; ///< Name of output stream for message facility.
  
  // --- END Configuration variables -------------------------------------------


  // --- BEGIN Setup variables -------------------------------------------------
  
  geo::GeometryCore const& fGeom;
  detinfo::DetectorClocksData fDetClocks;
  detinfo::DetectorTimings fDetTimings;
  
  // --- END Setup variables ---------------------------------------------------


  // --- BEGIN Internal variables ----------------------------------------------
  
  /// The beam gate (as a `icarus::trigger::OpticalTriggerGate` object).
  icarus::trigger::OpticalTriggerGate const fBeamGate;
  
  /// Beam gate start and stop time in simulation scale.
  std::pair<simulation_time, simulation_time> const fBeamGateSim;
  
  /// Main ROOT tree: event ID.
  details::EventIDTree fIDTree;
  
  // these trees are "optional" so that in the future they can be disabled;
  // but currently (v08_52_00) no option for disabling them is implemented
  
  /// ROOT tree: event information.
  std::optional<details::EventInfoTree> fEventTree;
  
  /// ROOT tree: trigger gates.
  std::optional<TriggerGateTree> fTriggerGateTree;
  
  /// Helper to fill a `EventInfo_t` from an _art_ event.
  details::EventInfoExtractor const eventInfoExtractor;
  
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
   *        is undefined.
   * 
   */
  TriggerGatesInfo extractTriggerInfo(art::Event const& event) const;
  
  
  /// Returns the centroid (middle point) of all the centers of optical detectors
  /// contributing to the specified `gate`, in world coordinates [cm]
  geo::Point_t gateChannelCentroid
    (icarus::trigger::OpticalTriggerGateData_t const& gate) const;
  
  
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
  
} // TriggerGateTree::checkSizes()


//------------------------------------------------------------------------------
void TriggerGateTree::assignTriggerGatesInfo(TriggerGatesInfo const& info) {
  
  fNChannels = info.TriggerGates.size();
  fOpDetPos.clear();
  fNOpenings.clear();
  fOpeningTime.clear();
  for
    (auto const& [ iChannel, channelInfo ]: util::enumerate(info.TriggerGates))
  {
    
    fOpDetPos.push_back(channelInfo.center);
    fNOpenings.push_back(channelInfo.nOpenings);
    
    // accepting the fallback value when there is no interaction
    // (that is `max()`)
    fOpeningTime.push_back(channelInfo.firstOpenTime.value());
    
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
  , fTriggerGatesTag(config().TriggerGatesTag())
  , fLogCategory(config().LogCategory())
  // setup
  , fGeom(*lar::providerFrom<geo::Geometry>())
  , fDetClocks(art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob())
  , fDetTimings(fDetClocks)
  // other data
  , fBeamGate(
      icarus::trigger::BeamGateMaker{ fDetTimings }
        .make(config().BeamGateDuration())
      )
  , fBeamGateSim(
    fDetTimings.toSimulationTime(fDetTimings.BeamGateTime()),
    fDetTimings.toSimulationTime(fDetTimings.BeamGateTime())
      + config().BeamGateDuration()
    )
  , fIDTree(*(art::ServiceHandle<art::TFileService>()
      ->make<TTree>(config().EventTreeName().c_str(), "Event information")
      ))
  , fEventTree(fIDTree.tree())
  , fTriggerGateTree(fIDTree.tree())
  , eventInfoExtractor(
    config().GeneratorTags(),              // truthTags
    config().EnergyDepositTags(),          // edepTags
    fBeamGateSim,                          // inSpillTimes
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
  
  
} // icarus::trigger::MakeTriggerSimulationTree::MakeTriggerSimulationTree()


//------------------------------------------------------------------------------
void icarus::trigger::MakeTriggerSimulationTree::analyze
  (art::Event const& event)
{
  
  details::EventInfo_t const eventInfo = eventInfoExtractor(event);
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
   *   4. fill the information
   * 
   */
  
  using detinfo::timescales::optical_tick;
  using icarus::trigger::OpticalTriggerGateData_t;
  
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
    auto const beamAndGate = OpticalTriggerGateData_t::Mul(gate, fBeamGate);
    
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
    // 2.4. fill the information
    //
    info.TriggerGates.push_back({
      gateChannelCentroid(gate), // 2.1. get centroid of associated detectors
      nOpenings,
        // users should ignore this if `nOpenings == 0`:
      (nOpenings == 0U)
        ? std::numeric_limits<simulation_time>::max()
        : fDetTimings.toSimulationTime
          (detinfo::timescales::optical_tick{ firstOpenTick })
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
DEFINE_ART_MODULE(icarus::trigger::MakeTriggerSimulationTree)


//------------------------------------------------------------------------------
