/**
 * @file   ExtractEnergyDepositionSummary_module.cc
 * @brief  Module producing a summary deposited energy data product.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateStruct.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.h"
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h"
#include "icaruscode/Utilities/DetectorClocksHelpers.h" // makeDetTimings()...
#include "icaruscode/IcarusObj/SimEnergyDepositSummary.h"
#include "icarusalg/Utilities/ChangeMonitor.h" // ThreadSafeChangeMonitor
#include "icarusalg/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h" // detinfo::DetectorTimings
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <cassert>
/*
#include <map>
*/

//------------------------------------------------------------------------------
namespace icarus::trigger { class ExtractEnergyDepositionSummary; }

/**
 * @brief Produces energy deposition summary data products.
 * 
 * The data product produced by this module is an alternative to storing the
 * whole energy deposition, which can be huge.
 * 
 * This module takes its information from `sim::SimEnergyDeposits` data product
 * in input, simply providing a sum of it within the specified time interval.
 * 
 * 
 * Output data products
 * =====================
 * 
 * * a single `icarus::SimEnergyDepositSummary`: sums of energies in a spill
 *   time, in a pre-spill time, and total, everywhere and in active volume only.
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<sim::SimEnergyDeposits>` (multiple supported) _or_
 *   `std::vector<sim::SimChannel>` to extract the deposited energy from.
 * 
 * 
 * Service requirements
 * ---------------------
 * 
 * The following services are _required_:
 * 
 * * `Geometry` for the determination of the active volume(s)
 * * `DetectorClocksService` for the determination of beam gate and for time 
 *     backtracking if energy deposition information is from `sim::SimChannel`
 * * `DetectorPropertiesService` for time backtracking if energy deposition
 * *   information is from `sim::SimChannel`
 * * _art_ message facility
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description ExtractEnergyDepositionSummary`.
 * 
 * * `EnergyDepositTags`
 *     (list of input tags, default: `[ "largeant:TPCActive" ]`): a list of
 *     data products with energy depositions. In alternative, `SimChannelTag`
 *     can be specified instead.
 * * `SimChannelTag` (input tag): tag of the data product with
 *     `sim::SimChannels` to extract the energy information from. If specified,
 *     it overrides `EnergyDepositTags` and extracts energy information from the
 *     `sim::SimChannel` data product identified by this tag. This option is
 *     for flawed input samples only, and `EnergyDepositTags` is recommended
 *     instead.
 * * `BeamGateDuration` (time, _mandatory_): the duration of the beam
 *     gate; _the time requires the unit to be explicitly specified_: use
 *     `"1.6 us"` for BNB, `9.5 us` for NuMI (also available as
 *     `BNB_settings.spill_duration` and `NuMI_settings.spill_duration` in
 *     `trigger_icarus.fcl`).
 * * `BeamGateStart` (time, default: 0 &micro;s): open the beam gate this long
 *      after the nominal beam gate time.
 * * `PreSpillWindow` (time, default: 10 &micro;s): duration of the pre-spill
 *      window.
 * * `PreSpillWindowGap` (time, default: 0 &micro;s): gap from the end of
 *      pre-spill window to the start of beam gate.
 * * `OutputCategory` (string, default: `"ExtractEnergyDepositionSummary"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * 
 */
class icarus::trigger::ExtractEnergyDepositionSummary: public art::EDProducer {
  
    public:
  
  using microseconds = util::quantities::intervals::microseconds;
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> EnergyDepositTags {
      Name("EnergyDepositTags"),
      Comment("label of energy deposition data product(s) in the detector"),
      std::vector<art::InputTag>{ "largeant:TPCActive" }
      };

    fhicl::OptionalAtom<art::InputTag> SimChannelTag {
      Name("SimChannelTag"),
      Comment("label of source sim::SimChannels data product in the detector"),
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

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("tag of the module output to console via message facility"),
      "ExtractEnergyDepositionSummary"
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit ExtractEnergyDepositionSummary(Parameters const& config);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  using EDepTags_t = details::EventInfoExtractorMaker::EDepTags_t; // alias
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Duration of the gate during with global optical triggers are accepted.
  microseconds fBeamGateDuration;
  
  /// Start of the beam gate with respect to `BeamGate()`.
  microseconds fBeamGateStart;
  
  microseconds fPreSpillWindow; ///< Duration of the pre-spill gate.
  
  microseconds fPreSpillStart; ///< Start of the pre-spill gate.
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  struct TimingPack_t {
    
    /// Detector clocks information.
    detinfo::DetectorClocksData const detClocks;
    
    /// Detector properties information.
    detinfo::DetectorPropertiesData const detProps;
    
    /// Detector timing conversion utility.
    detinfo::DetectorTimings const detTimings;
    
    TimingPack_t(
      detinfo::DetectorClocksData clocksData,
      detinfo::DetectorPropertiesData propsData
      )
      : detClocks{ std::move(clocksData) }
      , detProps{ std::move(propsData) }
      , detTimings{ detinfo::makeDetectorTimings(detClocks) }
      {}
    
  }; // TimingPack_t
  
  geo::GeometryCore const& fGeom; ///< Access to detector geometry information.
  
  /// All information about detector timings and properties.
  std::optional<TimingPack_t> fDetProps;
  
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  /// Functor returning whether a gate has changed.
  icarus::ns::util::ThreadSafeChangeMonitor<icarus::trigger::BeamGateStruct>
    fBeamGateChangeCheck;

  details::EventInfoExtractorMaker const fEventInfoExtractorMaker;
  
  // --- END Algorithms --------------------------------------------------------
  
  
  /// Fills a energy deposition tag object from the specified configuration.
  static EDepTags_t makeEnergyDepSourceTag(
    fhicl::Sequence<art::InputTag> const& energyDepositTags,
    fhicl::OptionalAtom<art::InputTag> const& simChannelTag
    );
  
}; // icarus::trigger::ExtractEnergyDepositionSummary



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::ExtractEnergyDepositionSummary
//------------------------------------------------------------------------------
icarus::trigger::ExtractEnergyDepositionSummary::ExtractEnergyDepositionSummary
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fBeamGateDuration     (config().BeamGateDuration())
  , fBeamGateStart        (config().BeamGateStart())
  , fPreSpillWindow       (config().PreSpillWindow())
  , fPreSpillStart
      (fBeamGateStart - config().PreSpillWindowGap() - fPreSpillWindow)
  , fLogCategory(config().OutputCategory())
  // services
  , fGeom(*(lar::providerFrom<geo::Geometry>()))
  , fDetProps{
      config().SimChannelTag.hasValue()
        ? std::make_optional<TimingPack_t>(
          art::ServiceHandle<detinfo::DetectorClocksService const>()
            ->DataForJob(),
          art::ServiceHandle<detinfo::DetectorPropertiesService const>()
            ->DataForJob()
          )
        : std::nullopt
      }
  // algorithms
  , fEventInfoExtractorMaker(
      {}                                               // truthTags: none
    , makeEnergyDepSourceTag(config().EnergyDepositTags, config().SimChannelTag)
                                                       // edepTags
    , fGeom                                            // geom
    , (fDetProps? &(fDetProps->detProps): nullptr)     // detProps
    , (fDetProps? &(fDetProps->detTimings): nullptr)   // detTimings
    , fLogCategory                                     // logCategory
    , consumesCollector()                              // consumesCollector
    )
{
  produces<icarus::SimEnergyDepositSummary>();
} // icarus::trigger::ExtractEnergyDepositionSummary::ExtractEnergyDepositionSummary()


//------------------------------------------------------------------------------
void icarus::trigger::ExtractEnergyDepositionSummary::produce(art::Event& event)
{
  
  //
  // prepare the information extractor and extract the information
  //
  
  // we need to convert the two relevant gates with the proper parameters
  auto const detTimings = icarus::ns::util::makeDetTimings(event);
  
  auto const beamGate = icarus::trigger::makeBeamGateStruct
    (detTimings, fBeamGateDuration, fBeamGateStart);
  
  if (auto oldGate = fBeamGateChangeCheck(beamGate); oldGate) {
    mf::LogWarning(fLogCategory)
      << "Beam gate has changed from " << oldGate->asOptTickRange()
      << " to " << beamGate.asOptTickRange() << " (optical tick)!";
  }
  
  details::EventInfo_t const info = fEventInfoExtractorMaker(
    beamGate.asSimulationRange(),
    icarus::trigger::makeBeamGateStruct
      (detTimings, fPreSpillStart, fPreSpillStart + fPreSpillWindow)
      .asSimulationRange()
    )(event);
  assert(info.hasDepEnergy());
  
  //
  // move the information into the data product
  //
  icarus::SimEnergyDepositSummary summary;
  summary.Total          = info.DepositedEnergy().value();
  summary.Spill          = info.DepositedEnergyInSpill().value();
  summary.PreSpill       = info.DepositedEnergyInPreSpill().value();
  summary.Active         = info.DepositedEnergyInActiveVolume().value();
  summary.SpillActive    = info.DepositedEnergyInSpillInActiveVolume().value();
  summary.PreSpillActive = info.DepositedEnergyInPreSpillInActiveVolume().value();

  //
  // finish
  //
  event.put
    (std::make_unique<icarus::SimEnergyDepositSummary>(std::move(summary)));
  
} // icarus::trigger::ExtractEnergyDepositionSummary::produce()


//------------------------------------------------------------------------------
auto icarus::trigger::ExtractEnergyDepositionSummary::makeEnergyDepSourceTag(
  fhicl::Sequence<art::InputTag> const& energyDepositTags,
  fhicl::OptionalAtom<art::InputTag> const& simChannelTag
  ) -> EDepTags_t
{
  EDepTags_t edepTags;
  
  if (auto const tag = util::fhicl::getOptionalValue(simChannelTag)) {
    edepTags
     = icarus::trigger::details::EventInfoExtractor::SimChannelsInputTag{ *tag }
     ;
  }
  else {
    edepTags = energyDepositTags();
  }
  
  return edepTags;
} // icarus::trigger::ExtractEnergyDepositionSummary::makeEnergyDepSourceTag()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::ExtractEnergyDepositionSummary)


//------------------------------------------------------------------------------
