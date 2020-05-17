/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.h
 * @brief  Functions dealing with `icarus::trigger::details::EventInfo_t`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 15, 2020
 * @see    icaruscode/PMT/Trigger/Algorithms/details/EventInfoUtils.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFOUTILS_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFOUTILS_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/details/EventInfo_t.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time
#include "lardataalg/Utilities/quantities/energy.h" // gigaelectronvolt, ...
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t
#include "nusimdata/SimulationBase/MCTruth.h"

// framework information
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h" // mf namespace

// C/C++ standard libraries
#include <vector>
#include <string>


//------------------------------------------------------------------------------
//---forward declarations
//---
namespace geo {
  class TPCGeo;
  class GeometryCore;
} // namespace geo

//------------------------------------------------------------------------------
namespace icarus::trigger::details { class EventInfoExtractor; }

/**
 * @brief Extracts from `event` the relevant information on the physics event.
 *
 * This returns a `EventInfo_t` object filled as completely as possible.
 *
 * This class can read any data product in the event, and is completely
 * configured at construction time via constructor arguments.
 * If information is needed from a data product that is not read yet,
 * the following actions are needed:
 *
 * 1. decide the input tag of the data product(s) to be read; this is usually
 *    delegated to the user via a configuration parameter;
 * 2. the required input tags need to be stored into data members;
 * 3. _art_ needs to be informed of the intention of reading the data products
 *    via `consumes()` or equivalent calls; this class can take care of that,
 *    if a `art::ConsumesCollector` object is provided;
 * 4. the data product can be read with the usual means.
 *
 */
class icarus::trigger::details::EventInfoExtractor {
  
    public:
  
  using simulation_time = detinfo::timescales::simulation_time;
  using TimeSpan_t = std::pair<simulation_time, simulation_time>;
  
  /**
    * @name Constructors
    * 
    * Constructors are available to take advantage of specific _art_ features,
    * and for the generic _art_/_gallery_ compatible interfaces.
    */
  /// @{
  
  /**
    * @brief Constructor: configures the object.
    * @param truthTags list of truth information data products to be read
    * @param edepTags list of energy deposition data products to be read
    * @param inSpillTimes start and end of spill, in simulation time
    * @param geom LArSoft geometry service provider
    * @param logCategory name of message facility stream to sent messages to
    * @see `EventInfoExtractor(std::vector<art::InputTag> const&, ConsumesColl&)`
    * 
    * In _art_ framework context, prefer the 
    * `EventInfoExtractor(std::vector<art::InputTag> const&, ConsumesColl&)`
    * constructor.
    */
  EventInfoExtractor(
    std::vector<art::InputTag> truthTags,
    std::vector<art::InputTag> edepTags,
    TimeSpan_t inSpillTimes,
    geo::GeometryCore const& geom,
    std::string logCategory = "EventInfoExtractor"
    );
  
  /**
    * @brief Constructor: configures, and declares consuming data product.
    * @tparam ConsumesColl type with `art::ConsumesCollector` interface
    * @param truthTags list of truth information data products to be read
    * @param edepTags list of energy deposition data products to be read
    * @param geom LArSoft geometry service provider
    * @param logCategory name of message facility stream to sent messages to
    * @param consumesCollector object to declare the consumed products to
    * 
    * Typical usage is within a _art_ module constructor:
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
    * icarus::trigger::details::EventInfoExtractor const eventInfoFrom {
    *   { art:InputTag{ "generator" } },
    *   consumesCollector()
    *   };
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * 
    */
  template <typename ConsumesColl>
  EventInfoExtractor(
    std::vector<art::InputTag> const& truthTags,
    std::vector<art::InputTag> const& edepTags,
    TimeSpan_t inSpillTimes,
    geo::GeometryCore const& geom,
    std::string const& logCategory,
    ConsumesColl& consumesCollector
    );
  
  /// @}
  
  
  // @{
  /**
    * @brief Returns an `EventInfo_t` object with information from `event`.
    * @tparam Event type of event to read data from (_art_ or _gallery_ event)
    * @param event event record to read the information from
    * @return a `EventInfo_t` with information extracted from `event`
    * 
    */
  template <typename Event>
  EventInfo_t extractInfo(Event const& event) const;
  
  template <typename Event>
  EventInfo_t operator() (Event const& event) const
    { return extractInfo(event); }
    
  // @}
  
    private:
  
  // --- BEGIN -- Configuration variables --------------------------------------
  
  /// List of truth data product tags (`std::vector<simb::MCTruth>`).
  std::vector<art::InputTag> const fGeneratorTags;
  
  /// List of energy deposition product tags
  /// (`std::vector<sim::SimEnergyDeposit>`).
  std::vector<art::InputTag> const fEnergyDepositTags;
  
  std::string const fLogCategory; ///< Stream name for message facility.

  // --- END -- Configuration variables ----------------------------------------
  
  
  // --- BEGIN -- Set up  ------------------------------------------------------
  
  geo::GeometryCore const& fGeom; ///< Geometry service provider.
  
  TimeSpan_t const fInSpillTimes; ///< Start and stop time for "in spill" label.
  
  // --- END -- Set up  --------------------------------------------------------
  
  
  /// Fills `info` record with generation information from `truth`.
  void fillGeneratorInfo(EventInfo_t& info, simb::MCTruth const& truth) const;
  
  /// Fills `info` record with generated neutrino information from `truth`.
  void fillGeneratorNeutrinoInfo
    (EventInfo_t& info, simb::MCTruth const& truth) const;

  /// Adds the energy depositions from `energyDeposits` into `info` record.
  void addEnergyDepositionInfo(
    EventInfo_t& info, std::vector<sim::SimEnergyDeposit> const& energyDeposits
    ) const;

  /// Returns in which TPC volume `point` falls in (`nullptr` if none).
  geo::TPCGeo const* pointInTPC(geo::Point_t const& point) const;

  /// Returns in which active TPC volume `point` falls in (`nullptr` if none).
  geo::TPCGeo const* pointInActiveTPC(geo::Point_t const& point) const;

}; // class icarus::trigger::details::EventInfoExtractor


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename ConsumesColl>
icarus::trigger::details::EventInfoExtractor::EventInfoExtractor(
  std::vector<art::InputTag> const& truthTags,
  std::vector<art::InputTag> const& edepTags,
  TimeSpan_t inSpillTimes,
  geo::GeometryCore const& geom,
  std::string const& logCategory,
  ConsumesColl& consumesCollector
  )
  : EventInfoExtractor(truthTags, edepTags, inSpillTimes, geom, logCategory)
{
  
  for (art::InputTag const& inputTag: fGeneratorTags)
    consumesCollector.template consumes<std::vector<simb::MCTruth>>(inputTag);
  
  for (art::InputTag const& inputTag: fEnergyDepositTags) {
    consumesCollector.template consumes<std::vector<sim::SimEnergyDeposit>>
      (inputTag);
  }
  
} // icarus::trigger::details::EventInfoExtractor::EventInfoExtractor(coll)


// -----------------------------------------------------------------------------
template <typename Event>
auto icarus::trigger::details::EventInfoExtractor::extractInfo
  (Event const& event) const -> EventInfo_t
{
  
  EventInfo_t info;
  
  //
  // generator information
  //
  for (art::InputTag const& inputTag: fGeneratorTags) {
  
    auto const& truthRecords
      = event.template getByLabel<std::vector<simb::MCTruth>>(inputTag);
    
    for (simb::MCTruth const& truth: truthRecords) {
      
      fillGeneratorInfo(info, truth);
      
    } // for truth records
    
  } // for generators
  
  //
  // propagation in the detector
  //
  for (art::InputTag const& edepTag: fEnergyDepositTags) {
    
    auto const& energyDeposits
      = event.template getByLabel<std::vector<sim::SimEnergyDeposit>>(edepTag);
    mf::LogTrace(fLogCategory)
      << "Event " << event.id() << " has " << energyDeposits.size()
      << " energy deposits recorded in '" << edepTag.encode() << "'";
    
    addEnergyDepositionInfo(info, energyDeposits);
    
  } // for all energy deposit tags
  
  mf::LogTrace(fLogCategory) << "Event " << event.id() << ": " << info;
  
  return info;
} // icarus::trigger::details::EventInfoExtractor::extractInfo()


//------------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHM_DETAILS_EVENTINFOUTILS_H
