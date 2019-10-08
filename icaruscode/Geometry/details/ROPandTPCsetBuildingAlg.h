/**
 * @file   icaruscode/Geometry/details/ROPandTPCsetBuildingAlg.h
 * @brief  Algorithm discovering TPC sets and readout planes for ICARUS.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 10, 2019
 * @see    `icaruscode/Geometry/details/ROPandTPCsetBuildingAlg.cxx`
 */

#ifndef ICARUSCODE_GEOMETRY_DETAILS_ROPANDTPCSETBUILDINGALG_H
#define ICARUSCODE_GEOMETRY_DETAILS_ROPANDTPCSETBUILDINGALG_H


// ICARUS libraries
#include "icaruscode/Geometry/details/GeometryObjectCollections.h"

// LArSoft libraries
#include "larcorealg/Geometry/GeometryData.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/ReadoutDataContainers.h"
#include "larcorealg/Geometry/GeometryDataContainers.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h" // readout::TPCsetID, ...
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// C/C++ standard library
#include <string>
#include <vector>
#include <tuple>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::details {
  
  //
  // forward declarations
  //
  class ROPandTPCsetBuildingAlg;
  
  //
  // data type definitions
  //
  
  /// Type of collection of TPCs (pointers to `geo::TPCGeo`).
  using TPCColl_t = std::vector<geo::TPCGeo const*>;
  
  /// Type of collection of planes (pointers to `geo::PlaneGeo`).
  using PlaneColl_t = std::vector<geo::PlaneGeo const*>;
  
  
} // namespace icarus::details


// -----------------------------------------------------------------------------
/**
 * @brief Extracts TPC sets and readout planes from a list of cryostats.
 * 
 * This algorithm parses the wire plane and TPC content of the cryostats and
 * groups planes which share the same drift coordinate.
 * Afterwards, it arranges them into LArSoft readout planes, defining at the
 * same time the TPC sets that contain them.
 * 
 * 
 * Algorithm workflow
 * ===================
 * 
 * An example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarus::details::ROPandTPCsetBuildingAlg builder("BuildingAlg");
 * auto results = builder.run(Cryostats);
 * 
 * // planes in each ROP
 * readout::ROPDataContainer<icarus::details::PlaneColl_t> ROPplanes
 *   = std::move(results.ROPplanes());
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will run the algorithm on `Cryostats` and extract one of the results, the
 * map of wire planes in each readout plane.
 * It will use the default configuration of the algorithm, except for the
 * message category being customized into `"BuildingAlg"`.
 * 
 * 
 * Configuration
 * ==============
 * 
 * See the
 * @ref ICARUS_ROPandTPCsetBuildingAlg_Config "constructor documentation".
 * 
 * 
 * 
 */
class icarus::details::ROPandTPCsetBuildingAlg {
  
    public:
  
  // import definitions
  using TPCColl_t = icarus::details::TPCColl_t;
  using PlaneColl_t = icarus::details::PlaneColl_t;
  
  
  /// Full group of results of the algorithm.
  using ResultsBase_t = std::tuple<
    std::vector<unsigned int>,                  // TPCsetCount
    readout::TPCsetDataContainer<TPCColl_t>,    // TPCsetTPCs
    readout::TPCsetDataContainer<unsigned int>, // ROPcount
    readout::ROPDataContainer<PlaneColl_t>,     // ROPplanes
    geo::TPCDataContainer<readout::TPCsetID>,   // TPCtoTPCset
    geo::PlaneDataContainer<readout::ROPID>     // PlaneToROP
    >;
  struct Results_t: public ResultsBase_t {
    
    Results_t(ResultsBase_t&& data): ResultsBase_t(std::move(data)) {}
    
    /// Output: number of TPC sets in each cryostat.
    std::vector<unsigned int> TPCsetCount() &&
      { return std::move(std::get<0U>(*this)); }
    
    /// Output: TPC's in each TPC set.
    readout::TPCsetDataContainer<TPCColl_t> TPCsetTPCs() &&
      { return std::move(std::get<1U>(*this)); }
    
    /// Output: number of readout planes in each TPC set.
    readout::TPCsetDataContainer<unsigned int> ROPcount() &&
      { return std::move(std::get<2U>(*this)); }
    
    /// Output: planes in each of the readout plane.
    readout::ROPDataContainer<PlaneColl_t> ROPplanes() &&
      { return std::move(std::get<3U>(*this)); }
    
    /// Output: number of readout planes in each TPC set.
    geo::TPCDataContainer<readout::TPCsetID> TPCtoTPCset() &&
      { return std::move(std::get<4U>(*this)); }
    
    /// Output: planes in each of the readout plane.
    geo::PlaneDataContainer<readout::ROPID> PlaneToROP() &&
      { return std::move(std::get<5U>(*this)); }
    
  }; // struct Results_t
  
  
  /// Construction: use default configuration.
  ROPandTPCsetBuildingAlg() = default;
  
  /**
   * @brief Construction: specify the algorithm configuration.
   * @param logCategory message facility category to write messages into
   * @param usePlaneIDforROPID whether to use `ROPnumberFromPlanes()` algorithm
   * 
   * 
   * Configuration
   * ==============
   * 
   * @anchor ICARUS_ROPandTPCsetBuildingAlg_Config
   * 
   * * *logCategory* (default: `"ROPandTPCsetBuildingAlg"`): name of the stream
   *     destination used for output messages via message facility
   * 
   */
  ROPandTPCsetBuildingAlg(std::string const& logCategory)
    : fLogCategory(logCategory)
    {}
  
  
  /**
   * @brief Runs the algorithm as configured from start to end.
   * @return the results of the algorithm
   *
   * Results are in the form of a `tuple` that can be used to construct
   * variables with structured binding declarations, or saved to extract
   * the results one by one (see `Results_t` for the interface).
   * 
   */
  Results_t run(geo::GeometryData_t::CryostatList_t const& Cryostats);
  
  
    private:
  
  // --- BEGIN -- Configuration ------------------------------------------------
  
  /// Category to write messages into.
  std::string fLogCategory = "ROPandTPCsetBuildingAlg";
  
  // --- END -- Configuration --------------------------------------------------
  
  
  geo::GeometryData_t::CryostatList_t const* fCryostats = nullptr;
  
  // --- BEGIN -- Output -------------------------------------------------------
  /// @name Output
  /// @{
  
  /// Output: number of TPC sets in each cryostat.
  std::vector<unsigned int> fTPCsetCount;
  
  /// Output: TPC's in each TPC set.
  readout::TPCsetDataContainer<TPCColl_t> fTPCsetTPCs;
  
  /// Output: number of readout planes in each TPC set.
  readout::TPCsetDataContainer<unsigned int> fROPcount;
  
  /// Output: planes in each of the readout plane.
  readout::ROPDataContainer<PlaneColl_t> fROPplanes;
  
  /// Output: the TPC set each TPC belongs to.
  geo::TPCDataContainer<readout::TPCsetID> fTPCtoTPCset;
  
  /// Output: the ROP each wire plane belongs to.
  geo::PlaneDataContainer<readout::ROPID> fPlaneToROP;
  
  /// @}
  // --- END -- Output ---------------------------------------------------------
  
  /// Highest number of TPC sets in any cryostat.
  unsigned int fMaxTPCsets = 0U;
  
  /// Highest number of ROPs in any TPC set.
  unsigned int fMaxROPs = 0U;
  
  
  /// Destroys all the result data members.
  void clear();
  
  /**
   * @brief Extracts composition of all readout planes.
   * @tparam Pred type of the predicate `pred` with one `geo::PlaneGeo` argument
   * @param standalonePlane predicate to determine if a plane is stand-alone
   * @return readout plane composition
   * 
   * This function returns a list: for each cryostat (first index),
   * a list of TPC sets (second index: TPC set within the cryostat)
   * each listing its TPC (third index: runs through all TPCs).
   * 
   * If `standalonePlane` is specified, it is a functor that returns whether its
   * argument, a `geo::PlaneGeo` object, should not be grouped at all (i.e.,
   * `standalonePlane(plane)` returns `false` if the `plane` can be grouped to
   * others).
   * Each plane for which `standalonePlane` evaluates `true` has their own
   * "group" of one plane.
   */
  template <typename Pred>
  std::vector<std::vector<PlaneColl_t>> groupPlanesAndTPCs
    (Pred standalonePlane);
  
  /**
   * @brief Extracts composition of all readout planes.
   * @return readout plane composition
   * 
   * This function returns a list: for each cryostat (first index),
   * a list of TPC sets (second index: TPC set within the cryostat)
   * each listing its TPC (third index: runs through all TPCs).
   */
  std::vector<std::vector<PlaneColl_t>> groupPlanesAndTPCs();
  
  
  /**
   * @brief Extracts all the TPC sets covered by any of the plane groups.
   * @param groupedPlanes wire planes, grouped by cryostat.
   * @return TPC sets as sets of TPC IDs
   * 
   * The return value is a container that for each cryostat (first index),
   * hosts a list of TPC sets (second index: TPC set within the cryostat),
   * each listing its TPC (third index: runs through all TPCs).
   * The order of the TPC sets within a cryostat is undefined.
   * This function also fills `fMaxTPCsets`.
   * 
   * The input parameter `groupedPlanes` follows the same structure as the
   * argument `AllTPCsOnROPs` of `groupPlanesIntoROPs()`.
   */
  std::vector<std::vector<std::vector<geo::TPCID>>> extractTPCsetsFromROPs
    (std::vector<std::vector<PlaneColl_t>> const& planes);
  
  
  /**
   * @brief Extracts the final set of TPC sets
   * @param AllTPCsOnROPs groups of TPC IDs planes in a readout planes belong to
   * 
   * The function collects sets of TPCs from all the ROP's.
   * It fills `fTPCsetCount` and `fTPCsetTPCs`.
   * It requires `fMaxTPCsets` to have been already set.
   * The input argument `AllTPCsOnROPs` is a vector with index the number of
   * cryostat. For each cryostat, a list of readout planes is presented
   * (unsorted and unlabelled), in each of which the list of TPCs the planes in
   * the readout planes belong to is stored.
   */
  void fillTPCsInSet
    (std::vector<std::vector<std::vector<geo::TPCID>>> const& AllTPCsOnROPs);
  
  
  /**
   * @brief Assigns each of the readout planes to a TPC set.
   * @param AllTPCsOnROPs groups of TPC IDs planes in a readout planes belong to
   * @param AllPlanesInROPs planes grouped by ROP
   * @return a map of readout planes in each TPC set, each with plane content
   * 
   * This function fills `MaxROPs`. In addition, it returns the list by TPC set
   * of the groups of sensitive planes into each readout plane in the set.
   * It requires `fTPCsetTPCs` to have been computed already.
   * 
   * The input argument `AllTPCsOnROPs` is a vector with index the number of
   * cryostat. For each cryostat, a list of readout planes is presented
   * (unsorted and unlabelled), in each of which the list of TPCs the planes in
   * the readout planes belong to is stored.
   * 
   */
  readout::TPCsetDataContainer<std::vector<PlaneColl_t>> groupPlanesIntoROPs(
    std::vector<std::vector<std::vector<geo::TPCID>>> const& AllTPCsOnROPs, std::vector<std::vector<PlaneColl_t>>&& AllPlanesInROPs
    );
  
  /**
   * @brief Builds final readout plane information.
   * @param PlanesInProtoROPs plane-to-ROP assignments in each TPC set
   * 
   * This function fills `fROPcount` and `fROPplanes`.
   * 
   * It requires `fMaxROPs`, `fTPCsetTPCs` and `fTPCsetCount` to have been
   * computed already.
   */
  void fillPlanesInROP(
    readout::TPCsetDataContainer<std::vector<PlaneColl_t>> const& PlanesInProtoROPs
    );
  
  /**
   * @brief Creates the map from each TPC to its TPC set.
   * 
   * This function fills `fTPCtoTPCset`.
   * 
   * It requires `fCryostats` and `fTPCsetTPCs` to have been computed already.
   */
  void fillTPCtoTPCsetMap();
  
  
  /**
   * @brief Creates the map from each wire plane to its readout plane.
   * 
   * This function fills `fPlaneToROP`.
   * 
   * It requires `fCryostats` and `fROPplanes` to have been computed
   * already.
   */
  void fillPlaneToROPmap();
  
  
  /**
   * @brief Returns the planes grouped by their drift coordinate.
   * @tparam PlaneColl a `geo::PlaneGeo` collection (forward-iterable)
   * @param planes collection of the wire planes to be grouped
   * @param tolerance maximum distance in drift for planes in a single group
   * @return a collection of groups of (pointers to) planes
   * 
   * The `planes` in the input collection are split into groups which share the
   * same drift coordinate. The `tolerance` parameter defines how different that
   * coordinate can be for planes to be grouped together.
   * The result is expressed as a collection of collections of planes: each
   * element of the outer collection is a group of constant pointers to
   * `geo::PlaneGeo` objects representing all the planes with the same drift
   * coordinate.
   * The sorting algorithm is guaranteed to be stable, i.e. if any two planes
   * _A_ and _B_ have the same drift coordinate and in the `planes` collection
   * _A_ is before _B_, in the resulting group _A_ will also be before _B_.
   * 
   * In the current ICARUS geometry description the drift direction is on _x_
   * axis, and _x_ is the drift coordinate.
   */
  template <typename PlaneColl>
  static std::vector<std::vector<geo::PlaneGeo const*>> groupPlanesByDriftCoord
    (PlaneColl const& planes, double tolerance = 0.1);
  
  
  /// Returns a collection with a TPC ID for each plane in the list of `planes`.
  static std::vector<geo::TPCID> extractTPCIDs
    (std::vector<geo::PlaneGeo const*> const& planes);
  
  
  /// Returns whether all the TPCs covered by a ROP are in a given TPC set.
  static bool isROPinTPCset(
    std::vector<geo::TPCID> const& ROPTPCIDs,
    std::vector<geo::TPCID> const& TPCsetTPCIDs
    );
  
  
  /**
   * @brief Returns ROP number matching the plane number shared by all `planes`.
   * @return a ROP number, or `readout::ROPID::getInvalidID()` on failure
   * 
   * The algorithms verifies that all elements in `planes` have an `ID()` with
   * the same plane number (`geo::PlaneID::Plane`), and returns it.
   * If the plane list is empty or if not all planes in the list have the same
   * plane number, `readout::ROPID::getInvalidID()` is returned.
   */
  static readout::ROPID::ROPID_t ROPnumberFromPlanes(PlaneColl_t const& planes);
  
  
}; // class icarus::details::ROPandTPCsetBuildingAlg


// ----------------------------------------------------------------------------

#endif // ICARUSCODE_GEOMETRY_DETAILS_ROPANDTPCSETBUILDINGALG_H
