/**
 * @file   icaruscode/Geometry/details/ROPandTPCsetBuildingAlg.cxx
 * @brief  Algorithm discovering TPC sets and readout planes for ICARUS.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 10, 2019
 * @see    `icaruscode/Geometry/details/ROPandTPCsetBuildingAlg.h`
 */

// library header
#include "icaruscode/Geometry/details/ROPandTPCsetBuildingAlg.h"

// LArSoft libraries
#include "larcorealg/Geometry/details/extractMaxGeometryElements.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::size()
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h" // readout::TPCsetID, ...
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h" // cet::exception

// C/C++ standard library
#include <string>
#include <vector>
#include <tuple>
#include <algorithm> // std::max(), std::transform()
#include <utility> // std::move(), std::pair, std::declval()
#include <type_traits> // std::decay_t
#include <cmath> // std::abs()
#include <cassert>


namespace {
  
  // --------------------------------------------------------------------------
  
  // Creates a STL vector with the result of the transformation of `coll`.
  template <typename Coll, typename Op>
  auto transformCollection(Coll const& coll, Op op) {
    
    using Result_t
      = std::decay_t<decltype(op(std::declval<typename Coll::value_type>()))>;
    
    std::vector<Result_t> transformed;
    transformed.reserve(util::size(coll));
    std::transform
      (coll.begin(), coll.end(), std::back_inserter(transformed), op);
    return transformed;
  } // transformCollection()


  // --------------------------------------------------------------------------
  
} // local namespace


// -----------------------------------------------------------------------------
// ---  icarus::details::StableMerger
// -----------------------------------------------------------------------------
namespace icarus::details {
  template <typename SetColl> class StableMerger;
  template <typename SetColl> SetColl stableMerge(SetColl const& sets);
}

template <typename SetColl>
class icarus::details::StableMerger {
  
    public:
  using SetColl_t = SetColl;
  
  static SetColl_t merge(SetColl_t const& sets);
  
    private:
  using Coll_t = typename SetColl_t::value_type;
  using Value_t = typename Coll_t::value_type;
  
  static bool addNewValue(Coll_t& coll, Value_t const& value);
  
  static Coll_t mergeColls(Coll_t const& a, Coll_t const& b);
  
  static bool overlap(Coll_t const& a, Coll_t const& b);
  
  static SetColl_t mergePass(SetColl_t const& sets);
  
  
}; // class icarus::details::StableMerger<>


template <typename SetColl>
SetColl icarus::details::stableMerge(SetColl const& sets)
  { return icarus::details::StableMerger<SetColl>::merge(sets); }


// -----------------------------------------------------------------------------
template <typename SetColl>
bool icarus::details::StableMerger<SetColl>::addNewValue
  (Coll_t& coll, Value_t const& value)
{
  if (std::find(coll.begin(), coll.end(), value) != coll.end()) return false;
  coll.push_back(value);
  return true;
} // icarus::details::StableMerger<>::mergeColls()


// -----------------------------------------------------------------------------
template <typename SetColl>
auto icarus::details::StableMerger<SetColl>::mergeColls
  (Coll_t const& a, Coll_t const& b) -> Coll_t
{
  Coll_t merged;
  for (Value_t const& value: a) addNewValue(merged, value);
  for (Value_t const& value: b) addNewValue(merged, value);
  return merged;
} // icarus::details::StableMerger<>::mergeColls()


// -----------------------------------------------------------------------------
template <typename SetColl>
bool icarus::details::StableMerger<SetColl>::overlap
  (Coll_t const& a, Coll_t const& b)
{
  // brute force since we can't assume ordering
  for (auto const& aElem: a) {
    if (std::find(b.begin(), b.end(), aElem) != b.end()) return true;
  } // for
  return false;
} // icarus::details::StableMerger<>::overlap()


// -----------------------------------------------------------------------------
template <typename SetColl>
auto icarus::details::StableMerger<SetColl>::mergePass(SetColl_t const& sets)
  -> SetColl_t
{
  
  SetColl_t mergedSets;
  for (Coll_t const& set: sets) {
    
    //
    // overlap check
    //
    Coll_t* mergeWith = nullptr;
    for (Coll_t& mergedSet: mergedSets) {
      if (!overlap(mergedSet, set)) continue;
      mergeWith = &mergedSet;
      break;
    }
    
    if (mergeWith) {
      *mergeWith = (mergeWith->size() < set.size())
        ? mergeColls(set, *mergeWith)
        : mergeColls(*mergeWith, set)
        ;
    }
    else {
      mergedSets.push_back(set);
    }
    
  } // for input sets
  return mergedSets;
} // icarus::details::StableMerger<>::mergePass()


// -----------------------------------------------------------------------------
template <typename SetColl>
auto icarus::details::StableMerger<SetColl>::merge(SetColl_t const& sets)
  -> SetColl_t
{
  //
  // do merging until no more overlap is detected (i.e. no merge is performed)
  // (for example, if the sets are {0} {1} {0, 1}, the first pass will merge
  // {0, 1} with {0} but the resulting sets {0, 1} {1} still overlap
  //
  SetColl_t mergedSets { sets };
  std::size_t nSets;
  do {
    
    nSets = mergedSets.size();
    
    mergedSets = mergePass(mergedSets);
    
  } while (mergedSets.size() != nSets);
  
  return mergedSets;
} // icarus::details::StableMerger<SetColl>::merge()


// -----------------------------------------------------------------------------
// ---  icarus::details::StableMerger
// -----------------------------------------------------------------------------
namespace icarus::details { class ROPnumberDispatcher; }

/**
 * @brief Algorithm assigning IDs to readout planes.
 * 
 * The workflow is the following:
 * 
 * 1. the algorithm is initialized with a TPC set to be taken as model
 *    (constructor);
 * 2. a new TPC set is declared (`setTPCset()`);
 * 3. for each readout plane, an ID is requested (`assignID()`);
 * 4. if more TPC sets need to follow the same reference, repeat steps 2 and 3.
 * 
 * The assignment of the ROP number is in blocks according to the view.
 * The distribution between different views is first come, first serve: the
 * first ROP within a certain view conquers the lowest number still available
 * for all the ROP's within that view.
 * 
 * The algorithm assumes that all TPC sets have the same view composition as
 * the model one used in the constructor. If this requirement is violated,
 * the behavior is undefined (but it's likely that a ROP ID will be assigned
 * multiple times).
 */
class icarus::details::ROPnumberDispatcher {
  
  /// Type of collection of `T` by view.
  template <typename T>
  using DataByView_t = std::map<geo::View_t, T>;
  
  std::string const fLogCategory; ///< Log category name for messages.
  
  /// The first ROP number for each view.
  DataByView_t<readout::ROPID::ROPID_t> const fFirstROPno;
  
  /// The first available ROP number for each view.
  DataByView_t<readout::ROPID::ROPID_t> fAvailableROPno;
  
  /// The TPC set being served.
  readout::TPCsetID fTPCset;
  
  /// Returns the number of planes with each view.
  DataByView_t<unsigned int> ViewROPcounts
    (std::vector<PlaneColl_t> const& ROPplanes) const;
  
  /// Returns the fist ROP number for each view in `ROP`.
  DataByView_t<readout::ROPID::ROPID_t> PreferredROPrangesByView
    (std::vector<PlaneColl_t> const& ROPs) const;
  
    public:
  /// Constructor: uses `refROP` as reference to assign number ranges to views.
  ROPnumberDispatcher
    (std::vector<PlaneColl_t> const& refROP, std::string const& logCategory)
    : fLogCategory(logCategory)
    , fFirstROPno(PreferredROPrangesByView(refROP))
    , fAvailableROPno(fFirstROPno)
    {}
  
  /// Rests the object to work with the specified TPC set next.
  void setTPCset(readout::TPCsetID const& TPCset)
    { fTPCset = TPCset; fAvailableROPno = fFirstROPno; }
  
  /// Returns the next available ID for the specified ROP.
  readout::ROPID assignID(PlaneColl_t const& ROP)
    { return { fTPCset, fAvailableROPno.at(ROPview(ROP))++ }; }
  
  /// Returns the view of the `planes` (`geo::kUnknown` if mixed or none).
  static geo::View_t ROPview(PlaneColl_t const& planes);
  
}; // class icarus::details::ROPnumberDispatcher


// -----------------------------------------------------------------------------
auto icarus::details::ROPnumberDispatcher::PreferredROPrangesByView
  (std::vector<PlaneColl_t> const& ROPs) const
  -> DataByView_t<readout::ROPID::ROPID_t>
{
  /*
   * ROP numbers are assigned in blocks of views, on a first come, first serve
   * basis.
   * This, the first ROP on a view books as many ROP numbers as there are ROPs
   * with that view.
   */
  DataByView_t<unsigned int> const viewROPcounts = ViewROPcounts(ROPs);
  
  DataByView_t<readout::ROPID::ROPID_t> FirstROPno;
  readout::ROPID::ROPID_t nextAvailable = 0;
  
  for (PlaneColl_t const& planes: ROPs) {
    
    geo::View_t const view = ROPview(planes);
    
    //
    // using the facts that:
    // * std::map::emplace() does not insert anything if the key already exists
    // * std::map::emplace() returns as `second` whether it insert or not
    //
    if (!FirstROPno.emplace(view, nextAvailable).second) continue;
    
    mf::LogTrace(fLogCategory)
      << "ROPs on view " << geo::PlaneGeo::ViewName(view)
      << " preferring numbers from " << nextAvailable
      << " (" << viewROPcounts.at(view) << " numbers reserved)";
    
    // so, insertion happened; let's reserve numbers
    nextAvailable += viewROPcounts.at(view);
    
  } // for ROPs
  
  return FirstROPno;
  
} // icarus::details::ROPnumberDispatcher::PreferredROPrangesByView()


// ----------------------------------------------------------------------------
auto icarus::details::ROPnumberDispatcher::ViewROPcounts
  (std::vector<PlaneColl_t> const& ROPplanes) const
  -> DataByView_t<unsigned int>
{
  DataByView_t<unsigned int> viewROPcounts;
  
  for (PlaneColl_t const& planes: ROPplanes) {
    if (planes.empty()) continue;
    
    geo::View_t const view = ROPview(planes);
    if (view == geo::kUnknown) {
      mf::LogWarning log(fLogCategory);
      log << "A ROP spans multiple views:";
      for (geo::PlaneGeo const* plane: planes) {
        log << " " << geo::PlaneGeo::ViewName(plane->View())
          << " (" << plane->ID() << ")";
      }
    } // if unknown
    
    ++(viewROPcounts[view]);
    
  } // for ROPs
  
  return viewROPcounts;
} // icarus::details::ROPnumberDispatcher::ViewROPcounts()


// -----------------------------------------------------------------------------
geo::View_t icarus::details::ROPnumberDispatcher::ROPview
  (PlaneColl_t const& planes)
{
  auto iPlane = planes.begin();
  auto const pend = planes.end();
  if (iPlane == pend) return geo::kUnknown;
  geo::View_t const view = (*iPlane)->View();
  while (++iPlane != pend) {
    if ((*iPlane)->View() !=view) return geo::kUnknown;
  } // while
  return view;
} // icarus::details::ROPnumberDispatcher::ROPview()


// -----------------------------------------------------------------------------
// ---  icarus::details::ROPandTPCsetBuildingAlg
// -----------------------------------------------------------------------------
auto icarus::details::ROPandTPCsetBuildingAlg::run
  (geo::GeometryData_t::CryostatList_t const& Cryostats) -> Results_t
{
  /*
   *  The goals:
   *  (G1) number of actual TPC sets per cryostat
   *  (G2) composition in TPC's of each TPC set
   *  (G3) number of actual ROP's in each TPC set
   *  (G4) composition in readout planes of each ROP
   *
   *  The plan:
   *  
   *  1. extract the composition of all the readout planes
   *     and the TPCs each of them spans
   *  2. extract the final set of TPC sets by collecting sets of TPCs
   *     from all the ROP's
   *  3. assign each of the readout planes to a TPC set
   *  4. sort out the readout plane assignments into the final sets
   */
  
  //
  // input and setup
  //
  mf::LogDebug(fLogCategory)
    << "Building TPC sets and readout planes from " << Cryostats.size()
    << " cryostats";
  
  fCryostats = &Cryostats;
  
  //
  // output
  //
  clear();
  
  //
  // 1. extract the composition of all the readout planes
  //    and the TPCs they span
  //
  auto standaloneHorizontalWires = [](geo::PlaneGeo const& plane)
    { return std::abs(plane.ThetaZ()) < 1e-3; };
  
  std::vector<std::vector<PlaneColl_t>> AllPlanesInROPs
    = groupPlanesAndTPCs(standaloneHorizontalWires);
  
  std::vector<std::vector<std::vector<geo::TPCID>>> const AllTPCsInTPCsets
    = extractTPCsetsFromROPs(AllPlanesInROPs);
  
  //
  // 2. extract the final set of TPC sets by collecting sets of TPCs
  //    from all the ROP's
  //
  fillTPCsInSet(AllTPCsInTPCsets);
  
  //
  // 3. assign each of the readout planes to a TPC set
  //
  readout::TPCsetDataContainer<std::vector<PlaneColl_t>> const PlanesInProtoROPs
    = groupPlanesIntoROPs(AllTPCsInTPCsets, std::move(AllPlanesInROPs));
  
  //
  // 4. sort out the readout plane assignments into the final sets
  //
  fillPlanesInROP(PlanesInProtoROPs);
  
  //
  // 5. invert the maps
  //
  fillTPCtoTPCsetMap();
  fillPlaneToROPmap();
  
  //
  // output
  //
  assert(!fTPCsetCount.empty());
  assert(!fTPCsetTPCs .empty());
  assert(!fROPcount   .empty());
  assert(!fTPCtoTPCset.empty());
  assert(!fPlaneToROP .empty());
  return ResultsBase_t{
    std::move(fTPCsetCount),
    std::move(fTPCsetTPCs ),
    std::move(fROPcount   ),
    std::move(fROPplanes  ),
    std::move(fTPCtoTPCset),
    std::move(fPlaneToROP )
    };
} // icarus::details::ROPandTPCsetBuildingAlg::run()


// -----------------------------------------------------------------------------
void icarus::details::ROPandTPCsetBuildingAlg::clear() {
  fTPCsetCount.clear();
  fTPCsetTPCs.reset();
  fROPcount.reset();
  fROPplanes.reset();
  fTPCtoTPCset.reset();
  fPlaneToROP.reset();
  fMaxTPCsets = 0U;
  fMaxROPs = 0U;
} // icarus::details::ROPandTPCsetBuildingAlg::clear()


// -----------------------------------------------------------------------------
template <typename Pred>
auto icarus::details::ROPandTPCsetBuildingAlg::groupPlanesAndTPCs
  (Pred standalonePlane)
  -> std::vector<std::vector<PlaneColl_t>>
{
  /*
   * extract the composition of all the readout planes
   * and the TPCs each of them spans
   */
  
  //
  // input check
  //
  assert(fCryostats);
  geo::GeometryData_t::CryostatList_t const& Cryostats = *fCryostats;
  
  //
  // output
  //
  
  // for each cryostat (first index),
  // a list of readout plane candidates
  //   (second index: readout plane within the cryostat)
  // each listing its wire planes (third index: runs through all wire planes):
  std::vector<std::vector<PlaneColl_t>> AllPlanesInROPs;
  
  // debug
  auto logPlanes = [this](PlaneColl_t const& planes)
    {
      mf::LogTrace log(fLogCategory);
      log << "New group of " << planes.size() << " planes:";
      for (geo::PlaneGeo const* plane: planes)
        log << " <" << plane->ID() << ">";
    };
  
  //
  // the algorithm
  //
  // extract the information one cryostat at a time
  for (geo::CryostatGeo const& cryo: Cryostats) {
    
    //
    // connect all wire planes with the same drift direction to build
    // readout plane sets
    //
    std::vector<geo::PlaneGeo const*> planes;
    for (geo::TPCGeo const& tpc: cryo.IterateTPCs()) {
      for (geo::PlaneGeo const& plane: tpc.IteratePlanes()) {
        planes.push_back(&plane);
      } // for planes
    } // for TPCs
    
    std::vector<PlaneColl_t> protoGroups = groupPlanesByDriftCoord(planes);
    
    //
    // Are all the planes in the groups expected to be grouped?
    // if not, add them as their own group before the group.
    // We could have filtered out the planes not to be grouped before grouping,
    // which makes a lot of sense but complicates the sorting afterwards.
    //
    std::vector<PlaneColl_t> groupedPlanes;
    for (PlaneColl_t const& protoGroup: protoGroups) {
      if (protoGroup.empty()) continue;
      PlaneColl_t group;
      for (geo::PlaneGeo const* plane: protoGroup) {
        if (standalonePlane(*plane)) {
          groupedPlanes.push_back({ plane });
          logPlanes(groupedPlanes.back());
        }
        else group.push_back(plane);
      }
      if (!group.empty()) {
        groupedPlanes.push_back(std::move(group));
        logPlanes(groupedPlanes.back());
      }
    } // for protogroups
    
    //
    // save the results of the cryostat to detector-level data collection
    //
    AllPlanesInROPs.push_back(std::move(groupedPlanes));
    
  } // for cryostats
  
  return AllPlanesInROPs;
  
} // icarus::details::ROPandTPCsetBuildingAlg::groupPlanesAndTPCs()


// -----------------------------------------------------------------------------
auto icarus::details::ROPandTPCsetBuildingAlg::groupPlanesAndTPCs()
  -> std::vector<std::vector<PlaneColl_t>>
{ return groupPlanesAndTPCs([](geo::PlaneGeo const&){ return false; }); }


// -----------------------------------------------------------------------------
std::vector<std::vector<std::vector<geo::TPCID>>>
icarus::details::ROPandTPCsetBuildingAlg::extractTPCsetsFromROPs
  (std::vector<std::vector<PlaneColl_t>> const& groupedPlanes)
{
  
  // output: for each cryostat (first index),
  // a list of TPC sets (second index: TPC set within the cryostat)
  // each listing its TPC (third index: runs through all TPCs):
  std::vector<std::vector<std::vector<geo::TPCID>>> AllTPCsOnROPs;
  
  //
  // output
  //
  unsigned int& MaxTPCsets = fMaxTPCsets;
  
  MaxTPCsets = 0U;
  
  // extract the information one cryostat at a time
  for (std::vector<PlaneColl_t> const& planeGroupsInCryostat: groupedPlanes) {
    
    //
    // for each readout plane, collect the set of TPCs that contain its planes;
    // keep them as IDs since they are easier to procure.
    // Should these be sorted?
    // (yes) because the set is the same regardless the order of TPCs => set
    // (no) because we want to preserve the original order => vector
    //
    std::vector<std::vector<geo::TPCID>> TPCsOnROPs
      = transformCollection(planeGroupsInCryostat, extractTPCIDs);
    
    // reduce the duplicate TPC sets and merge overlapping sets
//     auto const iUniqueTPCend = std::unique(TPCsOnROPs.begin(), TPCsOnROPs.end());
//     TPCsOnROPs.erase(iUniqueTPCend, TPCsOnROPs.end());
    
    TPCsOnROPs = icarus::details::stableMerge(TPCsOnROPs);
    
    //
    // save the results of the cryostat to detector-level data collection
    //
    MaxTPCsets
      = std::max(MaxTPCsets, static_cast<unsigned int>(TPCsOnROPs.size()));
    AllTPCsOnROPs.push_back(std::move(TPCsOnROPs));
    
  } // for cryostats
  mf::LogTrace("ICARUSChannelMapAlg")
    << "The maximum number of TPC sets in any cryostat is " << MaxTPCsets;
  
  return AllTPCsOnROPs;
  
} // icarus::details::ROPandTPCsetBuildingAlg::extractTPCsetsFromROPs()


// -----------------------------------------------------------------------------
void icarus::details::ROPandTPCsetBuildingAlg::fillTPCsInSet
  (std::vector<std::vector<std::vector<geo::TPCID>>> const& AllTPCsOnROPs)
{
  /*
   * extract the final set of TPC sets by collecting sets of TPCs
   * from all the ROP's
   * (deliver G1 and G2)
   */
  
  //
  // input check
  //
  assert(fCryostats);
  geo::GeometryData_t::CryostatList_t const& Cryostats = *fCryostats;
  
  unsigned int const MaxTPCsets = fMaxTPCsets;
  
  //
  // output
  //
  // actual size of each TPC set: goal (G1)
  fTPCsetCount.clear();
  fTPCsetCount.resize(AllTPCsOnROPs.size(), 0U);
  std::vector<unsigned int>& TPCsetCount = fTPCsetCount;
  
  fTPCsetTPCs.resize(AllTPCsOnROPs.size(), MaxTPCsets); // goal (G2)
  readout::TPCsetDataContainer<std::vector<geo::TPCGeo const*>>& TPCsetTPCs
    = fTPCsetTPCs;
  
  //
  // algorithm
  //
  for (auto&& [ c, TPCsets ]: util::enumerate(AllTPCsOnROPs)) {
    geo::CryostatGeo const& cryo = Cryostats[c];
    TPCsetCount[c] = TPCsets.size();
    mf::LogTrace("ICARUSChannelMapAlg")
      << "Cryostat " << cryo.ID() << " has " << TPCsetCount[c] << " TPC sets";
    for (auto&& [ s, TPCIDs ]: util::enumerate(TPCsets)) {
      
      // convert the IDs into pointers to the TPC objects
      readout::TPCsetID const tpcsetid
        { cryo.ID(), static_cast<readout::TPCsetID::TPCsetID_t>(s) };
      TPCsetTPCs[tpcsetid] = transformCollection(
        TPCIDs, [&cryo](geo::TPCID const& tpcid){ return &(cryo.TPC(tpcid)); }
        );
      
      { // local block for debug output
        auto const& TPCs = TPCsetTPCs[tpcsetid];
        mf::LogTrace log("ICARUSChannelMapAlg");
        log << tpcsetid << " has " << TPCs.size() << " TPCs:";
        for (geo::TPCGeo const* tpc: TPCs) log << " " << tpc->ID() << ";";
      } // local block for debug output
      
    } // for TPC sets
  } // for cryostats
} // icarus::details::ROPandTPCsetBuildingAlg::fillTPCsInSet()


// -----------------------------------------------------------------------------
auto icarus::details::ROPandTPCsetBuildingAlg::groupPlanesIntoROPs(
  std::vector<std::vector<std::vector<geo::TPCID>>> const& AllTPCsOnROPs,
  std::vector<std::vector<PlaneColl_t>>&& AllPlanesInROPs
)
  -> readout::TPCsetDataContainer<std::vector<PlaneColl_t>>
{
  /*
   * assigns each of the readout planes to a TPC set
   */
  
  //
  // input check
  //
  assert(fCryostats);
  geo::GeometryData_t::CryostatList_t const& Cryostats = *fCryostats;
  
  assert(!fTPCsetTPCs.empty());
  readout::TPCsetDataContainer<TPCColl_t> const& TPCsetTPCs = fTPCsetTPCs;
  
  //
  // output
  //
  readout::TPCsetDataContainer<std::vector<PlaneColl_t>> PlanesInProtoROPs
    (TPCsetTPCs.dimSize<0U>(), TPCsetTPCs.dimSize<1U>());
  unsigned int& MaxROPs = fMaxROPs;
  
  
  unsigned int nErrors = 0U; // count errors, then bail out at the end
  // we still don't know the maximum number of ROPs in a TPC set
  MaxROPs = 0U;
  for (auto&& [ c, ROPs ]: util::enumerate(AllPlanesInROPs)) {
    // find which TPC set number each ROP belongs to;
    // brute force approach: check the content of each TPC set until we find
    // one that matched; then we are happy.
    // We assume that each ROP has one plane from *each* of the TPC's in a set.
    geo::CryostatGeo const& cryo = Cryostats[c];
    std::vector<std::vector<geo::TPCID>> const& TPCsets = AllTPCsOnROPs[c];
    for (std::vector<geo::PlaneGeo const*>& ROPplanes: ROPs) {
      std::vector<geo::TPCID> const ROPTPCIDs = extractTPCIDs(ROPplanes);
      
      // find the TPC set
      readout::TPCsetID tpcsetid; // destination set, invalid by default
      for (auto&& [ s, TPCset ]: util::enumerate(TPCsets)) {
        // here is where we assume the ROP has one plane per TOC in the TPC set;
        // we are also gambling that the order is the same
        if (!isROPinTPCset(ROPTPCIDs, TPCset)) continue;
        tpcsetid = { cryo.ID(), static_cast<readout::TPCsetID::TPCsetID_t>(s) };
        break;
      } // for all sets
      
      if (!tpcsetid) { // long error message to help debugging
        mf::LogError log(fLogCategory);
        log << "Candidate ROP did not match any TPC set.";
        log << "\nROP planes:";
        for (geo::PlaneGeo const* plane: ROPplanes)
          log << " (" << plane->ID() << ")";
        log << "\nAvailable TPC sets (" << TPCsets.size() << "):";
        for (auto&& [ s, TPCset ]: util::enumerate(TPCsets)) {
          readout::TPCsetID const tpcsetid
            { cryo.ID(), static_cast<readout::TPCsetID::TPCsetID_t>(s) };
          log << "\n - " << tpcsetid << ", " << TPCset.size() << " TPC's:";
          for (geo::TPCID const& tpcid: TPCset) log << " (" << tpcid << ")";
        } // for TPC sets
        ++nErrors;
        continue;
      } // if no TPC set matched
      
      // we have found the (first) TPC set matching the set of TPCs in ROP;
      // we store (move) the information into the proper cell
      auto& TPCsetROPs = PlanesInProtoROPs[tpcsetid];
      TPCsetROPs.push_back(std::move(ROPplanes));
      MaxROPs = std::max(MaxROPs, static_cast<unsigned int>(TPCsetROPs.size()));
      
    } // for TPC sets
  } // for cryostats
  
  AllPlanesInROPs.clear(); // we have already depleted it anyway
  
  if (nErrors > 0) {
    throw cet::exception(fLogCategory)
      << "Encountered " << nErrors
      << " errors while assigning TPC sets to ROPs (see error messages above)\n"
      ;
  } // if errors
  
  mf::LogTrace(fLogCategory)
    << "The maximum number of readout planes in any TPC set is " << MaxROPs;
  
  return PlanesInProtoROPs;
  
} // icarus::details::ROPandTPCsetBuildingAlg::groupPlanesIntoROPs()


// -----------------------------------------------------------------------------
void icarus::details::ROPandTPCsetBuildingAlg::fillPlanesInROP(
  readout::TPCsetDataContainer<std::vector<PlaneColl_t>> const& PlanesInProtoROPs
) {
  //
  // sort out the readout plane assignments into the final sets
  // (deliver G3 and G4: see `buildReadoutPlanes()`)
  //
  
  assert(!fTPCsetTPCs.empty());
  assert(!fTPCsetCount.empty());
  assert(fMaxROPs > 0U);
  
  //
  // prepare the input
  //
  auto const& TPCsetTPCs = fTPCsetTPCs;
  auto const& TPCsetCount = fTPCsetCount;
  
  //
  // prepare the output
  //
  // this is the actual number of readout planes in each TPC set: goal (G3)
  fROPcount.resize(TPCsetTPCs.dimSize<0U>(), TPCsetTPCs.dimSize<1U>(), 0U);
  fROPplanes.resize
    (TPCsetTPCs.dimSize<0U>(), TPCsetTPCs.dimSize<1U>(), fMaxROPs); // goal (G4)
  
  auto& ROPcount = fROPcount;
  auto& ROPplanes = fROPplanes;
  
  // assign the preferred numbers according to the first TPC set:
  icarus::details::ROPnumberDispatcher ROPIDdispatcher
    (PlanesInProtoROPs[{ 0, 0 }], fLogCategory);
  
  // (double) loop through all TPC sets:
  for (auto [ c, nTPCsets ]: util::enumerate(TPCsetCount)) {
    readout::CryostatID const cryoid
      { static_cast<readout::CryostatID::CryostatID_t>(c) };
    for (readout::TPCsetID::TPCsetID_t s = 0; s < nTPCsets; ++s) {
      readout::TPCsetID const tpcsetid { cryoid, s };
      
      // all the ROPs in this TPC set, with their wire plane content:
      std::vector<PlaneColl_t> const& ROPs = PlanesInProtoROPs[tpcsetid];
      
      ROPcount[tpcsetid] = ROPs.size();
      
      ROPIDdispatcher.setTPCset(tpcsetid);
      
      for (PlaneColl_t const& planes: ROPs) {
        
        readout::ROPID ropid = ROPIDdispatcher.assignID(planes);
        {
          mf::LogTrace log(fLogCategory);
          log << "Readout plane " << ropid << " assigned with " << planes.size()
            << " planes:";
          for (geo::PlaneGeo const* plane: planes)
            log << " (" << plane->ID() << ")";
        }
        if (!ROPplanes[ropid].empty()) {
          //
          // If this happens, it may be that the geometry is not compatible
          // with the algorithm, or just a bug.
          // Enabling the debug stream will show which planes are assigned
          // each time, including the two conflicting assignments.
          //
          throw cet::exception(fLogCategory)
            << "Logic error: ROPID " << ropid
            << " has already been assigned!\n";
        }
        ROPplanes[ropid] = std::move(planes);
        
      } // for all ROPs in the TPC set
    } // for all TPC sets in cryostat
  } // for all cryostats
  
} // icarus::details::ROPandTPCsetBuildingAlg::fillPlanesInROP()


// ----------------------------------------------------------------------------
void icarus::details::ROPandTPCsetBuildingAlg::fillTPCtoTPCsetMap() {
  
  //
  // invert the TPC sets map content
  //
  
  assert(!fCryostats->empty());
  assert(!fTPCsetTPCs.empty());
  
  //
  // prepare the input
  //
  auto const& Cryostats = *fCryostats;
  auto const& TPCsetTPCs = fTPCsetTPCs;
  
  //
  // output
  //
  auto const [ NCryostats, MaxTPCs ]
    = geo::details::extractMaxGeometryElements<2U>(Cryostats);

  MF_LOG_TRACE(fLogCategory)
    << "Detected " << NCryostats << " cryostats."
    << "\nDetected at most " << MaxTPCs << " TPCs per cryostat."
    ;
  fTPCtoTPCset.resize(NCryostats, MaxTPCs, {});
  
  for (auto c
    : util::counter<geo::CryostatID::CryostatID_t>(TPCsetTPCs.dimSize<0>())
    )
  {
    geo::CryostatID const cid { c };
    for (auto s: util::counter<readout::TPCsetID::TPCsetID_t>(TPCsetTPCs.dimSize<1>()))
    {
      readout::TPCsetID const sid { cid, s };
      
      for (geo::TPCGeo const* TPC: TPCsetTPCs[sid]) {
        assert(TPC && TPC->ID());
        fTPCtoTPCset[TPC->ID()] = sid;
        MF_LOG_TRACE(fLogCategory) << TPC->ID() << " => " << sid;
      } // for TPCs in TPC set
      
    } // for TPC sets in cryostat
    
  } // for cryostats
  
} // icarus::details::ROPandTPCsetBuildingAlg::fillTPCtoTPCsetMap()


// ----------------------------------------------------------------------------
void icarus::details::ROPandTPCsetBuildingAlg::fillPlaneToROPmap() {
  
  //
  // invert the TPC sets map content
  //
  
  assert(!fCryostats->empty());
  assert(!fROPplanes.empty());
  
  //
  // prepare the input
  //
  auto const& Cryostats = *fCryostats;
  auto const& ROPplanes = fROPplanes;
  
  //
  // output
  //
  auto const [ NCryostats, MaxTPCs, MaxPlanes ]
    = geo::details::extractMaxGeometryElements<3U>(Cryostats);

  MF_LOG_TRACE(fLogCategory)
    << "Detected " << NCryostats << " cryostats."
    << "\nDetected at most " << MaxTPCs << " TPCs per cryostat."
    << "\nDetected at most " << MaxPlanes << " planes per TPC."
    ;
  fPlaneToROP.resize(NCryostats, MaxTPCs, MaxPlanes, {});
  
  for (auto c
    : util::counter<geo::CryostatID::CryostatID_t>(ROPplanes.dimSize<0>())
    )
  {
    geo::CryostatID const cid { c };
    for (auto s: util::counter<readout::TPCsetID::TPCsetID_t>(ROPplanes.dimSize<1>()))
    {
      readout::TPCsetID const sid { cid, s };
      for (auto r: util::counter<readout::ROPID::ROPID_t>(ROPplanes.dimSize<2>()))
      {
        readout::ROPID const rid { sid, r };
      
        for (geo::PlaneGeo const* plane: ROPplanes[rid]) {
          assert(plane && plane->ID());
          fPlaneToROP[plane->ID()] = rid;
          MF_LOG_TRACE(fLogCategory) << plane->ID() << " => " << rid;
        } // for planes in readout plane
      
      } // for readout planes in TPC set
      
    } // for TPC sets in cryostat
    
  } // for cryostats
  
} // icarus::details::ROPandTPCsetBuildingAlg::fillPlaneToROPmap()


// ----------------------------------------------------------------------------
template <typename PlaneColl>
auto icarus::details::ROPandTPCsetBuildingAlg::groupPlanesByDriftCoord
  (PlaneColl const& planes, double tolerance /* = 0.1 */)
  -> std::vector<PlaneColl_t>
{
  std::map<double, PlaneColl_t> groupedByDrift;
  
  geo::Vector_t const driftDir = geo::Xaxis();
  
  /*
   * we can't std::sort by drift coordinate because that would break the
   * stability of the sorting (we are guaranteeing that if plane _A_, drift
   * coordinate _x_, is before plane _B_, drift coordinate _x - epsilon_,
   * and _epsilon_ is smaller than `tolerance`, then in the group plane _A_
   * should still be of plane _B_.
   */
  for (geo::PlaneGeo const* plane: planes) {
    
    double const planeD = plane->GetCenter<geo::Point_t>().Dot(driftDir);
    // find the group before the plane
    
    auto iGroup = groupedByDrift.lower_bound(planeD);
    
    //
    // if the group we found is compatible with the plane, we add it in;
    // the key is the upper limit of the group coordinate range
    //
    if ((iGroup != groupedByDrift.end())
      && (iGroup->first - planeD <= tolerance))
    {
      iGroup->second.push_back(plane);
      continue;
    }
    //
    // the plane has a drift coordinate too small to join the group we found;
    // or we did not find any; so we create a new group centered on this plane;
    // the key is the right end of the allowed range, to simplify both search
    // and check
    //
    groupedByDrift.emplace_hint
      (iGroup, planeD + tolerance / 2.0, PlaneColl_t({ plane }));
    
  } // for
  
  //
  // moving the groups into the result structure; they'll be sorted by
  // increasing nominal drift coordinate
  //
  std::vector<PlaneColl_t> groups;
  groups.reserve(groupedByDrift.size());
  for (auto&& group: groupedByDrift)
    groups.push_back(std::move(std::get<1U>(group)));
  return groups;
  
} // icarus::details::ROPandTPCsetBuildingAlg::groupPlanesByDriftCoord()


// ----------------------------------------------------------------------------
std::vector<geo::TPCID> icarus::details::ROPandTPCsetBuildingAlg::extractTPCIDs
  (std::vector<geo::PlaneGeo const*> const& planes)
{
  return transformCollection(planes,
    [](geo::PlaneGeo const* plane){ return plane->ID().asTPCID(); }
    );
} // icarus::details::ROPandTPCsetBuildingAlg::extractTPCIDs()


// -----------------------------------------------------------------------------
readout::ROPID::ROPID_t
icarus::details::ROPandTPCsetBuildingAlg::ROPnumberFromPlanes
  (PlaneColl_t const& planes)
{
  readout::ROPID::ROPID_t r = readout::ROPID::getInvalidID();
  for (geo::PlaneGeo const* plane: planes) {
    if (!plane) continue;
    
    auto const fromPlane
      = static_cast<readout::ROPID::ROPID_t>(plane->ID().Plane);
    
    if (r == readout::ROPID::getInvalidID())
      r = fromPlane;
    else if (r != fromPlane)
      return readout::ROPID::getInvalidID();
    
  } // for
  return r;
} // icarus::details::ROPandTPCsetBuildingAlg::ROPnumberFromPlanes()


// -----------------------------------------------------------------------------
bool icarus::details::ROPandTPCsetBuildingAlg::isROPinTPCset(
  std::vector<geo::TPCID> const& ROPTPCIDs,
  std::vector<geo::TPCID> const& TPCsetTPCIDs
) {
  
  for (geo::TPCID const& ROPTPC: ROPTPCIDs) {
    bool found = false;
    for (geo::TPCID const& TPC: TPCsetTPCIDs) {
      if (TPC != ROPTPC) continue;
      found = true;
      break;
    } // for all TPCs in set
    if (found) continue;
    return false;
  } // for all TPCs in ROP
  return true;
} // icarus::details::ROPandTPCsetBuildingAlg::isROPinTPCset()


// ----------------------------------------------------------------------------
