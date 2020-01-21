/**
 * @file   icaruscode/Geometry/ICARUSChannelMapAlg.cxx
 * @brief  Channel mapping algorithms for ICARUS detector: implementation file.
 * @date   October 19, 2019
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/Geometry/ICARUSChannelMapAlg.h
 */

// library header
#include "icaruscode/Geometry/ICARUSChannelMapAlg.h"

// ICARUS libraries
#include "icaruscode/Geometry/details/ROPandTPCsetBuildingAlg.h"

// LArSoft libraries
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/ReadoutDataContainers.h"
#include "larcorealg/Geometry/details/extractMaxGeometryElements.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::printBacktrace()
#include "larcorealg/CoreUtils/StdUtils.h" // util::size()
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h" // cet::exception

// C/C++ libraries
#include <vector>
#include <set>
#include <algorithm> // std::transform(), std::find()
#include <utility> // std::move()
#include <iterator> // std::back_inserter()
#include <functional> // std::mem_fn()
#include <tuple>


namespace {
  
  // ---------------------------------------------------------------------------
  fhicl::ParameterSet getOptionalParameterSet
    (fhicl::OptionalDelegatedParameter const& param)
  {
    fhicl::ParameterSet pset; // empty by default
    param.get_if_present(pset);
    return pset;
  } // fhicl::ParameterSet getOptionalParameterSet()
  
  
  // ---------------------------------------------------------------------------
  
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


  // ---------------------------------------------------------------------------
  
} // local namespace


// -----------------------------------------------------------------------------
icarus::ICARUSChannelMapAlg::ICARUSChannelMapAlg(Config const& config)
  : fSorter(getOptionalParameterSet(config.Sorter))
  {}


// -----------------------------------------------------------------------------
void icarus::ICARUSChannelMapAlg::Initialize(geo::GeometryData_t const& geodata)
{
  // This is the only INFO level message we want this object to produce;
  // given the dynamic nature of the channel mapping choice,
  // it's better for the log to have some indication of chosen channel mapping.
  mf::LogInfo("ICARUSChannelMapAlg")
    << "Initializing ICARUSChannelMapAlg channel mapping algorithm.";
  
  buildReadoutPlanes(geodata.cryostats);
  
  fillChannelToWireMap(geodata.cryostats);
  
  MF_LOG_TRACE("ICARUSChannelMapAlg")
    << "ICARUSChannelMapAlg::Initialize() completed.";
  
} // icarus::ICARUSChannelMapAlg::Initialize()


// -----------------------------------------------------------------------------
void icarus::ICARUSChannelMapAlg::Uninitialize() {
  
  fReadoutMapInfo.clear();
  
  fChannelToWireMap.clear();
  
  fPlaneInfo.clear();
  
} // icarus::ICARUSChannelMapAlg::Uninitialize()


//------------------------------------------------------------------------------
std::vector<geo::WireID> icarus::ICARUSChannelMapAlg::ChannelToWire
  (raw::ChannelID_t channel) const
{
  //
  // input check
  //
  assert(!fPlaneInfo.empty());
  
  //
  // output
  //
  std::vector<geo::WireID> AllSegments;
  
  //
  // find the ROP with that channel
  //
  icarus::details::ChannelToWireMap::ChannelsInROPStruct const* channelInfo
    = fChannelToWireMap.find(channel);
  if (!channelInfo) {
    throw cet::exception("Geometry")
      << "icarus::ICARUSChannelMapAlg::ChannelToWire(" << channel
      << "): invalid channel requested (must be lower than "
      << Nchannels() << ")\n";
  }
  
  //
  // find the wire planes in that ROP
  //
  PlaneColl_t const& planes = ROPplanes(channelInfo->ropid);
  
  //
  // associate one wire for each of those wire planes to the channel
  //
  AllSegments.reserve(planes.size()); // this is sometimes (often?) too much
  for (geo::PlaneGeo const* plane: planes) {
    
    geo::PlaneID const& pid = plane->ID();
    ChannelRange_t const& channelRange = fPlaneInfo[pid].channelRange();
    
    if (!channelRange.contains(channel)) continue;
    AllSegments.emplace_back
      (pid, static_cast<geo::WireID::WireID_t>(channel - channelRange.begin()));
    
  } // for planes in ROP
  
  return AllSegments;
  
} // icarus::ICARUSChannelMapAlg::ChannelToWire()


//------------------------------------------------------------------------------
unsigned int icarus::ICARUSChannelMapAlg::Nchannels() const {
  
  return fChannelToWireMap.nChannels();
  
} // icarus::ICARUSChannelMapAlg::Nchannels()


//------------------------------------------------------------------------------
unsigned int icarus::ICARUSChannelMapAlg::Nchannels
  (readout::ROPID const& ropid) const 
{
  icarus::details::ChannelToWireMap::ChannelsInROPStruct const* ROPinfo
    = fChannelToWireMap.find(ropid);
  return ROPinfo? ROPinfo->nChannels: 0U;
} // icarus::ICARUSChannelMapAlg::Nchannels(ROPID)


//------------------------------------------------------------------------------
double icarus::ICARUSChannelMapAlg::WireCoordinate
  (double YPos, double ZPos, geo::PlaneID const& planeID) const
{
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ICARUSChannelMapAlg");
  e << "ICARUSChannelMapAlg does not support `WireCoordinate()` call."
    "\nPlease update calling software to use geo::PlaneGeo::WireCoordinate()`:"
    "\n";
  
  lar::debug::printBacktrace(e, 4U);
  
  throw e;
} // icarus::ICARUSChannelMapAlg::WireCoordinate()


//------------------------------------------------------------------------------
geo::WireID icarus::ICARUSChannelMapAlg::NearestWireID
  (const TVector3& worldPos, geo::PlaneID const& planeID) const
{
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ICARUSChannelMapAlg");
  e << "ICARUSChannelMapAlg does not support `NearestWireID()` call."
    "\nPlease update calling software to use geo::PlaneGeo::NearestWireID()`:"
    "\n";
  
  lar::debug::printBacktrace(e, 3U);
  
  throw e;
} // icarus::ICARUSChannelMapAlg::NearestWireID()


//------------------------------------------------------------------------------
raw::ChannelID_t icarus::ICARUSChannelMapAlg::PlaneWireToChannel
  (geo::WireID const& wireID) const
{
  return fPlaneInfo[wireID].firstChannel() + wireID.Wire;
} // icarus::ICARUSChannelMapAlg::PlaneWireToChannel()


//------------------------------------------------------------------------------
std::set<geo::PlaneID> const& icarus::ICARUSChannelMapAlg::PlaneIDs() const {
  
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ICARUSChannelMapAlg");
  e << "ICARUSChannelMapAlg does not support `PlaneIDs()` call."
    "\nPlease update calling software to use geo::GeometryCore::IteratePlanes()`"
    "\n";
  
  lar::debug::printBacktrace(e, 3U);
  
  throw e;
  
} // icarus::ICARUSChannelMapAlg::PlaneIDs()


//------------------------------------------------------------------------------
unsigned int icarus::ICARUSChannelMapAlg::NTPCsets
  (readout::CryostatID const& cryoid) const
{
  return HasCryostat(cryoid)? TPCsetCount(cryoid): 0U;
} // icarus::ICARUSChannelMapAlg::NTPCsets()


//------------------------------------------------------------------------------
/// Returns the largest number of TPC sets any cryostat in the detector has
unsigned int icarus::ICARUSChannelMapAlg::MaxTPCsets() const {
  assert(fReadoutMapInfo);
  return fReadoutMapInfo.MaxTPCsets();
} // icarus::ICARUSChannelMapAlg::MaxTPCsets()


//------------------------------------------------------------------------------
/// Returns whether we have the specified TPC set
/// @return whether the TPC set is valid and exists
bool icarus::ICARUSChannelMapAlg::HasTPCset
  (readout::TPCsetID const& tpcsetid) const
{
  return
    HasCryostat(tpcsetid)? (tpcsetid.TPCset < TPCsetCount(tpcsetid)): false;
} // icarus::ICARUSChannelMapAlg::HasTPCset()


//------------------------------------------------------------------------------
readout::TPCsetID icarus::ICARUSChannelMapAlg::TPCtoTPCset
  (geo::TPCID const& tpcid) const
{
  return tpcid? TPCtoTPCset()[tpcid]: readout::TPCsetID{};
} // icarus::ICARUSChannelMapAlg::TPCtoTPCset()


//------------------------------------------------------------------------------
std::vector<geo::TPCID> icarus::ICARUSChannelMapAlg::TPCsetToTPCs
  (readout::TPCsetID const& tpcsetid) const
{
  std::vector<geo::TPCID> TPCs;
  if (!tpcsetid) return TPCs;
  
  auto const& TPClist = TPCsetTPCs(tpcsetid);
  TPCs.reserve(TPClist.size());
  std::transform(TPClist.begin(), TPClist.end(), std::back_inserter(TPCs),
    std::mem_fn(&geo::TPCGeo::ID)
    );
  return TPCs;
} // icarus::ICARUSChannelMapAlg::TPCsetToTPCs()


//------------------------------------------------------------------------------
geo::TPCID icarus::ICARUSChannelMapAlg::FirstTPCinTPCset
  (readout::TPCsetID const& tpcsetid) const
{
  if (!tpcsetid) return {};
  
  auto const& TPClist = TPCsetTPCs(tpcsetid);
  return TPClist.empty()? geo::TPCID{}: TPClist.front()->ID();
  
} // icarus::ICARUSChannelMapAlg::FirstTPCinTPCset()


//------------------------------------------------------------------------------
unsigned int icarus::ICARUSChannelMapAlg::NROPs
  (readout::TPCsetID const& tpcsetid) const
{
  return HasTPCset(tpcsetid)? ROPcount(tpcsetid): 0U;
} // icarus::ICARUSChannelMapAlg::NROPs()


//------------------------------------------------------------------------------
unsigned int icarus::ICARUSChannelMapAlg::MaxROPs() const {
  assert(fReadoutMapInfo);
  return fReadoutMapInfo.MaxROPs();
} // icarus::ICARUSChannelMapAlg::MaxROPs()

//------------------------------------------------------------------------------
bool icarus::ICARUSChannelMapAlg::HasROP(readout::ROPID const& ropid) const {
  return HasTPCset(ropid)? (ropid.ROP < ROPcount(ropid)): false;
} // icarus::ICARUSChannelMapAlg::HasROP()


//------------------------------------------------------------------------------
  /**
   * @brief Returns the ID of the ROP planeid belongs to, or invalid if none
   * @param planeid ID of the plane
   * @return the ID of the corresponding ROP, or invalid ID when planeid is
   *
   * In this mapping, readout planes and wire planes are mapped one-to-one.
   * The returned value mirrors the plane ID in the readout space.
   * If the plane ID is not valid, an invalid readout plane ID is returned.
   * Note that this check is performed on the validity of the plane ID, that
   * does not necessarily imply that the plane specified by the ID actually
   * exists.
   */
readout::ROPID icarus::ICARUSChannelMapAlg::WirePlaneToROP
  (geo::PlaneID const& planeid) const
{
  return planeid? PlaneToROP(planeid): readout::ROPID{};
} // icarus::ICARUSChannelMapAlg::WirePlaneToROP()


//------------------------------------------------------------------------------
std::vector<geo::PlaneID> icarus::ICARUSChannelMapAlg::ROPtoWirePlanes
  (readout::ROPID const& ropid) const
{
  std::vector<geo::PlaneID> Planes;
  if (!ropid) return Planes;
  
  auto const& PlaneList = ROPplanes(ropid);
  Planes.reserve(PlaneList.size());
  std::transform(PlaneList.begin(), PlaneList.end(), std::back_inserter(Planes),
    std::mem_fn(&geo::PlaneGeo::ID)
    );
  return Planes;
} // icarus::ICARUSChannelMapAlg::ROPtoWirePlanes()


//------------------------------------------------------------------------------
std::vector<geo::TPCID> icarus::ICARUSChannelMapAlg::ROPtoTPCs
  (readout::ROPID const& ropid) const
{
  std::vector<geo::TPCID> TPCs;
  if (!ropid) return TPCs;
  
  /*
   * We use the same algorithm as for extracting the plane IDs
   * (they implicitly convert to TPC ID... kind of them).
   * The algorithm does not test for duplication, i.e. in theory it could
   * produce lists with the same TPC ID being present multiple times
   * from different planes.
   * But this is not expected in this mapping, where each TPC holds at most
   * one wire plane for each view, and the planes in a ROP are all on the same
   * view. It might be matter of an assertion, but it's too complex to fit in an
   * `assert()` call.
   */
  auto const& PlaneList = ROPplanes(ropid);
  TPCs.reserve(PlaneList.size());
  std::transform(PlaneList.begin(), PlaneList.end(), std::back_inserter(TPCs),
    std::mem_fn(&geo::PlaneGeo::ID)
    );
  return TPCs;
} // icarus::ICARUSChannelMapAlg::ROPtoTPCs()


//------------------------------------------------------------------------------
readout::ROPID icarus::ICARUSChannelMapAlg::ChannelToROP
  (raw::ChannelID_t channel) const
{
  if (!raw::isValidChannelID(channel)) return {};
  icarus::details::ChannelToWireMap::ChannelsInROPStruct const* info
    = fChannelToWireMap.find(channel);
  return info? info->ropid: readout::ROPID{};
} // icarus::ICARUSChannelMapAlg::ChannelToROP()


//------------------------------------------------------------------------------
raw::ChannelID_t icarus::ICARUSChannelMapAlg::FirstChannelInROP
  (readout::ROPID const& ropid) const
{
  if (!ropid) return raw::InvalidChannelID;
  icarus::details::ChannelToWireMap::ChannelsInROPStruct const* info
    = fChannelToWireMap.find(ropid);
  return info? info->firstChannel: raw::InvalidChannelID;
} // icarus::ICARUSChannelMapAlg::FirstChannelInROP()


//------------------------------------------------------------------------------
geo::PlaneID icarus::ICARUSChannelMapAlg::FirstWirePlaneInROP
  (readout::ROPID const& ropid) const
{
  if (!ropid) return {};
  PlaneColl_t const& planes = ROPplanes(ropid);
  return planes.empty()? geo::PlaneID{}: planes.front()->ID();
} // icarus::ICARUSChannelMapAlg::FirstWirePlaneInROP()


//------------------------------------------------------------------------------
bool icarus::ICARUSChannelMapAlg::HasCryostat
  (readout::CryostatID const& cryoid) const
{
  assert(fReadoutMapInfo);
  return cryoid.Cryostat < fReadoutMapInfo.NCryostats();
} // icarus::ICARUSChannelMapAlg::HasCryostat()


//------------------------------------------------------------------------------
void icarus::ICARUSChannelMapAlg::fillChannelToWireMap
  (geo::GeometryData_t::CryostatList_t const& Cryostats)
{
  
  //
  // input check
  //
  assert(fReadoutMapInfo);
  assert(!Cryostats.empty());
  
  //
  // output setup
  //
  assert(fPlaneInfo.empty());
  std::array<unsigned int, 3U> maxSizes
    = geo::details::extractMaxGeometryElements<3U>(Cryostats);

  fPlaneInfo.resize(maxSizes[0U], maxSizes[1U], maxSizes[2U]);
  
  
  raw::ChannelID_t nextChannel = 0; // next available channel
  
  // once again we do not have iteration facilities from `geo::GeometryCore`
  // available yet, so we go the nested loop way and bite it
  for (geo::CryostatGeo const& cryo: Cryostats) {
    
    readout::CryostatID const cid { cryo.ID() };
    
    auto const nTPCsets 
      = static_cast<readout::TPCsetID::TPCsetID_t>(TPCsetCount(cid));
    
    for (readout::TPCsetID::TPCsetID_t s: util::counter(nTPCsets)) {
      
      readout::TPCsetID const sid { cid, s };
      
      auto const nROPs = static_cast<readout::ROPID::ROPID_t>(ROPcount(sid));
      
      for (readout::ROPID::ROPID_t r: util::counter(nROPs)) {
        
        mf::LogTrace log("ICARUSChannelMapAlg");
        
        readout::ROPID const rid { sid, r };
        log << "ROP: " << rid;
        
        PlaneColl_t const& planes = ROPplanes(rid);
        log << " (" << planes.size() << " planes):";
        assert(!planes.empty());
        
        raw::ChannelID_t const firstROPchannel = nextChannel;
        
        auto iPlane = util::begin(planes);
        auto const pend = util::end(planes);
        
        // assign available channels to all wires of the first plane
        nextChannel += (*iPlane)->Nwires();
        fPlaneInfo[(*iPlane)->ID()] = {
          ChannelRange_t{ firstROPchannel, nextChannel },
          rid
          };
        log << " [" << (*iPlane)->ID() << "] "
          << fPlaneInfo[(*iPlane)->ID()].firstChannel()
          << " -- " << fPlaneInfo[(*iPlane)->ID()].lastChannel() << ";";
        
        geo::Point_t lastWirePos = (*iPlane)->LastWire().GetCenter<geo::Point_t>();
        
        while (++iPlane != pend) {
          
          geo::PlaneGeo const& plane = **iPlane;
          
          // find out which wire matches the last wire from the previous plane;
          // if there is no such wire, an exception will be thrown,
          // which is ok to us since we believe it should not happen.
          geo::WireID const lastMatchedWireID
            = plane.NearestWireID(lastWirePos);
          
          /*
          mf::LogTrace("ICARUSChannelMapAlg")
            << (*std::prev(iPlane))->ID() << " W:" << ((*std::prev(iPlane))->Nwires() - 1)
            << " ending at " << (*std::prev(iPlane))->LastWire().GetEnd<geo::Point_t>()
            << " matched " << lastMatchedWireID
            << " which starts at " << plane.Wire(lastMatchedWireID).GetStart<geo::Point_t>()
            ;
          */
          
          //
          // the first channel in this plane (`iPlane`) is the one associated
          // to the first wire in the plane, which has local wire number `0`;
          // the last channel from the previous plane (`nextChannel - 1`)
          // is associated to the matched wire (`lastMatchedWireID`),
          // which has some wire number (`lastMatchedWireID.Wire`).
          //
          auto const nWires = plane.Nwires();
          raw::ChannelID_t const firstChannel
            = (nextChannel - 1) - lastMatchedWireID.Wire;
          nextChannel = firstChannel + nWires;
          
          fPlaneInfo[plane.ID()] = { { firstChannel, nextChannel }, rid };
          log << " [" << plane.ID() << "] "
            << fPlaneInfo[plane.ID()].firstChannel() << " -- "
            << fPlaneInfo[plane.ID()].lastChannel() << ";";
          
          // update for the next iteration
          lastWirePos = plane.LastWire().GetCenter<geo::Point_t>();
          
        } // while
        
        unsigned int const nChannels = nextChannel - firstROPchannel;
        fChannelToWireMap.addROP(rid, firstROPchannel, nChannels);
        log
          << " => " << nChannels << " channels starting at " << firstROPchannel;
        
      } // for readout plane
      
    } // for TPC set
    
  } // for cryostat
  
  fChannelToWireMap.setEndChannel(nextChannel);
  mf::LogTrace("ICARUSChannelMapAlg")
    << "Counted " << fChannelToWireMap.nChannels() << " channels.";
  
} // icarus::ICARUSChannelMapAlg::fillChannelToWireMap()


// -----------------------------------------------------------------------------
void icarus::ICARUSChannelMapAlg::buildReadoutPlanes
  (geo::GeometryData_t::CryostatList_t const& Cryostats)
{
  // the algorithm is delegated:
  icarus::details::ROPandTPCsetBuildingAlg builder("ICARUSChannelMapAlg");
  
  auto results = builder.run(Cryostats);
  
  fReadoutMapInfo.set(
    std::move(results).TPCsetCount(), std::move(results).TPCsetTPCs(),
    std::move(results).ROPcount(), std::move(results).ROPplanes(),
    std::move(results).TPCtoTPCset(), std::move(results).PlaneToROP()
    );
  
} // icarus::ICARUSChannelMapAlg::buildReadoutPlanes()


// ----------------------------------------------------------------------------
geo::SigType_t icarus::ICARUSChannelMapAlg::SignalTypeForChannelImpl
  (raw::ChannelID_t const channel) const
{
  /*
   * The current mapping algorithm groups the wire planes by _x_, starting from
   * the lowest ones.
   * In a cathode-in-the-middle geometry, this makes the wire planes with higher
   * number, which are farther from the cathode, to be picked first.
   * 
   * The same algorithm assigns the first readout plane to the lowest number.
   * Since the first readout plane contains the first wire plane group, and that
   * is the farthest from the cathode and has the largest number, the first
   * readout plane, with number C:0 S:0 R:0, is assigned to a group pf planes
   * including the likes of C:0 T:0 P:2.
   * 
   * Without bothering explicitly with coordinates, those assumptions mean that
   * the ROPs with the lowest number (R:0) and the planes with the highest (P:2)
   * are the farthest from the cathode, and they are therefore the collection
   * planes, while all others are induction planes.
   * 
   * We have the choice of which of the two facts to rely upon.
   * The readout one is somehow more fragile because it relies on one more
   * convention (the assignment of readout plane numbers to plane groups) that
   * may change without the rest being affected.
   * On the other end, it's faster to implement because the access to the
   * ROP from the channel number is more direct, and the number to compare to
   * is always `0` instead of the actual number of planes minus one.
   * Here we go the latter way.
   */
  
  icarus::details::ChannelToWireMap::ChannelsInROPStruct const* channelInfo
    = fChannelToWireMap.find(channel);
  
  return channelInfo?
    ((channelInfo->ropid.ROP == 0)
      ? geo::kCollection
      : geo::kInduction
    )
    : geo::kMysteryType
    ;
  
} // icarus::ICARUSChannelMapAlg::SignalTypeForChannelImpl()


// ----------------------------------------------------------------------------

