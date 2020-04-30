/**
 * @file   icaruscode/PMT/LibraryMappingTools/ICARUSPhotonMappingTransformations.cxx
 * @brief  A photon mapping identity transformation: implementation.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 3, 2019
 * @see    `icaruscode/PMT/LibraryMappingTools/ICARUSPhotonMappingTransformations.h`
 * @see    `icaruscode/PMT/LibraryMappingTools/ICARUSPhotonMappingTransformations.cc`
 * 
 */


// ICARUS libraries
#include "icaruscode/PMT/LibraryMappingTools/ICARUSPhotonMappingTransformations.h"

// LArSoft libraries
#include "lardataalg/Utilities/StatCollector.h" // lar::util::MinMaxCollector<>
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h" // MF_LOG_DEBUG()

// C/C++ standard libraries
#include <algorithm> // std::iota()
#include <string> // std::to_string()
#include <cassert>


//------------------------------------------------------------------------------
phot::ICARUSPhotonMappingTransformations::ICARUSPhotonMappingTransformations
  (Config const& config)
  : fDumpMapping(config.DumpMapping())
  , fGeom(lar::providerFrom<geo::Geometry>())
  , fNOpDetChannels(fGeom? fGeom->NOpDets(): 0)
{
  LibraryIndexToOpDetMap libraryIndices;
  if (!config.CryostatChannelRemap(libraryIndices)) {
    assert(fNOpDetChannels > 0U); // if no mapping is specified, we need to know
    unsigned int nCryoChannels = fGeom->Cryostat(0U).NOpDet();
    libraryIndices.resize(nCryoChannels);
    std::iota(libraryIndices.begin(), libraryIndices.end(), 0U);
  }
  prepareMappings(libraryIndices);
  
} // phot::ICARUSPhotonMappingTransformations::ICARUSPhotonMappingTransformations()


//------------------------------------------------------------------------------
geo::Point_t phot::ICARUSPhotonMappingTransformations::detectorToLibrary
  (geo::Point_t const& location) const
{
  geo::CryostatID onCryo = whichCryostat(location);
  return onCryo? location + fTranslations[onCryo.Cryostat]: location;
} // phot::ICARUSPhotonMappingTransformations::detectorToLibrary()


//------------------------------------------------------------------------------
auto phot::ICARUSPhotonMappingTransformations::opDetsToLibraryIndicesImpl
  (geo::Point_t const& location) const -> OpDetToLibraryIndexMap const&
{
  geo::CryostatID onCryo = whichCryostat(location);
  return onCryo
    ? fOpDetToLibraryIndexMaps[onCryo.Cryostat]: fInvalidOpDetToLibraryIndexMap;
} // phot::ICARUSPhotonMappingTransformations::opDetsToLibraryIndicesImpl()


//------------------------------------------------------------------------------
auto phot::ICARUSPhotonMappingTransformations::libraryIndicesToOpDetsImpl
  (geo::Point_t const& location) const -> LibraryIndexToOpDetMap const& 
{
  geo::CryostatID onCryo = whichCryostat(location);
  return onCryo
    ? fLibraryIndexToOpDetMaps[onCryo.Cryostat]: fInvalidLibraryIndexToOpDetMap;
} // phot::ICARUSPhotonMappingTransformations::libraryIndicesToOpDetsImpl()


//------------------------------------------------------------------------------
void phot::ICARUSPhotonMappingTransformations::prepareGeometryMapping() {
  
  /*
   * geometry transformation:
   * 1. determine the switch coordinate (on x axis) between cryostats
   * 2. determine the amount of shift required
   */
  
  //
  // (1) geometry transformation
  //
  // (1.1) determine the switch coordinate (on x axis) between cryostats
  
  geo::Length_t const C0endX = fGeom->Cryostat(0U).MaxX();
  geo::Length_t const C1startX = fGeom->Cryostat(1U).MinX();
  if (C0endX > C1startX) {
    throw std::runtime_error(
      "phot::ICARUSPhotonMappingTransformations::prepareMappings(): "
      "C:0 ends at x=" + std::to_string(C0endX)
      + ", C:1 starts at x=" + std::to_string(C1startX)
      + "... this algorithm does not understand this geometry."
      );
  }
  fSwitchPoint = (C0endX + C1startX) / 2.0; // pick the middle point
  MF_LOG_DEBUG("ICARUSPhotonMappingTransformations")
    << "Switching from cryostat 0 to 1 when x > " << fSwitchPoint << " cm";
  
  // (1.2) determine the amount of shift required
  geo::Point_t const refPoint = fGeom->Cryostat(0).BoundingBox().Min();
  for (geo::CryostatGeo const& cryo: fGeom->IterateCryostats()) {
    fTranslations.push_back(refPoint - cryo.BoundingBox().Min());
  } // for all cryostats
  
} // phot::ICARUSPhotonMappingTransformations::prepareGeometryMapping()

  
//------------------------------------------------------------------------------
void phot::ICARUSPhotonMappingTransformations::prepareLibraryMappings
  (LibraryIndexToOpDetMap const& libraryIndices)
{
  
  /*
   * 2. library transformation
   *    1. determine the range of PMT channels in the two cryostats
   *    2. determine the shift amount (hint: 180)
   *    3. fill the maps with the shift
   */
  
  //
  // (2) library transformation
  //
  // (2.1) determine the range of PMT channels in the two cryostats
  /*
   * This algorithm is a bit cumbersome, but with the complications of optical
   * detector geometry in LArSoft (which are usually irrelevant for ICARUS)
   * it's a safer way.
   * In short: we test for each possible channel in which cryostat its optical
   * detector is, and extract ranges of channel numbers for each cryostat.
   * In the end, the result is pretty much known:
   * 0-179 on C:0 and 180-359 for C:1.
   */
  std::vector<lar::util::MinMaxCollector<OpDetID_t>> opDetChannelRangeByCryostat
    (fGeom->Ncryostats());
  lar::util::MinMaxCollector<OpDetID_t> opDetChannelRange;
  
  unsigned int const maxOpChannel = fGeom->MaxOpChannel();
  for (unsigned int channel = 0; channel < maxOpChannel; ++channel) {
    
    geo::OpDetGeo const& opDet = fGeom->OpDetGeoFromOpChannel(channel);
    auto const& opDetID = opDet.ID();
    
    if (opDetID.isValid) {
      opDetChannelRange.add(OpDetID_t(channel));
      opDetChannelRangeByCryostat[opDetID.Cryostat].add(OpDetID_t(channel));
    }
    
  } // for all channels
  
  // (2.2) determine the shift amount (hint: 180)
  fChannelShifts.clear();
  for (geo::CryostatGeo const& cryo: fGeom->IterateCryostats()) {
    
    // sanity checks
    geo::CryostatID const cid = cryo.ID();
    
    auto const& channelRange = opDetChannelRangeByCryostat[cid.Cryostat];
    if (!channelRange.has_data()) {
      throw std::runtime_error(
        "phot::ICARUSPhotonMappingTransformations::prepareMappings(): "
        + cid.toString() + " ends up with no optical channels??"
        );
    }
    
    auto const nChannels = channelRange.max() + 1 - channelRange.min();
    if ((unsigned int) nChannels != cryo.NOpDet()) {
      throw std::runtime_error(
        "phot::ICARUSPhotonMappingTransformations::prepareMappings(): "
        + cid.toString() + " expected to have "
        + std::to_string(cryo.NOpDet()) + " optical channels, we end up with "
        + std::to_string(nChannels) + " ("
        + std::to_string(channelRange.min()) + " - "
        + std::to_string(channelRange.max()) + ")"
        );
    } // if
    
    // relative shift with respect to the first cryostat:
    fChannelShifts.push_back
      (channelRange.min() - opDetChannelRangeByCryostat.front().min());
    
  } // for cryostats
  
  
  // (2.3) fill the maps with the shift
  
  // (2.3.1) std::vector<LibraryIndexToOpDetMap> fLibraryIndexToOpDetMaps
  //         mapping library index => optical detector ID
  //         indexed by cryostat number of the source
  // 
  // For ICARUS, this mapping is expected to be straightforward:
  //  * [ 0, 179 ] => [ 0, 179 ] if `location` is in the first cryostat (C:0)
  //  * [ 0, 179 ] => [ 180, 359 ] if `location` is in the other one (C:1)
  //
  fLibraryIndexToOpDetMaps.clear();
  for (geo::CryostatGeo const& cryo: fGeom->IterateCryostats()) {
    
    // this is the library for sources in cryostat `cryo`;
    // it is positioned at `cryo.ID().Cryostat` in `fLibOpDetIDmaps`;
    // the library is expected to have only 180 entries,
    // one per optical detector within `cryo`
    
    auto const nChannels = cryo.NOpDet();
    if (libraryIndices.size() != nChannels) {
      throw std::runtime_error(
        "phot::ICARUSPhotonMappingTransformations::prepareLibraryMappings(): "
        "the internal mapping should cover " + std::to_string(nChannels)
        + " channels but it covers " + std::to_string(libraryIndices.size())
        + " instead."
        );
    }
    
    auto const nFirst = fChannelShifts[cryo.ID().Cryostat];
    LibraryIndexToOpDetMap libraryIndexToOpDetMap { libraryIndices }; // copy
    for (auto& c: libraryIndexToOpDetMap) c += nFirst;
    
    fLibraryIndexToOpDetMaps.push_back(std::move(libraryIndexToOpDetMap));
    
  } // for channel range
  
  // if the position is invalid, no library data is actually present,
  // and the mapping is moot; we choose the rank of this mootness to be the same
  // as the first valid mapping
  fInvalidLibraryIndexToOpDetMap.clear();
  fInvalidLibraryIndexToOpDetMap.resize
    (fLibraryIndexToOpDetMaps.front().size(), InvalidOpDetID);
  
  // (2.3.2) std::vector<OpDetToLibraryIndexMap> fOpDetToLibraryIndexMaps
  //         mapping detector optical ID => library index
  //         
  // For ICARUS, this mapping is expected to be:
  //  * if `location` is in the first cryostat (C:0):
  //      * [ 0, 179 ]   => [ 0, 179 ]
  //      * [ 180, 359 ] => invalid indices
  //  * if `location` is in the second cryostat (C:1):
  //      * [ 0, 179 ]   => invalid indices
  //      * [ 180, 359 ] => [ 0, 179 ]
  // 
  fOpDetToLibraryIndexMaps.clear();
  for (auto const& libToOpDetMap: fLibraryIndexToOpDetMaps) {
    
    fOpDetToLibraryIndexMaps.push_back
      (invertMapping(libToOpDetMap, maxOpChannel, InvalidLibraryIndex));
    
  } // for each library mapping
  
  fInvalidOpDetToLibraryIndexMap.clear();
  fInvalidOpDetToLibraryIndexMap.resize(maxOpChannel, InvalidLibraryIndex);
  
} // phot::ICARUSPhotonMappingTransformations::prepareMappings()

  
//------------------------------------------------------------------------------
void phot::ICARUSPhotonMappingTransformations::prepareMappings
  (LibraryIndexToOpDetMap const& libraryIndices)
{
  
  /*
   * 1. geometry transformation:
   *    1. determine the switch coordinate (on x axis) between cryostats
   *    2. determine the amount of shift required
   * 2. library transformation
   *    1. determine the range of PMT channels in the two cryostats
   *    2. determine the shift amount (hint: 180)
   *    3. fill the maps with the shift
   */
  
  mf::LogInfo("ICARUSPhotonMappingTransformations")
    << "Photon visibility mapping tool: 'ICARUSPhotonMappingTransformations'";
    
  //
  // (1) geometry transformation
  //
  prepareGeometryMapping();
  
  //
  // (2) library transformation
  //
  prepareLibraryMappings(libraryIndices);
  
  // debug:
  if (fDumpMapping) dumpMapping();
  
} // phot::ICARUSPhotonMappingTransformations::prepareMappings()

  
//------------------------------------------------------------------------------
namespace {
  
  std::string channelIndex
    (phot::IPhotonMappingTransformations::OpDetID_t channel)
  {
    return (channel == phot::IPhotonMappingTransformations::InvalidOpDetID)
      ? "<invalid>": std::to_string(channel)
      ;
  } // channelIndex()
  
  std::string libraryIndex
    (phot::IPhotonMappingTransformations::LibraryIndex_t libIndex)
  {
    return
      (libIndex == phot::IPhotonMappingTransformations::InvalidLibraryIndex)
      ? "<invalid>": std::to_string(libIndex)
      ;
  } // libraryIndex()
  
} // local namespace


void phot::ICARUSPhotonMappingTransformations::dumpMapping() const {
  mf::LogInfo log("ICARUSPhotonMappingTransformations");
  
  log << "ICARUSPhotonMappingTransformations mapping";
  
  log << "\nMapping of geometry: '" << fGeom->DetectorName() << "':"
    << "\n  - " << fGeom->Ncryostats() << " cryostats"
    << "\n  - optical channels: " << fNOpDetChannels
    << "\n  - maximum optical detector channel number: " << fGeom->MaxOpChannel()
    << "\n"
    ;
  
  log << "\nThe switch point is at x=" << fSwitchPoint << " cm.";
  
  constexpr unsigned int PageBreak = 12;
  
  for (auto const& cryo: fGeom->IterateCryostats()) {
    
    auto const c = cryo.ID().Cryostat;
    
    log << "\n" << cryo.ID() << ":"
      << "\n  * optical detectors:        " << cryo.NOpDet()
      << "\n  * channel -> library shift: " << fChannelShifts[c]
      << "\n  * translation:              " << fTranslations[c];
    
    unsigned int pager;
    
    // print library to detector mapping
    auto const& libOpDetIDmap = fLibraryIndexToOpDetMaps[c];
    log
      << "\n  * mapping [library index @" << cryo.ID()
        << "] => [optical detector] (" << libOpDetIDmap.size()
        << " entries):";
      ;
    pager = PageBreak;
    for (auto [ i, opDetID ]: util::enumerate(libOpDetIDmap)) {
      if (++pager >= PageBreak) {
        pager = 0;
        log << "\n    ";
      }
      log << "  [" << i << "] => [" << channelIndex(opDetID) << "];";
    } // for
    
    // print detector to library mapping
    auto const& opDetLibMap = fOpDetToLibraryIndexMaps[c];
    pager = PageBreak;
    log
      << "\n  * mapping [optical detector] => [library index @"
        << cryo.ID() << "] (" << opDetLibMap.size() << " entries):";
      ;
    pager = PageBreak;
    for (std::size_t i = 0; i < opDetLibMap.size(); ++i) {
      if (++pager >= PageBreak) {
        pager = 0;
        log << "\n    ";
      }
      log << "  [" << i << "] => [" << libraryIndex(opDetLibMap[i]) << "];";
    }
    
  } // for cryostats
  
  {
    log << "\n"
      << "\nMapping [library index] => [no optical detector] ("
      << fInvalidLibraryIndexToOpDetMap.size() << " entries):";
    unsigned int pager = PageBreak;
    for (std::size_t i = 0; i < fInvalidLibraryIndexToOpDetMap.size(); ++i) {
      if (++pager >= PageBreak) {
        pager = 0;
        log << "\n    ";
      }
      log << "  [" << i << "] => ["
        << libraryIndex(fInvalidLibraryIndexToOpDetMap[i]) << "];";
    } // for
    log << "\n";
  } // anonymous scope block
  
  {
    log << "\n"
      << "\nMapping [optical detector] => [no library] ("
      << fInvalidOpDetToLibraryIndexMap.size() << " entries):";
    unsigned int pager = PageBreak;
    for (std::size_t i = 0; i < fInvalidOpDetToLibraryIndexMap.size(); ++i) {
      if (++pager >= PageBreak) {
        pager = 0;
        log << "\n    ";
      }
      log << "  [" << i << "] => ["
        << channelIndex(fInvalidOpDetToLibraryIndexMap[i]) << "];";
    } // for
    log << "\n";
  } // // anonymous scope block
  
} // phot::ICARUSPhotonMappingTransformations::dumpMapping()


//------------------------------------------------------------------------------

