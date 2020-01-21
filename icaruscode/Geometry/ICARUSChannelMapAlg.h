/**
 * @file   icaruscode/Geometry/ICARUSChannelMapAlg.h
 * @brief  Channel mapping algorithms for ICARUS detector.
 * @date   October 19, 2019
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    `icaruscode/Geometry/ICARUSChannelMapAlg.cxx`
 */

#ifndef ICARUSCODE_GEOMETRY_ICARUSCHANNELMAPALG_H
#define ICARUSCODE_GEOMETRY_ICARUSCHANNELMAPALG_H

// ICARUS libraries
#include "icaruscode/Geometry/details/ChannelToWireMap.h"
#include "icaruscode/Geometry/details/GeometryObjectCollections.h"

// LArSoft libraries
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/GeometryData.h"
#include "larcorealg/Geometry/GeometryDataContainers.h"
#include "larcorealg/Geometry/ReadoutDataContainers.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// framework libraries
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <vector>
#include <cassert>


// -----------------------------------------------------------------------------
// forward declarations
namespace icarus {
  
  class ICARUSChannelMapAlg;
  
} // namespace icarus


// -----------------------------------------------------------------------------
/**
 * @brief Channel mapping for ICARUS detector with split wires.
 * 
 * This channel mapping is designed for a detector description including two
 * TPC volumes in each drift volume.
 * 
 * Each of the four drift volumes in the ICARUS detector (one on each of the two
 * side of the cathode, on each of the two cryostats/modules) is described by
 * two "logical" TPC volumes. The division of the drift volume in two TPCs is
 * not physical, but it is a convenient description to accommodate the fact that
 * the first induction plane wires are split in two 9-meter halves which are
 * read independently.
 * 
 * In this scheme, the physical wires from the second induction and the
 * collection planes may be described by two `geo::WireGeo` objects associated
 * to the same readout channel, with the two `geo::WireGeo` belonging to a
 * different `geo::TPCGeo` and `geo::PlaneGeo` (but the planes are of the same
 * type and on the same view).
 * Therefore, physical wires on the second induction and collection planes are
 * associated with a single `geo::WireGeo` unless they cross the virtual border
 * between TPCs half way on the beam direction, in which case they are
 * associated with two `geo::WireGeo`. The first induction wires are always
 * associated with a single `geo::WireGeo`.
 * Each physical wire is associated and identified with a single TPC readout
 * channel.
 * 
 * The numbering of channels is as follows:
 * 
 * * the numbering is driven by physical wire planes
 * * plane `2` goes first, then planes `1` and `0` follow; plane `0` is the
 *   plane closest to the drift volume, that is the first induction plane;
 * * the internal order of the channel is defined by `geo::PlaneGeo` internal
 *   sorting algorithms, but it is guaranteed to be such that contiguous wires
 *   within the same logical wire plane have contiguous wire numbers, with no
 *   jumps and holes allowed;
 * * TPC `0` and `1`, which are "logical" TPCs describing the same drift volume,
 *   go first, then TPCs `2` and `3`
 * * cryostat `0` goes first, then cryostat `1`
 * 
 * This implies that the first wire plane being enumerated is `C:0 T:0 P:2`,
 * and then `C:0 T:1 P:2` to complete the collection plane of the first drift
 * volume, then the same thing with `C:0 T:0 P:1` and `C:0 T:1 P:1` for the
 * second induction plane, and finally `C:0 T:0 P:0` and `C:0 T:1 P:0` for the
 * first induction plane; then, `C:0 T:0 P:1` and `C:0 T:1 P:1` to complete the
 * second induction plane.
 * Thereafter, the enumeration proceeds with `C:0 T:2 P:2` (covering `C:0 T:2`
 * and `C:0 T:3`) and then with the other cryostat in the same TPC and plane
 * order.
 * 
 */
class icarus::ICARUSChannelMapAlg: public geo::ChannelMapAlg {
  
  // import definitions
  using TPCColl_t = icarus::details::TPCColl_t;
  using PlaneColl_t = icarus::details::PlaneColl_t;
  
    public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::OptionalDelegatedParameter Sorter {
      Name("Sorter"),
      Comment("configuration of the geometry object sorter")
      };
    
  }; // struct Config
  

  /// Constructor.
  ICARUSChannelMapAlg(Config const& config);
  
  /// Prepares the algorithm extracting information from the geometry.
  virtual void Initialize(geo::GeometryData_t const& geodata) override;
  
  
  /// Frees the resources of this algorithm.
  virtual void Uninitialize() override;
  
  
  // --- BEGIN -- Channel mapping ----------------------------------------------
  /// @name Channel mapping
  /// @{
  /**
   * @brief Returns a collection of ID of wires connected to the `channel`.
   * @param channel TPC readout channel number
   * @return collection of the wire IDs associated with `channel`
   * @throws cet::exception (category: "Geometry") if non-existent channel
   * 
   * If the TPC readout `channel` is invalid or non-existing, an exception is
   * thrown.
   * In ICARUS valid channels are expected to be associated with at least one
   * wire.
   */
  virtual std::vector<geo::WireID> ChannelToWire(raw::ChannelID_t channel) const
    override;
  
  /// Returns the number of readout channels (ID's go `0` to `Nchannels()`).
  virtual unsigned int Nchannels() const override;
  
  /// @brief Returns the number of channels in the specified ROP
  /// @return number of channels in the specified ROP, 0 if non-existent
  virtual unsigned int Nchannels(readout::ROPID const& ropid) const override;
  
  //@{
  virtual raw::ChannelID_t PlaneWireToChannel
    (geo::WireID const& wireID) const override;
  virtual raw::ChannelID_t PlaneWireToChannel(unsigned int plane,
                                              unsigned int wire,
                                              unsigned int tpc,
                                              unsigned int cstat) const override
    { return PlaneWireToChannel(geo::WireID(cstat, tpc, plane, wire)); }
  //@}

  /// @}
  // --- END -- Channel mapping ------------------------------------------------
  
  
  /**
   * @name Deprecated functions.
   *
   * These methods are legacy and might be replaced by `geo::GeometryCore`
   * calls.
   */
  /// @{
  //@{
  virtual double WireCoordinate
    (double YPos, double ZPos, geo::PlaneID const& planeID) const override;
  virtual double WireCoordinate(double YPos, double ZPos,
                               unsigned int PlaneNo,
                               unsigned int TPCNo,
                               unsigned int cstat) const override
    { return WireCoordinate(YPos, ZPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
  //@}
  
  //@{
  virtual geo::WireID NearestWireID
    (const TVector3& worldPos, geo::PlaneID const& planeID) const override;
  virtual geo::WireID NearestWireID(const TVector3& worldPos,
                               unsigned int    PlaneNo,
                               unsigned int    TPCNo,
                               unsigned int    cstat) const override
    { return NearestWireID(worldPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
  //@}
  
  virtual std::set<geo::PlaneID> const& PlaneIDs() const override;
  
  /// @}
  
  
  //
  // TPC set interface
  //
  /// @name TPC set mapping
  /// @{
  /**
   * @brief Returns the total number of TPC sets in the specified cryostat.
   * @param cryoid cryostat ID
   * @return number of TPC sets in the cryostat, or 0 if no cryostat found
   */
  virtual unsigned int NTPCsets
    (readout::CryostatID const& cryoid) const override;

  /// Returns the largest number of TPC sets any cryostat in the detector has.
  virtual unsigned int MaxTPCsets() const override;

  /// Returns whether we have the specified TPC set.
  /// @return whether the TPC set is valid and exists
  virtual bool HasTPCset(readout::TPCsetID const& tpcsetid) const override;

  /**
   * @brief Returns the ID of the TPC set the specified TPC belongs to
   * @param tpcid ID of the TPC
   * @return the ID of the corresponding TPC set, or invalid ID when tpcid is
   *
   * If the TPC ID is not valid, an invalid TPC set ID is returned.
   * Note that this check is performed on the validity of the TPC ID, that
   * does not necessarily imply that the TPC specified by the ID actually
   * exists.
   */
  virtual readout::TPCsetID TPCtoTPCset
    (geo::TPCID const& tpcid) const override;

  /**
   * @brief Returns a list of ID of TPCs belonging to the specified TPC set
   * @param tpcsetid ID of the TPC set to convert into TPC IDs
   * @return the list of TPCs, empty if TPC set is invalid
   *
   * In this mapping, TPC sets and TPCs are mapped one-to-one.
   * The returned list contains always one entry, unless the specified TPC
   * set ID is invalid, in which case the list is empty.
   * Note that the check is performed on the validity of the TPC set ID, that
   * does not necessarily imply that the TPC set specified by the ID actually
   * exists. Check the existence of the TPC set first (HasTPCset()).
   * Behaviour on valid, non-existent TPC set IDs is undefined.
   */
  virtual std::vector<geo::TPCID> TPCsetToTPCs
    (readout::TPCsetID const& tpcsetid) const override;

  /// Returns the ID of the first TPC belonging to the specified TPC set
  virtual geo::TPCID FirstTPCinTPCset
    (readout::TPCsetID const& tpcsetid) const override;

  /// @} TPC set mapping


  
  // --- BEGIN -- Readout plane interface --------------------------------------
  /// @name Readout plane mapping
  /// @{

  /**
   * @brief Returns the total number of readout planes in the specified TPC set.
   * @param tpcsetid TPC set ID
   * @return number of readout planes in the TPC sets, or `0` if ID is invalid
   *
   * The validity of the requested TPC set is performed like in `HasTPCset()`.
   */
  virtual unsigned int NROPs
    (readout::TPCsetID const& tpcsetid) const override;

  /// Returns the largest number of ROPs a TPC set in the detector has.
  virtual unsigned int MaxROPs() const override;

  /// Returns whether we have the specified ROP
  /// @return whether the readout plane is valid and exists
  virtual bool HasROP(readout::ROPID const& ropid) const override;

  /**
   * @brief Returns the ID of the ROP planeid belongs to, or invalid if none.
   * @param planeid ID of the plane
   * @return the ID of the corresponding ROP, or invalid ID when `planeid` is
   *
   * If the plane ID is not valid, an invalid readout plane ID is returned.
   * Note that this check is performed on the validity of the plane ID, that
   * does not necessarily imply that the plane specified by the ID actually
   * exists.
   */
  virtual readout::ROPID WirePlaneToROP
    (geo::PlaneID const& planeid) const override;

  /**
   * @brief Returns a list of ID of wire planes belonging to the specified ROP.
   * @param ropid ID of the readout plane to convert into wire planes
   * @return the list of wire plane IDs, empty if readout plane ID is invalid
   *
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists.
   * 
   * In this mapping, readout planes contain one or two wire planes each,
   * depending on the view.
   */
  virtual std::vector<geo::PlaneID> ROPtoWirePlanes
    (readout::ROPID const& ropid) const override;

  /**
   * @brief Returns a list of ID of TPCs the specified ROP spans
   * @param ropid ID of the readout plane
   * @return the list of TPC IDs, empty if readout plane ID is invalid
   *
   * In this mapping, readout planes and wire planes are mapped one-to-one.
   * The returned list contains always one entry, unless the specified readout
   * plane ID is invalid, in which case the list is empty.
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists. Check if the ROP exists with HasROP().
   * The behaviour on non-existing readout planes is undefined.
   */
  virtual std::vector<geo::TPCID> ROPtoTPCs
    (readout::ROPID const& ropid) const override;

  /// Returns the ID of the ROP the channel belongs to (invalid if none).
  virtual readout::ROPID ChannelToROP
    (raw::ChannelID_t channel) const override;

  /**
   * @brief Returns the ID of the first channel in the specified readout plane.
   * @param ropid ID of the readout plane
   * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
   *
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists. Check if the ROP exists with HasROP().
   * The behaviour for non-existing readout planes is undefined.
   */
  virtual raw::ChannelID_t FirstChannelInROP
    (readout::ROPID const& ropid) const override;

  /**
   * @brief Returns the ID of the first plane belonging to the specified ROP.
   * 
   * The wire planes within a readout plane are supposed to be sorted by beam
   * (_z_) coordinate, so that the first one should be the closest to the
   * beam entrance point.
   */
  virtual geo::PlaneID FirstWirePlaneInROP
    (readout::ROPID const& ropid) const override;

  /// @}
  // --- END -- Readout plane interface ----------------------------------------

  
  /// Return the sorter.
  virtual geo::GeoObjectSorter const& Sorter() const override
    { return fSorter; }

    private:
  
  // --- BEGIN -- Readout element information ----------------------------------
  /**
   * @name Readout element information
   * 
   * The geometry and readout data containers have currently no support for
   * resizing and their size is assigned on construction.
   * 
   * Access should happen via the corresponding member functions.
   * 
   */
  /// @{
  
  using ChannelRange_t = icarus::details::ChannelRange_t; // import type
  
  /// Collected information about TPC sets and readout planes in the geometry.
  struct ReadoutMappingInfo_t {
    /// Number of TPC sets in each cryostat.
    std::vector<unsigned int> fTPCsetCount;
    
    /// All `geo::TPCGeo` objects in each TPC set, sorted by increasing _z_.
    readout::TPCsetDataContainer<TPCColl_t> fTPCsetTPCs;
    
    /// Number of readout planes in each TPC set.
    readout::TPCsetDataContainer<unsigned int> fROPcount;
    
    /// All `geo::PlaneGeo` objects in each readout plane, sorted by _z_.
    readout::ROPDataContainer<PlaneColl_t> fROPplanes;
    
    /// The TPC set each TPC belongs to.
    geo::TPCDataContainer<readout::TPCsetID> fTPCtoTPCset;
    
    /// The ROP each wire plane belongs to.
    geo::PlaneDataContainer<readout::ROPID> fPlaneToROP;
    
    ReadoutMappingInfo_t() = default;
    
    void set(
      std::vector<unsigned int>&& TPCsetCount,
      readout::TPCsetDataContainer<TPCColl_t>&& TPCsetTPCs,
      readout::TPCsetDataContainer<unsigned int>&& ROPcount,
      readout::ROPDataContainer<PlaneColl_t>&& ROPplanes,
      geo::TPCDataContainer<readout::TPCsetID>&& TPCtoTPCset,
      geo::PlaneDataContainer<readout::ROPID>&& PlaneToROP
      )
      {
        fTPCsetCount = std::move(TPCsetCount);
        fTPCsetTPCs  = std::move(TPCsetTPCs );
        fROPcount    = std::move(ROPcount   );
        fROPplanes   = std::move(ROPplanes  );
        fTPCtoTPCset = std::move(TPCtoTPCset);
        fPlaneToROP  = std::move(PlaneToROP );
        assert(fTPCsetCount.size() == fTPCsetTPCs.dimSize<0U>());
        assert(fTPCsetCount.size() == fROPcount.dimSize<0U>());
        assert(fTPCsetCount.size() == fROPplanes.dimSize<0U>());
        assert(fTPCsetCount.size() == fTPCtoTPCset.dimSize<0U>());
        assert(fTPCsetCount.size() == fPlaneToROP.dimSize<0U>());
        assert(fTPCsetTPCs.dimSize<1U>() == fROPcount.dimSize<1U>());
        assert(fTPCsetTPCs.dimSize<1U>() == fROPplanes.dimSize<1U>());
        assert(fTPCtoTPCset.dimSize<1U>() == fPlaneToROP.dimSize<1U>());
      } // set()
    
    unsigned int NCryostats() const
      { return fROPplanes.dimSize<0U>(); }
    unsigned int MaxTPCsets() const { return fROPplanes.dimSize<1U>(); }
    unsigned int MaxROPs() const { return fROPplanes.dimSize<2U>(); }
    
    /// Frees the memory and leaves the object unusable until next `set()`.
    void clear()
      {
        fTPCsetCount.clear(); fTPCsetTPCs.clear();
        fROPcount.clear(); fROPplanes.clear();
        fTPCtoTPCset.clear(); fPlaneToROP.clear();
      }
    
    /// Returns whether all the data containers are initialized.
    operator bool() const
      {
        return !fTPCsetCount.empty() && !fTPCsetTPCs.empty()
          && !fROPcount.empty() && !fROPplanes.empty()
          && !fTPCtoTPCset.empty() && !fPlaneToROP.empty();
      }
    
  }; // ReadoutMappingInfo_t
  
  /// Collection of information on one plane.
  struct PlaneInfo_t {
    
    ChannelRange_t fChannelRange; ///< Range of channels covered by the plane.
    readout::ROPID fROPID; ///< Which readout plane this wire plane belongs to.
    
    /// Returns the range of channels covered by the wire plane.
    constexpr ChannelRange_t const& channelRange() const
      { return fChannelRange; }
    
    /// Returns the ID of the last channel in the range.
    constexpr raw::ChannelID_t firstChannel() const
      { return fChannelRange.begin(); }
    
    /// Returns the ID of the last channel in the range.
    constexpr raw::ChannelID_t lastChannel() const
      { return fChannelRange.end() - 1; }
    
    /// Returns the ID of the channel after the last in the range.
    constexpr raw::ChannelID_t endChannel() const
      { return fChannelRange.end(); }
    
    /// Returns the ID of the readout plane this wire plane belongs to.
    constexpr readout::ROPID ROP() const { return fROPID; }
    
  }; // struct PlaneInfo_t
  
  
  /// Information about TPC sets and readout planes in the geometry.
  ReadoutMappingInfo_t fReadoutMapInfo;
  
  /// Mapping of channels to wire planes and ROP's.
  icarus::details::ChannelToWireMap fChannelToWireMap;
  
  /// Range of channels covered by each of the wire planes.
  geo::PlaneDataContainer<PlaneInfo_t> fPlaneInfo;
  
  
  /// @}
  // --- END -- Readout element information ------------------------------------
  
  
  // --- BEGIN -- Sorting ------------------------------------------------------
  /// Algorithms to sort geometry elements.
  geo::GeoObjectSorterStandard fSorter;
  
  // --- END -- Sorting --------------------------------------------------------
  
  
  // --- BEGIN -- Readout element information access ---------------------------
  /// @name Readout element information access
  /// @{
  
  /// Returns the number of TPC sets in each cryostat.
  std::vector<unsigned int> const& TPCsetCount() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fTPCsetCount; }
  
  /// Returns the number of TPC sets in the specified cryostat `cid`.
  unsigned int TPCsetCount(readout::CryostatID const& cid) const
    { return TPCsetCount()[cid.Cryostat]; }
  
  /// All `geo::TPCGeo` objects in each TPC set, sorted by increasing _z_.
  readout::TPCsetDataContainer<TPCColl_t> const& TPCsetTPCs() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fTPCsetTPCs; }
  
  /// All `geo::TPCGeo` objects in the specified TPC set `sid`.
  TPCColl_t const& TPCsetTPCs(readout::TPCsetID const& sid) const
    { return TPCsetTPCs()[sid]; }
  
  /// Number of readout planes in each TPC set.
  readout::TPCsetDataContainer<unsigned int> const& ROPcount() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fROPcount; }
  
  /// Number of readout planes in the specified TPC set `sid`.
  unsigned int ROPcount(readout::TPCsetID const& sid) const
    { return ROPcount()[sid]; }
  
  /// All `geo::PlaneGeo` objects in each readout plane, sorted by _z_.
  readout::ROPDataContainer<PlaneColl_t> const& ROPplanes() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fROPplanes; }
  
  /// All `geo::PlaneGeo` objects in the specified readout plane `rid`.
  PlaneColl_t const& ROPplanes(readout::ROPID const& rid) const
    { return ROPplanes()[rid]; }
  
  /// The TPC set including each TPC.
  geo::TPCDataContainer<readout::TPCsetID> const& TPCtoTPCset() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fTPCtoTPCset; }
  
  /* // something similar to this already belongs to the interface
  /// The TPC set the specified TPC `tid` belongs to.
  readout::TPCsetID const& TPCtoTPCset(geo::TPCID const& tid) const
    { return TPCtoTPCset()[tid]; }
  */
  
  /// The readout plane including each wire plane.
  geo::PlaneDataContainer<readout::ROPID> const& PlaneToROP() const
    { assert(fReadoutMapInfo); return fReadoutMapInfo.fPlaneToROP; }
  
  /// The readout plane the specified wire plane `pid` belongs to.
  readout::ROPID const& PlaneToROP(geo::PlaneID const& pid) const
    { return PlaneToROP()[pid]; }
  
  /// @}
  // --- END -- Readout element information access -----------------------------
  
  
  /// Returns whether the specified cryostat is known to the mapping.
  bool HasCryostat(readout::CryostatID const& cryoid) const;
  
  
  /**
   * @brief Fills the information about readout channel mapping.
   * @param Cryostats the sorted list of cryostats in the detector
   * 
   * 
   * 
   * The readout information must have been already filled
   * (`buildReadoutPlanes()`).
   * 
   */
  void fillChannelToWireMap
    (geo::GeometryData_t::CryostatList_t const& Cryostats);
  
  
  /**
   * @brief Fills information about the TPC set and readout plane structure.
   * @param Cryostats the sorted list of cryostats in the detector
   * 
   * This method extracts and fills the following information:
   * * the number of TPC sets in each cryostat (`fTPCsetCount`);
   * * all `geo::TPCGeo` objects in each TPC set, sorted by increasing _z_
   *     (`fTPCsetTPCs`);
   * * the number of readout planes in each TPC set (`fROPcount`);
   * * all `geo::PlaneGeo` objects in each readout plane, sorted by increasing
   *     _z_ (`fROPplanes`).
   * 
   * Cryostats and its components are expected to be already in the final order
   * and with all their ID's set.
   * 
   */
  void buildReadoutPlanes(geo::GeometryData_t::CryostatList_t const& Cryostats);
  
  
  /// Returns the type of signal on the specified `channel`.
  virtual geo::SigType_t SignalTypeForChannelImpl
    (raw::ChannelID_t const channel) const override;
  
  
}; // class icarus::ICARUSChannelMapAlg


#endif // ICARUSCODE_GEOMETRY_ICARUSCHANNELMAPALG_H

