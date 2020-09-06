/**
 * @file   icaruscode/PMT/LibraryMappingTools/ICARUSPhotonMappingTransformations.h
 * @brief  Photon library mapping for ICARUS geometry.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 3, 2019
 * @see    `icaruscode/PMT/LibraryMappingTools/ICARUSPhotonMappingTransformations_tool.cc`
 * 
 */

#ifndef ICARUSCODE_LIGHT_LIBRARYMAPPINGTOOLS_ICARUSPHOTONMAPPINGTRANSFORMATIONS_H
#define ICARUSCODE_LIGHT_LIBRARYMAPPINGTOOLS_ICARUSPHOTONMAPPINGTRANSFORMATIONS_H

// LArSoft libraries
#include "larsim/PhotonPropagation/LibraryMappingTools/IPhotonMappingTransformations.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Vector_t

// framework libraries
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>


namespace phot {
  
  /**
   * @brief Photon library mapping for ICARUS geometry.
   * 
   * This is an implementation of `phot::IPhotonMappingTransformation` interface
   * to exploit ICARUS detector symmetries.
   * 
   * The current implementation only exploits the fact that the two cryostats
   * are identical, and ignores the symmetry respect to the cathode within each
   * of the cryostats.
   * 
   * The required library is expected to cover the first 180 PMT channels (the
   * ones pertaining the first cryostat, `C:0`), and the full volume of the TPC.
   * 
   * When requested a point in the first cryostat (`C:0`), the point is used
   * directly in the visibility library to get visibility on PMT channels 0 to
   * 179 (the lower half). The mapping then maps the remaining PMT channels
   * (180 to 359) to be invalid and with default value of `0`.
   * 
   * When requested a point in the other cryostat (`C:1`), the point is
   * translated into the first one to get visibility of PMT channels 0 to 179
   * (the lower half). The mapping then translates those channels into the range
   * 180 to 359, and maps the remaining PMT channels (0 to 179) to be invalid
   * and with default value of `0`.
   * 
   * Note that the content of the library for channels 180 to 359 is always
   * ignored, and it may well be absent.
   */
  class ICARUSPhotonMappingTransformations
    : public phot::IPhotonMappingTransformations
  {
    
      public:
    struct Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::OptionalSequence<OpDetID_t> CryostatChannelRemap {
        Name("CryostatChannelRemap"),
        Comment("internal library mapping (new library ID for each old one)")
        };
      
      // no configuration required
      fhicl::Atom<bool> DumpMapping{
        Name("DumpMapping"),
        Comment("dump the mapping into console"),
        false // default
        };
      
    }; // struct Config
    
    using Parameters = art::ToolConfigTable<Config>;
    
    
    /// Constructor.
    ICARUSPhotonMappingTransformations(Config const& config);
    
    /// Constructor: ignores the configuration.
    ICARUSPhotonMappingTransformations(Parameters const& config)
      : ICARUSPhotonMappingTransformations(config())
      {}
    
    
    // --- BEGIN Geometry mapping interface ------------------------------------
    /// @name Geometry mapping interface
    /// @{
    
    /**
     * @brief Returns the representation within the library of a detector
     *        location.
     * @param location position in world coordinates [cm]
     * @return a vector expressing `location` in the library space
     * 
     * The returned vector is shifted as it were on the first cryostat.
     * 
     * No exception is ever thrown.
     */
    virtual geo::Point_t detectorToLibrary
      (geo::Point_t const& location) const override;
    
    /// @}
    // --- END Geometry mapping interface --------------------------------------
    
    
    // --- BEGIN Optical detector mapping interface ----------------------------
    /// @name Optical detector mapping interface
    /// @{
    
    /**
     * @brief Returns the library index for the specified optical detector.
     * @param location where the scintillation source is in world frame [cm]
     * @param opDetID optical detector identifier
     *              (as used, e.g., in `geo::GeometryCore::OpDetGeoFromOpDet()`)
     * @return index corresponding to `opDetID`, or `InvalidLibraryIndex`
     * @throw std::out_of_range if input optical detector ID can't be handled TODO not implemented this way (yet?)
     * @see `opDetFromLibrary()`
     * 
     * The mapping reduces the optical detector ID to the corresponding one in
     * the first cryostat.
     */
    virtual LibraryIndex_t opDetToLibraryIndex
      (geo::Point_t const& location, OpDetID_t opDetID) const override
      { return opDetsToLibraryIndicesImpl(location)[opDetID]; }
    
    /**
     * @brief Returns the optical detector ID for the specified library index.
     * @param location where the scintillation source is in world frame [cm]
     * @param libIndex library index to be mapped
     * @return optical detector corresponding to `libIndex` (as used, e.g., in
     *         `geo::GeometryCore::OpDetGeoFromOpDet()`), or `InvalidOpDetID`
     * @throw std::out_of_range if input library index can't be handled TODO not implemented this way (yet?)
     * @see `opDetToLibraryIndex()`, `libraryIndicesToOpDets()`
     */
    virtual OpDetID_t libraryIndexToOpDet
      (geo::Point_t const& location, LibraryIndex_t libIndex) const override
      { return libraryIndicesToOpDetsImpl(location)[libIndex]; }
    
    
    /**
     * @brief Returns a map of library indices as function of optical detectors.
     * @param location where the scintillation source is in world frame [cm]
     * @return library indices for all optical detectors
     * @see `opDetToLibraryIndex()`, `opDetsToLibraryIndices()`
     */
    virtual OpDetToLibraryIndexMap const& opDetsToLibraryIndices
      (geo::Point_t const& location) const override
      { return opDetsToLibraryIndicesImpl(location); }
    
    /**
     * @brief Expected number of mappings of optical detector into library
     *        index.
     * @return the expected size of the mapping of optical detectors
     * @see `opDetsToLibraryIndices()`
     * 
     * This is effectively the number of available optical detectors, as well
     * as the size of the mapping as returned by `opDetsToLibraryIndices()`.
     */
    virtual std::size_t opDetMappingSize() const override
      { return opDetsToLibraryIndicesImpl(geo::origin()).size(); }
    
    
    /**
     * @brief Returns a map of identifiers of optical detectors for each library
     *        index, for the library appropriate around `location`
     * @param location where the scintillation source is in world frame [cm]
     * @return optical detector identifiers for all library indices
     * @see `opDetsToLibraryIndices()`, `libraryIndexToOpDet()`
     * 
     * We use this mapping when we have information from the visibility library
     * appropriate for a source at the specified `location`: this information
     * is addressed via a "library index" with the reduced range (i.e. 0 to 179)
     * and we need to know which are the optical detector channels that
     * visibility covers.
     * 
     * For ICARUS, this mapping is expected to be straightforward:
     *  * [ 0, 179 ] => [ 0, 179 ] if `location` is in the first cryostat (C:0)
     *  * [ 0, 179 ] => [ 180, 359 ] if `location` is in the other one (C:1)
     * 
     * The specified `location` is used to provide context in a similar
     * fashion as `detectorToLibrary()` does. It can be used to choose the
     * correct mapping among the available ones.
     * 
     * The returned value is a mapping object (see `LibraryIndexToOpDetMap`
     * documentation for the interface). If a library index does not map to any
     * optical detector in the library at `location` (which is unusual!),
     * the optical detector ID corresponding to it is `InvalidOpDetID`.
     */
    virtual LibraryIndexToOpDetMap const& libraryIndicesToOpDets
      (geo::Point_t const& location) const override
      { return libraryIndicesToOpDetsImpl(location); }
    
    /**
     * @brief Expected size of the mapping from library to optical detectors.
     * @param location where the scintillation source is in world frame [cm]
     * @return the expected size of the mapping from library indices to
     *         optical detectors
     * @see `libraryIndicesToOpDets()`
     * 
     * This is effectively the size of the mapping returned by
     * `libraryIndicesToOpDets()`. It represents how many library indices are
     * provided by the library for the specified `location`, that is `180`.
     */
    virtual std::size_t libraryMappingSize
      (geo::Point_t const& location) const override
      { return libraryIndicesToOpDets(location).size();}
    
    
    /// @}
    // --- END Optical detector identifier mapping interface -------------------
    
    
      protected:
    //
    // configuration parameters
    //
    bool fDumpMapping = false; ///< Whether to dump mapping on initialization.
    
    //
    // tool setup
    //
    
    /// Detector geometry service provider. Not really used.
    geo::GeometryCore const* fGeom = nullptr;
    
    
    //
    // for geometry transformation
    //
    std::vector<geo::Vector_t> fTranslations; ///< Translation of the point.
    geo::Length_t fSwitchPoint; ///< Switch coordinate on x axis [cm]
    
    //
    // for optical channel transformation
    //
    unsigned int fNOpDetChannels; /// Total number of optical detector channels.
    /// Amount of channel number shifting indexed by cryostat. Not really used.
    std::vector<int> fChannelShifts;
    /// Library to detector channel mappings, indexed by cryostat number.
    std::vector<LibraryIndexToOpDetMap> fLibraryIndexToOpDetMaps;
    /// A library-to-detector mapping for invalid points.
    LibraryIndexToOpDetMap fInvalidLibraryIndexToOpDetMap; // currently not used
    /// Detector channel to library mappings, indexed by cryostat number.
    std::vector<OpDetToLibraryIndexMap> fOpDetToLibraryIndexMaps;
    /// A detector-to-library mapping for invalid points.
    OpDetToLibraryIndexMap fInvalidOpDetToLibraryIndexMap; // currently not used
    
    //
    // mapping helper functions
    //
    
    /// Returns which cryostat better contain `point`. Never invalid so far.
    geo::CryostatID whichCryostat(geo::Point_t const& point) const
      { return geo::CryostatID{ (point.X() > fSwitchPoint)? 1U: 0U }; }
    
    //
    // implementation methods (non-virtual)
    //
    OpDetToLibraryIndexMap const& opDetsToLibraryIndicesImpl
      (geo::Point_t const& location) const;
      
    LibraryIndexToOpDetMap const& libraryIndicesToOpDetsImpl
      (geo::Point_t const& location) const;


    void prepareGeometryMapping();
    void prepareLibraryMappings(LibraryIndexToOpDetMap const& libraryIndices);
    
    /// Extracts the necessary information for mapping from the geometry.
    void prepareMappings(LibraryIndexToOpDetMap const& libraryIndices);
    
    /// Writes the current mapping information into the console. Debug stuff.
    void dumpMapping() const;
    
    /**
     * @brief Inverts a given mapping.
     * @tparam OutputIndex the output index of the mapping to invert
     * @tparam IndexIndex the input index of the mapping to invert
     * @tparam Container container used for the mapping (be it `std::vector`)
     * @param directMap the mapping (input index to output index) to be inverted
     * @param size the total number of output indices (some may be not mapped)
     * @param invalidIndex input index to associate to unmapped output indices
     * @return a mapping (output index to input index) reversing `directMap`
     * 
     * The resulting, "inverse" mapping has exactly `size` entries.
     * If the input mapping maps out beyond that, those indices are not mapped
     * back.
     */
    template<
      typename OutputIndex, typename InputIndex,
      template <typename...> typename Container
      >
    static Container<InputIndex> invertMapping(
      Container<OutputIndex> const& directMap,
      std::size_t size,
      InputIndex invalidIndex
      );
    
  }; // class ICARUSPhotonMappingTransformations
  
  
} // namespace phot


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
/**
 * @brief Inverts a given mapping.
 * @tparam OutputIndex the output index of the mapping to invert
 * @tparam IndexIndex the input index of the mapping to invert
 * @tparam Container type of container used for the mapping (be it `std::vector`)
 * @param directMap the mapping (input index to output index) to be inverted
 * @param size the total number of output indices (some may be not mapped)
 * @param invalidIndex the input index to associate to unmapped output indices
 * @return a mapping (output index to input index) reversing `directMap`
 * 
 * The resulting, "inverse" mapping has exactly `size` entries.
 * If the input mapping maps out beyond that, those indices are not mapped back.
 */
template
  <typename OutputIndex, typename InputIndex, template <typename...> typename Container>
auto phot::ICARUSPhotonMappingTransformations::invertMapping (
  Container<OutputIndex> const& directMap,
  std::size_t size,
  InputIndex invalidIndex
) -> Container<InputIndex> {
  
  OutputIndex const endOutputIndex = size;
  
  Container<InputIndex> inverseMap(size, invalidIndex);
  for (std::size_t i = 0U; i < directMap.size(); ++i) {
    OutputIndex o = directMap[i];
    if (o < endOutputIndex) inverseMap[o] = i;
  }
  return inverseMap;
}


//------------------------------------------------------------------------------


#endif // ICARUSCODE_LIGHT_LIBRARYMAPPINGTOOLS_ICARUSPHOTONMAPPINGTRANSFORMATIONS_H
