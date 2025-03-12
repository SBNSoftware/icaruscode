/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapProviderBase.h
 * @author T. Usher (factorised by G. Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapProviderBase.cxx
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPROVIDERBASE_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPROVIDERBASE_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMapProvider.h"
#include "icarusalg/Utilities/mfLoggingClass.h"

// framework libraries
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <string>
#include <memory> // std::unique_ptr<>


// -----------------------------------------------------------------------------
namespace icarusDB {
  template <typename ChMapAlg> class ICARUSChannelMapProviderBase;
}
/**
 * @brief Base implementation of ICARUS channel mapping database interface.
 * @tparam ChMapAlg type of channel mapping helper to be used
 * 
 * 
 * Caches
 * -------
 * 
 * This implementation relies on a database backend reading the full information
 * from the database and caching the result, every time a new run period is
 * requested. To help users who in turn cache information from this object to
 * track whether _their_ caches are invalidated by such occurrences, this object
 * uses the `util::CacheCounter` facility to version every cache.
 * Users will want to use the `util::CacheGuard` tool to monitor the cache.
 * Three caches are tracked, with tags `"TPC"`, `"PMT"` and `"CRT"`, plus
 * the overall cache tracking (tagged with an empty name `""`).
 * At the moment of writing, the three caches are actually updated all at the
 * same times.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * * `ChannelMappingTool` (algorithm configuration): see the configuration
 *     of the database backend access helper (`ChMapAlg`)
 * * `DiagnosticOutput` (flag, default: `false`): verbosely load parameters
 *     and in general prints more messages to console
 * * `LogCategory` (string, default: same as `ChannelMappingTool.LogCategory`):
 *     name of the messagefacility category used to send messages to console,
 *     useful for filtering messages.
 * 
 * 
 * Implementation details
 * =======================
 * 
 * The only reason for the existence of this base class is that due to their
 * original form the two service providers implementing access to SQLite and to
 * PostgreSQL database share most of the code.
 * So that code is collected in a base class.
 * 
 * In fact, the main difference between the previous factor form and the current
 * one is that then the database backend was selected via _art_ tool within a
 * single provider, while now each provider has its backend statically coded in
 * and the choice happens at service level (the service level choice was dummy
 * since there was a single service and a single provider).
 * 
 * The service provider interface is inherited with virtual inheritance, because
 * it is foreseen that the service interface also inherits from the service
 * provider interface, and virtual inheritance guarantees that there is only
 * one underlying interface class (and ultimately it is this one which provides
 * the actual implementation of the interface).
 * 
 */
template <typename ChMapAlg>
class icarusDB::ICARUSChannelMapProviderBase
  : virtual public IICARUSChannelMapProvider
  , private icarus::ns::util::mfLoggingClass
{
  
    public:
  /// The helper used to extract data from the database.
  using ChannelMappingAlg_t = ChMapAlg;
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Table<typename ChannelMappingAlg_t::Config> ChannelMappingTool {
      Name{ "ChannelMappingTool" },
      Comment{ "parameters for the channel mapping algorithm" }
      };
    
    fhicl::Atom<bool> DiagnosticOutput{
      Name{ "DiagnosticOutput" },
      Comment{ "enables verbose output to console" },
      false // default
      };
    
    fhicl::OptionalAtom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "name of the console stream to send messages to" }
      };
    
  }; // Config
  
  using Parameters = fhicl::Table<Config>;
  
  
  /// Constructor: configures the service provider via FHiCL configuration.
  ICARUSChannelMapProviderBase(Config const& config);
  
  /// Constructor: configures the service provider via FHiCL configuration.
  ICARUSChannelMapProviderBase(Parameters const& params)
    : ICARUSChannelMapProviderBase{ params() }
    {}
  
  
  /// --- BEGIN --- Data period selection --------------------------------------
  /// @name Data period selection
  /// @{
  
  /// Loads the mapping for `run`, returns whether a new mapping was loaded.
  virtual bool forRun(int run) override;
  
  /// Loads the mapping for `period`, returns whether a new mapping was loaded.
  virtual bool forPeriod(icarusDB::RunPeriod period) override;
  
  /// @}
  /// --- END ----- Data period selection --------------------------------------
  
  
  /// --- BEGIN --- TPC information --------------------------------------------
  /// @name TPC information
  /// @{
  
  /// --- BEGIN - - TPC fragment information - - - - - - - - - - - - - - - - -
  /// @name TPC fragment information
  /// @{
  
  /// Returns whether the specified `ID` is a known TPC fragment ID.
  virtual bool hasFragmentID(unsigned int ID) const override;
  /// Returns the number of TPC fragment IDs known to the service.

  /// Returns the number of known TPC fragments.
  virtual unsigned int nTPCfragmentIDs() const override;
  
  /// Returns the name of the crate served by the specified `fragmentID`.
  virtual std::string const& getCrateName (unsigned int fragmentID) const
    override;
  
  /// Returns the list of board IDs included in the specified `fragmentID`.
  virtual ReadoutIDVec const& getReadoutBoardVec(unsigned int fragmentID) const
    override;
  
  /// @}
  /// --- END - - - TPC fragment information - - - - - - - - - - - - - - - - -
  
  /**
   * @brief Returns the full TPC channel mapping.
   * 
   * The returned mapping is an associative container, associating each board
   * ID (`boardID`) to the following information (in a `std::tuple`):
   *  * `[0]` slot number the board is on (like `getBoardSlot(boardID)`)
   *  * `[1]` list of channels served by the board, and their plane
   *          (like `getChannelPlanePair(boardID)`)
   * 
   */
  virtual TPCReadoutBoardToChannelMap const& getReadoutBoardToChannelMap() const
    override;
  
  /// --- BEGIN - - TPC board information  - - - - - - - - - - - - - - - - - -
  /// @name TPC board information
  /// @{
  
  /// Returns whether there is a board with the specified `ID`.
  virtual bool hasBoardID(unsigned int ID) const override;
  
  /// Returns the number of TPC board IDs known to the service.
  virtual unsigned int nTPCboardIDs() const override;
  
  /// Returns the number of slot the `boardID` is on.
  virtual unsigned int getBoardSlot(unsigned int boardID) const override;
  
  /// Returns a list of channels served by the `boardID` and for each the plane
  /// it is on (`0` to `2`).
  virtual ChannelPlanePairVec const& getChannelPlanePair
    (unsigned int boardID) const override;

  /// @}
  /// --- END - - - TPC board information  - - - - - - - - - - - - - - - - - -
  
  
  /// --- BEGIN --- PMT information --------------------------------------------
  /// @name PMT information
  /// @{
  
  /// Returns whether the specified fragment `ID` is known to the mapping.
  virtual bool hasPMTDigitizerID(unsigned int ID) const override;
  
  /// Returns the number of PMT fragment IDs known to the mapping.
  virtual unsigned int nPMTfragmentIDs() const override;
  
  /// Returns records on all the PMT channels covered by the fragment `ID`.
  virtual PMTdigitizerInfoVec const& getPMTchannelInfo (unsigned int ID) const
    override;
  
  /// @}
  /// --- END ----- PMT information --------------------------------------------


  /// --- BEGIN --- CRT information --------------------------------------------
  /// @name CRT information
  /// @{
  
  /// Returns the sim Mac address corresponding to the specified side CRT
  /// hardware address.
  virtual unsigned int getSimMacAddress(unsigned int hwmacaddress) const
    override;
  
  /// Returns the sim Mac address corresponding to the specified top CRT
  /// hardware address.
  virtual unsigned int gettopSimMacAddress(unsigned int) const override;

  /// Returns the Gain and Pedestal for Side CRT.
  virtual std::pair<double, double> getSideCRTCalibrationMap
    (int mac5, int chan) const override;
  
  /// @}
  /// --- END ----- CRT information --------------------------------------------

  /// Returns the channel mapping database key for the specified PMT fragment ID.
  static constexpr unsigned int PMTfragmentIDtoDBkey(unsigned int fragmentID);
  
  /// Returns the PMT fragment ID for the specified channel mapping database
  /// key.
  static constexpr unsigned int DBkeyToPMTfragmentID(unsigned int DBkey);


    protected:
  
  // --- BEGIN --- Configuration parameters ------------------------------------
      
  bool const fDiagnosticOutput;
  
  std::string const fLogCategory;
  
  // --- END ----- Configuration parameters ------------------------------------
  
  
  /// The helper extracting information from the database.
  ChannelMappingAlg_t fChannelMappingAlg;
  
  
  // --- BEGIN --- Cache -------------------------------------------------------
  
  icarusDB::TPCFragmentIDToReadoutIDMap fTPCFragmentToReadoutMap;
  
  icarusDB::TPCReadoutBoardToChannelMap fTPCReadoutBoardToChannelMap;

  icarusDB::PMTFragmentToDigitizerChannelMap fPMTFragmentToDigitizerMap; 

  icarusDB::CRTChannelIDToHWtoSimMacAddressPairMap
    fCRTChannelIDToHWtoSimMacAddressPairMap;

  icarusDB::TopCRTHWtoSimMacAddressPairMap fTopCRTHWtoSimMacAddressPairMap;

  icarusDB::SideCRTChannelToCalibrationMap fSideCRTChannelToCalibrationMap;
  
  // --- END ----- Cache -------------------------------------------------------
  

  /// Has the channel mapping tool fill the mapping caches.
  void readFromDatabase();
  
  /// Returns the list of records of all channels in the PMT readout board with
  /// the specified fragment.
  /// @returns a pointer to the list, or `nullptr` if invalid fragment
  icarusDB::PMTdigitizerInfoVec const* findPMTfragmentEntry
    (unsigned int fragmentID) const;
  
  /// Returns an exception signed by this object.
  cet::exception myException() const
    { return cet::exception{ "ICARUSChannelMapProviderBase" }; }
  
}; // icarusDB::ICARUSChannelMapProviderBase


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------

#include "icaruscode/Decode/ChannelMapping/ICARUSChannelMapProviderBase.tcc"

// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPROVIDERBASE_H
