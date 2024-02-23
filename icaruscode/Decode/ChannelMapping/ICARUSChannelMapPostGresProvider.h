/**
 * @file   icaruscode/Decode/ChannelMapping/ICARUSChannelMapPostGresProvider.h
 * @author T. Usher (factorised by G. Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ICARUSChannelMapPostGresProvider.cxx
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPOSTGRESPROVIDER_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPOSTGRESPROVIDER_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"
#include "icaruscode/Decode/ChannelMapping/ChannelMapPostGres.h"

// framework libraries
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <string>
#include <memory> // std::unique_ptr<>


// -----------------------------------------------------------------------------
namespace icarusDB { class ICARUSChannelMapPostGresProvider; }
/**
 * @brief Interface to the PostgreSQL ICARUS channel mapping database.
 * @see icarusDB::ICARUSChannelMapPostGres
 * 
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * * `ChannelMappingTool` (algorithm configuration): see
 *     `icarusDB::ChannelMapPostGres` configuration.
 * * `DiagnosticOutput` (flag, default: `false`): verbosely load parameters
 *     and in general prints more messages to console
 * * `LogCategory` (string, default: `"ICARUSChannelMapPostGresProvider"`):
 *     name of the messagefacility category used to send messages to console,
 *     useful for filtering messages.
 * 
 * 
 */
class icarusDB::ICARUSChannelMapPostGresProvider: public IICARUSChannelMap {
  
    public:
  /// The helper used to extract data from the database.
  using ChannelMappingAlg_t = icarusDB::ChannelMapPostGres;
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Table<ChannelMappingAlg_t::Config> ChannelMappingTool {
      Name{ "ChannelMappingTool" },
      Comment{ "parameters for the channel mapping algorithm" }
      };
    
    fhicl::Atom<bool> DiagnosticOutput{
      Name{ "DiagnosticOutput" },
      Comment{ "enables verbose output to console" },
      false // default
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "name of the console stream to send messages to" },
      "ICARUSChannelMapPostGresProvider" // default
      };
    
  }; // Config
  
  using Parameters = fhicl::Table<Config>;
  
  
  /// Constructor: configures the service provider via FHiCL configuration.
  ICARUSChannelMapPostGresProvider(Config const& config);
  
  /// Constructor: configures the service provider via FHiCL configuration.
  ICARUSChannelMapPostGresProvider(Parameters const& params)
    : ICARUSChannelMapPostGresProvider{ params() }
    {}
  
  
  /// --- BEGIN --- Data period selection ------------------------------------
  /// @name Data period selection
  /// @{
  
  /// Loads the mapping for `run`, returns whether a new mapping was loaded.
  virtual bool forRun(int run) override;
  
  /// Loads the mapping for `period`, returns whether a new mapping was loaded.
  virtual bool forPeriod(icarusDB::RunPeriod period) override;
  
  /// @}
  /// --- END ----- Data period selection ------------------------------------
  
  
  /// --- BEGIN --- TPC information ------------------------------------------
  /// @name TPC information
  /// @{
  
  /// --- BEGIN - - TPC fragment information - - - - - - - - - - - - - - - - -
  /// @name TPC fragment information
  /// @{
  
  /// Returns whether the specified `ID` is a known TPC fragment ID.
  virtual bool hasFragmentID(const unsigned int ID) const override;
  /// Returns the number of TPC fragment IDs known to the service.

  /// Returns the number of known TPC fragments.
  virtual unsigned int nTPCfragmentIDs() const override;
  
  /// Returns the name of the crate served by the specified `fragmentID`.
  virtual std::string const& getCrateName
    (const unsigned int fragmentID) const override;
  
  /// Returns the list of board IDs included in the specified `fragmentID`.
  virtual ReadoutIDVec const& getReadoutBoardVec
    (const unsigned int fragmentID) const override;
  
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
  virtual bool hasBoardID(const unsigned int ID) const override;
  /// Returns the number of TPC board IDs known to the service.
  virtual unsigned int nTPCboardIDs() const override;
  /// Returns the number of slot the `boardID` is on.
  virtual unsigned int getBoardSlot(const unsigned int boardID) const override;
  /// Returns a list of channels served by the `boardID` and for each the plane
  /// it is on (`0` to `2`).
  virtual ChannelPlanePairVec const& getChannelPlanePair
    (const unsigned int boardID) const override;

  /// @}
  /// --- END - - - TPC board information  - - - - - - - - - - - - - - - - - -
  
  
  /// --- BEGIN --- PMT information ------------------------------------------
  /// @name PMT information
  /// @{
  
  /// Returns whether the specified fragment `ID` is known to the mapping.
  virtual bool hasPMTDigitizerID(const unsigned int ID) const override;
  
  /// Returns the number of PMT fragment IDs known to the mapping.
  virtual unsigned int nPMTfragmentIDs() const override;
  
  /// Returns a list of triplets: 
  virtual DigitizerChannelChannelIDPairVec const& getChannelIDPairVec
    (const unsigned int) const override;
  
  /// @}
  /// --- END ----- PMT information ------------------------------------------


  /// --- BEGIN --- CRT information ------------------------------------------
  /// @name CRT information
  /// @{
  
  /// Returns the sim Mac address corresponding to the specified side CRT hardware address.
  virtual unsigned int getSimMacAddress(const unsigned int hwmacaddress) const override;
  /// Returns the sim Mac address corresponding to the specified top CRT hardware address.
  virtual unsigned int gettopSimMacAddress(const unsigned int) const override;

  /// Returns the Gain and Pedestal for Side CRT.
  virtual std::pair<double, double> getSideCRTCalibrationMap
    (int mac5, int chan) const override;
  
  /// @}
  /// --- END ----- CRT information ------------------------------------------

  /// Returns the channel mapping database key for the specified PMT fragment ID.
  static constexpr unsigned int PMTfragmentIDtoDBkey(unsigned int fragmentID);
  
  /// Returns the PMT fragment ID for the specified channel mapping database key.
  static constexpr unsigned int DBkeyToPMTfragmentID(unsigned int DBkey);


    private:
  
  // --- BEGIN --- Configuration parameters ------------------------------------
      
  bool const fDiagnosticOutput;
  
  std::string const fLogCategory;
  
  // --- END ----- Configuration parameters ------------------------------------
  
  
  // --- BEGIN --- Cache -------------------------------------------------------
  IChannelMapping::TPCFragmentIDToReadoutIDMap   fFragmentToReadoutMap;
    
  IChannelMapping::TPCReadoutBoardToChannelMap   fReadoutBoardToChannelMap;

  IChannelMapping::FragmentToDigitizerChannelMap fFragmentToDigitizerMap; 

  IChannelMapping::CRTChannelIDToHWtoSimMacAddressPairMap fCRTChannelIDToHWtoSimMacAddressPairMap;

  IChannelMapping::TopCRTHWtoSimMacAddressPairMap fTopCRTHWtoSimMacAddressPairMap;

  IChannelMapping::SideCRTChannelToCalibrationMap fSideCRTChannelToCalibrationMap;

  // --- END ----- Cache -------------------------------------------------------
  
  /// The helper extracting information from the database.
  std::unique_ptr<IChannelMapping> fChannelMappingAlg;


  /// Has the channel mapping tool fill the mapping caches.
  virtual void readFromDatabase();

  
  /// Returns the list of board channel-to-PMT channel ID mapping within the
  /// specified fragment.
  /// 
  /// @returns a pointer to the mapping list, or `nullptr` if invalid fragment
  DigitizerChannelChannelIDPairVec const* findPMTfragmentEntry
    (unsigned int fragmentID) const;
  
  
}; // icarusDB::ICARUSChannelMapPostGresProvider


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_ICARUSCHANNELMAPPOSTGRESPROVIDER_H
