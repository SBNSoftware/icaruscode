/**
 * @file   icaruscode/Decode/ChannelMapping/ChannelMapSQLite.h
 * @author T. Usher (factorised by Gianluca Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ChannelMapSQLite.cxx
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_CHANNELMAPSQLITE_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_CHANNELMAPSQLITE_H

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"
#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"
#include "icarusalg/Utilities/mfLoggingClass.h"

// framework libraries
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// C++ standard libraries
#include <set>
#include <string>


namespace icarusDB { class ChannelMapSQLite; };

/**
 * @brief Interface with ICARUS channel mapping SQLite database.
 * 
 * This interface fills "standard" data structures from a SQLite database.
 * These data structures are then managed by a ICARUS service provider
 * (`icarusDB::IICARUSChannelMap`).
 * 
 */
class icarusDB::ChannelMapSQLite
  : public IChannelMapping, private icarus::ns::util::mfLoggingClass
{
  /// Tells `fhicl::Table` to ignore the parameter for _art_ tool configuration.
  struct IgnoreArtToolConfig
    { std::set<std::string> operator()() { return { "tool_type" }; } };
  
    public:
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<std::string> DBFileName {
      Name{ "DBFileName" },
      Comment{ "name of the channel mapping database SQLite file" }
      };
    
    fhicl::Atom<std::string> CalibDBFileName {
      Name{ "CalibDBFileName" },
      Comment{ "name of the CRT calibration database SQLite file" }
      };
    
    fhicl::Atom<std::string> Tag {
      Name{ "Tag" },
      Comment{ "tag/version of the database to load" }
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of the streams to send console messages to" },
      "ChannelMapSQLite" // default
      };
    
  }; // Config
  
  
  /// Constructor: configures the object.
  explicit ChannelMapSQLite(Config const& config);
  
  
  /**
   * @brief   Prepares the object for queries pertaining the specified period.
   * @param   period the period to be prepared for
   * @return  whether values cached from the previous period are invalidated
   * 
   * Periods are defined in `icarusDB::RunPeriods` class.
   * 
   * If the return value is `true`, caller should follow by rebuilding all the
   * maps they need to use.
   * If the return value is `false`, all maps built after the previous call to
   * `SelectPeriod()` are still current for `period`: building new ones is not
   * necessary, but it's not harmful either (at most, wasteful).
   */
  virtual bool SelectPeriod(RunPeriod period) override;
  
  
  /**
   *  @brief Define the returned data structures for a mapping between TPC Fragment IDs
   *         and the related crate and readout information. 
   *         Then define the function interface to fill these data structures 
   */
  virtual int BuildTPCFragmentIDToReadoutIDMap(TPCFragmentIDToReadoutIDMap&) const override;
  
  /**
   *  @brief Define the returned data structures for a mapping between TPC readout boards
   *         and the channel information 
   *         Then define the function interface to fill these data structures 
   */
  virtual int BuildTPCReadoutBoardToChannelMap(TPCReadoutBoardToChannelMap&) const override;
  
  /**
   *  @brief Define the returned data structures for a mapping between PMT Fragment IDs
   *         and the related crate and readout information. 
   *         Then define the function interface to fill these data structures 
   */
  virtual int BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap&) const override;
  
  
  /**
   *  @brief Define the returned data structures for a mapping between CRT hardware mac_address
   *         to the simulated mac_address.
   *         Then define the function interface to fill these data structures
   */
  virtual int BuildCRTChannelIDToHWtoSimMacAddressPairMap(CRTChannelIDToHWtoSimMacAddressPairMap&) const override;
  virtual int BuildTopCRTHWtoSimMacAddressPairMap(TopCRTHWtoSimMacAddressPairMap&) const override;

  virtual int BuildSideCRTCalibrationMap(SideCRTChannelToCalibrationMap&) const override;
  
    private:

  /// Record of all relevant table names.
  struct TableNames_t {
    std::string TPCfragmentMap;
    std::string TPCreadoutBoardMap;
    std::string PMTfragmentMap;
    std::string CRTsideMap;
    std::string CRTtopMap;
    // NOTE: CRT side calibration is is a different type of database which uses a timestamp
  };
  
  /// SQLite callback function: SQLite will call this including `argc` arguments
  /// each from a column named from `colNames` with the value in `argv`.
  using SQLiteCallbackFunc_t = int(*)(void* data, int argc, char** argv, char** colNames);


  // --- BEGIN --- Configuration parameters ------------------------------------
  
  std::string const fDBFileName;      ///< File name of our SQLite database.
  std::string const fCalibDBFileName; ///< File name of our side CRT calibration SQLite database.
  std::string const fTag;             ///< Tag for conditioned database.
  
  std::string const fLogCategory;     ///< Message facility category name.
  
  // --- END ----- Configuration parameters ------------------------------------

  /// The set of tables being served. Chosen by `SelectRun()`.
  TableNames_t const* fCurrentTable = nullptr;

  /**
   * @brief Reads a full `table` from the database and process data with `func`.
   * @param table name of the table to read
   * @param func callback function to be called to process the fetched data
   * @param data pointer to a memory data region
   * @return an error code, `SQLITE_OK` on success
   * 
   * The full content of the table named `table` is read from the configured
   * database. A connection to the database is opened, used and finally closed.
   * 
   * The callback function `func` is executed with the newly loaded data.
   * The pointer `data` is passed to that function, which can reinterpret it
   * and modify it in the usual untyped C-style data transfer.
   */
  int GetDataset
    (std::string const& table, SQLiteCallbackFunc_t func, void* data) const;
  
  
  /// The list of names of tables.
  static std::array<TableNames_t, RunPeriods::NPeriods> const TableNameSets;
  
  
  // --- BEGIN --- SQLite callback functions -----------------------------------
  /**
   * @name SQLite callback functions
   * 
   * These (static) functions are executed during SQLite calls (in
   * `GetDataset()`) to translate the data from a table into something we can
   * use.
   * They all adhere to the prototype `SQLiteCallbackFunc_t`, with:
   * 
   * @param[out] dataOut a pointer to a data area; this pointer can be casted
   *                     to the actual data type and then filled
   * @param argc the number of columns in the table
   * @param argv the values on columns in one database row (`argc` in count)
   * @param azColName the names of the columns (`argc` in count)
   */
  /// @{
  
  /**
   * @brief Builds a map between the TPC Fragment IDs and the readout board IDs.
   * 
   * Here we expect there will be some number of boards per Fragment.
   */
  static int buildTPCFragmentIDToReadoutIDMap_callback
    (void* dataOut, int argc, char**argv, char** azColName);
  
  /**
   * @brief Builds a map between the TPC readout board IDs and the associated
   *        channels.
   *
   * For each readout board ID we expect a number of channels back from the data
   * base. So the returned data structure will be a map of readout ID to a
   * vector of channels.
   */
  static int buildTPCReadoutBoardToChannelMap_callback
    (void* dataOut, int argc, char**argv, char** azColName);
  
  /// Fills the PMT fragment map with the information from one channel.
  static int buildFragmentToDigitizerChannelMap_callback
    (void* dataOut, int argc, char**argv, char** azColName);
  
  /// Fills the channel mapping for side CRT.
  static int buildCRTChannelIDToHWtoSimMacAddressPairMap_callback
    (void* dataOut, int argc, char**argv, char** azColName);
  
  /// Fills the channel mapping for top CRT.
  static int buildTopCRTHWtoSimMacAddressPairMap_callback
    (void* dataOut, int argc, char**argv, char** azColName);
  
  /// @}
  // --- END ----- SQLite callback functions -----------------------------------
  
  
}; // icarusDB::ChannelMapSQLite


#endif // ICARUSCODE_DECODE_CHANNELMAPPING_CHANNELMAPSQLITE_H
