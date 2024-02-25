/**
 * @file   icaruscode/Decode/ChannelMapping/ChannelMapSQLite.cxx
 * @author T. Usher (factorised by Gianluca Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ChannelMapSQLite.h
 */

// library header
#include "icaruscode/Decode/ChannelMapping/ChannelMapSQLite.h"
#include "icaruscode/Decode/ChannelMapping/PositionFinder.h"

// ICARUS libraries

// LArSoft libraries
#include "larevt/CalibrationDBI/Providers/DBFolder.h"
#include "larcorealg/CoreUtils/counter.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/search_path.h"

// SQLite
#include "sqlite3.h"

// Guildelines library
#include "gsl/span"

// C++ standard libraries
#include <algorithm> // std::transform()
#include <string>
#include <memory>
#include <cctype> // std::tolower()


// -----------------------------------------------------------------------------
namespace {
  
  /// Tracks a SQLite resource object and frees it on destruction (RAII).
  template <typename PTR = void*>
  class SQLiteResource {
    
    PTR fResource = nullptr; ///< Pointer to the resource data.
    
      public:
    
    using ResourcePtr_t = PTR; ///< Type of the resource data pointer.
    
    /// Constructor: acquires and manages the specified `resource`.
    SQLiteResource(ResourcePtr_t resource = nullptr): fResource{ resource } {}
    
    // do not copy
    SQLiteResource(SQLiteResource const&) = delete;
    SQLiteResource& operator= (SQLiteResource&) = delete;
    
    /// Move constructor: `from` loses the ownership of its resource.
    SQLiteResource(SQLiteResource&& from): fResource{ from.fResource }
      { from.fResource = nullptr; }
    
    /// Move assignment: the current resource is released.
    SQLiteResource& operator= (SQLiteResource&& from)
      {
        if (from.fResource != fResource)
          { free(); std::swap(fResource, from.fResource); }
        return *this;
      }
    
    /// Destructor: releases the resource.
    ~SQLiteResource() { free(); }
    
    //@{
    /// Returns the pointer to the resource data.
    ResourcePtr_t resource() const { return fResource; }
    ResourcePtr_t const* constReference() const { return &fResource; }
    ResourcePtr_t operator() () const { return fResource; }
    operator ResourcePtr_t() const { return resource(); }
    //@}
    
    //@{
    /// Returns _the_ pointer to the address of the resource.
    /// @note The value pointed to is the actual address that this object
    ///       tracks; this allows using this object with functions which set
    ///       the address, but also provides a route to circumvent RAII,
    ///       since the resource pointed by the overwritten value won't be
    ///       freed.
    ResourcePtr_t* reference() { return &fResource; }
    ResourcePtr_t const* reference() const { return constReference(); }
    ResourcePtr_t* operator& () const { return &fResource; }
    ResourcePtr_t* operator& () { return &fResource; }
    //@}
    
    /// Acquires the specified `resource` after freeing the current one.
    SQLiteResource& operator= (ResourcePtr_t resource)
      {
        if (resource != fResource) { free(); fResource = resource; }
        return *this;
      }
    
    /// Releases the currently tracked resource, and stops the tracking.
    void free()
      { if (fResource) { sqlite3_free(fResource); fResource = nullptr; } }
    
  }; // SQLiteResource
  
  
  /// Tracks a SQLite database object and closes it on destruction (RAII).
  class SQLiteConnection {
    
    sqlite3* fDatabase = nullptr; ///< Pointer to the database object.
    
      public:
    SQLiteConnection(sqlite3* database): fDatabase{ database } {}
    
    // this resource handle can't be copied
    SQLiteConnection(SQLiteConnection const&) = delete;
    SQLiteConnection& operator= (SQLiteConnection&) = delete;
    
    /// Move constructor: steals the database `from` the specified handle.
    SQLiteConnection(SQLiteConnection&& from): fDatabase{ from.fDatabase }
      { from.fDatabase = nullptr; }
    /// Move assignment: closes the current database and steals the `from` one.
    SQLiteConnection& operator= (SQLiteConnection&& from)
      {
        if (from.fDatabase != fDatabase)
          { close(); std::swap(fDatabase, from.fDatabase); }
        return *this;
      }
    
    /// Destructor: closes the tracked database, if any. No reaction to failure.
    ~SQLiteConnection() { close(); }
    
    //@{
    /// Returns a pointer to the database object.
    sqlite3* database() const { return fDatabase; }
    operator sqlite3*() const { return database(); }
    //@}
    
    /// Closes the tracked database, if any. No reaction to failure.
    void close()
      { if (fDatabase) { sqlite3_close(fDatabase); fDatabase = nullptr; } }
    
  }; // SQLiteConnection


  /// Debug SQLite callback: prints the content of the row to console.
  int datadump_callback [[maybe_unused]]
    (void*, int argc, char **argv, char **azColName)
  {
    for(int const i: util::counter(argc)) {
      std::cout << "column: '" << azColName[i] << "' - value: "
        << (argv[i] ? argv[i] : "<NULL>") << " ";
    }
    std::cout << "\n" << std::endl;
    return 0;
  } // datadump_callback()

} // local namespace


// -----------------------------------------------------------------------------
std::array
  <icarusDB::ChannelMapSQLite::TableNames_t, icarusDB::RunPeriods::NPeriods>
const icarusDB::ChannelMapSQLite::TableNameSets{
  TableNames_t{
      "readout_boards"  // TPCfragmentMap
    , "daq_channels"    // TPCreadoutBoardMap
    , "pmt_placements"  // PMTfragmentMap
    , "feb_channels"    // CRTsideMap
    , "crtfeb"          // CRTtopMap
  },
  TableNames_t{
      "readout_boards"            // TPCfragmentMap
    , "daq_channels"              // TPCreadoutBoardMap
    , "pmt_placements_23Aug2023"  // PMTfragmentMap
    , "feb_channels"              // CRTsideMap
    , "crtfeb"                    // CRTtopMap
  },
  TableNames_t{
      "readout_boards"            // TPCfragmentMap
    , "daq_channels"              // TPCreadoutBoardMap
    , "pmt_placements_29Aug2023"  // PMTfragmentMap
    , "feb_channels"              // CRTsideMap
    , "crtfeb"                    // CRTtopMap
  }
}; // icarusDB::ChannelMapSQLite::TableNameSets


// -----------------------------------------------------------------------------
icarusDB::ChannelMapSQLite::ChannelMapSQLite(Config const& config)
  : icarus::ns::util::mfLoggingClass{ config.LogCategory() }
  , fDBFileName     { config.DBFileName() }
  , fCalibDBFileName{ config.CalibDBFileName() }
  , fTag            { config.Tag() }
{
}


// -----------------------------------------------------------------------------
bool icarusDB::ChannelMapSQLite::SelectPeriod(RunPeriod period) {
  
  auto const iPeriod = static_cast<unsigned int>(period);
  TableNames_t const* newSet = &(TableNameSets.at(iPeriod));
  
  if (fCurrentTable == newSet) {
    mfLogDebug() << "Period #" << iPeriod << " already selected";
  }
  else if (fCurrentTable) {
    mfLogDebug() << "Switched from period #"
      << (fCurrentTable - &(TableNameSets[0])) << " to #" << iPeriod;
  }
  else {
    mfLogDebug() << "Switching to period #" << iPeriod;
  }
  
  if (fCurrentTable == newSet) return false;
  
  fCurrentTable = newSet;
  return true;
  
} // icarusDB::ChannelMapSQLite::SelectPeriod()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::GetDataset
  (std::string const& table, SQLiteCallbackFunc_t func, void* data) const
{
  
  // find the SQLite database file
  std::string fullFileName;
  cet::search_path searchPath("FW_SEARCH_PATH");
  if (!searchPath.find_file(fDBFileName, fullFileName)) {
    throw cet::exception{ "ChannelMapSQLite" }
      << "GetDataset(): can't find input file: '" << fDBFileName << "'\n";
  }
  mfLogDebug() << "Database file: '" << fullFileName << "'.";
  
  // set up to open the database
  sqlite3* database = nullptr;
  int rc = sqlite3_open(fullFileName.c_str(), &database);
  if (rc) {
    throw cet::exception{ "ChannelMapSQLite" }
      << "GetDataset(): can't open the database, return code:"
      << sqlite3_errmsg(database) << "'\n";
  }
  SQLiteConnection const DBguard{ database };
  
  // read the full table (SQLite will also call func after fetching the data)
  std::string const select = "SELECT * FROM " + table;
  SQLiteResource<char*> zErrMsg;
  rc = sqlite3_exec(database, select.c_str(), func, data, &zErrMsg);
  if (rc != SQLITE_OK) {
    mfLogError() << "GetDataset(): SQL error: " << zErrMsg;
  }
  else {
    mfLogDebug() << "GetDataset(): successfully read SQLite database table '"
      << table << "'";
  }
  
  return rc;
} // icarusDB::ChannelMapSQLite::GetDataset()


// -----------------------------------------------------------------------------
// ---  TPC
// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::buildTPCFragmentIDToReadoutIDMap_callback
  (void* dataOut, int argc, char**argv, char** azColName)
{
  constexpr unsigned int tpcIdentifier(0x00001000);
  // find the position of the columns we need
  /*
   * [20240224] 9 columns:
   * [0] "readout_board_id" [1] "flange_id"   [2] "chimney_number"
   * [3] "tpc_id"           [4] "create_time" [5] "create_user"
   * [6] "update_time"      [7] "update_user" [8] "fragement_id"
   */
  // technical detail: `PositionFinder` uses `operator==` to compare, and the
  // input is C-strings; we force the other term of comparison to C++ strings
  using namespace std::string_literals;
  auto const [ ReadoutBoardIDcolumn, FlangeIDcolumn, FragmentIDcolumn ]
    = icarus::ns::util::PositionFinder{ gsl::span{azColName, (unsigned) argc} }
             ("readout_board_id"s,  "flange_id"s,   "fragement_id"s   );
  
  auto& fragmentBoardMap
    = *static_cast<icarusDB::TPCFragmentIDToReadoutIDMap*>(dataOut);

  // Include a by hand mapping of fragment ID to crate
  // Note that we now know we can get this from the "flanges" table... so an upgrade coming soon...
  static std::map<std::size_t, std::string> const flangeIDToCrateMap = {
    {  19, "WW01T" },
    {  68, "WW01M" },
    {  41, "WW01B" },
    {  11, "WW02"  },
    {  17, "WW03"  },
    {  36, "WW04"  },
    {  18, "WW05"  },
    {  58, "WW06"  },
    {  71, "WW07"  },
    {  14, "WW08"  },
    {  25, "WW09"  },
    {  34, "WW10"  },
    {  67, "WW11"  },
    {  33, "WW12"  },
    {  87, "WW13"  },
    {  10, "WW14"  },
    {  59, "WW15"  },
    {  95, "WW16"  },
    {  22, "WW17"  },
    {  91, "WW18"  },
    {  61, "WW19"  },
    {  55, "WW20T" },
    {  97, "WW20M" },
    { 100, "WW20B" },
    {  83, "WE01T" },
    {  85, "WE01M" },
    {   7, "WE01B" },
    {  80, "WE02"  },
    {  52, "WE03"  },
    {  32, "WE04"  },
    {  70, "WE05"  },
    {  74, "WE06"  },
    {  46, "WE07"  },
    {  81, "WE08"  },
    {  63, "WE09"  },
    {  30, "WE10"  },
    {  51, "WE11"  },
    {  90, "WE12"  },
    {  23, "WE13"  },
    {  93, "WE14"  },
    {  92, "WE15"  },
    {  88, "WE16"  },
    {  73, "WE17"  },
    {   1, "WE18"  },
    {  66, "WE19"  },
    {  48, "WE20T" },
    {  13, "WE20M" },
    {  56, "WE20B" },
    {  94, "EW01T" },
    {  77, "EW01M" },
    {  72, "EW01B" },
    {  65, "EW02"  },
    {   4, "EW03"  },
    {  89, "EW04"  },
    {  37, "EW05"  },
    {  76, "EW06"  },
    {  49, "EW07"  },
    {  60, "EW08"  },
    {  21, "EW09"  },
    {   6, "EW10"  },
    {  62, "EW11"  },
    {   2, "EW12"  },
    {  29, "EW13"  },
    {  44, "EW14"  },
    {   9, "EW15"  },
    {  31, "EW16"  },
    {  98, "EW17"  },
    {  38, "EW18"  },
    {  99, "EW19"  },
    {  53, "EW20T" },
    {  82, "EW20M" },
    {  35, "EW20B" },
    {  96, "EE01T" },
    {  28, "EE01M" },
    {  16, "EE01T" },
    {  69, "EE02"  },
    {  20, "EE03"  },
    {  79, "EE04"  },
    {  50, "EE05"  },
    {  45, "EE06"  },
    {  84, "EE07"  },
    {  42, "EE08"  },
    {  39, "EE09"  },
    {  26, "EE10"  },
    {  64, "EE11"  },
    {  43, "EE12"  },
    {  47, "EE13"  },
    {  15, "EE14"  },
    {   3, "EE15"  },
    {  27, "EE16"  },
    {  24, "EE17"  },
    {  40, "EE18"  },
    {  75, "EE19"  },
    {  86, "EE20T" },
    {  54, "EE20M" },
    {   8, "EE20B" }
  }; // flangeIDToCrateMap[]

  unsigned int const fragmentID
    = std::stol(argv[FragmentIDcolumn], nullptr, 16);
  
  if (fragmentID & tpcIdentifier) {
    if (fragmentBoardMap.find(fragmentID) == fragmentBoardMap.end()) {
      unsigned int const flangeID = std::stol(argv[FlangeIDcolumn]);
      fragmentBoardMap[fragmentID].first = flangeIDToCrateMap.at(flangeID);
    }

    unsigned int const readoutID = std::stol(argv[ReadoutBoardIDcolumn]);
    fragmentBoardMap[fragmentID].second.emplace_back(readoutID);
  }

  return 0;
} // icarusDB::buildTPCReadoutBoardToChannelMap_callback()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::BuildTPCFragmentIDToReadoutIDMap
  (TPCFragmentIDToReadoutIDMap& fragmentBoardMap) const
{
  assert(fCurrentTable);
  const std::string dataType(fCurrentTable->TPCfragmentMap);
  
  // Recover the data from the database
  int error = GetDataset
    (dataType, buildTPCFragmentIDToReadoutIDMap_callback, &fragmentBoardMap);
  
  // If there was an error the function above would have printed a message
  // so bail out
  if (error) {
    throw cet::exception("ChannelMapSQLite")
      << "BuildTPCFragmentIDToReadoutIDMap(): error " << error
      << " in reading the database.\n";
  }
  
  return error;
} // icarusDB::ChannelMapSQLite::BuildTPCFragmentIDToReadoutIDMap()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::buildTPCReadoutBoardToChannelMap_callback
  (void* dataOut, int argc, char**argv, char** azColName)
{
  
  // find the position of the columns we need
  /*
   * [20240224] 19 columns:
   *  [0] "channel_id"     [1] "wire_number"         [2] "readout_board_id"
   *  [3] "chimney_number" [4] "readout_board_slot"  [5] "channel_number"
   *  [6] "create_time"    [7] "create_user"         [8] "update_time"
   *  [9] "update_user"   [10] "plane"              [11] "cable_label_number"
   * [12] "channel_type"
   */
  // technical detail: `PositionFinder` uses `operator==` to compare, and the
  // input is C-strings; we force the other term of comparison to C++ strings
  using namespace std::string_literals;
  auto const [
     ChannelIDcolumn,     ReadoutBoardIDcolumn, ReadoutBoardSlotColumn,
     ChannelNumberColumn, PlaneIdentifierColumn
  ] = icarus::ns::util::PositionFinder{ gsl::span{azColName, (unsigned) argc} }(
    "channel_id"s,       "readout_board_id"s,  "readout_board_slot"s,
    "channel_number"s,   "plane"s
    );
  
  auto& rbChanMap= *static_cast<TPCReadoutBoardToChannelMap*>(dataOut);
  
  unsigned int readoutBoardID = std::stol(argv[ReadoutBoardIDcolumn]);
  
  if (rbChanMap.find(readoutBoardID) == rbChanMap.end()) {
    unsigned int readoutBoardSlot = std::stol(argv[ReadoutBoardSlotColumn]);

    rbChanMap[readoutBoardID].first = readoutBoardSlot;
    rbChanMap[readoutBoardID].second.resize(ChannelsPerTPCreadoutBoard);
  }
  
  unsigned int channelNum = std::stol(argv[ChannelNumberColumn]);
  unsigned int channelID  = std::stol(argv[ChannelIDcolumn]);
  
  std::string const fragmentBuffer = argv[PlaneIdentifierColumn];
  unsigned int const plane = TPCplaneIdentifierToPlane(fragmentBuffer);
  if (plane >= 3) {
    mf::LogError{ "ChannelMapSQLite" } << "YIKES!!! Plane is " << plane
      << " for channel " << channelID << " with type " << fragmentBuffer;
  }
  
  rbChanMap[readoutBoardID].second[channelNum] = { channelID, plane };
  
  return 0;
} // icarusDB::ChannelMapSQLite::buildTPCReadoutBoardToChannelMap_callback()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::BuildTPCReadoutBoardToChannelMap
  (TPCReadoutBoardToChannelMap& rbChanMap) const
{
  assert(fCurrentTable);
  const std::string  dataType(fCurrentTable->TPCreadoutBoardMap);
  
  // Recover the data from the database
  int error = GetDataset
    (dataType, buildTPCReadoutBoardToChannelMap_callback, &rbChanMap);
  
  // If there was an error the function above would have printed a message
  // so bail out
  if (error) {
    throw cet::exception{ "ChannelMapSQLite" }
      << "BuildTPCReadoutBoardToChannelMap(): error " << error
      << " in reading the database!\n";
  }
  
  return error;
} // icarusDB::ChannelMapSQLite::BuildTPCReadoutBoardToChannelMap()


// -----------------------------------------------------------------------------
// --- PMT
// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::buildPMTFragmentToDigitizerChannelMap_callback
  (void* dataOut, int argc, char**argv, char** azColName)
{
  // find the position of the columns we need
  /*
   * [20240224] 19 columns:
   * [0] "pmt_id"               [1] "pmt_sn"            [2] "sector_label"
   * [3] "ch_number"            [4] "pmt_position_code" [5] "hv_cable_label"
   * [6] "signal_cable_label"   [7] "light_fiber_label" [8] "digitizer_label"
   * [9] "digitizer_ch_number" [10] "hv_supply_label" [11] "hv_supply_ch_number"
   * [12] "create_time"        [13] "update_user"      [14] "update_time"
   * [15] "create_user"        [16] "pmt_in_tpc_plane" [17] "channel_id"
   * [18] "fragment_id"
   * Look for the ones we care of:
   */
  // technical detail: `PositionFinder` uses `operator==` to compare, and the
  // input is C-strings; we force the other term of comparison to C++ strings
  using namespace std::string_literals;
  auto const [
     LaserChannelColumn,   DigitizerChannelColumn, ChannelIDcolumn,
     FragmentIDcolumn
  ] = icarus::ns::util::PositionFinder{ gsl::span{azColName, (unsigned) argc} }(
    "light_fiber_label"s, "digitizer_ch_number"s, "channel_id"s,
    "fragment_id"s
    );
  
  auto& fragmentToDigitizerChannelMap
    = *static_cast<FragmentToDigitizerChannelMap*>(dataOut);
  
  // Start extracting info
  unsigned int const fragmentID         = std::stol(argv[FragmentIDcolumn]);
  unsigned int const digitizerChannelNo = std::stol(argv[DigitizerChannelColumn]);
  unsigned int const channelID          = std::stol(argv[ChannelIDcolumn]);

  // Read the laser channel; format is `L-<number>`. <number> is int from [1-41]
  std::string const laserChannelLabel = argv[LaserChannelColumn];
  unsigned int laserChannel;
  try {
    laserChannel = std::stol(laserChannelLabel.substr(2));
  }
  catch(...) {
    throw cet::exception{ "ChannelMapSQLite" }
      << "Failed to convert laser channel '" << laserChannelLabel
      << "' into a channel number for channel " << channelID << "!\n";
  }

  // Fill the map
  fragmentToDigitizerChannelMap[fragmentID].emplace_back
    (digitizerChannelNo, channelID, laserChannel);
  
  return 0;
} // ...::ChannelMapSQLite::buildPMTFragmentToDigitizerChannelMap_callback()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::BuildPMTFragmentToDigitizerChannelMap
  (FragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap) const
{
  // clearing is cleansing
  fragmentToDigitizerChannelMap.clear();
  
  // Recover the information from the database on the mapping
  assert(fCurrentTable);
  const std::string  dataType(fCurrentTable->PMTfragmentMap);
  
  // Recover the data from the database
  int error = GetDataset(
    dataType, buildPMTFragmentToDigitizerChannelMap_callback,
    &fragmentToDigitizerChannelMap
  );
  
  // If there was an error the function above would have printed a message
  // so bail out
  if (error) {
    throw cet::exception{ "ChannelMapSQLite" }
      << "BuildFragmentToDigitizerChannelMap(): encountered error (code: "
      << error << ") in reading the database\n";
  }
  
  return error;
} // icarusDB::ChannelMapSQLite::BuildPMTFragmentToDigitizerChannelMap()


// -----------------------------------------------------------------------------
// --- CRT
// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::buildCRTChannelIDToHWtoSimMacAddressPairMap_callback
  (void* dataOut, int argc, char**argv, char** azColName)
{
  /*
   * [20240224] 13 columns:
   *  [0] "feb_id"            [1] "feb_channel"  [2] "pedestal"
   *  [3] "threshold_adjust"  [4] "bias"         [5] "hg"
   *  [6] "create_time"       [7] "update_user"  [8] "update_time"
   *  [9] "create_user"      [10] "channel_id"  [11] "feb_index"
   * [12] "mac_address"
   */
  // technical detail: `PositionFinder` uses `operator==` to compare, and the
  // input is C-strings; we force the other term of comparison to C++ strings
  using namespace std::string_literals;
  auto const [ ChannelIDcolumn, SimMACaddressColumn, HWaddressColumn ]
    = icarus::ns::util::PositionFinder{ gsl::span{azColName, (unsigned) argc} }
             ("channel_id"s,   "feb_index"s,        "mac_address"s   );
  
  auto& crtChannelIDToHWtoSimMacAddressPairMap
    = *static_cast<CRTChannelIDToHWtoSimMacAddressPairMap*>
    (dataOut);

  // extracting info
  auto const valueOrNone
    = [](const char* s){ return std::strcmp(s, "None") == 0? 0: std::stol(s); };
  
  unsigned int const channelID     = valueOrNone(argv[ChannelIDcolumn]);
  unsigned int const simmacaddress = valueOrNone(argv[SimMACaddressColumn]);
  unsigned int const hwmacaddress  = valueOrNone(argv[HWaddressColumn]); 
  
  // fill the map
  crtChannelIDToHWtoSimMacAddressPairMap[channelID]
    = std::make_pair(hwmacaddress, simmacaddress);

  return 0;
} // ...::ChannelMapSQLite::buildCRTChannelIDToHWtoSimMacAddressPairMap_callback


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::BuildCRTChannelIDToHWtoSimMacAddressPairMap
  (CRTChannelIDToHWtoSimMacAddressPairMap& crtChannelIDToHWtoSimMacAddressPairMap)
  const
{
  // clearing is cleansing
  crtChannelIDToHWtoSimMacAddressPairMap.clear();
  
  // Recover the information from the database on the mapping
  assert(fCurrentTable);
  const std::string  dataType(fCurrentTable->CRTsideMap);
  
  // Recover the data from the database
  int error = GetDataset(
    dataType,
    buildCRTChannelIDToHWtoSimMacAddressPairMap_callback,
    &crtChannelIDToHWtoSimMacAddressPairMap
  );
  
  // If there was an error the function above would have printed a message
  // so bail out
  if (error) {
    throw cet::exception{ "ChannelMapSQLite" }
      << "BuildCRTChannelIDToHWtoSimMacAddressPairMap(): error " << error
      << " reading the database!\n";
  }
  
  return error;
} // icarusDB::ChannelMapSQLite::BuildCRTChannelIDToHWtoSimMacAddressPairMap()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::buildTopCRTHWtoSimMacAddressPairMap_callback
  (void* dataOut, int argc, char**argv, char** azColName)
{
  /*
   * [20240224] 42 columns:
   *  [0] "feb_barcode"  [1] "serialnum"    [2] "mac_add8b"
   *  [3] "mac_add"      [4] "voltage"      [5] "ch0"
   *  [6] "ch1"          [7] "ch2"          [8] "ch3"
   *  [9] "ch4"         [10] "ch5"         [11] "ch6"
   * [12] "ch7"         [13] "ch8"         [14] "ch9"
   * [15] "ch10"        [16] "ch11"        [17] "ch12"
   * [18] "ch13"        [19] "ch14"        [20] "ch15"
   * [21] "ch16"        [22] "ch17"        [23] "ch18"
   * [24] "ch19"        [25] "ch20"        [26] "ch21" 
   * [27] "ch22"        [28] "ch23"        [29] "ch24" 
   * [30] "ch25"        [31] "ch26"        [32] "ch27" 
   * [33] "ch28"        [34] "ch29"        [35] "ch30" 
   * [36] "ch31"        [37] "create_time" [38] "update_user" 
   * [39] "update_time" [40] "create_user" [41] "feb_index"
   */
  // technical detail: `PositionFinder` uses `operator==` to compare, and the
  // input is C-strings; we force the other term of comparison to C++ strings
  using namespace std::string_literals;
  auto const [ SimMACaddressColumn, HWaddressColumn      ]
    = icarus::ns::util::PositionFinder{ gsl::span{azColName, (unsigned) argc} }
             ("feb_index"s,        "mac_add"s            );
  
  auto& topcrtHWtoSimMacAddressPairMap
    = *static_cast<TopCRTHWtoSimMacAddressPairMap*>(dataOut);
  
  auto const valueOrNone
    = [](const char* s){ return std::strcmp(s, "None") == 0? 0: std::stol(s); };
  
  // start extracting info
  unsigned int const simmacaddress = valueOrNone(argv[SimMACaddressColumn]);
  unsigned int const hwmacaddress  = valueOrNone(argv[HWaddressColumn]);
      
  // fill the map
  topcrtHWtoSimMacAddressPairMap[hwmacaddress] = simmacaddress;

  return 0;
} // icarusDB::ChannelMapSQLite::buildTopCRTHWtoSimMacAddressPairMap_callback()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::BuildTopCRTHWtoSimMacAddressPairMap
  (TopCRTHWtoSimMacAddressPairMap& topcrtHWtoSimMacAddressPairMap) const
{
  // clearing is cleansing
  topcrtHWtoSimMacAddressPairMap.clear();
  
  // Recover the information from the database on the mapping
  assert(fCurrentTable);
  const std::string  dataType(fCurrentTable->CRTtopMap);
  
  // Recover the data from the database
  int error = GetDataset(
    dataType,
    buildTopCRTHWtoSimMacAddressPairMap_callback,
    &topcrtHWtoSimMacAddressPairMap
    );
  
  // If there was an error the function above would have printed a message
  // so bail out
  if (error) {
    throw cet::exception{ "ChannelMapSQLite" }
      << "BuildTopCRTHWtoSimMacAddressPairMap(): error " << error
      << " reading the database!\n";
  }
  
  return error;
} // icarusDB::ChannelMapSQLite::BuildTopCRTHWtoSimMacAddressPairMap()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::BuildSideCRTCalibrationMap
  (SideCRTChannelToCalibrationMap& sideCRTChannelToCalibrationMap) const
{
  // clearing is cleansing
  sideCRTChannelToCalibrationMap.clear();
  
  std::string fullFileName;
  cet::search_path searchPath("FW_SEARCH_PATH");
  if (!searchPath.find_file(fCalibDBFileName + ".db", fullFileName)) {
    throw cet::exception{ "ChannelMapSQLite" }
      << "BuildSideCRTCalibrationMap(): can't find calibration input file: '"
      << fCalibDBFileName << "'\n";
  }
  mfLogDebug() << "Found calibration input file: '" << fCalibDBFileName << "'";
  
  lariov::DBFolder database(fCalibDBFileName, "", "", fTag, true, false);
  
  database.UpdateData(1'638'918'271'000'000'000ULL);
  
  std::vector<unsigned int> channels;
  database.GetChannelList(channels);
  
  for (unsigned int const channel: channels) {
    
    long mac5, chan;
    double  gain, ped;
    
    // extracting info
    database.GetNamedChannelData(channel, "mac5", mac5);
    database.GetNamedChannelData(channel, "localchannel", chan);
    database.GetNamedChannelData(channel, "gain", gain);
    database.GetNamedChannelData(channel, "pedestal", ped);
    
    // fill the map
    sideCRTChannelToCalibrationMap.emplace
      (std::make_pair((int) mac5, (int) chan), std::make_pair(gain, ped));
  } // for

  return 0;
} // icarusDB::ChannelMapSQLite::BuildSideCRTCalibrationMap()


// -----------------------------------------------------------------------------
