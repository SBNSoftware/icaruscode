/**
 * @file   icaruscode/Decode/ChannelMapping/ChannelMapSQLite.cxx
 * @author T. Usher (factorised by Gianluca Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ChannelMapSQLite.h
 */

// library header
#include "icaruscode/Decode/ChannelMapping/ChannelMapSQLite.h"

// ICARUS libraries

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
// #include "wda.h"

// LArSoft libraries
#include "larevt/CalibrationDBI/Providers/DBFolder.h"
#include "larcorealg/CoreUtils/counter.h"
#include "cetlib/search_path.h"

// SQLite
#include <sqlite3.h>

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
// The aim of this function is to build a map between the
// TPC Fragment IDs and the readout board IDs. Here we 
// expect there will be some number of boards per Fragment
//-----------------------------------------------------
int icarusDB::ChannelMapSQLite::buildTPCFragmentIDToReadoutIDMap_callback
  (void* dataOut, int argc, char**argv, char** azColName)
{
  const unsigned int tpcIdentifier(0x00001000);

  auto& fragmentBoardMap = *static_cast
    <icarusDB::IChannelMapping::TPCFragmentIDToReadoutIDMap*>(dataOut);

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

  constexpr std::size_t ReadoutBoardIDcolumn = 0;
  constexpr std::size_t FlangeIDcolumn = 1;
  constexpr std::size_t FragmentIDcolumn = 8;
  
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
  constexpr std::size_t ChannelIDcolumn        =  0;
  constexpr std::size_t ReadoutBoardIDcolumn   =  2;
  constexpr std::size_t ReadoutBoardSlotColumn =  4;
  constexpr std::size_t ChannelNumberColumn    =  5;
  constexpr std::size_t FragmentBufferColumn   = 10;
  
  auto& rbChanMap
    = *static_cast<IChannelMapping::TPCReadoutBoardToChannelMap*>(dataOut);
  
  unsigned int readoutBoardID = std::stol(argv[ReadoutBoardIDcolumn]);
  
  if (rbChanMap.find(readoutBoardID) == rbChanMap.end()) {
    unsigned int readoutBoardSlot = std::stol(argv[ReadoutBoardSlotColumn]);

    rbChanMap[readoutBoardID].first = readoutBoardSlot;
    rbChanMap[readoutBoardID].second.resize(ChannelsPerTPCreadoutBoard);
  }
  
  unsigned int channelNum = std::stol(argv[ChannelNumberColumn]);
  unsigned int channelID  = std::stol(argv[ChannelIDcolumn]);
  
  std::string const fragmentBuffer = argv[FragmentBufferColumn];
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
int icarusDB::ChannelMapSQLite::buildFragmentToDigitizerChannelMap_callback
  (void* dataOut, int argc, char**argv, char** azColName)
{
  constexpr std::size_t LaserChannelColumn     =  7;
  constexpr std::size_t DigitizerChannelColumn =  9;
  constexpr std::size_t ChannelIDcolumn        = 17;
  constexpr std::size_t FragmentIDcolumn       = 18;
  
  auto& fragmentToDigitizerChannelMap
    = *static_cast<IChannelMapping::FragmentToDigitizerChannelMap*>(dataOut);
  
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
} // icarusDB::ChannelMapSQLite::buildFragmentToDigitizerChannelMap_callback()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::BuildFragmentToDigitizerChannelMap
  (FragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap) const
{
  // clearing is cleansing
  fragmentToDigitizerChannelMap.clear();
  
  // Recover the information from the database on the mapping
  assert(fCurrentTable);
  const std::string  dataType(fCurrentTable->PMTfragmentMap);
  
  // Recover the data from the database
  int error = GetDataset(
    dataType, buildFragmentToDigitizerChannelMap_callback,
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
} // icarusDB::ChannelMapSQLite::BuildFragmentToDigitizerChannelMap()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapSQLite::buildCRTChannelIDToHWtoSimMacAddressPairMap_callback
  (void* dataOut, int argc, char**argv, char** azColName)
{
  constexpr std::size_t ChannelIDcolumn     = 10;
  constexpr std::size_t SimMACaddressColumn = 11;
  constexpr std::size_t HWaddressColumn     = 12;
  
  auto& crtChannelIDToHWtoSimMacAddressPairMap
    = *static_cast<IChannelMapping::CRTChannelIDToHWtoSimMacAddressPairMap*>
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
  constexpr std::size_t SimMACaddressColumn = 41;
  constexpr std::size_t HWaddressColumn     =  3;
  
  auto& topcrtHWtoSimMacAddressPairMap
    = *static_cast<IChannelMapping::TopCRTHWtoSimMacAddressPairMap*>(dataOut);
  
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
