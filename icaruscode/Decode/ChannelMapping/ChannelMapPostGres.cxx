/**
 * @file   icaruscode/Decode/ChannelMapping/ChannelMapPostGres.cxx
 * @brief  Interface with ICARUS channel mapping PostGres database.
 * @author T. Usher (factorised by G. Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ChannelMapPostGres.h
 */

// library header
#include "icaruscode/Decode/ChannelMapping/ChannelMapPostGres.h"

// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"
#include "icaruscode/Decode/ChannelMapping/PositionFinder.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "wda.h" // getDataWithTimeout(), getTuple(), getStringValue(), ...

// C/C++ standard libraries
#include <string>
#include <array>
#include <map>
#include <memory>
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarusDB::details {
  
  struct WDADatasetDeleter
    { void operator() (Dataset dataset) { releaseDataset(dataset); } };
  struct WDATupleDeleter
    { void operator() (Tuple tuple) { releaseTuple(tuple); } };
  
  struct WDADataset: std::unique_ptr<void, WDADatasetDeleter> {
    WDADataset(Dataset dataset)
      : std::unique_ptr<void, WDADatasetDeleter>{ dataset } {}
    operator Dataset() const { return get(); }
  };
  
  struct WDATuple: std::unique_ptr<void, WDATupleDeleter> {
    WDATuple(Tuple tuple): std::unique_ptr<void, WDATupleDeleter>{ tuple } {}
    operator Tuple() const { return get(); }
  };
  
  
  /// Wrapper of `PositionFinder` to extract the names from a `WDATuple`.
  class WDAPositionFinder
    : public icarus::ns::util::PositionFinder<std::vector<std::string>>
  {
    
    std::vector<std::string> fNames; ///< Storage for the column names
    
      public:
    WDAPositionFinder(WDATuple&& namesTuple)
      : icarus::ns::util::PositionFinder<std::vector<std::string>>{ fNames }
      , fNames{ tupleToVector(namesTuple) }
      {}
    
    static std::vector<std::string> tupleToVector(WDATuple& namesTuple)
      {
        unsigned int const n = getNfields(namesTuple);
        std::vector<std::string> names{ n };
        for (unsigned int const i: util::counter(n)) {
          std::string& name = names[i]; // util::enumerate() fails on l-values
          int error = 0;
          name.resize(128, '\0');
          std::size_t const length
            = getStringValue(namesTuple, i, name.data(), name.size(), &error);
          if (error) { // we know very little of the context... sorry, use gdb
            throw cet::exception{ "ChannelMapPostGres" }
              << "Failed (code: " << error << ") to extract column #" << i
              << " name from a channel mapping database table";
          }
          name.erase(length);
        }
        return names;
      } // tupleToVector()
    
  }; // WDAPositionFinder
  
} // namespace icarusDB::details


// -----------------------------------------------------------------------------
std::array<std::string, icarusDB::RunPeriods::NPeriods> const
icarusDB::ChannelMapPostGres::PMTTimestampSet{
   "start"
  ,"23aug2023"
  ,"29aug2023"
};


// -----------------------------------------------------------------------------
icarusDB::ChannelMapPostGres::ChannelMapPostGres(Config const& config)
  : icarus::ns::util::mfLoggingClass{ config.LogCategory() }
  , fDBURL                { config.DatabaseURL() }
  , fCRTcalibrationDBURL  { config.CRTcalibrationDatabaseURL() }
  , fDatabaseAccessTimeout{ config.DatabaseAccessTimeout() }
  {}


// -----------------------------------------------------------------------------
bool icarusDB::ChannelMapPostGres::SelectPeriod(RunPeriod period) {
  
  auto const iPeriod = static_cast<unsigned int>(period);
  std::string newPMTTimestamp = PMTTimestampSet.at(iPeriod);
  
  if (fCurrentPMTTimestamp == newPMTTimestamp) {
    mfLogDebug() << "Period " << newPMTTimestamp << " already selected";
  }
  else if ( !fCurrentPMTTimestamp.empty() ) {
    mfLogDebug() << "Switched from period " << fCurrentPMTTimestamp
      << " to " << newPMTTimestamp;
  }
  else {
    mfLogDebug() << "Switching to period " << newPMTTimestamp;
  }
  
  if (fCurrentPMTTimestamp == newPMTTimestamp) return false;
  
  fCurrentPMTTimestamp = newPMTTimestamp;
  return true;
  
} // icarusDB::ChannelMapPostGres::SelectPeriod()


// -----------------------------------------------------------------------------
icarusDB::details::WDADataset icarusDB::ChannelMapPostGres::GetDataset
  (std::string const& name, std::string url, std::string const& dataType) const
{
  if (!dataType.empty()) url += "&t=" + dataType;
  int error{ 0 };
  details::WDADataset dataPtr{
    getDataWithTimeout
      (url.c_str(), name.c_str(), fDatabaseAccessTimeout, &error)
    };
  int const status = getHTTPstatus(dataPtr);
  if (status != 200) {
    mfLogError() << "libwda error: HTTP code=" << status << ": '"
      << getHTTPmessage(dataPtr) << "'";
  }
  if (error) {
    auto e = myException()
      << "GetDataset(): database access  failed with error " << error
      << "\nDatabase URL: '" << url << "', table: '" << name << "', type: '"
      << dataType << "'\n";
    if (status != 200) {
      e << "\nlibwda error: HTTP code=" << status << ": '"
        << getHTTPmessage(dataPtr) << "'";
    }
    throw e;
  }
  return dataPtr;
} // icarusDB::ChannelMapPostGres::GetDataset()


// -----------------------------------------------------------------------------
icarusDB::details::WDADataset icarusDB::ChannelMapPostGres::GetCRTCaldata
  (std::string const& name, std::string url) const
{
  
  return GetDataset(name, std::move(url), "");
  
}


// -----------------------------------------------------------------------------
// ---  TPC
// -----------------------------------------------------------------------------
int icarusDB::ChannelMapPostGres::BuildTPCFragmentIDToReadoutIDMap
  (TPCFragmentIDToReadoutIDMap& fragmentBoardMap) const
{
  using namespace std::string_literals;
  constexpr unsigned int tpcIdentifier(0x00001000);
  
  // Recover the data from the database
  std::string const name { "icarus_hw_readoutboard" };
  std::string const dataType { "readout_boards" };
  details::WDADataset dataset = GetDataset(name, fDBURL, dataType);
  
  // Include a by hand mapping of fragment ID to position in "corner" chimneys
  // (this information is not present in the database)
  static std::map<unsigned int, std::string> const flangeIDsuffixes = {
    {  19, "T"s }, {  68, "M"s }, {  41, "B"s }, // WW01
    {  55, "T"s }, {  97, "M"s }, { 100, "B"s }, // WW20
    {  83, "T"s }, {  85, "M"s }, {   7, "B"s }, // WE01
    {  48, "T"s }, {  13, "M"s }, {  56, "B"s }, // WE20
    {  94, "T"s }, {  77, "M"s }, {  72, "B"s }, // EW01
    {  53, "T"s }, {  82, "M"s }, {  35, "B"s }, // EW20
    {  96, "T"s }, {  28, "M"s }, {  16, "B"s }, // EE01
    {  86, "T"s }, {  54, "M"s }, {   8, "B"s }  // EE20
  }; // flangeIDsuffixes[]
  
  auto const flangeSuffix = [](int flangeID) -> std::string
    {
      auto const it = flangeIDsuffixes.find(flangeID);
      return (it == flangeIDsuffixes.end())? std::string{}: it->second;
    };

  // Loop through the data to recover the channels;
  // find the position of the columns we need from the first row
  /*
   * [20240224] 9 columns:
   *  [0] "readout_board_id" [1] "flange_id"   [2] "chimney_number"
   *  [3] "tpc_id"           [4] "create_time" [5] "create_user"
   *  [6] "update_time"      [7] "update_user" [8] "fragement_id"
   */
  // technical detail: `PositionFinder` uses `operator==` to compare, and the
  // input is C-strings; we force the other term of comparison to C++ strings
  auto const [
    ReadoutBoardIDcolumn, FlangeIDcolumn,   FragmentIDcolumn,
    TPCIDcolumn,          ChimneyNoColumn
  ] = details::WDAPositionFinder{ getTuple(dataset, 0) }(
   "readout_board_id"s,  "flange_id"s,     "fragement_id"s,
   "tpc_id"s,            "chimney_number"s
    );
  
  for (int const row: util::counter(1, getNtuples(dataset))) {
    // Recover the row
    details::WDATuple tuple { getTuple(dataset, row) };
    if (!tuple) continue;
    
    int error = 0;
    
    // Note that the fragment ID is stored in the database as a string which
    // reads as a hex number, meaning we have to read back as a string and
    // decode to get the numerical value
    Expected const fragmentIDString
      = getStringFromTuple(tuple, FragmentIDcolumn);
    if (!fragmentIDString) {
      throw myException() << "Error (code: " << fragmentIDString.code()
        << " on row " << row
        << ") retrieving TPC fragment ID from channel mapping database\n";
    }
    
    unsigned int const fragmentID = std::stol(fragmentIDString, nullptr, 16);
    if (!(fragmentID & tpcIdentifier)) continue;
    
    if (fragmentBoardMap.find(fragmentID) == fragmentBoardMap.end()) {
      
      unsigned int const flangeID
        = getLongValue(tuple, FlangeIDcolumn, &error);
      if (error) {
        throw myException() << "Error (code: " << error << " on row " << row
          << ") retrieving TPC flange ID from channel mapping database\n";
      }
      
      // build the flange name
      Expected const chimneyNoString = getStringFromTuple(tuple, ChimneyNoColumn);
      if (!chimneyNoString) {
        throw myException() << "Error (code: " << chimneyNoString.code()
          << " on row " << row
          << ") retrieving chimney number from channel mapping database\n";
      }
      
      Expected const TPCIDstring = getStringFromTuple(tuple, TPCIDcolumn);
      if (!TPCIDstring) {
        throw myException() << "Error (code: " << TPCIDstring.code()
          << " on row " << row
          << ") retrieving TPC tag from channel mapping database\n";
      }
      
      std::string const flangeTag = TPCIDstring.value()
        + std::string(std::max(0, 2-(int)chimneyNoString.value().length()), '0')
        + chimneyNoString.value() + flangeSuffix(flangeID);
      
      fragmentBoardMap[fragmentID].first = flangeTag;
    } // if first time with this flange
    
    unsigned int const readoutID
      = getLongValue(tuple, ReadoutBoardIDcolumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") retrieving TPC readout board ID from channel mapping database\n";
    }
    
    fragmentBoardMap[fragmentID].second.emplace_back(readoutID);
    
  } // for rows

  return 0;
} // icarusDB::ChannelMapPostGres::BuildTPCFragmentIDToReadoutIDMap()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapPostGres::BuildTPCReadoutBoardToChannelMap
  (TPCReadoutBoardToChannelMap& rbChanMap) const
{
  // Recover the data from the database
  std::string const name { "icarus_hardware_prd" };
  std::string const dataType { "daq_channels" };
  details::WDADataset dataset = GetDataset(name, fDBURL, dataType);
  
  // find the position of the columns we need from the first row
  /*
   * [20240224] 13 columns:
   *  [0] "channel_id"      [1] "wire_number"         [2] "readout_board_id" 
   *  [3] "chimney_number"  [4] "readout_board_slot"  [5] "channel_number"
   *  [6] "create_time"     [7] "create_user"         [8] "update_time"
   *  [9] "update_user"    [10] "plane"              [11] "cable_label_number"
   * [12] "channel_type"
   */
  // technical detail: `PositionFinder` uses `operator==` to compare, and the
  // input is C-strings; we force the other term of comparison to C++ strings
  using namespace std::string_literals;
  auto const [
     ChannelIDcolumn,     ReadoutBoardIDcolumn, ReadoutBoardSlotColumn,
     ChannelNumberColumn, PlaneIdentifierColumn
  ]= details::WDAPositionFinder{ getTuple(dataset, 0) }(
    "channel_id"s,       "readout_board_id"s,  "readout_board_slot"s,
    "channel_number"s,   "plane"s
    );
  
  // Loop through the data to recover the channels,
  // making sure to skip the first (header) row
  for(int const row: util::counter(1, getNtuples(dataset))) {
    // Recover the row
    details::WDATuple tuple { getTuple(dataset, row) };
    if (!tuple) continue;
    int error = 0;
    
    unsigned int const readoutBoardID
      = getLongValue(tuple, ReadoutBoardIDcolumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading TPC readout board ID\n";
    }
    
    if (rbChanMap.find(readoutBoardID) == rbChanMap.end()) {
      unsigned int const readoutBoardSlot
        = getLongValue(tuple, ReadoutBoardSlotColumn, &error);
      
      if (error) {
        throw myException() << "Error (code: " << error << " on row " << row
          << ") reading TPC readout board slot\n";
      }
      
      rbChanMap[readoutBoardID].first = readoutBoardSlot;
      rbChanMap[readoutBoardID].second.resize(ChannelsPerTPCreadoutBoard);
    } // if new board ID
    
    unsigned int const channelNum
      = getLongValue(tuple, ChannelNumberColumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading TPC channel number\n";
    }
    
    unsigned int const channelID = getLongValue(tuple, ChannelIDcolumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading TPC channel ID\n";
    }
    
    // Recover the plane identifier
    Expected const fragmentBuffer
      = getStringFromTuple(tuple, PlaneIdentifierColumn);
    if (!fragmentBuffer) {
      throw myException() << "Error (code: " << fragmentBuffer.code()
        << " on row " << row << ") reading plane type\n";
    }
    // Make sure lower case... (sigh...)
    unsigned int const plane = TPCplaneIdentifierToPlane(fragmentBuffer);
    if (plane >= 3) {
      mf::LogError{ "ChannelMapSQLite" } << "YIKES!!! Plane is " << plane
        << " for channel " << channelID << " with type "
        << fragmentBuffer.value();
    }
    
    rbChanMap[readoutBoardID].second[channelNum]
      = std::make_pair(channelID, plane);
  } // for rows

  return 0; // on error, an exception was thrown
  
} // icarusDB::ChannelMapPostGres::BuildTPCReadoutBoardToChannelMap()


// -----------------------------------------------------------------------------
// ---  PMT
// -----------------------------------------------------------------------------
int icarusDB::ChannelMapPostGres::BuildPMTFragmentToDigitizerChannelMap
  (PMTFragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap) const
{
  assert( !fCurrentPMTTimestamp.empty() );
  
  // clearing is cleansing
  fragmentToDigitizerChannelMap.clear();
  
  // Recover the information from the database on the mapping
  std::string const name { "Pmt_placement" };
  std::string const dataType { "pmt_placements" };
  std::string const period_query
    { "&w=period_active:eq:" + fCurrentPMTTimestamp };
  details::WDADataset dataset
    = GetDataset(name, fDBURL + period_query, dataType);
  
  // find the position of the columns we need from the first row
  /*
   * [20240224] 20 columns:
   *  [0] "pmt_id"               [1] "period_active"      [2] "pmt_in_tpc_plane"
   *  [3] "channel_id"           [4] "pmt_sn"             [5] "sector_label" 
   *  [6] "ch_number"            [7] "pmt_position_code"  [8] "hv_cable_label"
   *  [9] "signal_cable_label"  [10] "light_fiber_label" [11] "digitizer_label"
   * [12] "digitizer_ch_number" [13] "hv_supply_label" [14]"hv_supply_ch_number"
   * [15] "fragment_id"         [16] "create_time"       [17] "update_user"
   * [18] "update_time"         [19] "create_user"
   */
  // technical detail: `PositionFinder` uses `operator==` to compare, and the
  // input is C-strings; we force the other term of comparison to C++ strings
  using namespace std::string_literals;
  auto const [
     LaserChannelColumn,     DigitizerChannelColumn, ChannelIDcolumn,
     FragmentIDcolumn,       DigitizerLabelColumn,   LVDSconnectorColumn,
     AdderConnectorColumn
  ]= details::WDAPositionFinder{ getTuple(dataset, 0) }(
    "light_fiber_label"s,   "digitizer_ch_number"s, "channel_id"s,
    "fragment_id"s,         "digitizer_label"s,     "FPGA_connector_DIO"s,
    "adder_connector_DIO"s
    );
  
  // Ok, now we can start extracting the information
  // We do this by looping through the database and building the map from that
  // NOTE that we skip the first row because that is just the labels
  for (int const row: util::counter(1, getNtuples(dataset))) {
    // Recover the row
    details::WDATuple tuple { getTuple(dataset, row) };
    if (!tuple) continue;
    
    int error = 0;
    
    PMTChannelInfo_t chInfo;
    auto const nFields = static_cast<std::size_t>(getNfields(tuple));
    
    // digitizer label
    Expected const digitizerLabel
      = getStringFromTuple(tuple, DigitizerLabelColumn);
    if (!digitizerLabel) {
      throw myException() << "Error (code: " << digitizerLabel.code()
        << " on row " << row
        << ") retrieving PMT digitizer from channel mapping database\n";
    }
    chInfo.digitizerLabel = digitizerLabel;
    
    // fragment id
    unsigned long const fragmentID
      = getLongValue(tuple, FragmentIDcolumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading PMT fragment ID\n";
    }
    
    // digitizer channel number
    chInfo.digitizerChannelNo
      = getLongValue(tuple, DigitizerChannelColumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading PMT readout board channel number\n";
    }
    
    // LArsoft channel ID
    chInfo.channelID = getLongValue(tuple, ChannelIDcolumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading PMT channel ID\n";
    }
    // laser channel number
    Expected laserChannelLabel = getStringFromTuple(tuple, LaserChannelColumn);
    if (!laserChannelLabel) {
      throw myException() << "Error (code: " << laserChannelLabel.code()
        << " on row " << row
        << ") retrieving PMT laser channel from channel mapping database\n";
    }
    // will throw on error:
    chInfo.laserChannelNo = std::stol(laserChannelLabel.value().substr(2));
    
    // LVDS connector and bit
    if (LVDSconnectorColumn < nFields) {
      Expected LVDSconnectorLabel
        = getStringFromTuple(tuple, LVDSconnectorColumn);
      if (!LVDSconnectorLabel) {
        throw myException() << "Error (code: " << LVDSconnectorLabel.code()
          << " on row " << row
          << ") retrieving LVDS connector from channel mapping database\n";
      }
      auto const [ LVDSconnector, LVDSbit ]
        = splitIntegers<2, unsigned short int>(LVDSconnectorLabel, "-");
      chInfo.LVDSconnector = LVDSconnector;
      chInfo.LVDSbit = LVDSbit;
    }
    
    // adder connector and bit
    if (AdderConnectorColumn < nFields) {
      Expected adderConnectorLabel
        = getStringFromTuple(tuple, AdderConnectorColumn);
      if (!adderConnectorLabel) {
        throw myException() << "Error (code: " << adderConnectorLabel.code()
          << " on row " << row
          << ") retrieving adder connector from channel mapping database\n";
      }
      auto const [ adderConnector, adderBit ]
        = splitIntegers<2, unsigned short int>(adderConnectorLabel, "-");
      chInfo.adderConnector = adderConnector;
      chInfo.adderBit = adderBit;
    }
    
    // fill the map
    fragmentToDigitizerChannelMap[fragmentID].push_back(std::move(chInfo));
    
  } // for
  
  return 0;
} // icarusDB::ChannelMapPostGres::BuildPMTFragmentToDigitizerChannelMap()


// -----------------------------------------------------------------------------
// ---  CRT
// -----------------------------------------------------------------------------
int icarusDB::ChannelMapPostGres::BuildCRTChannelIDToHWtoSimMacAddressPairMap
  (CRTChannelIDToHWtoSimMacAddressPairMap& crtChannelIDToHWtoSimMacAddressPairMap)
  const
{
  // clearing is cleansing
  crtChannelIDToHWtoSimMacAddressPairMap.clear();
  
  // Recover the information from the database on the mapping
  std::string const name { "Feb_channels" };
  std::string const dataType { "feb_channels" };
  details::WDADataset dataset = GetDataset(name, fDBURL, dataType);
  
  // find the position of the columns we need from the first row
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
  auto const [ ChannelIDcolumn, SimMacAddressColumn, HWMacAddressColumn ]
    = details::WDAPositionFinder{ getTuple(dataset, 0) }
             ("channel_id"s,   "feb_index"s,        "mac_address"s      );
  
  // Ok, now we can start extracting the information
  // We do this by looping through the database and building the map from that
  for(int const row: util::counter(1, getNtuples(dataset))) {
    // Recover the row
    details::WDATuple tuple { getTuple(dataset, row) };
    if (!tuple) continue;
      
    int error = 0;
    
	  // Recover the simmacaddress
    unsigned int const simmacaddress
      = getLongValue(tuple, SimMacAddressColumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading side CRT SimMac address\n";
    }

    // Now recover the hwmacaddress
    unsigned int const hwmacaddress
      = getLongValue(tuple, HWMacAddressColumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading side CRT hardware Mac address\n";
    }
    
    // Finally, get the LArsoft channel ID
    unsigned int const channelID = getLongValue(tuple, ChannelIDcolumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading side CRT channel ID\n";
    }
    
    // Fill the map
    crtChannelIDToHWtoSimMacAddressPairMap[channelID]
      = std::make_pair(hwmacaddress,simmacaddress);
    
  } // for

  return 0;
} // icarusDB::ChannelMapPostGres::BuildCRTChannelIDToHWtoSimMacAddressPairMap()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapPostGres::BuildTopCRTHWtoSimMacAddressPairMap
  (TopCRTHWtoSimMacAddressPairMap& topcrtHWtoSimMacAddressPairMap) const
{
  // clearing is cleansing
  topcrtHWtoSimMacAddressPairMap.clear();
  
  // Recover the information from the database on the mapping
  std::string const name { "topcrt_febs" };
  std::string const dataType { "crtfeb" };

  details::WDADataset dataset = GetDataset(name, fDBURL, dataType);
  
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
  auto const [ SimMacAddressColumn, HWMacAddressColumn ]
    = details::WDAPositionFinder{ getTuple(dataset, 0) }
             ("feb_index"s,        "mac_add"s      );
  
  // Ok, now we can start extracting the information
  // We do this by looping through the database and building the map from that
  for(int const row: util::counter(1, getNtuples(dataset))) {
    // Recover the row
    details::WDATuple tuple { getTuple(dataset, row) };
    if (!tuple) continue;
    int error = 0;

    // Recover the simmacaddress
    unsigned int const simmacaddress
      = getLongValue(tuple, SimMacAddressColumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading top CRT SimMac address\n";
    }
    
    // Now recover the hwmacaddress
    unsigned int const hwmacaddress
      = getLongValue(tuple, HWMacAddressColumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading top CRT hardware Mac address\n";
    }

    // Fill the map
    topcrtHWtoSimMacAddressPairMap[hwmacaddress] = simmacaddress;
  } // for

  return 0;
} // icarusDB::ChannelMapPostGres::BuildTopCRTHWtoSimMacAddressPairMap()


// -----------------------------------------------------------------------------
int icarusDB::ChannelMapPostGres::BuildSideCRTCalibrationMap
  (SideCRTChannelToCalibrationMap& sideCRTChannelToCalibrationMap) const
{

  std::string const name { "SideCRT_calibration_data" };
  
  details::WDADataset dataset = GetCRTCaldata(name, fCRTcalibrationDBURL);
  printDatasetError(dataset); // probably never prints since on error throws
  
  // the shape of this database is different from the others,
  // there is no header row and there is some "introductory" data
  // (which we'll skip) in the first rows

  for (int const row: util::counter(getNtuples(dataset))) {
    
    // Get the row with an array of double
    details::WDATuple tuple { getTuple(dataset, row) };
    if (!tuple) continue;
    
    int const ncols = getNfields(tuple);//check number of columns
    
    // first few rows aren't actual data and have ncols==1: exclude those
    if (ncols < 5) continue;
    
    // assign values from the database to variables so we can use them
    int error = 0;
    
    int const mac5 = static_cast<int>(getDoubleValue(tuple, 1, &error));
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading CRT calibration Mac5 address\n";
    }
    
    int const channel = static_cast<int>(getDoubleValue(tuple, 2, &error));
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading CRT calibration channel\n";
    }
    
    double const gain = getDoubleValue(tuple, 3, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading CRT calibration gain for channel " << channel << "\n";
    }
    
    double const ped = getDoubleValue(tuple, 4, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading CRT calibration pedestal for channel " << channel << "\n";
    }
    
    // add the association between the two pairs to the map object
    sideCRTChannelToCalibrationMap.emplace
      (std::make_pair(mac5, channel), std::make_pair(gain, ped));

  } // for rows

  return 0;
  
} // icarusDB::ChannelMapPostGres::BuildSideCRTCalibrationMap()


// -----------------------------------------------------------------------------
bool icarusDB::ChannelMapPostGres::printDatasetError
  (details::WDADataset const& dataset) const
{
  int const status = getHTTPstatus(dataset);
  if (status == 200) return false;
  
  mfLogError() << "libwda error: HTTP code=" << status << ": '"
    << getHTTPmessage(dataset) << "'";
  return true;
}


// -----------------------------------------------------------------------------
template <std::size_t BufferSize /* = 32 */>
auto icarusDB::ChannelMapPostGres::getStringFromTuple
  (details::WDATuple const& tuple, std::size_t column) -> Expected<std::string>
{
  int error = 0;
  std::string buffer(BufferSize, '\0');
  std::size_t const length
    = getStringValue(tuple, column, buffer.data(), buffer.size(), &error);
  if (error) return error;
  buffer.erase(length);
  return buffer;
} // icarusDB::ChannelMapPostGres::getStringFromTuple<>()


// -----------------------------------------------------------------------------
