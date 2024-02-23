/**
 * @file   icaruscode/Decode/ChannelMapping/ChannelMapPostGres.cxx
 * @brief  Interface with ICARUS channel mapping PostGres database.
 * @author T. Usher (factorised by Gianluca Petrillo, petrillo@slac.stanford.edu)
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
  
} // local namespace


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
int icarusDB::ChannelMapPostGres::BuildTPCFragmentIDToReadoutIDMap
  (TPCFragmentIDToReadoutIDMap& fragmentBoardMap) const
{
  const unsigned int tpcIdentifier(0x00001000);
  
  // Recover the data from the database
  std::string const name { "icarus_hw_readoutboard" };
  std::string const dataType { "readout_boards" };
  details::WDADataset dataset = GetDataset(name, fDBURL, dataType);
  
  // Include a by hand mapping of fragment ID to crate
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
  
  // Loop through the data to recover the channels
  // NOTE that we skip the first row because that is just the labels
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
      fragmentBoardMap[fragmentID].first = flangeIDToCrateMap.at(flangeID);
    }
    
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
  
  constexpr std::size_t ChannelIDcolumn        =  0;
  constexpr std::size_t ReadoutBoardIDcolumn   =  2;
  constexpr std::size_t ReadoutBoardSlotColumn =  4;
  constexpr std::size_t ChannelNumberColumn    =  5;
  constexpr std::size_t PlaneIdentifierColumn  = 10;
  
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
int icarusDB::ChannelMapPostGres::BuildFragmentToDigitizerChannelMap
  (FragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap) const
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

  // Ok, now we can start extracting the information
  // We do this by looping through the database and building the map from that
  // NOTE that we skip the first row because that is just the labels
  for (int const row: util::counter(1, getNtuples(dataset))) {
    // Recover the row
    details::WDATuple tuple { getTuple(dataset, row) };
    if (!tuple) continue;
    
    constexpr std::size_t ChannelIDcolumn          =  3;
    constexpr std::size_t LaserChannelColumn       = 10;
    constexpr std::size_t DigitizerColumn          = 11;
    constexpr std::size_t DigitizerChannelNoColumn = 12;
    constexpr std::size_t FragmentIDcolumn         = 15;
    
    int error = 0;
    
    // nice... and currently unused
    Expected const digitizerLabel = getStringFromTuple(tuple, DigitizerColumn);
    if (!digitizerLabel) {
      throw myException() << "Error (code: " << digitizerLabel.code()
        << " on row " << row
        << ") retrieving PMT digitizer from channel mapping database\n";
    }
    
    // fragment id
    unsigned long const fragmentID
      = getLongValue(tuple, FragmentIDcolumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading PMT fragment ID\n";
    }
    
    // digitizer channel number
    unsigned int const digitizerChannelNo
      = getLongValue(tuple, DigitizerChannelNoColumn, &error);
    if (error) {
      throw myException() << "Error (code: " << error << " on row " << row
        << ") reading PMT readout board channel number\n";
    }
    
    // LArsoft channel ID
    unsigned int const channelID = getLongValue(tuple, ChannelIDcolumn, &error);
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
    unsigned int const laserChannel
      = std::stol(laserChannelLabel.value().substr(2));
    
    // fill the map
    fragmentToDigitizerChannelMap[fragmentID].emplace_back
      (digitizerChannelNo, channelID, laserChannel);
    
  } // for
  
  return 0;
} // icarusDB::ChannelMapPostGres::BuildFragmentToDigitizerChannelMap()


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
  
  constexpr std::size_t ChannelIDcolumn     = 10;
  constexpr std::size_t SimMacAddressColumn = 11;
  constexpr std::size_t HWMacAddressColumn  = 12;
  
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
  
  constexpr std::size_t SimMacAddressColumn = 41;
  constexpr std::size_t HWMacAddressColumn  =  3;

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
