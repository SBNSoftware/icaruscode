/**
 * @file   icaruscode/Decode/ChannelMapping/ChannelMapPostGres.h
 * @brief  Interface with ICARUS channel mapping PostGres database.
 * @author T. Usher (factorised by Gianluca Petrillo, petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/ChannelMapping/ChannelMapPostGres.cxx
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_CHANNELMAPPOSTGRES_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_CHANNELMAPPOSTGRES_H


// ICARUS libraries
#include "icaruscode/Decode/ChannelMapping/IChannelMapping.h"
#include "icaruscode/Decode/ChannelMapping/RunPeriods.h"
#include "icarusalg/Utilities/mfLoggingClass.h"

// framework libraries
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <variant>
#include <string>
#include <cstdint> // std::uint32_t


namespace icarusDB {
  class ChannelMapPostGres;
  
  namespace details {
    struct WDADataset;
    struct WDATuple;
  }
  
} // namespace icarusDB

/**
 * @brief Helper to extract mapping information from ICARUS PostGres database.
 * 
 * This helper fills "standard" data structures from a PostGres database.
 * These data structures are then managed by a ICARUS service provider
 * (`icarusDB::IICARUSChannelMap`).
 */
class icarusDB::ChannelMapPostGres
  : public IChannelMapping, private icarus::ns::util::mfLoggingClass
{
    public:
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<std::string> DatabaseURL {
      Name{ "DatabaseURL" },
      Comment{ "full URL of the database to be accessed" },
      R"(https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=icarus_hardware_prd)"
      };
    
    fhicl::Atom<std::string> CRTcalibrationDatabaseURL {
      Name{ "CRTcalibrationDatabaseURL" },
      Comment{ "full URL of the CRT calibration database to be accessed" },
      R"(https://dbdata0vm.fnal.gov:9443/icarus_con_prod/app/data?f=crt_gain_reco_data&t=1638918270)"
      };
    
    fhicl::Atom<unsigned int> DatabaseAccessTimeout {
      Name{ "DatabaseAccessTimeout" },
      Comment{ "timeout for database accesses [s]" },
      200 // default
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of the streams to send console messages to" },
      "ChannelMapPostGres" // default
      };
    
  }; // Config
  
  
  /// Constructor: configuration via FHiCL.
  explicit ChannelMapPostGres(Config const& config);

  /**
   * @brief  Prepares the object for queries pertaining the specified period.
   * @param  period the period to be prepared for
   * @return whether values cached from the previous period are invalidated
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
  virtual int BuildPMTFragmentToDigitizerChannelMap
    (PMTFragmentToDigitizerChannelMap&) const override;


  /**
   *  @brief Define the returned data structures for a mapping between CRT hardware mac_address
   *         to the simulated mac_address.
   *         Then define the function interface to fill these data structures
   */
  virtual int BuildCRTChannelIDToHWtoSimMacAddressPairMap(CRTChannelIDToHWtoSimMacAddressPairMap&) const override;
  virtual int BuildTopCRTHWtoSimMacAddressPairMap(TopCRTHWtoSimMacAddressPairMap&) const override;

  /**
   *  @brief Define the returned data structures for a mapping between Side CRT Channels and
   *         their calibration values.
   *         Then define the function interface to fill these data structures
   */
  virtual int BuildSideCRTCalibrationMap(SideCRTChannelToCalibrationMap&) const override;

    private:

  /// Loosely based on C++-20 `std::expected`.
  template <typename T, typename E = int>
  class Expected {
    std::variant<T, E> fValue;
      public:
    constexpr Expected(E e = E{}): fValue{ std::move(e) } {}
    constexpr Expected(T t): fValue{ std::move(t) } {}
    
    Expected<T, E>& operator= (E e) { fValue = std::move(e); return *this; }
    Expected<T, E>& operator= (T t) { fValue = std::move(t); return *this; }
    template <typename... Args>
    T& emplace(Args&&... args)
      { return fValue.emplace(std::forward<Args>(args)...); }
    
    bool has_value() const { return std::holds_alternative<T>(fValue); }
    
    E code() const { return std::get<E>(fValue); }
    T const& value() const { return std::get<T>(fValue); }
    
    operator T const& () const { return value(); }
    T const& operator*() const { return value(); }
    
    explicit operator bool() const { return has_value(); }
    
  }; // Expected
  
  
  // --- BEGIN --- Configuration parameters ------------------------------------
  
  std::string const fDBURL; ///< Complete URL for hardware database access.
  
  /// Complete URL for CRT calibration database access.
  std::string const fCRTcalibrationDBURL;
  
  unsigned int const fDatabaseAccessTimeout; ///< In seconds.
  
  // --- END ----- Configuration parameters ------------------------------------
  
  
  // --- BEGIN --- Cache -------------------------------------------------------
  
  /// The PMT timestamp being served. Chosen by `SelectPeriod()`.
  std::string fCurrentPMTTimestamp = "";
  
  // --- END ----- Cache -------------------------------------------------------

  /**
   * @brief Reads a full table from the specified database table.
   * @param name name of the database table to read data from
   * @param url URL of the database to be accessed requested
   * @param dataType type of data being requested
   * @param[out] dataSet where the fetch data will be saved
   * @return a "smart" pointer to the data
   * @throw cet::exception (category: `"ChannelMapPostGres"`) on DB access error
   * 
   * Database access is performed via Fermilab's libwda.
   */
  details::WDADataset GetDataset
    (std::string const& name, std::string url, std::string const& dataType)
    const;

  
  /**
   * @brief Reads a full table from the specified database table.
   * @param name name of the database table to read data from
   * @param url URL of the database to be accessed requested
   * @param[out] dataSet where the fetch data will be saved
   * @return a "smart" pointer to the data
   * @throw cet::exception (category: `"ChannelMapPostGres"`) on DB access error
   * @see `GetDataset()`
   * 
   * Function used specifically for accessing CRT calibration data.
   * 
   * @note [Tyler Boone] Adding this so I can control how this workflow goes
   */
  /// Fetches CRT calibration data from the database.
  details::WDADataset GetCRTCaldata
    (std::string const& name, std::string url) const;
  
  
  /// Emits an error message if libwda dataset contains one.
  /// @return whether the dataset was in an error state
  bool printDatasetError(details::WDADataset const& dataset) const;
  
  
  /// The list of available PMT timestamps.
  static std::array<std::string, RunPeriods::NPeriods> const PMTTimestampSet;
  
  /// Returns the string on the specified `column` of the `tuple`.
  /// @return the requested string, or an error code from libwda on error
  template <std::size_t BufferSize = 32>
  static Expected<std::string> getStringFromTuple
    (details::WDATuple const& tuple, std::size_t column);

  /// Returns an exception object pre-approved with our signature.
  static cet::exception myException()
    { return cet::exception{ "ChannelMapPostGres" }; }
  
}; // icarusDB::ChannelMapPostGres


#endif // ICARUSCODE_DECODE_CHANNELMAPPING_CHANNELMAPPOSTGRES_H

