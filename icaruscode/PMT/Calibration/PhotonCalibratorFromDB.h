/**
 * @file   icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h
 * @brief  Implementation of optical hit photoelectron calibration from database.
 * @author Michael Carrigan, Matteo Vicenzi
 * @see    icaruscode/PMT/Calibration/PhotonCalibratorFromDB.cxx
 */

#ifndef ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H
#define ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H

// LArSoft libraries
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larevt/CalibrationDBI/Providers/DBFolder.h"

// ART includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TableAs.h"

// C++ standard libraries
#include <cstdint> // std::uint64_t
#include <map>
#include <string>


namespace icarusDB::details {
    
  /// Structure for single channel calibration.
  struct PhotonCalibratorInfo {

    /// Area [positive, ADC count sum] of response to single photoelectron.
    double speArea = -1.0;
    
    double speAreaErr = -1.0;  /// Uncertainty on `speArea` [ADC count sum]
    double speFitWidth = -1.0;
    double speFitWidthErr = -1.0;
    
  };
  
} // icarusDB::details

// -----------------------------------------------------------------------------
namespace icarusDB { class PhotonCalibratorFromDB; }
/**
 * @brief Optical hit photoelectron calibration service with data from database.
 * 
 * This service implements the `IPhotonCalibrator` interface.
 * 
 * The calibration factors are loaded from a database, per PMT channel and per
 * calibration period.
 * They represent the "area" of a single photoelectron response, that is the
 * integral of the response signal in time, divided by the sampling time.
 * The typical way to use this provider to get the hit photoelectrons is
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * double hitPhotoelectrons
 *   (calib::IPhotonCalibrator const& calibrator, recob::OpHit const& hit)
 * {
 *   assert(calibrator.UseArea());
 *   return calibrator.PE(hit.Area(), hit.OpChannel());
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * The example highlights that this service provider operates on the hit area.
 * Note however that typically the `recob::OpHit` are already calibrated with
 * an implementation of `calib::IPhotonCalibrator` and directly using
 * `recob::OpHit::PE()` is better.
 * 
 * The database interface is accessed only on `readCalibrationFromDB()` calls,
 * and the relevant information is cached.
 * 
 * The service provider is thread-safe in the standard C++ way: `const` methods
 * are safe to call in any thread, non-const methods are not.
 * 
 */
class icarusDB::PhotonCalibratorFromDB: public calib::IPhotonCalibrator {
  
    
  public:

    /// FHiCL configuration of the calibrator.
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      struct DefaultCalib {
        fhicl::Atom<double> SPEArea {
          Name{ "SPEArea" },
          Comment{ "area (ADC sum) of the response to a single photoelectron" }
          };
      }; // DefaultCalib
      
      fhicl::Atom<std::string> DBname {
        Name{ "DBname" },
        Comment{ "the SPE area database name" },
        "pmt_speareas_data"
        };
      
      fhicl::Atom<bool> Verbose {
        Name{ "Verbose" },
        Comment{ "enable additional messages for debugging" },
        false
        };
      
      fhicl::Atom<std::string> LogCategory {
        Name{ "LogCategory" },
        Comment{ "name of the message stream where to send console output" },
        "PhotonCalibratorFromDB"
        };
      
      fhicl::Atom<std::string> AreaTag {
        Name{ "AreaTag" },
        Comment{ "the database version (tag) to use" },
        ""
        };
      
      fhicl::TableAs<details::PhotonCalibratorInfo, DefaultCalib> Defaults {
        Name{ "Defaults" },
        Comment{ "values used for channels not present in the database" }
        };
      
    }; // Config
  
  
    /// Constructor: reads the FHiCL configuration (no access to database yet).
    PhotonCalibratorFromDB(const Config& config);

    /**
     * @brief Converts the specified value in ADC into photoelectrons.
     * @param adcs area of the hit (see `UseArea()`)
     * @return photoelectrons corresponding to the specified `adcs` area
     */
    double PE(double adcs, int channel) const override;

    /**
     * @brief Whether calibration parameter is area or peak amplitude.
     * @return `false` (this calibration is amplitude-based)
     */
    bool UseArea() const override;

    /// Prepares the calibration for data from the specified `run`.
    void readCalibrationFromDB(unsigned int run);

    /**
     * @brief Returns the database timestamp-like tag appropriate to the `run`.
     * 
     * The backend (`lariov::DBFolder`) only takes 19-digit numbers
     * (timestamp in nanoseconds), but our database tables are currently using
     * run numbers, baked so that a run number `XXXXX` results into the
     * timestamp `1'000'0XX'XXX'000'000'000`.
     */
    std::uint64_t RunToDatabaseTimestamp( unsigned int run ) const;

  private:

    using PhotonCalibratorInfo = details::PhotonCalibratorInfo;

    PhotonCalibratorInfo const fCalibDefaults; ///< Default calibration values.
    bool const fVerbose;
    std::string const fLogCategory;

    lariov::DBFolder fDB; ///< Cached database interface.
    
    /// Map: channel to calibration information.
    std::map<int, PhotonCalibratorInfo> fDatabaseSPECalibrations;

    /// Internal access to the channel calibration record; returns defaults if not present.
    PhotonCalibratorInfo const& getChannelCalibOrDefault(int channelID) const{
      auto const it = fDatabaseSPECalibrations.find(channelID);
      return (it == fDatabaseSPECalibrations.end())? fCalibDefaults: it->second;
    }


}; // class icarus::PhotonCalibratorStandard


#endif // ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H
