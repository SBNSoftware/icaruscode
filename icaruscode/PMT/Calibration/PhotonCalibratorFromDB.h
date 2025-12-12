/**
 * @file   icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h
 * @brief  Implementation of optical hit photoelectron calibration from database.
 * @author Michael Carrigan, Matteo Vicenzi
 */

#ifndef ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H
#define ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H

// LArSoft libraries
#include "larreco/Calibrator/IPhotonCalibrator.h"

// Database interface helpers
#include "larevt/CalibrationDBI/Providers/DBFolder.h"
#include "larevt/CalibrationDBI/IOVData/TimeStampDecoder.h"

// Message facility
#include "messagefacility/MessageLogger/MessageLogger.h"

// ART includes
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"

namespace icarusDB::details {
    
  /// Structure for single channel corrections
  struct PhotonCalibratorInfo { 

    double speArea = 256.658;
    double speAreaErr = -1.0;
    double speFitWidth = -1.0;
    double speFitWidthErr = -1.0;
    
  };
  
} // icarusDB::details

// -----------------------------------------------------------------------------
namespace icarusDB { class PhotonCalibratorFromDB; }
/**
* @brief Optical hit photoelectron calibration service with data from database.
* 
* More description here.
* 
*/
class icarusDB::PhotonCalibratorFromDB: public calib::IPhotonCalibrator {
  
    
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      // FHiCL configuration here
      
    };
  
    using Parameters = fhicl::Table<Config>;
  
    PhotonCalibratorFromDB(const fhicl::ParameterSet& pset);

    /// Convert the specified value in ADC into photoelectrons.
    double PE(double adcs, int channel) const override; 

    bool UseArea() const override;

    void readCalibrationFromDB(unsigned int run);

    uint64_t RunToDatabaseTimestamp( uint32_t run ) const;

  private:

    bool fVerbose;
    std::string fLogCategory;
    std::string fAreaTag;

    using PhotonCalibratorInfo = details::PhotonCalibratorInfo;
    static constexpr PhotonCalibratorInfo CorrectionDefaults {}; ///< Default values

    std::map<unsigned int, PhotonCalibratorInfo> fDatabaseSPECalibrations;

    /// Internal access to the channel correction record; returns defaults if not present.
    PhotonCalibratorInfo const& getChannelCorrOrDefault(unsigned int channelID) const{
      auto const it = fDatabaseSPECalibrations.find(channelID);
      return (it == fDatabaseSPECalibrations.end())? CorrectionDefaults: it->second;
    }


}; // class icarus::PhotonCalibratorStandard


#endif // ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H
