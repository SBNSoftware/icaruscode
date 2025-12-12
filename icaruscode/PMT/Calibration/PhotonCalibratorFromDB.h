/**
 * @file   icaruscode/PMT/Calibration/PhotonCalibratorFromDB.h
 * @brief  Implementation of optical hit photoelectron calibration from database.
 * @author Michael Carrigan, Matteo Vicenzi
 */

#ifndef ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H
#define ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H

// LArSoft libraries
#include "larreco/Calibrator/IPhotonCalibrator.h"

// ART includes
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/ParameterSet.h"

// -----------------------------------------------------------------------------
namespace calib {
  /**
  * @brief Optical hit photoelectron calibration service with data from database.
  * 
  * More description here.
  * 
  */
  class PhotonCalibratorFromDB: public calib::IPhotonCalibrator {
    
      
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

    private:

      bool fVerbose;
      std::string fLogCategory;

  }; // class icarus::PhotonCalibratorStandard

}
#endif // ICARUSCODE_PMT_CALIBRATION_PHOTONCALIBRATORFROMDB_H
