/**
 * @file   icaruscode/Timing/ICARUSPMTTimingCorrections_service.cc
 * @brief  Service for the PMT timing corrections.
 * @author Andrea Scarpelli (ascarpell@bnl.gov)
 */

#ifndef ICARUSCODE_TIMING_PMTIMINGCORRECTIONSPROVIDER_H
#define ICARUSCODE_TIMING_PMTIMINGCORRECTIONSPROVIDER_H

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "icaruscode/Timing/PMTTimingCorrections.h"

// Database interface helpers
#include "wda.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>

namespace icarusDB{ class PMTTimingCorrectionsProvider; }
/**
 * @brief 
 * 
 * This module reads 
 * 
 * Input
 * ------
 * 
 * 
 * Output
 * -------
 * 
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * 
 */
class icarusDB::PMTTimingCorrectionsProvider : public PMTTimingCorrections {

    public: 

        PMTTimingCorrectionsProvider(const fhicl::ParameterSet& pset);

        void readTimeCorrectionDatabase(const art::Run& run);

        double getTriggerCableDelay( const unsigned int & channelID ) {
            return fDatabaseTimingCorrections[channelID].triggerCableDelay;
        };

        double getResetCableDelay( const unsigned int & channelID ) {
            return fDatabaseTimingCorrections[channelID].resetCableDelay;
        };

        double getLaserCorrections( const unsigned int & channelID ) {
            return fDatabaseTimingCorrections[channelID].laserCableDelay;
        };

        double getCosmicsCorrections( const unsigned int & channelID ) {
            return fDatabaseTimingCorrections[channelID].cosmicsCorrections;
        };

    private:

        std::string fUrl;

        unsigned int fTimeout;

        bool fCorrectCablesDelay;

        bool fVerbose = false; ///< Whether to print the configuration we read.
  
        std::string fLogCategory; ///< Category tag for messages.

        /// Interface to LArSoft configuration for detector timing.
        detinfo::DetectorClocksData const fClocksData;

        struct PMTTimeCorrectionsDB {

            double triggerCableDelay=0;

            double resetCableDelay=0;

            double laserCableDelay=0;

            double cosmicsCorrections=0;
        };

        std::map<unsigned int, PMTTimeCorrectionsDB> fDatabaseTimingCorrections;

        int ConnectToDataset(const std::string& name, 
            const uint32_t &run, Dataset& dataset ) const;

        void ReadPMTCablesCorrections(const uint32_t & run);

        void ReadLaserCorrections(const uint32_t & run);

        void ReadCosmicsCorrections(const uint32_t & run);

}; // services class

#endif 