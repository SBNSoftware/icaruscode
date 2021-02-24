/**
 * @file   icaruscode/Utilities/ReadArtConfiguration.h
 * @brief  Utilities to extract _art_ FHiCL configuration from different sources.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 18, 2021
 * @see    icaruscode/Utilities/ReadArtConfiguration.cxx
 * 
 */

#ifndef ICARUSCODE_UTILITIES_READARTCONFIGURATION_H
#define ICARUSCODE_UTILITIES_READARTCONFIGURATION_H


// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetID.h"
#include "fhiclcpp/ParameterSetRegistry.h" // also defines ParameterSetID hash
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TFile.h"

// C/C++ standard libraries
#include <map>


// -----------------------------------------------------------------------------
namespace util {
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Reads and returns the _art_ configuration stored in `sourceDir`.
   * @param file ROOT file where the configuration is stored
   * @return the full configuration
   * @throw cet::exception (category: `"readConfigurationFromArtFile"`) on error
   * 
   * The configuration is expected to be stored by _art_ in the way it does
   * for _art_ ROOT files.
   * 
   * The configuration is returned as a map of parameter set ID to parameter
   * set.
   */
  std::map<fhicl::ParameterSetID, fhicl::ParameterSet>
    readConfigurationFromArtFile(TFile& file);
  
  
  // ---------------------------------------------------------------------------
  
} // namespace util

// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_READARTCONFIGURATION_H
