/**
 * @file   icaruscode/Utilities/ReadArtConfiguration.cxx
 * @brief  Utilities to extract _art_ FHiCL configuration from different sources.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 18, 2021
 * @see    icaruscode/Utilities/ReadArtConfiguration.h
 * 
 */


// library header
#include "icaruscode/Utilities/ReadArtConfiguration.h"

// framework libraries
// #include "art_root_io/GetFileFormatEra.h"
#include "art_root_io/RootDB/SQLite3Wrapper.h"
// #include "art_root_io/RootDB/tkeyvfs.h"
// #include "canvas/Persistency/Provenance/FileFormatVersion.h"
// #include "canvas/Persistency/Provenance/ParameterSetBlob.h"
// #include "canvas/Persistency/Provenance/ParameterSetMap.h"
// #include "canvas/Persistency/Provenance/rootNames.h"
// #include "cetlib_except/exception.h"
// #include "cetlib/container_algorithms.h"
// #include "cetlib/exempt_ptr.h"
#include "fhiclcpp/ParameterSetRegistry.h"
// #include "fhiclcpp/make_ParameterSet.h"

// ROOT libraries
#include "TTree.h"

// C/C++ standard libraries
// #include <unordered_map>
// #include <vector>
// #include <mutex>
// #include <utility> // std::pair<>
#include <string>
// #include <functional> // std::hash<>
#include <memory> // std::unique_ptr<>
// #include <optional>
// #include <limits> // std::numeric_limits<>
// #include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
std::map<fhicl::ParameterSetID, fhicl::ParameterSet>
util::readConfigurationFromArtFile(TFile& file)
{
  // OMGOMGOMG this is SO copied from art `config_dumper` source!!
  
  // Open the DB
  art::SQLite3Wrapper sqliteDB(&file, "RootFileDB");
  fhicl::ParameterSetRegistry::importFrom(sqliteDB);
  fhicl::ParameterSetRegistry::stageIn();
  
  std::map<fhicl::ParameterSetID, fhicl::ParameterSet> config;
  for (auto const& idAndPSset: fhicl::ParameterSetRegistry::get())
    config.emplace(idAndPSset);
  
  return config;
} // util::readConfigurationFromArtFile()


// -----------------------------------------------------------------------------
