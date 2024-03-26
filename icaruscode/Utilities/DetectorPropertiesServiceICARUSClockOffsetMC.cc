////////////////////////////////////////////////////////////////////////
//
//  \file DetectorProperties_service.cc
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// icaruscode includes
#include "icaruscode/Utilities/DetectorPropertiesServiceICARUSClockOffsetMC.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/ServicePack.h" // lar::extractProviders()
#include "messagefacility/MessageLogger/MessageLogger.h"

// Art includes
#include "art_root_io/RootDB/SQLite3Wrapper.h"

#include "TFile.h"
#include "TTree.h"

namespace detinfo {

  //--------------------------------------------------------------------
  DetectorPropertiesServiceICARUSClockOffsetMC::DetectorPropertiesServiceICARUSClockOffsetMC(
    fhicl::ParameterSet const& pset,
    art::ActivityRegistry& reg)
    : fProp{pset,
            lar::providerFrom<geo::Geometry>(),
            lar::providerFrom<detinfo::LArPropertiesService>(),
            std::set<std::string>({"InheritNumberTimeSamples"})}
    , fPS{pset}
    , fInheritNumberTimeSamples{pset.get<bool>("InheritNumberTimeSamples", false)}
  {
    reg.sPostOpenFile.watch(this, &DetectorPropertiesServiceICARUSClockOffsetMC::postOpenFile);
  }

  //--------------------------------------------------------------------
  //  Callback called after input file is opened.

  void DetectorPropertiesServiceICARUSClockOffsetMC::postOpenFile(const std::string& filename)
  {
    // Use this method to figure out whether to inherit configuration
    // parameters from previous jobs.

    // There is no way currently to correlate parameter sets saved in
    // sqlite RootFileDB with process history (from MetaData tree).
    // Therefore, we use the approach of scanning every historical
    // parameter set in RootFileDB, and finding all parameter sets
    // that appear to be DetectorPropertiesService configurations.  If all
    // historical parameter sets are in agreement about the value of
    // an inherited parameter, then we accept the historical value,
    // print a message, and override the configuration parameter.  In
    // cases where the historical configurations are not in agreement
    // about the value of an inherited parameter, we ignore any
    // historical parameter values that are the same as the current
    // configured value of the parameter (that is, we resolve the
    // conflict in favor of parameters values that are different than
    // the current configuration).  If two or more historical values
    // differ from the current configuration, throw an exception.
    // Note that it is possible to give precendence to the current
    // configuration by disabling inheritance for that configuration
    // parameter.

    // Don't do anything if no parameters are supposed to be inherited.

    if (!fInheritNumberTimeSamples) return;

    // The only way to access art service metadata from the input file
    // is to open it as a separate TFile object.  Do that now.

    if (filename.empty()) { return; }

    std::unique_ptr<TFile> file{TFile::Open(filename.c_str(), "READ")};
    if (!file) { return; }

    if (file->IsZombie() || !file->IsOpen()) { return; }

    // Open the sqlite datatabase.

    art::SQLite3Wrapper sqliteDB(file.get(), "RootFileDB");

    // Loop over all stored ParameterSets.

    unsigned int iNumberTimeSamples = 0; // Combined value of NumberTimeSamples.
    unsigned int nNumberTimeSamples = 0; // Number of NumberTimeSamples parameters seen.

    sqlite3_stmt* stmt = nullptr;
    sqlite3_prepare_v2(sqliteDB, "SELECT PSetBlob from ParameterSets;", -1, &stmt, nullptr);
    while (sqlite3_step(stmt) == SQLITE_ROW) {
      fhicl::ParameterSet ps;
      ps = fhicl::ParameterSet::make(reinterpret_cast<char const*>(sqlite3_column_text(stmt, 0)));
      // Is this a DetectorPropertiesService parameter set?

      if (isDetectorPropertiesServiceICARUSClockOffsetMC(ps)) {

        // Check NumberTimeSamples

        auto const newNumberTimeSamples = ps.get<unsigned int>("NumberTimeSamples");

        // Ignore parameter values that match the current configuration.

        if (newNumberTimeSamples != fPS.get<unsigned int>("NumberTimeSamples")) {
          if (nNumberTimeSamples == 0)
            iNumberTimeSamples = newNumberTimeSamples;
          else if (newNumberTimeSamples != iNumberTimeSamples) {
            throw cet::exception(__FUNCTION__)
              << "Historical values of NumberTimeSamples do not agree: " << iNumberTimeSamples
              << " " << newNumberTimeSamples << "\n";
          }
          ++nNumberTimeSamples;
        }
      }
    }

    // Done looping over parameter sets.
    // Now decide which parameters we will actually override.

    if (nNumberTimeSamples != 0 && iNumberTimeSamples != fProp.NumberTimeSamples()) {
      mf::LogInfo("DetectorPropertiesServiceStandard")
        << "Overriding configuration parameter NumberTimeSamples using "
           "historical value.\n"
        << "  Configured value:        " << fProp.NumberTimeSamples() << "\n"
        << "  Historical (used) value: " << iNumberTimeSamples << "\n";
      fProp.SetNumberTimeSamples(iNumberTimeSamples);
    }
  }

  //--------------------------------------------------------------------
  //  Determine whether a parameter set is a DetectorPropertiesService configuration.

  bool DetectorPropertiesServiceICARUSClockOffsetMC::isDetectorPropertiesServiceICARUSClockOffsetMC(
    const fhicl::ParameterSet& ps) const
  {
    // This method uses heuristics to determine whether the parameter
    // set passed as argument is a DetectorPropertiesService configuration
    // parameter set.

    return (ps.get<std::string>("service_type", "") == "DetectorPropertiesService") &&
           (ps.get<std::string>("service_provider", "") == "DetectorPropertiesServiceICARUSClockOffsetMC");
  }

} // namespace detinfo
