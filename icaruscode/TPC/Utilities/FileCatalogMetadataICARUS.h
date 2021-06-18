////////////////////////////////////////////////////////////////////////
//
// Name:  FileCatalogMetadataICARUS.h.  
//
// Purpose:  Art service adds icarus-specific per-job sam metadata.
//
//           FCL parameters:
//
//           FCLName        - FCL file name.
//           ProjectName    - Project name.
//           ProjectStage   - Project stage.
//           ProjectVersion - Project version.
//           ProductionName - Production name.
//           ProductionType - Production type.
//
//           Above values will be added in internal metadata of artroot
//           output files whenever this service is included in job
//           configuration (service does not need to be called).  The
//           public interface of this service consists of accessors that
//           allow other code to discover above metadata parameters.
//
// Created:  30-Dec-2019,  M. Wospakrik
//  based on the MicroBooNE version by H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#ifndef FILECATALOGMETADATAICARUS_H
#define FILECATALOGMETADATAICARUS_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace icarusutil {

  class FileCatalogMetadataICARUS
  {
  public:

    // Constructor, destructor.

    FileCatalogMetadataICARUS(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~FileCatalogMetadataICARUS() = default;

    // Accessors.

    const std::string& GetFCLName() const {return fFCLName;}
    const std::string& GetProjectName() const {return fProjectName;}
    const std::string& GetProjectStage() const {return fProjectStage;}
    const std::string& GetProjectVersion() const {return fProjectVersion;}
    const std::string& GetProjectSoftware() const {return fProjectSoftware;}
    const std::string& GetProductionName() const {return fProductionName;}
    const std::string& GetProductionType() const {return fProductionType;}
    int Merge() const {return fMerge;}
    const std::vector<std::string>& Parameters() const {return fParameters;}


  private:

    // Callbacks.

    void postBeginJob();

    // Data members.

    std::string fFCLName;
    std::string fProjectName;
    std::string fProjectStage;
    std::string fProjectVersion;
    std::string fProjectSoftware;
    std::string fProductionName; //Production parameter, do not use if not running a production
    std::string fProductionType; //Production parameter, do not use if not running a production
    int fMerge;
    std::vector<std::string> fParameters;
  };

} // namespace icarusutil

DECLARE_ART_SERVICE(icarusutil::FileCatalogMetadataICARUS, LEGACY)

#endif
