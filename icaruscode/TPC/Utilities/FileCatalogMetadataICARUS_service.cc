////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataICARUS_service.cc.  
//
// Purpose:  Implementation for FileCatalogMetadataICARUS.
//
// Created:  30-Dec-2019,  M. Wospakrik
//  based on the MicroBooNE version by H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#include "icaruscode/TPC/Utilities/FileCatalogMetadataICARUS.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"

//--------------------------------------------------------------------
// Constructor.

icarusutil::FileCatalogMetadataICARUS::
FileCatalogMetadataICARUS(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
{
  // Get parameters.
  std::cout << "FileCatalogMetadataICARUS::parameterset begin" << std::endl;
  fFCLName = pset.get<std::string>("FCLName");
  fProjectName = pset.get<std::string>("ProjectName");
  fProjectStage = pset.get<std::string>("ProjectStage");
  fProjectVersion = pset.get<std::string>("ProjectVersion");    
  fProjectSoftware = pset.get<std::string>("ProjectSoftware","");    
  fProductionName = pset.get<std::string>("ProductionName","");  //Leave as default value if not running a production   
  fProductionType = pset.get<std::string>("ProductionType",""); //Leave as default value if not running a production

  fMerge = pset.get<int>("Merge", 0);
  fParameters = pset.get<std::vector<std::string> >("Parameters", std::vector<std::string>());
  // It doesn't make sense for parameter vector to have an odd number of elements.
  if(fParameters.size() % 2 != 0) {
    throw cet::exception("FileCatalogMetadataICARUS") 
      << "Parameter vector has odd number of elements.\n";
  }


  // Register for callbacks.

  reg.sPostBeginJob.watch(this, &FileCatalogMetadataICARUS::postBeginJob);
  std::cout << "FileCatalogMetadataICARUS::parameterset ends" << std::endl;

}

//--------------------------------------------------------------------
// PostBeginJob callback.
// Insert per-job metadata via FileCatalogMetadata service.
void icarusutil::FileCatalogMetadataICARUS::postBeginJob()
{
  // Get art metadata service.

  std::cout << "FileCatalogMetadataICARUS::postBeginJob() begin" << std::endl;
  art::ServiceHandle<art::FileCatalogMetadata> mds;

  // Add metadata.
  mds->addMetadata("fcl.name", fFCLName);
  mds->addMetadata("icarus_project.name", fProjectName);
  mds->addMetadata("icarus_project.stage", fProjectStage);
  mds->addMetadata("icarus_project.version", fProjectVersion);
  mds->addMetadata("icarus_project.software", fProjectSoftware);
  mds->addMetadata("production.name", fProductionName);
  mds->addMetadata("production.type", fProductionType);
  std::ostringstream ostr;
  ostr << fMerge;
  mds->addMetadata("merge.merge", ostr.str());
  mds->addMetadata("merge.merged", "0");
  for(unsigned int i=0; i<fParameters.size(); i += 2)
    mds->addMetadata(fParameters[i], fParameters[i+1]);
  std::cout << "FileCatalogMetadataICARUS::postBeginJob() end" << std::endl;
}

DEFINE_ART_SERVICE(icarusutil::FileCatalogMetadataICARUS)
