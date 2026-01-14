// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larevt/CalibrationDBI/Providers/DBFolder.h"

// Tool include
#include "larreco/Calorimetry/INormalizeCharge.h"

// Services
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Lab helpers
#include "wda.h"

// C++
#include <string>

namespace icarus {
  namespace calo {

class NormalizeTPCSQL : public INormalizeCharge
{
public:
  NormalizeTPCSQL(fhicl::ParameterSet const &pset);

  void configure(const fhicl::ParameterSet& pset) override;
  double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) override;

private:
  // Configuration
  std::string fDBFileName;
  std::string fDBTag;
  bool fVerbose;
  int fMC;

  lariov::DBFolder fDB;

  // Class to hold data from DB
  class ScaleInfo {
  public:
    std::map<unsigned, double> scale;
  };

  // Helpers
  ScaleInfo GetScaleInfo(uint64_t run);

  // Cache run requests
  std::map<uint64_t, ScaleInfo> fScaleInfos;
};

DEFINE_ART_CLASS_TOOL(NormalizeTPCSQL)

  } // end namespace calo
} // end namespace icarus


icarus::calo::NormalizeTPCSQL::NormalizeTPCSQL(fhicl::ParameterSet const &pset):
  fDBFileName(pset.get<std::string>("DBFileName")),
  fDBTag(pset.get<std::string>("DBTag")),
  fVerbose(pset.get<bool>("Verbose", false)),
  fMC(pset.get<int>("MC")),
  fDB(fDBFileName, "", "", fDBTag, true, false) {}

void icarus::calo::NormalizeTPCSQL::configure(const fhicl::ParameterSet& pset) {}

icarus::calo::NormalizeTPCSQL::ScaleInfo icarus::calo::NormalizeTPCSQL::GetScaleInfo(uint64_t run) {

  std::cout << "NormalizeTPCSQL Tool -- Getting scale info for run: " << run << std::endl;

  // check the cache
  if (fScaleInfos.count(run)) {
    return fScaleInfos.at(run);
  }

  // Look up the run
  //
  // Translate the run into a fake "timestamp"
  fDB.UpdateData((run+1000000000)*1000000000);

  // Collect the run info
  ScaleInfo thisscale;

  // Iterate over the rows
  for (unsigned ch = 0; ch < 4; ch++) {
    double scale;
    fDB.GetNamedChannelData(ch, "scale", scale);

    thisscale.scale[ch] = scale;
  }
  // Set the cache
  fScaleInfos[run] = thisscale;

  return thisscale;
}

double icarus::calo::NormalizeTPCSQL::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {
  
  std::cout << "NormalizeTPCSQL Tool -- MC Flag: " << fMC << " Run: " << e.id().runID().run() << ", Subrun: " << e.id().subRunID().run() << std::endl;

  // Get the info
  uint64_t runID = -1;
  switch (fMC) {
    case 1:
      runID = 1;
      break;
    case 2:
      runID = 9400;
      break;
    case 3:
      runID = 3;
      break;
    case 4:
      runID = 4;
      break;
    case 5:
      runID = 5;
      break;
    default:
      runID = e.id().runID().run();
      break;
  }

  ScaleInfo const& i = GetScaleInfo(runID);

  // Lookup the TPC, cryo
  unsigned tpc = hit.WireID().TPC;
  unsigned cryo = hit.WireID().Cryostat;
  // Get the TPC index
  unsigned itpc = 2*cryo + tpc/2;

  double scale = 1;

  // TODO: what to do if no scale is found? throw an exception??
  if (i.scale.count(itpc)) scale = i.scale.at(itpc);

  if (fVerbose) std::cout << "NormalizeTPCSQL Tool -- Data at itpc: " << itpc << " scale: " << scale << std::endl;

  return dQdx * scale;
}

