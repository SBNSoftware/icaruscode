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

class NormalizeTPCPerPlaneSQL : public INormalizeCharge
{
public:
  NormalizeTPCPerPlaneSQL(fhicl::ParameterSet const &pset);

  void configure(const fhicl::ParameterSet& pset) override;
  double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) override;

private:
  // Configuration
  std::string fDBFileName;
  std::string fDBTag;
  bool fVerbose;

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

DEFINE_ART_CLASS_TOOL(NormalizeTPCPerPlaneSQL)

  } // end namespace calo
} // end namespace icarus


icarus::calo::NormalizeTPCPerPlaneSQL::NormalizeTPCPerPlaneSQL(fhicl::ParameterSet const &pset):
  fDBFileName(pset.get<std::string>("DBFileName")),
  fDBTag(pset.get<std::string>("DBTag")),
  fVerbose(pset.get<bool>("Verbose", false)),
  fDB(fDBFileName, "", "", fDBTag, true, false) {}

void icarus::calo::NormalizeTPCPerPlaneSQL::configure(const fhicl::ParameterSet& pset) {}

icarus::calo::NormalizeTPCPerPlaneSQL::ScaleInfo icarus::calo::NormalizeTPCPerPlaneSQL::GetScaleInfo(uint64_t run) {
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
  for (unsigned ch = 0; ch < 12; ch++) {
    double scale;
    fDB.GetNamedChannelData(ch, "scale", scale);

    thisscale.scale[ch] = scale;
  }
  // Set the cache
  fScaleInfos[run] = thisscale;

  return thisscale;
}

double icarus::calo::NormalizeTPCPerPlaneSQL::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {
  // Get the info
  ScaleInfo i = GetScaleInfo(e.id().runID().run());

  // Lookup the TPC, cryo
  unsigned tpc = hit.WireID().TPC;
  unsigned cryo = hit.WireID().Cryostat;
  unsigned plane = hit.WireID().Plane;

  // Get the TPC-Plane index
  unsigned itpc_plane = 2*cryo + tpc/2 + plane*4;

  double scale = 1;

  // TODO: what to do if no scale is found? throw an exception??
  if (i.scale.count(itpc_plane)) scale = i.scale.at(itpc_plane);

  if (fVerbose) std::cout << "NormalizeTPCPerPlaneSQL Tool -- Data at Cryo: " << cryo << " TPC: " << tpc << " Plane: " << plane << " itpc_plane: " << itpc_plane << " scale: " << scale << std::endl;

  return dQdx * scale;
}

