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

// C++
#include <string>

namespace icarus {
  namespace calo {

class NormalizeTPCLocal : public INormalizeCharge
{
public:
  NormalizeTPCLocal(fhicl::ParameterSet const &pset);

  void configure(const fhicl::ParameterSet& pset) override;
  double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) override;

private:
  // Configuration
  std::vector<double> fCalConstants;
  bool fVerbose;
};

DEFINE_ART_CLASS_TOOL(NormalizeTPCLocal)

  } // end namespace calo
} // end namespace icarus


icarus::calo::NormalizeTPCLocal::NormalizeTPCLocal(fhicl::ParameterSet const &pset):
  fCalConstants(pset.get<std::vector<double>>("CalConstants")),
  fVerbose(pset.get<bool>("Verbose", false)) {}

void icarus::calo::NormalizeTPCLocal::configure(const fhicl::ParameterSet& pset) {}

double icarus::calo::NormalizeTPCLocal::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {
  // Lookup the TPC, cryo
  unsigned tpc = hit.WireID().TPC;
  unsigned cryo = hit.WireID().Cryostat;
  // Get the TPC index
  unsigned itpc = 2*cryo + tpc/2;

  double scale = fCalConstants.at(itpc);

  if (fVerbose) std::cout << "NormalizeTPCLocal Tool -- Data at itpc: " << itpc << " scale: " << scale << std::endl;

  return dQdx * scale;
}

