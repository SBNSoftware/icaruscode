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

class NormalizeDrift : public INormalizeCharge
{
public:
  NormalizeDrift(fhicl::ParameterSet const &pset);

  void configure(const fhicl::ParameterSet& pset) override;
  double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) override;

private:
  // Configuration
  int fTimeout;
  std::string fURL;
  bool fVerbose;

  // Class to hold data from DB
  class RunInfo {
  public:
    double tau_EE;
    double tau_EW;
    double tau_WE;
    double tau_WW;
  };

  // Helpers
  RunInfo GetRunInfo(uint32_t run);
  std::string URL(uint32_t run);

  // Cache run requests
  std::map<uint32_t, RunInfo> fRunInfos;
};

DEFINE_ART_CLASS_TOOL(NormalizeDrift)

  } // end namespace calo
} // end namespace icarus


icarus::calo::NormalizeDrift::NormalizeDrift(fhicl::ParameterSet const &pset) {
  this->configure(pset);
}

void icarus::calo::NormalizeDrift::configure(const fhicl::ParameterSet& pset) {
  fURL = pset.get<std::string>("URL");
  fTimeout = pset.get<unsigned>("Timeout");
  fVerbose = pset.get<bool>("Verbose", false);
}

std::string icarus::calo::NormalizeDrift::URL(uint32_t run) {
  return fURL + std::to_string(run);
}

icarus::calo::NormalizeDrift::RunInfo icarus::calo::NormalizeDrift::GetRunInfo(uint32_t run) {
  // check the cache
  if (fRunInfos.count(run)) {
    return fRunInfos.at(run);
  }

  // Otherwise, look it up
  int error = 0;
  std::string url = URL(run);

  if (fVerbose) std::cout << "NormalizeDrift Tool -- New Run info, requesting data from url:\n" << url << std::endl;

  Dataset d = getDataWithTimeout(url.c_str(), "", fTimeout, &error);

  if (error) {
    throw cet::exception("NormalizeDrift") << "Calibration Database access failed. URL: (" << url << ") Error Code: " << error;
  }

  if (fVerbose) std::cout << "NormalizeDrift Tool -- Received HTTP response:\n" << getHTTPmessage(d) << std::endl;

  if (getHTTPstatus(d) != 200) {
    throw cet::exception("NormalizeDrift") 
      << "Calibration Database access failed. URL: (" << url
      << "). HTTP error status: " << getHTTPstatus(d) << ". HTTP error message: " << getHTTPmessage(d);
  }


  // Check all the TPC's are set
  std::vector<bool> tpc_set(4, false);
  RunInfo thisrun;

  // Iterate over the rows
  // Should be 4: one for each TPC
  // The first 4 rows are metadata
  for (unsigned row = 4; row < 8; row++) {
    Tuple tup = getTuple(d, row);

    int err = 0;
    // Get the channel number
    int ch = getLongValue(tup, 0, &err);
    if (error) {
      throw cet::exception("NormalizeDrift") << "Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 0. Error Code: " << error;
    }

    // .. and the purity
    double tau = getDoubleValue(tup, 1, &err);
    if (error) {
      throw cet::exception("NormalizeDrift") << "Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 1. Error Code: " << error;
    }

    // Check the channel number
    if (ch < 0 || ch > 3) {
      throw cet::exception("NormalizeDrift") << "Calibration Database access failed. URL: (" << url << ") Bad channel number: " << ch;
    }

    tpc_set.at(ch) = true;

    // Map channel to TPC
    if (ch == 0) thisrun.tau_EE = tau;
    if (ch == 1) thisrun.tau_EW = tau;
    if (ch == 2) thisrun.tau_WE = tau;
    if (ch == 3) thisrun.tau_WW = tau;
  }

  if (fVerbose) std::cout << "NormalizeDrift Tool -- Lifetime Data:" << "\nTPC EE: " << thisrun.tau_EE << "\nTPC EW: " << thisrun.tau_EW << "\nTPC WE: " << thisrun.tau_WE << "\nTPC WW: " << thisrun.tau_WW << std::endl;
 
  // Make sure all the channels are set
  for (unsigned tpc = 0; tpc < 4; tpc++) {
    if (!tpc_set[tpc]) {
      throw cet::exception("NormalizeDrift") << "Calibration Database access failed. URL: (" << url << ") TPC not set: " << tpc;
    }
  }

  // Set the cache
  fRunInfos[run] = thisrun;

  return thisrun;
}

double icarus::calo::NormalizeDrift::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {
  // Services
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  // Get the info
  RunInfo runelifetime = GetRunInfo(e.id().runID().run());

  // lookup the TPC
  double thiselifetime = -1;
  unsigned tpc = hit.WireID().TPC;
  unsigned cryo = hit.WireID().Cryostat;

  // EE
  if (cryo == 0 && (tpc == 0 || tpc == 1)) thiselifetime = runelifetime.tau_EE;
  // EW
  if (cryo == 0 && (tpc == 2 || tpc == 3)) thiselifetime = runelifetime.tau_EW;
  // WE
  if (cryo == 1 && (tpc == 0 || tpc == 1)) thiselifetime = runelifetime.tau_WE;
  // WW
  if (cryo == 1 && (tpc == 2 || tpc == 3)) thiselifetime = runelifetime.tau_WW;
  
  // Get the hit time
  double thit = clock_data.TPCTick2TrigTime(hit.PeakTime()) - t0;

  if (fVerbose) std::cout << "NormalizeDrift Tool -- Norm factor: " << exp(thit / thiselifetime) << " at TPC: " << tpc << " Cryo: " << cryo << " Time: " << thit << " Track T0: " << t0 << std::endl;

  // Scale
  if (thiselifetime > 0) {
    dQdx = dQdx*exp(thit / thiselifetime);
  }
  // TODO: what to do if no lifetime is found? throw an exception??
  else {}

  return dQdx;
}

