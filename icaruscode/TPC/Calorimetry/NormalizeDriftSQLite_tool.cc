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
#include <optional>
#include <cassert>

namespace icarus {
  namespace calo {

class NormalizeDriftSQLite : public INormalizeCharge
{
public:
  NormalizeDriftSQLite(fhicl::ParameterSet const &pset);

  void configure(const fhicl::ParameterSet& pset) override;
  void setup(const art::Event& e) override;
  double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) override;

private:
  // Configuration
  std::string fDBFileName;
  std::string fDBTag;
  bool fVerbose;

  lariov::DBFolder fDB;

  std::optional<detinfo::DetectorClocksData> fClockData; // need delayed construction

  // Class to hold data from DB
  class RunInfo {
  public:
    double tau_EE;
    double tau_EW;
    double tau_WE;
    double tau_WW;
  };

  // Helpers
  RunInfo GetRunInfo(uint64_t run);

  // Cache run requests
  std::map<uint32_t, RunInfo> fRunInfos;
};

DEFINE_ART_CLASS_TOOL(NormalizeDriftSQLite)

  } // end namespace calo
} // end namespace icarus


icarus::calo::NormalizeDriftSQLite::NormalizeDriftSQLite(fhicl::ParameterSet const &pset):
  fDBFileName(pset.get<std::string>("DBFileName")),
  fDBTag(pset.get<std::string>("DBTag")),
  fVerbose(pset.get<bool>("Verbose", false)),
  fDB(fDBFileName, "", "", fDBTag, true, false)
{}

void icarus::calo::NormalizeDriftSQLite::configure(const fhicl::ParameterSet& pset) {}

void icarus::calo::NormalizeDriftSQLite::setup(const art::Event& e) {
  fClockData.emplace(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e));
}

icarus::calo::NormalizeDriftSQLite::RunInfo icarus::calo::NormalizeDriftSQLite::GetRunInfo(uint64_t run) {
  // check the cache
  if (fRunInfos.count(run)) {
    return fRunInfos.at(run);
  }

  // Look up the run
  //
  // Translate the run into a fake "timestamp"
  fDB.UpdateData((run+1000000000)*1000000000);

  RunInfo thisrun;

  // Iterate over the rows
  // Should be 4: one for each TPC
  for (unsigned ch = 0; ch < 4; ch++) {
    double tau;

    fDB.GetNamedChannelData(ch, "elifetime", tau);

    // Map channel to TPC
    if (ch == 0) thisrun.tau_EE = tau;
    if (ch == 1) thisrun.tau_EW = tau;
    if (ch == 2) thisrun.tau_WE = tau;
    if (ch == 3) thisrun.tau_WW = tau;
  }

  if (fVerbose) std::cout << "NormalizeDriftSQLite Tool -- Lifetime Data:" << "\nTPC EE: " << thisrun.tau_EE << "\nTPC EW: " << thisrun.tau_EW << "\nTPC WE: " << thisrun.tau_WE << "\nTPC WW: " << thisrun.tau_WW << std::endl;

  // Set the cache
  fRunInfos[run] = thisrun;

  return thisrun;
}

double icarus::calo::NormalizeDriftSQLite::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {

  assert(fClockData);

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
  double thit = fClockData->TPCTick2TrigTime(hit.PeakTime()) - t0;

  if (fVerbose) std::cout << "NormalizeDriftSQLite Tool -- Norm factor: " << exp(thit / thiselifetime) << " at TPC: " << tpc << " Cryo: " << cryo << " Time: " << thit << " Track T0: " << t0 << std::endl;

  // Scale
  if (thiselifetime > 0) {
    dQdx = dQdx*exp(thit / thiselifetime);
  }
  // TODO: what to do if no lifetime is found? throw an exception??
  else {}

  return dQdx;
}

