#include "ICARUSDrifter.h"
#include "art/Framework/Principal/Event.h"
#include "larevt/CalibrationDBI/Providers/DBFolder.h"
#include "larcore/Geometry/Geometry.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include <memory>

WIRECELL_FACTORY(wclsICARUSDrifter, wcls::ICARUSDrifter,
                 wcls::IArtEventVisitor, WireCell::IDrifter)

using namespace WireCell;

wcls::ICARUSDrifter::ICARUSDrifter()
    : Drifter(),
      fDB(NULL) 
{
}

wcls::ICARUSDrifter::~ICARUSDrifter()
{
  if (fDB) delete fDB;
}

double wcls::ICARUSDrifter::GetLifetime(uint64_t run, uint64_t itpc) {
  // check the cache
  if (fLifetimes.count({run, itpc})) {
    if (fVerbose) std::cout << "Run: " << run << " ITPC: " << itpc << " SEEN. Returning: " << fLifetimes.at({run, itpc}) << ".\n";
    return fLifetimes.at({run, itpc});
  }

  // Look up the run
  //
  // Translate the run into a fake "timestamp"
  fDB->UpdateData((run+1000000000)*1000000000);

  double tau;
  fDB->GetNamedChannelData(itpc, "elifetime", tau);

  // Set the cache
  std::pair<uint64_t, uint64_t> toset {run, itpc};
  fLifetimes[toset] = tau;

  if (fVerbose) std::cout << "Run: " << run << " ITPC: " << itpc << " NEW. Returning: " << tau << ".\n";

  return tau;
}

void wcls::ICARUSDrifter::visit(art::Event & event)
{
    /// this function will be executed post WCT configuration !!!
    if (fVerbose) std::cout << "wcls ICARUSDrifter: electron lifetime read in!\n"; 
    
    if (fELifetimeCorrection) {
      double elifetime = GetLifetime(event.id().runID().run(), fTPC);
      Drifter::set_lifetime(elifetime*units::us);

      if (fVerbose) std::cout << "www: lifetime from database = " << elifetime << std::endl;
      if (fVerbose) std::cout << "www: lifetime from database = " << elifetime*units::us << std::endl;
    }
    else {
      Drifter::set_lifetime(1000.0*units::ms); // set to infinite
    }
}

void wcls::ICARUSDrifter::configure(const WireCell::Configuration& cfg)
{
    Drifter::configure(cfg);

    fDBFileName = cfg["DBFileName"].asString();
    fDBTag = cfg["DBTag"].asString();
    fVerbose = cfg["Verbose"].asBool();
    fDB = new lariov::DBFolder(fDBFileName, "", "", fDBTag, true, false);
    fTPC = cfg["TPC"].asInt();
    fELifetimeCorrection = cfg["ELifetimeCorrection"].asBool();
}
