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

class NormalizeWire : public INormalizeCharge
{
public:
  NormalizeWire(fhicl::ParameterSet const &pset);

  void configure(const fhicl::ParameterSet& pset) override;
  double Normalize(double dQdx, const art::Event &e, const recob::Hit &h, const geo::Point_t &location, const geo::Vector_t &direction, double t0) override;

private:
  // Configuration
  int fTimeout;
  std::string fURL;
  bool fVerbose;

  // Class to hold data from DB
  class ScaleInfo {
  public:
    std::map<unsigned, double> scale;
  };

  // Helpers
  ScaleInfo GetScaleInfo(uint64_t timestamp);
  std::string URL(uint64_t timestamp);

  // Cache timestamp requests
  std::map<uint64_t, ScaleInfo> fScaleInfos;
};

DEFINE_ART_CLASS_TOOL(NormalizeWire)

  } // end namespace calo
} // end namespace icarus


icarus::calo::NormalizeWire::NormalizeWire(fhicl::ParameterSet const &pset) {
  this->configure(pset);
}

void icarus::calo::NormalizeWire::configure(const fhicl::ParameterSet& pset) {
  fURL = pset.get<std::string>("URL");
  fTimeout = pset.get<unsigned>("Timeout");
  fVerbose = pset.get<bool>("Verbose", false);
}

std::string icarus::calo::NormalizeWire::URL(uint64_t timestamp) {
  return fURL + std::to_string(timestamp);
}

icarus::calo::NormalizeWire::ScaleInfo icarus::calo::NormalizeWire::GetScaleInfo(uint64_t timestamp) {
  // check the cache
  if (fScaleInfos.count(timestamp)) {
    return fScaleInfos.at(timestamp);
  }

  // Otherwise, look it up
  int error = 0;
  std::string url = URL(timestamp);

  if (fVerbose) std::cout << "NormalizeWire Tool -- New Scale info, requesting data from url:\n" << url << std::endl;

  Dataset d = getDataWithTimeout(url.c_str(), "", fTimeout, &error);
  if (error) {
    throw cet::exception("NormalizeWire") << "Calibration Database access failed. URL: (" << url << ") Error Code: " << error;
  }

  if (fVerbose) std::cout << "NormalizeWire Tool -- Received HTTP response:\n" << getHTTPmessage(d) << std::endl;

  if (getHTTPstatus(d) != 200) {
    throw cet::exception("NormalizeWire") 
      << "Calibration Database access failed. URL: (" << url
      << "). HTTP error status: " << getHTTPstatus(d) << ". HTTP error message: " << getHTTPmessage(d);
  }

  // Collect the timestamp info
  ScaleInfo thisscale;

  // Number of rows
  int n_tuple = getNtuples(d);
  if (n_tuple < 0) {
    throw cet::exception("NormalizeWire") << "Calibration Database access failed. URL: (" << url << ") Bad Tuple Number: " << n_tuple;
  }

  // Iterate over the rows
  // The first 4 are metadata
  for (unsigned row = 4; row < (unsigned)n_tuple; row++) {
    Tuple tup = getTuple(d, row);

    int err = 0;
    // Get the channel number
    int ch = getLongValue(tup, 0, &err);
    if (error) {
      throw cet::exception("NormalizeWire") << "Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 0. Error Code: " << error;
    }

    // and the scale
    double scale = getDoubleValue(tup, 1, &err);
    if (error) {
      throw cet::exception("NormalizeWire") << "Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 1. Error Code: " << error;
    }

    thisscale.scale[ch] = scale;
  }

  // Set the cache
  fScaleInfos[timestamp] = thisscale;

  return thisscale;
}

double icarus::calo::NormalizeWire::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {
  // Get the info
  ScaleInfo i = GetScaleInfo(e.time().timeHigh());

  // Lookup the channel
  unsigned channel = hit.Channel();

  double scale = 1;

  // TODO: what to do if no lifetime is found? throw an exception??
  if (i.scale.count(channel)) scale = i.scale.at(channel);

  if (fVerbose) std::cout << "NormalizeWire Tool -- Data at channel: " << channel << " scale: " << scale << std::endl;

  return dQdx / scale;
}

