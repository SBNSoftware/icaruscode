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

class NormalizeYZ : public INormalizeCharge
{
public:
  NormalizeYZ(fhicl::ParameterSet const &pset);

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
    class ScaleBin {
      public:
      int itpc;
      double ylo;
      double yhi;
      double zlo;
      double zhi;

      double scale;
    };

    float tzero; // Earliest time that this scale info is valid
    std::vector<ScaleBin> bins;
  };

  // Helpers
  const ScaleInfo& GetScaleInfo(uint64_t timestamp);
  std::string URL(uint64_t timestamp);

  // Cache timestamp requests
  std::map<uint64_t, ScaleInfo> fScaleInfos;
};

DEFINE_ART_CLASS_TOOL(NormalizeYZ)

  } // end namespace calo
} // end namespace icarus


icarus::calo::NormalizeYZ::NormalizeYZ(fhicl::ParameterSet const &pset) {
  this->configure(pset);
}

void icarus::calo::NormalizeYZ::configure(const fhicl::ParameterSet& pset) {
  fURL = pset.get<std::string>("URL");
  fTimeout = pset.get<unsigned>("Timeout");
  fVerbose = pset.get<bool>("Verbose", false);
}

std::string icarus::calo::NormalizeYZ::URL(uint64_t timestamp) {
  return fURL + std::to_string(timestamp);
}

const icarus::calo::NormalizeYZ::ScaleInfo& icarus::calo::NormalizeYZ::GetScaleInfo(uint64_t timestamp) {
  // check the cache
  if (fScaleInfos.count(timestamp)) {
    return fScaleInfos.at(timestamp);
  }

  // Otherwise, look it up
  int error = 0;
  std::string url = URL(timestamp);

  if (fVerbose) std::cout << "NormalizeYZ Tool -- New Scale info, requesting data from url:\n" << url << std::endl;

  Dataset d = getDataWithTimeout(url.c_str(), "", fTimeout, &error);
  if (error) {
    throw cet::exception("NormalizeYZ") << "Calibration Database access failed. URL: (" << url << ") Error Code: " << error;
  }

  if (fVerbose) std::cout << "NormalizeYZ Tool -- Received HTTP response:\n" << getHTTPmessage(d) << std::endl;

  if (getHTTPstatus(d) != 200) {
    throw cet::exception("NormalizeYZ") 
      << "Calibration Database access failed. URL: (" << url
      << "). HTTP error status: " << getHTTPstatus(d) << ". HTTP error message: " << getHTTPmessage(d);
  }

  // Collect the timestamp info
  ScaleInfo thisscale;

  // Get the First row to get tzero
  error = 0;
  Tuple tup = getTuple(d, 0);
  float tzero = getDoubleValue(tup, 0, &error);
  if (error) {
    throw cet::exception("NormalizeYZ")
      << "Calibration Database access failed. URL: (" << url
      << "). Failed on tuple access, row 0, col 0. Error Code: " << error;
  }
  
  if (fVerbose) std::cout << "NormalizeYZ Tool -- Obtained T0: " << tzero << std::endl;

  // Check if we've seen this t0 before
  bool found_scale_t0 = false;
  for (auto const &scale_pair: fScaleInfos) {
    const ScaleInfo &scale = scale_pair.second;
    if (scale.tzero == tzero) {
      thisscale = scale;
      found_scale_t0 = true;

      if (fVerbose) std::cout << "NormalizeYZ Tool -- Found prior matching T0 from timestamp: " << scale_pair.first << std::endl;

      break;
    }
  }

  if (found_scale_t0) {
    fScaleInfos[timestamp] = thisscale;
    return fScaleInfos.at(timestamp);
  }

  // We haven't seen this timestamp before and we haven't seen the valid t0 before.
  //
  // Process the HTTP response
  thisscale.tzero = tzero;

  // Number of rows
  int n_tuple = getNtuples(d);
  if (n_tuple < 0) {
    throw cet::exception("NormalizeYZ") << "NormalizeYZ Tool -- Calibration Database access failed. URL: (" << url << ") Bad Tuple Number: " << n_tuple;
  }

  // Iterate over the rows
  // The first 4 are metadata
  for (unsigned row = 4; row < (unsigned)n_tuple; row++) {
    Tuple tup = getTuple(d, row);

    int err = 0;
    // Get the TPC value
    char tpcbuf[10];
    int strl = getStringValue(tup, 1, tpcbuf, 10, &err);
    (void) strl;
    if (err) {
      throw cet::exception("NormalizeYZ") << "NormalizeYZ Tool -- Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 1. Error Code: " << err;
    }
    int itpc = -1;
    std::string tpcname(tpcbuf);
    if (tpcname == "EE") itpc = 0;
    else if (tpcname == "EW") itpc = 1;
    else if (tpcname == "WE") itpc = 2;
    else if (tpcname == "WW") itpc = 3;
    else {
      throw cet::exception("NormalizeYZ") << "NormalizeYZ Tool -- Bad TPC name (" << tpcname << ").";
    }

    // Get the bin limits
    double ylo = getDoubleValue(tup, 8, &err);
    if (err) {
      throw cet::exception("NormalizeYZ") << "NormalizeYZ Tool -- Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 8. Error Code: " << err;
    }
    double yhi = getDoubleValue(tup, 9, &err);
    if (err) {
      throw cet::exception("NormalizeYZ") << "NormalizeYZ Tool -- Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 9. Error Code: " << err;
    }
    double zlo = getDoubleValue(tup, 10, &err);
    if (err) {
      throw cet::exception("NormalizeYZ") << "NormalizeYZ Tool -- Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 10. Error Code: " << err;
    }
    double zhi = getDoubleValue(tup, 11, &err);
    if (err) {
      throw cet::exception("NormalizeYZ") << "NormalizeYZ Tool -- Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 11. Error Code: " << err;
    }

    // Get the scale
    double scale = getDoubleValue(tup, 4, &err);
    if (err) {
      throw cet::exception("NormalizeYZ") << "NormalizeYZ Tool -- Calibration Database access failed. URL: (" << url << ") Failed on tuple access, row: " << row << ", col 4. Error Code: " << err;
    }

    ScaleInfo::ScaleBin bin;
    bin.ylo = ylo;
    bin.yhi = yhi;
    bin.zlo = zlo;
    bin.zhi = zhi;
    bin.itpc = itpc;
    bin.scale = scale;

    thisscale.bins.push_back(bin);
  }

  // Set the cache
  fScaleInfos[timestamp] = thisscale;
  return fScaleInfos.at(timestamp);
}

double icarus::calo::NormalizeYZ::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {
  // Get the info
  ScaleInfo i = GetScaleInfo(e.time().timeHigh());

  double scale = 1;
  bool found_bin = false;;

  // compute itpc
  int cryo = hit.WireID().Cryostat;
  int tpc = hit.WireID().TPC;
  int itpc = cryo*2 + tpc/2;
  // position
  double y = location.y();
  double z = location.z();

  for (const ScaleInfo::ScaleBin &b: i.bins) {
    if (itpc == b.itpc &&
        (y >= b.ylo) && (y < b.yhi) &&
        (z >= b.zlo) && (z < b.zhi)) {
      found_bin = true;
      scale = b.scale;
      break;
    }
  }
  // TODO: what to do if no lifetime is found? throw an exception??
  (void) found_bin;

  if (fVerbose) std::cout << "NormalizeYZ Tool -- Data Cryo: " << cryo << " TPC: " << tpc << " iTPC: " << itpc << " Y: " << y << " Z: " << z << " scale: " << scale << std::endl;

  return dQdx / scale;
}

