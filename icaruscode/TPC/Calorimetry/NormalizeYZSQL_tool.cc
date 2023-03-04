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

class NormalizeYZSQL : public INormalizeCharge
{
public:
  NormalizeYZSQL(fhicl::ParameterSet const &pset);

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
    struct Point {
      int itpc;
      double y, z;
    };

    class ScaleBin {
      public:
      int itpc;
      double ylo;
      double yhi;
      double zlo;
      double zhi;

      double scale;

      bool contains(const Point& point) const noexcept;
      constexpr bool operator< (const ScaleBin& other) const noexcept;
    };

    struct BinComp {
      /// Returns if the point `b` is strictly after the bin `a`.
      bool operator() (const ScaleBin& a, const Point& b) const noexcept;
      /// Returns if the point `a` is strictly after the bin `b`.
      bool operator() (const Point& a, const ScaleBin& b) const noexcept;
    };

    float tzero; // Earliest time that this scale info is valid
    std::vector<ScaleBin> bins;

    /// Returns the bin containing the `point`, `nullptr` if none.
    ScaleBin const* findBin(const Point& point) const noexcept;

  };
  // Cache timestamp requests
  std::map<uint64_t, ScaleInfo> fScaleInfos;

  // Helpers
  const ScaleInfo GetScaleInfo(uint64_t timestamp);
};

DEFINE_ART_CLASS_TOOL(NormalizeYZSQL)

  } // end namespace calo
} // end namespace icarus


constexpr bool icarus::calo::NormalizeYZSQL::ScaleInfo::ScaleBin::operator<
  (const ScaleBin& other) const noexcept
{
  if (itpc != other.itpc) return itpc < other.itpc;
  if (yhi != other.yhi) return yhi < other.yhi;
  return zhi < other.zhi;
}

bool icarus::calo::NormalizeYZSQL::ScaleInfo::BinComp::operator()
  (const ScaleBin& a, const Point& b) const noexcept
{
  // the bin `a` must be strictly before the point `b`
  if (a.itpc != b.itpc) return a.itpc < b.itpc;
  if (a.yhi <= b.y) return true;
  if (a.ylo > b.y) return false;
  return a.zhi <= b.z;
}

bool icarus::calo::NormalizeYZSQL::ScaleInfo::BinComp::operator()
  (const Point& a, const ScaleBin& b) const noexcept
{
  // the point `a` must be strictly before the bin `b`
  if (a.itpc != b.itpc) return a.itpc < b.itpc;
  if (a.y < b.ylo) return true;
  if (a.y >= b.yhi) return false;
  return a.z < b.zlo;
}

bool icarus::calo::NormalizeYZSQL::ScaleInfo::ScaleBin::contains
  (const Point& point) const noexcept
{
  if (point.itpc != itpc) return false;
  if ((point.y < ylo) || (point.y >= yhi)) return false;
  return ((point.z >= zlo) && (point.z < zhi));
}


auto icarus::calo::NormalizeYZSQL::ScaleInfo::findBin
  (const Point& point) const noexcept -> ScaleBin const*
{
  auto itBin = std::upper_bound(bins.begin(), bins.end(), point, ScaleInfo::BinComp{});
  /*
   * because of how upper_bound() works:
   *  * if the point is before the first bin, `begin()` is returned
   *  * if the point is after the last bin, `end()` is returned
   *  * if the point is within the last bin, `end()` is returned
   *  * if the point is within the domain, the returned bin is the one after the desired one
   *  * if the point is in a TPC but beyond its last bin, first bin of next TPC is returned
   */
  return ((itBin != bins.cbegin()) && (--itBin)->contains(point))? &*itBin: nullptr;
}



icarus::calo::NormalizeYZSQL::NormalizeYZSQL(fhicl::ParameterSet const &pset):
  fDBFileName(pset.get<std::string>("DBFileName")),
  fDBTag(pset.get<std::string>("DBTag")),
  fVerbose(pset.get<bool>("Verbose", false)),
  fDB(fDBFileName, "", "", fDBTag, true, false) {}

void icarus::calo::NormalizeYZSQL::configure(const fhicl::ParameterSet& pset) {}

const icarus::calo::NormalizeYZSQL::ScaleInfo icarus::calo::NormalizeYZSQL::GetScaleInfo(uint64_t timestamp) {
  // check the cache
  if (fScaleInfos.count(timestamp)) {
    return fScaleInfos.at(timestamp);
  }

  // Prep data
  fDB.UpdateData(timestamp*1e9);

  // Collect the timestamp info
  ScaleInfo thisscale;

  // Lookup the channels
  std::vector<lariov::DBChannelID_t> channels;
  fDB.GetChannelList(channels);

  // Iterate over the channels
  thisscale.bins.reserve(channels.size());
  for (unsigned ch = 0; ch < channels.size(); ch++) {
    std::string tpcname;
    fDB.GetNamedChannelData(ch, "tpc", tpcname);
    int itpc = -1;
    if (tpcname == "EE") itpc = 0;
    else if (tpcname == "EW") itpc = 1;
    else if (tpcname == "WE") itpc = 2;
    else if (tpcname == "WW") itpc = 3;
    else {
      throw cet::exception("NormalizeYZSQL") << "NormalizeYZSQL Tool -- Bad TPC name (" << tpcname << ").";
    }

    // Bin limits
    double ylo, yhi, zlo, zhi;
    fDB.GetNamedChannelData(ch, "ylow", ylo);
    fDB.GetNamedChannelData(ch, "yhigh", yhi);
    fDB.GetNamedChannelData(ch, "zlow", zlo);
    fDB.GetNamedChannelData(ch, "zhigh", zhi);
    
    // Get the scale
    double scale;
    fDB.GetNamedChannelData(ch, "scale", scale);

    ScaleInfo::ScaleBin bin;
    bin.ylo = ylo;
    bin.yhi = yhi;
    bin.zlo = zlo;
    bin.zhi = zhi;
    bin.itpc = itpc;
    bin.scale = scale;

    thisscale.bins.push_back(bin);
  }
  std::sort(thisscale.bins.begin(), thisscale.bins.end());

  // Set the cache
  fScaleInfos[timestamp] = thisscale;

  return thisscale;
}

double icarus::calo::NormalizeYZSQL::Normalize(double dQdx, const art::Event &e, 
    const recob::Hit &hit, const geo::Point_t &location, const geo::Vector_t &direction, double t0) {
  // Get the info
  ScaleInfo i = GetScaleInfo(e.time().timeHigh());


  // compute itpc
  int cryo = hit.WireID().Cryostat;
  int tpc = hit.WireID().TPC;
  int itpc = cryo*2 + tpc/2;
  // position
  double y = location.y();
  double z = location.z();

  ScaleInfo::Point const point { itpc, y, z };
  ScaleInfo::ScaleBin const* b = i.findBin(point);

  double const scale = b? b->scale: 1;
  if (!b) {
    // TODO: what to do if no lifetime is found? throw an exception??
  }

  if (fVerbose) std::cout << "NormalizeYZSQL Tool -- Data Cryo: " << cryo << " TPC: " << tpc << " iTPC: " << itpc << " Y: " << y << " Z: " << z << " scale: " << scale << std::endl;

  return dQdx / scale;
}

