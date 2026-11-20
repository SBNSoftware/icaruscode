// std includes
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TGraph2D.h"
#include "TMath.h"

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

namespace wiremod {

class WireModScaleDebug : public art::EDAnalyzer {
public:
  explicit WireModScaleDebug(fhicl::ParameterSet const& pset);
  void beginJob() override;
  void analyze(art::Event const&) override {}

private:
  std::string fRatioFileName_XXW;
  std::vector<double> fThetaEdges;

  static void shapeGraphAngle(TGraph2D& graph)
  {
    double* xPtr = graph.GetX();
    double* yPtr = graph.GetY();
    double* zPtr = graph.GetZ();
    for (int idx = 0; idx < graph.GetN(); ++idx) {
      double x = *(xPtr + idx);
      double y = *(yPtr + idx);
      double z = *(zPtr + idx);
      graph.SetPoint(idx, x, TMath::DegToRad() * y, z);
    }
  }

  static void shapeGraphPosNegate(TGraph2D& graph, bool negate)
  {
    if (!negate) return;
    double* xPtr = graph.GetX();
    double* yPtr = graph.GetY();
    double* zPtr = graph.GetZ();
    for (int idx = 0; idx < graph.GetN(); ++idx) {
      double x = *(xPtr + idx);
      double y = *(yPtr + idx);
      double z = *(zPtr + idx);
      graph.SetPoint(idx, -x, y, z);
    }
  }
};

WireModScaleDebug::WireModScaleDebug(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}
  , fRatioFileName_XXW{pset.get<std::string>("RatioFileName_XXW", "Run4_Coarse_XW.root")}
  , fThetaEdges{pset.get<std::vector<double>>("ThetaEdges",
      {0, 6, 12, 18, 24, 30, 38, 46, 56, 90})}
{
}

void WireModScaleDebug::beginJob()
{
  const geo::GeometryCore* geometry = lar::providerFrom<geo::Geometry>();
  const geo::WireReadoutGeom* wireReadout = &(art::ServiceHandle<geo::WireReadout const>()->Get());

  TFile* ratioFile = new TFile(("./" + fRatioFileName_XXW).c_str(), "READ");
  if (!ratioFile || ratioFile->IsZombie()) {
    std::cout << "ERROR: could not open " << fRatioFileName_XXW
              << " in current directory. Run lar from the directory containing it.\n";
    return;
  }

  struct TPCSetInfo { unsigned cryo; unsigned tpcset; unsigned globalTPC; std::string name; };
  std::vector<TPCSetInfo> tpcsets = {
    {0, 0, 0, "EE"}, {0, 1, 1, "EW"}, {1, 0, 2, "WE"}, {1, 1, 3, "WW"}
  };

  std::cout
    << "TPC,plane,cathode_X,test_X,theta_xw_deg,theta_bin,"
    << "interp_raw_q,found_scale_q,interp_raw_sigma,found_scale_sigma,"
    << "nearest_node_theta_deg,expected_scale_q,expected_scale_sigma\n";

  for (auto const& tpcinfo : tpcsets) {
    readout::TPCsetID tpcsetid(tpcinfo.cryo, tpcinfo.tpcset);
    geo::TPCID tpcid = wireReadout->TPCsetToTPCs(tpcsetid).front();
    const geo::TPCGeo* tpcGeom = geometry->TPCPtr(tpcid);
    double cathodeX = tpcGeom->GetCathodeCenter().X();
    bool negate = (cathodeX <= 0);

    for (int plane = 0; plane < 3; ++plane) {
      std::string nameQ = "TPC" + std::to_string(tpcinfo.globalTPC) + "_plane" +
                           std::to_string(plane) + "_ratio_integral";
      std::string nameS = "TPC" + std::to_string(tpcinfo.globalTPC) + "_plane" +
                           std::to_string(plane) + "_ratio_width";

      TGraph2D* gQ_orig = dynamic_cast<TGraph2D*>(ratioFile->Get(nameQ.c_str()));
      TGraph2D* gS_orig = dynamic_cast<TGraph2D*>(ratioFile->Get(nameS.c_str()));
      if (!gQ_orig || !gS_orig) {
        std::cout << "MISSING graph " << nameQ << " or " << nameS << "\n";
        continue;
      }
      TGraph2D* gQ = static_cast<TGraph2D*>(gQ_orig->Clone());
      TGraph2D* gS = static_cast<TGraph2D*>(gS_orig->Clone());

      // mirrors WireModifierXXW_module.cc reconfigure() (InRadians:false, XAbs:true)
      shapeGraphAngle(*gQ);
      shapeGraphAngle(*gS);
      shapeGraphPosNegate(*gQ, negate);
      shapeGraphPosNegate(*gS, negate);

      double xmin = gQ->GetXmin(), xmax = gQ->GetXmax();
      double xmag = 0.5 * (std::abs(xmin) + std::abs(xmax));
      double testX = negate ? -xmag : xmag;

      double* xPtrQ = gQ->GetX();
      double* yPtrQ = gQ->GetY();
      double* zPtrQ = gQ->GetZ();
      double* zPtrS = gS->GetZ();
      int nPts = gQ->GetN();

      for (double theta_deg = 0.0; theta_deg <= 90.0; theta_deg += 2.0) {
        double theta_rad = theta_deg * TMath::DegToRad();

        double raw_q = gQ->Interpolate(testX, theta_rad);
        double raw_s = gS->Interpolate(testX, theta_rad);
        // production threshold check (WireModUtility.cc): > 0.001 else default 1.0
        double found_q = (raw_q > 0.001) ? raw_q : 1.0;
        double found_s = (raw_s > 0.001) ? raw_s : 1.0;

        // nearest real calibration grid node (2D nearest neighbor)
        double bestD = 1e18;
        int bestIdx = -1;
        for (int i = 0; i < nPts; ++i) {
          double dx = xPtrQ[i] - testX;
          double dy = yPtrQ[i] - theta_rad;
          double d = dx * dx * 1e-4 + dy * dy;  // theta-dominated distance (X fixed at midpoint)
          if (d < bestD) { bestD = d; bestIdx = i; }
        }
        double nearest_theta_deg = yPtrQ[bestIdx] * TMath::RadToDeg();
        double expected_q = zPtrQ[bestIdx];
        double expected_s = zPtrS[bestIdx];

        int thetaBin = -1;
        for (size_t b = 0; b + 1 < fThetaEdges.size(); ++b)
          if (theta_deg >= fThetaEdges[b] && theta_deg < fThetaEdges[b + 1]) { thetaBin = (int)b; break; }

        std::cout << std::fixed
          << tpcinfo.name << "," << plane << ","
          << std::setprecision(1) << cathodeX << ","
          << std::setprecision(1) << testX << ","
          << std::setprecision(1) << theta_deg << "," << thetaBin << ","
          << std::setprecision(4) << raw_q << "," << found_q << ","
          << raw_s << "," << found_s << ","
          << std::setprecision(1) << nearest_theta_deg << ","
          << std::setprecision(4) << expected_q << "," << expected_s
          << "\n";
      }
    }
  }

  ratioFile->Close();
}

DEFINE_ART_MODULE(WireModScaleDebug)

}  // namespace wiremod
