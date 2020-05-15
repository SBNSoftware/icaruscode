/**
 * @file   FastAndPoorGauss_test.cc
 * @brief  Unit test for `util::FastAndPoorGauss`.
 * @date   February 15, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    `icaruscode/Utilities/FastAndPoorGauss.h
 */

// Boost libraries
#define BOOST_TEST_MODULE FastAndPoorGauss
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()
#include <boost/test/tools/floating_point_comparison.hpp> // BOOST_CHECK_CLOSE()

// ICARUS libraries
#include "icaruscode/Utilities/FastAndPoorGauss.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"


// ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFitResult.h"


//------------------------------------------------------------------------------
template <std::size_t NSamples>
void Test(TFile* pFile = nullptr, bool statChecks = true) {
  
  constexpr unsigned int NPoints = 1'000'000U;
  
  std::string const NSamplesTag = std::to_string(NSamples);
  
  util::FastAndPoorGauss<NSamples> gauss;
  BOOST_TEST_MESSAGE("Testing sampling " << NSamples << " points.");
  
  TH1D HGaus(
    ("HGaus" + NSamplesTag).c_str(),
    ("Distribution from util::FastAndPoorGauss<" + NSamplesTag + ">").c_str(),
    1024, -5.0, +5.0
    );
  HGaus.SetDirectory(pFile);
  
  util::UniformSequence<> extract { NPoints };
  for (auto _ [[gnu::unused]]: util::counter(NPoints)) {
    double const u = extract();
    double const z = gauss(u);
    HGaus.Fill(z);
  } // for
  
  TF1* pGaus = new TF1(("FGaus" + NSamplesTag).c_str(), "gaus", -5.0, +5.0);
  TFitResultPtr pFitRes = HGaus.Fit(pGaus, "QS"); // Quiet, reSult
  BOOST_CHECK(pFitRes->IsValid());
  if (statChecks) BOOST_CHECK_LT(pFitRes->Chi2()/pFitRes->Ndf(), 5.0);
  pFitRes->Print();
  
  // mean should be smaller than twice its uncertainty
  if (statChecks) BOOST_CHECK_SMALL(pFitRes->Parameter(1), 2.*pFitRes->Error(1));
  if (statChecks) BOOST_CHECK_CLOSE(pFitRes->Parameter(2), 1.0, 1.0); // stddev (1% tolerance)
  
  if (pFile) {
    HGaus.Write();
    HGaus.SetDirectory(nullptr);
  }
  
} // void Test()


//------------------------------------------------------------------------------
//---  The tests
//---
BOOST_AUTO_TEST_CASE( TestCase ) {
  
  TFile F("FastAndPoorGauss_test.root", "RECREATE");
  
  Test<   1024U>(&F, false); // with low points the statistics is not very good;
  Test<   4096U>(&F, false); // so we skip the statistics checks
  Test<   8192U>(&F, false);
  Test<  16384U>(&F);
  Test<  32768U>(&F);
  Test<  65536U>(&F);
  Test< 131072U>(&F);
  Test< 262144U>(&F);
  Test< 524288U>(&F);
//   Test<1048576U>(&F);
  
  F.Write();
  
} // BOOST_AUTO_TEST_CASE( TestCase )

