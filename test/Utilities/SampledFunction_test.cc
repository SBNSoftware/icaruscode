/**
 * @file   SampledFunction_test.cc
 * @brief  Unit test for `util::SampledFunction`.
 * @date   February 14, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    `icaruscode/Utilities/SampledFunction.h
 */

// Boost libraries
#define BOOST_TEST_MODULE SampledFunction
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()
#include <boost/test/tools/floating_point_comparison.hpp> // BOOST_CHECK_CLOSE()

// ICARUS libraries
#include "icaruscode/Utilities/SampledFunction.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"


//------------------------------------------------------------------------------
void IdentityTest() {
  
  auto identity = [](double x){ return x; };
  
  constexpr gsl::index nSamples = 16;
  constexpr gsl::index nSubsamples = 4;
  constexpr double min = 0.0;
  constexpr double max = static_cast<double>(nSamples * nSubsamples);
  constexpr double step = (max - min) / nSamples;
  constexpr double substep = step / nSubsamples;
  
  
  // function with 10 samples, sampled 4 times
  util::SampledFunction sampled { identity, min, max, nSamples, nSubsamples };
  
  // 
  // Query
  // 
  BOOST_CHECK_EQUAL(sampled.size(), nSamples);
  BOOST_CHECK_EQUAL(sampled.nSubsamples(), nSubsamples);
  BOOST_CHECK_EQUAL(sampled.lower(), min);
  BOOST_CHECK_EQUAL(sampled.upper(), max);
  BOOST_CHECK_CLOSE(sampled.rangeSize(), max - min, 1e-6);
  BOOST_CHECK_CLOSE(sampled.stepSize(), step, 1e-6);
  BOOST_CHECK_CLOSE(sampled.substepSize(), substep, 1e-6);
  
  for (auto const iSub: util::counter(nSubsamples)) BOOST_TEST_CONTEXT("Subsample: " << iSub)
  {
    
    double subsampleOffset = iSub * substep;
    
    auto const& subSample = sampled.subsample(iSub);
    BOOST_TEST_MESSAGE
      ("Subsample #" << iSub << ": " << subSample.size() << " samples");
    auto itSample = subSample.begin();
    
    for (auto const iSample: util::counter(nSamples)) BOOST_TEST_CONTEXT("Sample: " << iSample)
    {
      double const expected_x
        = static_cast<double>(iSub + iSample * nSubsamples);
      double const expected_value = identity(expected_x); // I wonder how much
      
      gsl::index const stepIndex
        = sampled.stepIndex(expected_x + substep / 2.0, iSub);
      
      BOOST_CHECK(sampled.isValidStepIndex(stepIndex));
      BOOST_CHECK_EQUAL(stepIndex, iSample);
      
      BOOST_CHECK_EQUAL(sampled.value(iSample, iSub), expected_value);
      
      BOOST_CHECK_EQUAL(*itSample, expected_value);
      BOOST_TEST_MESSAGE("[" << iSample << "] " << *itSample);
      
      ++itSample;
    } // for all samples in the subsample
    
    BOOST_CHECK(itSample == subSample.end());
    
    BOOST_CHECK(!sampled.isValidStepIndex
      (sampled.stepIndex(subsampleOffset + (max - min), iSub))
      );
    
  } // for all subsamples
  
  
} // void IdentityTest()


//------------------------------------------------------------------------------
//---  The tests
//---
BOOST_AUTO_TEST_CASE( TestCase ) {
  
  IdentityTest();
  
} // BOOST_AUTO_TEST_CASE( TestCase )

