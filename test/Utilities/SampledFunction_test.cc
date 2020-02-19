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
  
  // [ -2.0 , 6.0 ] with step size 0.5 and 4 subsamples (0.125 substep size)
  constexpr gsl::index nSamples = 16;
  constexpr gsl::index nSubsamples = 4;
  constexpr double min = -2.0;
  constexpr double max = min + static_cast<double>(nSamples / 2);
  constexpr double range = max - min;
  constexpr double step = range / nSamples;
  constexpr double substep = step / nSubsamples;
  
  
  BOOST_TEST_MESSAGE("Test settings:"
    << "\nRange:      " << min << " -- " << max << " (range: " << range << ")"
    << "\nSamples:    " << nSamples << " (size: " << step << ")"
    << "\nSubsamples: " << nSubsamples << " (size: " << substep << ")"
    );
  
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
    
    double const subsampleStart = min + iSub * substep;
    
    auto const& subSample = sampled.subsample(iSub);
    BOOST_TEST_MESSAGE
      ("Subsample #" << iSub << ": " << subSample.size() << " samples");
    auto itSample = subSample.begin();
    
    for (auto const iSample: util::counter(-nSamples, 2*nSamples+1)) BOOST_TEST_CONTEXT("Sample: " << iSample)
    {
      bool const bInRange = (iSample >= 0) && (iSample < nSamples);
      
      double const expected_x
        = static_cast<double>(subsampleStart + iSample * step);
      double const expected_value = identity(expected_x); // I wonder how much
      
      if (bInRange) {
        BOOST_CHECK_EQUAL(sampled.value(iSample, iSub), expected_value);
        BOOST_CHECK_EQUAL(*itSample, expected_value);
        BOOST_TEST_MESSAGE("[" << iSample << "] " << *itSample);
        ++itSample;
      }
      
      // check lookup from within the substep
      for (double const shift: { 0.0, 0.25, 0.50, 0.75 }) BOOST_TEST_CONTEXT("Shift: " << shift) {
        double const expected_x_in_the_middle = expected_x + shift * substep;
        
        gsl::index const stepIndex
          = sampled.stepIndex(expected_x_in_the_middle, iSub);
        
        BOOST_CHECK_EQUAL(sampled.isValidStepIndex(stepIndex), bInRange);
        BOOST_CHECK_EQUAL(stepIndex, iSample);
        
        BOOST_CHECK_EQUAL
          (sampled.closestSubsampleIndex(expected_x_in_the_middle), iSub);
        
      } // for shift
      
    } // for all samples in the subsample
    
    BOOST_CHECK(itSample == subSample.end());
    
    BOOST_CHECK(!sampled.isValidStepIndex
      (sampled.stepIndex(subsampleStart + max - min, iSub))
      );
    
  } // for all subsamples
  
} // void IdentityTest()


//------------------------------------------------------------------------------
//---  The tests
//---
BOOST_AUTO_TEST_CASE( TestCase ) {
  
  IdentityTest();
  
} // BOOST_AUTO_TEST_CASE( TestCase )

