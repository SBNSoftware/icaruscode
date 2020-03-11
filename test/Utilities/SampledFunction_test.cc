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
void ExtendedRangeTest() {

  auto identity = [](double x){ return x; };

  // the following checks work for a monotonic function;
  // the value at x = stopAt should not be included in the function range.

  constexpr gsl::index nSubsamples = 4;
  constexpr double min = -2.0;
  constexpr double atLeast = +1.0;
  constexpr double stopBefore = +8.2; // deliberately avoid border effects
  constexpr double step = 0.5;
  constexpr double substep = step / nSubsamples;

  // this stop function does *not* include the stop value in the range
  // when that value would be the right limit of the range; for that to
  // be included, (y > stopValue) should be used.
  constexpr double stopValue = identity(stopBefore);
  auto const stopIf
    = [stopValue](double, double y){ return (y < 0) || (y >= stopValue); };

  // function with 10 samples, sampled 4 times
  util::SampledFunction sampled
    { identity, min, step, stopIf, nSubsamples, atLeast };

  constexpr gsl::index expected_nSamples
    = static_cast<gsl::index>(std::floor((stopBefore - min) / step));
  constexpr gsl::index expected_max = min + step * expected_nSamples;
  constexpr gsl::index expected_range = expected_max - min;

  //
  // Query
  //
  BOOST_CHECK_EQUAL(sampled.nSubsamples(), nSubsamples);
  BOOST_CHECK_EQUAL(sampled.lower(), min);
  BOOST_CHECK_CLOSE(sampled.stepSize(), step, 1e-6);
  BOOST_CHECK_CLOSE(sampled.substepSize(), substep, 1e-6);

  BOOST_CHECK_CLOSE(sampled.upper(), expected_max, 1e-6);
  BOOST_REQUIRE_EQUAL(sampled.size(), expected_nSamples);
  BOOST_CHECK_CLOSE(sampled.rangeSize(), expected_range, 1e-6);


  auto const nSamples = sampled.size();

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

    } // for all samples in the subsample

    BOOST_CHECK(itSample == subSample.end());

    BOOST_CHECK(!sampled.isValidStepIndex
      (sampled.stepIndex(subsampleStart + expected_max - min, iSub))
      );

  } // for all subsamples

} // void ExtendedRangeTest()


//------------------------------------------------------------------------------
//---  The tests
//---
BOOST_AUTO_TEST_CASE( TestCase ) {

  IdentityTest();
  ExtendedRangeTest();

} // BOOST_AUTO_TEST_CASE( TestCase )

