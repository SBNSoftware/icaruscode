/**
 * @file   test/PMT/Data/WaveformBaseline_test.cc
 * @brief  Unit test for `WaveformBaseline.h` header.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 11, 2020
 * @see    `icaruscode/PMT/Data/WaveformBaseline.h`
 *
 * The main purpose of this test is to make sure `TriggerGateData` code is
 * compiled, since it is just a header.
 */

// ICARUS libraries
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"
#include "icaruscode/Utilities/WaveformOperations.h"

// LArSoft libraries
#include "lardataalg/Utilities/quantities/electronics.h" // counts

// Boost libraries
#define BOOST_TEST_MODULE ( WaveformBaseline_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK_EQUAL()

// C/C++ standard library
#include <sstream>
#include <vector>
#include <cmath> // std::round()


// -----------------------------------------------------------------------------
// --- WaveformBaseline tests
// -----------------------------------------------------------------------------

void WaveformBaseline_value_test() {
  
  constexpr float baselineValue { 6.6f };
  
  icarus::WaveformBaseline const baseline { baselineValue };
  
  BOOST_CHECK_EQUAL(baseline.fBaseline, baselineValue);
  BOOST_CHECK_EQUAL(baseline.baseline(), baselineValue);
  BOOST_CHECK_EQUAL(baseline(), baselineValue);
  
  std::ostringstream sstr;
  sstr << baseline;
  std::string const baselineStr { sstr.str() };
  sstr.str("");
  sstr << baselineValue;
  std::string const baselineValueStr { sstr.str() };
  
  BOOST_CHECK_EQUAL(baselineStr, baselineValueStr);
  
} // WaveformBaseline_value_test()


void WaveformBaseline_documentation1_test() {
  
  std::ostringstream sstr;
  
  /*
   * icarus::WaveformBaseline const baseline { 1.2f };
   * 
   * std::cout << "Baseline: " << baseline << " ADC" << std::endl;
   */
  
  icarus::WaveformBaseline const baseline { 1.2f };
  
  sstr << "Baseline: " << baseline << " ADC";
  
  
  BOOST_CHECK_EQUAL(sstr.str(), "Baseline: 1.2 ADC");
  
} // WaveformBaseline_documentation1_test()


void WaveformBaseline_documentation2_test() {
  /*
   *
   */
  
  //
  // initialization
  //
  using ADCCount_t = util::quantities::counts;
  std::vector<ADCCount_t> const data
    { ADCCount_t{ -6 },  ADCCount_t{ -8 }, ADCCount_t{ -2 } };
  
  //
  // test code
  //
  icarus::WaveformBaseline const baseline { -1.8 };
  
  icarus::waveform_operations::NegativePolarityOperations<float> const
    waveOps { baseline() };
  
  std::vector<ADCCount_t> subtracted;
  subtracted.reserve(data.size());
  for (ADCCount_t sample: data) {
    subtracted.emplace_back
      (std::round(waveOps.subtractBaseline(sample.value())));
  }
  
  //
  // checks
  //
  BOOST_CHECK_EQUAL(subtracted.size(), data.size());
  
  BOOST_CHECK_EQUAL(subtracted[0U], ADCCount_t{ 4 });
  BOOST_CHECK_EQUAL(subtracted[1U], ADCCount_t{ 6 });
  BOOST_CHECK_EQUAL(subtracted[2U], ADCCount_t{ 0 });
  
} // WaveformBaseline_documentation2_test()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(WaveformBaseline_documentation_testcase) {
  
  WaveformBaseline_value_test();
  
  WaveformBaseline_documentation1_test();
  WaveformBaseline_documentation2_test();
  
} // BOOST_AUTO_TEST_CASE(WaveformBaseline_documentation_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
