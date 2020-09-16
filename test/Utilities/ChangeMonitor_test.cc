/**
 * @file   ChangeMonitor_test.cc
 * @brief  Unit test for utilities from `ChangeMonitor.h`.
 * @date   September 16, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    `icaruscode/Utilities/ChangeMonitor.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE ChangeMonitorTest
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// ICARUS libraries
#include "icaruscode/Utilities/ChangeMonitor.h"


//------------------------------------------------------------------------------
void documentationTest() {
  
  /*
   * // starts with no reference by default
   * icarus::ns::util::ChangeMonitor<int> monitor;
   * 
   * // first check just establishes the reference
   * int var = 0;
   * monitor(var); // this is also a `update()`, which returns no value
   * 
   * // reference is 0, new value is 1: a change is detected
   * if (monitor(1)) {
   *   std::cout << "Value has changed!" << std::endl;
   * }
   * 
   * var = 5; // this does not change the monitoring
   * // reference is now 1, new value is 1: no change is detected
   * if (monitor(1)) {
   *   std::cout << "Value has changed again!" << std::endl;
   * }
   * 
   * // more complex syntax (C++17) allows accessing the old reference value;
   * // reference is now 1, new value is 2: change is detected
   * if (auto newVal = monitor(2); newVal) {
   *   std::cout << "Value has changed from " << *newVal << " to 2!"
   *     << std::endl;
   * }
   */

  icarus::ns::util::ChangeMonitor<int> monitor;
  BOOST_CHECK(!monitor.hasReference());
  
  // first check just establishes the reference
  int var = 0;
  auto&& res1 = monitor(var); // this is also a `update()`, which returns no value
  BOOST_CHECK(!res1);
  BOOST_CHECK(monitor.hasReference());
  BOOST_CHECK_EQUAL(monitor.reference(), var);
  
  // reference is 0, new value is 1: a change is detected
  auto&& res2 = monitor(1);
  BOOST_CHECK(res2);
  BOOST_CHECK_EQUAL(res2.value(), 0);
  BOOST_CHECK(monitor.hasReference());
  BOOST_CHECK_EQUAL(monitor.reference(), 1);
  
  var = 5; // this does not change the monitoring
  // reference is now 1, new value is 1: no change is detected
  auto&& res3 = monitor(1);
  BOOST_CHECK(!res3);
  BOOST_CHECK(monitor.hasReference());
  BOOST_CHECK_EQUAL(monitor.reference(), 1);
  
  bool detected = false;
  if (auto prevVal = monitor(2); prevVal) {
    detected = true;
    BOOST_CHECK(prevVal);
    BOOST_CHECK_EQUAL(prevVal.value(), 1);
    BOOST_CHECK(monitor.hasReference());
    BOOST_CHECK_EQUAL(monitor.reference(), 2);
  }
  BOOST_CHECK(detected);
 
} // documentationTest()


//------------------------------------------------------------------------------
//---  The tests
//---
BOOST_AUTO_TEST_CASE( ChangeMonitorTestCase ) {
  
  documentationTest();
  
} // BOOST_AUTO_TEST_CASE( ChangeMonitorTestCase )


