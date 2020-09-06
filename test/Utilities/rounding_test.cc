/**
 * @file   rounding_test.cc
 * @brief  Unit test for utilities from `rounding.h`.
 * @date   April 17, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    `icaruscode/Utilities/rounding.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE RoundingTest
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// ICARUS libraries
#include "icaruscode/Utilities/rounding.h"


//------------------------------------------------------------------------------
void roundupTest() {
  
  BOOST_CHECK_EQUAL(util::roundup(23.0, 2.5),      25.0);
  BOOST_CHECK_EQUAL(util::roundup(22.0, 2.5, 1.0), 23.5); // 21.0 + 1.0 => 22.5 + 1.0
  BOOST_CHECK_EQUAL(util::roundup(23.0, 2.5, 1.0), 23.5); // 22.0 + 1.0 => 22.5 + 1.0
  BOOST_CHECK_EQUAL(util::roundup(24.0, 2.5, 1.0), 26.0); // 23.0 + 1.0 => 25.0 + 1.0
  
  BOOST_CHECK_EQUAL(util::roundup(23, 5),    25);
  BOOST_CHECK_EQUAL(util::roundup(22, 5, 3), 23); // 19 + 3 => 20 + 3
  BOOST_CHECK_EQUAL(util::roundup(23, 5, 3), 23); // 20 + 3 => 20 + 3
  BOOST_CHECK_EQUAL(util::roundup(24, 5, 3), 28); // 21 + 3 => 25 + 3
  
} // roundupTest()


//------------------------------------------------------------------------------
void rounddownTest() {
  
  BOOST_CHECK_EQUAL(util::rounddown(23.0, 2.5),      22.5);
  BOOST_CHECK_EQUAL(util::rounddown(23.0, 2.5, 1.0), 21.0); // 21.0 + 1.0 => 20.0 + 1.0
  BOOST_CHECK_EQUAL(util::rounddown(23.0, 2.5, 1.0), 21.0); // 22.0 + 1.0 => 20.0 + 1.0
  BOOST_CHECK_EQUAL(util::rounddown(24.0, 2.5, 1.0), 23.5); // 23.0 + 1.0 => 22.5 + 1.0
  
  BOOST_CHECK_EQUAL(util::rounddown(23, 5),    20);
  BOOST_CHECK_EQUAL(util::rounddown(22, 5, 3), 18); // 19 + 3 => 15 + 3
  BOOST_CHECK_EQUAL(util::rounddown(23, 5, 3), 23); // 20 + 3 => 20 + 3
  BOOST_CHECK_EQUAL(util::rounddown(24, 5, 3), 23); // 21 + 3 => 20 + 3
  
} // rounddownTest()


//------------------------------------------------------------------------------
//---  The tests
//---
BOOST_AUTO_TEST_CASE( RoundUpTestCase ) {
  
  roundupTest();
  
} // BOOST_AUTO_TEST_CASE( RoundUpTestCase )


BOOST_AUTO_TEST_CASE( RoundDownTestCase ) {
  
  rounddownTest();
  
} // BOOST_AUTO_TEST_CASE( RoundDownTestCase )

