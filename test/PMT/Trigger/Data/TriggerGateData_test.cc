/**
 * @file   test/PMT/Trigger/Data/TriggerGateData_test.cc
 * @brief  Unit test for `TriggerGateData.h` header
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 29, 2019
 * @see    `icaruscode/PMT/Trigger/Data/TriggerGateData.h`
 *
 * The main purpose of this test is to make sure `TriggerGateData` code is
 * compiled, since it is just a header.
 */

// Boost libraries
#define BOOST_TEST_MODULE ( TriggerGateData_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/TriggerGateData.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h"


// -----------------------------------------------------------------------------
// --- implementation detail tests
// -----------------------------------------------------------------------------

template<
  typename Tick,
  typename TickInterval = util::quantities::concepts::interval_of<Tick>
  >
void instantiate() {
  
  icarus::trigger::TriggerGateData<Tick, TickInterval> gate;
  
  BOOST_TEST_MESSAGE("Created: " << gate);
  
  BOOST_CHECK(gate.alwaysClosed());
  
} // instantiate()

// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(opticaltick_testcase) {
  
  instantiate<detinfo::timescales::optical_tick>();
  
} // BOOST_AUTO_TEST_CASE(opticaltick_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
