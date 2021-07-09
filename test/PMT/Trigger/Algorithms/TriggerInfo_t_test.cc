/**
 * @file TriggerInfo_t_test.cc
 * @brief Unit test for utilities in `TriggerInfo_t.h`
 * @date July 8, 2021
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h
 * 
 * Largely incomplete.
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/details/TriggerInfo_t.h"
#include "sbnobj/ICARUS/PMT/Trigger/Data/TriggerGateData.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h"

// Boost libraries
#define BOOST_TEST_MODULE ( TriggerGateData_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()


// -----------------------------------------------------------------------------
// --- implementation detail tests
// -----------------------------------------------------------------------------

icarus::trigger::TriggerGateData<int, int> TestedInputGate() {
  
  // prepare the input:
  icarus::trigger::TriggerGateData<int, int> gate;
  
  BOOST_CHECK(gate.alwaysClosed());
  
  /*
   *   ^                     12                    34
   * 4 |                       ===   20              ===
   * 3 |     -4       4    10     ==   === 26
   * 2 |       ===  2  ===   ==  15 ===      ===         40
   * 1 |             ==  7===      17   23=== 29==     37  ===
   * 0-*=======---=|=---,----,----,----,----,----,===-,-===,--===============
   *    -10  -5    0    5   10        20        30        40        50
   */
  
  gate.openBetween(-4, -1, 2); // -> 2
  gate.openAt ( 2);    // -> 1
  gate.openAt ( 4);    // -> 2
  gate.closeAt( 7);    // -> 1
  gate.openAt (10);    // -> 2
  gate.openAt (12, 2); // -> 4
  gate.closeAt(15);    // -> 3
  gate.closeAt(17);    // -> 2
  gate.openAt (20);    // -> 3
  gate.closeAt(23, 2); // -> 1
  gate.openAt (26);    // -> 2
  gate.closeAt(29);    // -> 1
  gate.closeAt(31);    // -> 0
  gate.openAt (34, 4); // -> 4
  gate.closeAt(37, 4); // -> 0
  gate.openAt (40);    // -> 1
  gate.closeAt(43);    // -> 0
  
  BOOST_CHECK(!gate.alwaysClosed());
  BOOST_CHECK_EQUAL(gate.openingCount(-9), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(-8), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(-7), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(-6), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(-5), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(-4), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(-3), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(-2), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(-1), 0);
  BOOST_CHECK_EQUAL(gate.openingCount( 0), 0);
  BOOST_CHECK_EQUAL(gate.openingCount( 1), 0);
  BOOST_CHECK_EQUAL(gate.openingCount( 2), 1);
  BOOST_CHECK_EQUAL(gate.openingCount( 3), 1);
  BOOST_CHECK_EQUAL(gate.openingCount( 4), 2);
  BOOST_CHECK_EQUAL(gate.openingCount( 5), 2);
  BOOST_CHECK_EQUAL(gate.openingCount( 6), 2);
  BOOST_CHECK_EQUAL(gate.openingCount( 7), 1);
  BOOST_CHECK_EQUAL(gate.openingCount( 8), 1);
  BOOST_CHECK_EQUAL(gate.openingCount( 9), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(10), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(11), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(12), 4);
  BOOST_CHECK_EQUAL(gate.openingCount(13), 4);
  BOOST_CHECK_EQUAL(gate.openingCount(14), 4);
  BOOST_CHECK_EQUAL(gate.openingCount(15), 3);
  BOOST_CHECK_EQUAL(gate.openingCount(16), 3);
  BOOST_CHECK_EQUAL(gate.openingCount(17), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(18), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(19), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(20), 3);
  BOOST_CHECK_EQUAL(gate.openingCount(21), 3);
  BOOST_CHECK_EQUAL(gate.openingCount(22), 3);
  BOOST_CHECK_EQUAL(gate.openingCount(23), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(24), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(25), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(26), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(27), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(28), 2);
  BOOST_CHECK_EQUAL(gate.openingCount(29), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(30), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(31), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(32), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(33), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(34), 4);
  BOOST_CHECK_EQUAL(gate.openingCount(35), 4);
  BOOST_CHECK_EQUAL(gate.openingCount(36), 4);
  BOOST_CHECK_EQUAL(gate.openingCount(37), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(38), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(39), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(40), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(41), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(42), 1);
  BOOST_CHECK_EQUAL(gate.openingCount(43), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(44), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(45), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(46), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(47), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(48), 0);
  BOOST_CHECK_EQUAL(gate.openingCount(49), 0);
 
  BOOST_CHECK_EQUAL(gate.findMaxOpen(  ), 12);
  BOOST_CHECK_EQUAL(gate.findMaxOpen(12), 12);
  BOOST_CHECK_EQUAL(gate.findMaxOpen(13), 13);
  BOOST_CHECK_EQUAL(gate.findMaxOpen(14), 14);
  BOOST_CHECK_EQUAL(gate.findMaxOpen(15), 34);
  
  auto [ lower2031, upper2031 ] = gate.openingRange(20, 30);
  BOOST_CHECK_EQUAL(lower2031, 1);
  BOOST_CHECK_EQUAL(upper2031, 4);
  
  return gate;
} // TestedInputGate()


void GateOpeningInfoExtractor_test() {
  
  // prepare the input:
  icarus::trigger::TriggerGateData<int, int> const gate { TestedInputGate() };
  
  /*
   *   ^                     12                    34
   * 4 |                       ===   20              ===
   * 3 |     -4       4    10     ==   === 26
   * 2 |       ===  2  ===   ==  15 ===      ===         40
   * 1 |             ==  7===      17   23=== 29==     37  ===
   * 0-*=======---=|=---,----,----,----,----,----,===-,-===,--===============
   *    -10  -5    0    5   10        20        30        40        50
   */
  
  icarus::trigger::details::GateOpeningInfoExtractor extract { gate, 2U };
  
  BOOST_CHECK_EQUAL(extract.openThreshold(),  2U);
  BOOST_CHECK_EQUAL(extract.closeThreshold(), 1U);
  BOOST_CHECK_EQUAL(extract.minGap(),         0U);
  BOOST_CHECK_EQUAL(extract.minWidth(),       1U);
  
  BOOST_CHECK(!extract.atEnd());
  
  auto opening = extract.findNextOpening();
  BOOST_CHECK(opening);
  BOOST_CHECK_EQUAL(opening.value().tick.value(), -4);
  BOOST_CHECK_EQUAL(opening.value().level,         2);
  BOOST_CHECK(!extract.atEnd());
  
  opening = extract.findNextOpening();
  BOOST_CHECK(opening);
  BOOST_CHECK_EQUAL(opening.value().tick.value(),  4);
  BOOST_CHECK_EQUAL(opening.value().level,         2);
  BOOST_CHECK(!extract.atEnd());
  
  opening = extract.findNextOpening();
  BOOST_CHECK(opening);
  BOOST_CHECK_EQUAL(opening.value().tick.value(), 10);
  BOOST_CHECK_EQUAL(opening.value().level,         4);
  BOOST_CHECK(!extract.atEnd());
  
  opening = extract.findNextOpening();
  BOOST_CHECK(opening);
  BOOST_CHECK_EQUAL(opening.value().tick.value(), 26);
  BOOST_CHECK_EQUAL(opening.value().level,         2);
  BOOST_CHECK(!extract.atEnd());
  
  opening = extract.findNextOpening();
  BOOST_CHECK(opening);
  BOOST_CHECK_EQUAL(opening.value().tick.value(), 34);
  BOOST_CHECK_EQUAL(opening.value().level,         4);
  BOOST_CHECK(extract.atEnd());
  
  opening = extract.findNextOpening();
  BOOST_CHECK(!opening);
  BOOST_CHECK(extract.atEnd());
  
} // GateOpeningInfoExtractor_test()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(GateOpeningInfoExtractor_testcase) {
  
  GateOpeningInfoExtractor_test();
  
} // BOOST_AUTO_TEST_CASE(GateOpeningInfoExtractor_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
