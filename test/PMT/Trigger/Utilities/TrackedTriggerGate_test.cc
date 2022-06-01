/**
 * @file TrackedTriggerGate_test.cc
 * @brief Unit test for some utilities in `TrackedTriggerGate.h`
 * @date January 13, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see icaruscode/PMT/Trigger/Utilities/TrackedTriggerGate.h
 * 
 * Largely incomplete.
 */

// ICARUS libraries
#include "icaruscode/PMT/Trigger/Utilities/TrackedOpticalTriggerGate.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h" // sbn::OpDetWaveformMeta


// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h"
#include "larcorealg/CoreUtils/DebugUtils.h"

// Boost libraries
#define BOOST_TEST_MODULE ( TrackedTriggerGate_test )
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <utility> // std::as_const(), std::move()
#include <type_traits> // std::is_same_v


// -----------------------------------------------------------------------------
// --- static tests
// -----------------------------------------------------------------------------

void static_tests() {
  
  using GateData_t = icarus::trigger::OpticalTriggerGateData_t::GateData_t;
  
  icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta> trackedGate;
  icarus::trigger::OpticalTriggerGateData_t gate;
  GateData_t gateData;
  
  //
  // gateIn()
  //
  static_assert(std::is_same_v<
    decltype(gateIn(trackedGate)),
    icarus::trigger::OpticalTriggerGateData_t&
    >);
  static_assert(std::is_same_v<
    decltype(gateIn(std::as_const(trackedGate))),
    icarus::trigger::OpticalTriggerGateData_t const&
    >);
  static_assert(std::is_same_v<
    decltype(gateIn(std::move(trackedGate))),
    icarus::trigger::OpticalTriggerGateData_t&&
    >);
  
  static_assert(std::is_same_v<
    decltype(gateIn(gate)),
    icarus::trigger::OpticalTriggerGateData_t&
    >);
  static_assert(std::is_same_v<
    decltype(gateIn(std::as_const(gate))),
    icarus::trigger::OpticalTriggerGateData_t const&
    >);
  static_assert(std::is_same_v<
    decltype(gateIn(std::move(gate))),
    icarus::trigger::OpticalTriggerGateData_t&&
    >);
  
  //
  // gateDataIn()
  //
  static_assert(std::is_same_v<
    decltype(gateDataIn(trackedGate)),
    GateData_t&
    >);
  static_assert(std::is_same_v<
    decltype(gateDataIn(std::as_const(trackedGate))),
    GateData_t const&
    >);
  static_assert(std::is_same_v<
    decltype(gateDataIn(std::move(trackedGate))),
    GateData_t&&
    >);
  static_assert(std::is_same_v<
    decltype(gateDataIn(gate)),
    GateData_t&
    >);
  static_assert(std::is_same_v<
    decltype(gateDataIn(std::as_const(gate))),
    GateData_t const&
    >);
  static_assert(std::is_same_v<
    decltype(gateDataIn(std::move(gate))),
    GateData_t&&
    >);
  
  static_assert(std::is_same_v<
    decltype(gateDataIn(gateData)),
    GateData_t&
    >);
  static_assert(std::is_same_v<
    decltype(gateDataIn(std::as_const(gateData))),
    GateData_t const&
    >);
  static_assert(std::is_same_v<
    decltype(gateDataIn(std::move(gateData))),
    GateData_t&&
    >);
  
  //
  // ok, now some non-static tests just for safety
  //
  BOOST_TEST(gateIn(trackedGate) == trackedGate.gate());
  BOOST_TEST(gateIn(gate) == gate);
  BOOST_TEST(gateDataIn(trackedGate) == trackedGate.gate().gateLevels());
  BOOST_TEST(gateDataIn(gate) == gate.gateLevels());
  BOOST_TEST(gateDataIn(gateData) == gateData);
  
} // static_tests()


// -----------------------------------------------------------------------------
// --- other tests
// -----------------------------------------------------------------------------
void TrackedTriggerGate_test() {
  
  /*
   * This test is a copy from TriggerGateData unit test.
   * It's basically a placeholder, since the real tests are the static ones.
   */
  
  icarus::trigger::TrackedOpticalTriggerGate<sbn::OpDetWaveformMeta> trackedGate;
  icarus::trigger::OpticalTriggerGateData_t& gate = trackedGate.gate();
  
  BOOST_CHECK((gate.alwaysClosed()));
  
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
  
  BOOST_CHECK((!gate.alwaysClosed()));
  BOOST_TEST((gate.openingCount(-9) ==  0));
  BOOST_TEST((gate.openingCount(-8) ==  0));
  BOOST_TEST((gate.openingCount(-7) ==  0));
  BOOST_TEST((gate.openingCount(-6) ==  0));
  BOOST_TEST((gate.openingCount(-5) ==  0));
  BOOST_TEST((gate.openingCount(-4) ==  2));
  BOOST_TEST((gate.openingCount(-3) ==  2));
  BOOST_TEST((gate.openingCount(-2) ==  2));
  BOOST_TEST((gate.openingCount(-1) ==  0));
  BOOST_TEST((gate.openingCount( 0) ==  0));
  BOOST_TEST((gate.openingCount( 1) ==  0));
  BOOST_TEST((gate.openingCount( 2) ==  1));
  BOOST_TEST((gate.openingCount( 3) ==  1));
  BOOST_TEST((gate.openingCount( 4) ==  2));
  BOOST_TEST((gate.openingCount( 5) ==  2));
  BOOST_TEST((gate.openingCount( 6) ==  2));
  BOOST_TEST((gate.openingCount( 7) ==  1));
  BOOST_TEST((gate.openingCount( 8) ==  1));
  BOOST_TEST((gate.openingCount( 9) ==  1));
  BOOST_TEST((gate.openingCount(10) ==  2));
  BOOST_TEST((gate.openingCount(11) ==  2));
  BOOST_TEST((gate.openingCount(12) ==  4));
  BOOST_TEST((gate.openingCount(13) ==  4));
  BOOST_TEST((gate.openingCount(14) ==  4));
  BOOST_TEST((gate.openingCount(15) ==  3));
  BOOST_TEST((gate.openingCount(16) ==  3));
  BOOST_TEST((gate.openingCount(17) ==  2));
  BOOST_TEST((gate.openingCount(18) ==  2));
  BOOST_TEST((gate.openingCount(19) ==  2));
  BOOST_TEST((gate.openingCount(20) ==  3));
  BOOST_TEST((gate.openingCount(21) ==  3));
  BOOST_TEST((gate.openingCount(22) ==  3));
  BOOST_TEST((gate.openingCount(23) ==  1));
  BOOST_TEST((gate.openingCount(24) ==  1));
  BOOST_TEST((gate.openingCount(25) ==  1));
  BOOST_TEST((gate.openingCount(26) ==  2));
  BOOST_TEST((gate.openingCount(27) ==  2));
  BOOST_TEST((gate.openingCount(28) ==  2));
  BOOST_TEST((gate.openingCount(29) ==  1));
  BOOST_TEST((gate.openingCount(30) ==  1));
  BOOST_TEST((gate.openingCount(31) ==  0));
  BOOST_TEST((gate.openingCount(32) ==  0));
  BOOST_TEST((gate.openingCount(33) ==  0));
  BOOST_TEST((gate.openingCount(34) ==  4));
  BOOST_TEST((gate.openingCount(35) ==  4));
  BOOST_TEST((gate.openingCount(36) ==  4));
  BOOST_TEST((gate.openingCount(37) ==  0));
  BOOST_TEST((gate.openingCount(38) ==  0));
  BOOST_TEST((gate.openingCount(39) ==  0));
  BOOST_TEST((gate.openingCount(40) ==  1));
  BOOST_TEST((gate.openingCount(41) ==  1));
  BOOST_TEST((gate.openingCount(42) ==  1));
  BOOST_TEST((gate.openingCount(43) ==  0));
  BOOST_TEST((gate.openingCount(44) ==  0));
  BOOST_TEST((gate.openingCount(45) ==  0));
  BOOST_TEST((gate.openingCount(46) ==  0));
  BOOST_TEST((gate.openingCount(47) ==  0));
  BOOST_TEST((gate.openingCount(48) ==  0));
  BOOST_TEST((gate.openingCount(49) ==  0));
 
  BOOST_TEST((gate.findMaxOpen(  ) ==  12));
  BOOST_TEST((gate.findMaxOpen(12) ==  12));
  BOOST_TEST((gate.findMaxOpen(13) ==  13));
  BOOST_TEST((gate.findMaxOpen(14) ==  14));
  BOOST_TEST((gate.findMaxOpen(15) ==  34));
  
  auto [ lower2031, upper2031 ] = gate.openingRange(20, 30);
  BOOST_TEST((lower2031 ==  1));
  BOOST_TEST((upper2031 ==  4));
  
} // TrackedTriggerGate_test()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(TrackedTriggerGate_testcase) {
  
  static_tests();
  
  TrackedTriggerGate_test();
  
} // BOOST_AUTO_TEST_CASE(TrackedTriggerGate_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
