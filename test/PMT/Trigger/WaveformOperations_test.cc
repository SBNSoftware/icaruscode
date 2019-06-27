/**
 * @file   WaveformOperations_test.cc
 * @brief  Unit test for waveform operations.
 * @date   June 27, 2019
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    `icaruscode/PMT/Trigger/WaveformOperations.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE WaveformOperationsTest
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// ICARUS libraries
#include "icaruscode/PMT/Trigger/WaveformOperations.h"


//------------------------------------------------------------------------------
void NegativePolarityTest() {
  
  using Sample_t = signed short int;
  
  using NegOp = icarus::waveform_operations::NegativePolarityOperations<Sample_t>;
  
  // negative polarity waveform with baseline +8000
  constexpr NegOp op(+8000);
  
  //
  // "absolute" operations
  //
  
  static_assert(NegOp::distance(+8000, +7000) == +1000);
  static_assert(NegOp::distance(+8000, +9000) == -1000);
  
  static_assert(NegOp::shiftBy(+8000, +1000) == +7000);
  static_assert(NegOp::shiftBy(+8000, -1000) == +9000);
  
  static_assert(NegOp::subtractBaseline(+7000, +8000) == +1000);
  static_assert(NegOp::subtractBaseline(+9000, +8000) == -1000);
  
  //
  // comparisons
  //
  
  static_assert( (NegOp::lessThan     (+8000, +7000)));
  static_assert(!(NegOp::lessThan     (+8000, +8000)));
  static_assert(!(NegOp::lessThan     (+8000, +9000)));
  static_assert(!(NegOp::greaterThan  (+8000, +7000)));
  static_assert(!(NegOp::greaterThan  (+8000, +8000)));
  static_assert( (NegOp::greaterThan  (+8000, +9000)));
  static_assert(!(NegOp::noLessThan   (+8000, +7000)));
  static_assert( (NegOp::noLessThan   (+8000, +8000)));
  static_assert( (NegOp::noLessThan   (+8000, +9000)));
  static_assert( (NegOp::noGreaterThan(+8000, +7000)));
  static_assert( (NegOp::noGreaterThan(+8000, +8000)));
  static_assert(!(NegOp::noGreaterThan(+8000, +9000)));
  
  //
  // operations relative to the baseline
  //
  
  static_assert(op.shiftFromBaseline(+1000) == +7000);
  static_assert(op.shiftFromBaseline(-1000) == +9000);
  
  static_assert(op.subtractBaseline(+7000) == +1000);
  static_assert(op.subtractBaseline(+9000) == -1000);
  
  
} // void NegativePolarityTest()


//------------------------------------------------------------------------------
void PositivePolarityTest() {
  
  using Sample_t = signed short int;
  
  using PosOp = icarus::waveform_operations::PositivePolarityOperations<Sample_t>;
  
  // negative polarity waveform with baseline +8000
  constexpr PosOp op(+8000);
  
  //
  // "absolute" operations
  //
  
  static_assert(PosOp::distance(+8000, +9000) == +1000);
  static_assert(PosOp::distance(+8000, +7000) == -1000);
  
  static_assert(PosOp::shiftBy(+8000, +1000) == +9000);
  static_assert(PosOp::shiftBy(+8000, -1000) == +7000);
  
  static_assert(PosOp::subtractBaseline(+9000, +8000) == +1000);
  static_assert(PosOp::subtractBaseline(+7000, +8000) == -1000);
  
  //
  // comparisons
  //
  
  static_assert(!(PosOp::lessThan     (+8000, +7000)));
  static_assert(!(PosOp::lessThan     (+8000, +8000)));
  static_assert( (PosOp::lessThan     (+8000, +9000)));
  static_assert( (PosOp::greaterThan  (+8000, +7000)));
  static_assert(!(PosOp::greaterThan  (+8000, +8000)));
  static_assert(!(PosOp::greaterThan  (+8000, +9000)));
  static_assert( (PosOp::noLessThan   (+8000, +7000)));
  static_assert( (PosOp::noLessThan   (+8000, +8000)));
  static_assert(!(PosOp::noLessThan   (+8000, +9000)));
  static_assert(!(PosOp::noGreaterThan(+8000, +7000)));
  static_assert( (PosOp::noGreaterThan(+8000, +8000)));
  static_assert( (PosOp::noGreaterThan(+8000, +9000)));
  
  //
  // operations relative to the baseline
  //
  
  static_assert(op.shiftFromBaseline(+1000) == +9000);
  static_assert(op.shiftFromBaseline(-1000) == +7000);
  
  static_assert(op.subtractBaseline(+9000) == +1000);
  static_assert(op.subtractBaseline(+7000) == -1000);
  
  
} // void PositivePolarityTest()


//------------------------------------------------------------------------------
//---  The tests
//---
BOOST_AUTO_TEST_CASE( AllTests ) {
  
  NegativePolarityTest();
  PositivePolarityTest();
  
} // BOOST_AUTO_TEST_CASE( AllTests )

