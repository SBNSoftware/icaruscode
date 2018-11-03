/**
 * @file   test/Utilities/frequency_test.cc
 * @brief  Unit test for `icaruscode/Utilities/quantities/frequency.h` header.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 2, 2018
 * @see    `icaruscode/Utilities/quantities.h`
 *
 */

// Boost libraries
#define BOOST_TEST_MODULE ( frequency_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "icaruscode/Utilities/quantities/frequency.h"

// C/C++ standard libraries
#include <string>
#include <type_traits> // std::decay_t<>


// -----------------------------------------------------------------------------
void test_frequency_literals() {
  
  using namespace util::quantities::frequency_literals;
  
  constexpr auto f_Hz = 2_Hz;
  static_assert(std::is_same<decltype(f_Hz), util::quantities::hertz const>());
  static_assert(f_Hz == 2.0);
  static_assert(f_Hz == 2.0_Hz);
  
  constexpr auto f_kHz = 20_kHz;
  static_assert
    (std::is_same<decltype(f_kHz), util::quantities::kilohertz const>());
  static_assert(f_kHz == 20.0);
  static_assert(f_kHz == 20.0_kHz);
  
  constexpr auto f_MHz = 200_MHz;
  static_assert
    (std::is_same<decltype(f_MHz), util::quantities::megahertz const>());
  static_assert(f_MHz == 200.0);
  static_assert(f_MHz == 200.0_MHz);
  
  constexpr auto f_GHz = 5_GHz;
  static_assert
    (std::is_same<decltype(f_GHz), util::quantities::gigahertz const>());
  static_assert(f_GHz == 5.0);
  static_assert(f_GHz == 5.0_GHz);
  
} // test_frequency_literals()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_frequency_special_operations() {
  
  using namespace util::quantities::time_literals;
  using namespace util::quantities::frequency_literals;
  
  // ---------------------------------------------------------------------------
  
  static_assert(std::is_same<decltype(5_us * 2_MHz), double>());
  static_assert(5_us * 2_MHz == 10.0);
  static_assert(5_ms * 2_MHz == 10'000.0);
  
  static_assert(std::is_same<decltype(2_MHz * 5_us), double>());
  static_assert(2_MHz * 5_us == 10.0);
  static_assert(2_MHz * 5_ms == 10'000.0);
  
  // ---------------------------------------------------------------------------
  
  static_assert
    (std::is_same<decltype(10 / 5_us), util::quantities::megahertz>());
  static_assert(10 / 5_us == 2_MHz);
  static_assert(10'000.0 / 5_ms == 2_MHz);
  
  // ---------------------------------------------------------------------------
  
  static_assert
    (std::is_same<decltype(10 / 2_MHz), util::quantities::microsecond>());
  static_assert(10 / 2_MHz == 5_us);
  static_assert(10'000.0 / 2_MHz == 5_ms);
  
  // ---------------------------------------------------------------------------
  
} // test_constexpr_operations()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(frequency_testcase) {
  
  test_frequency_literals();
  test_frequency_special_operations();
  
} // BOOST_AUTO_TEST_CASE(frequency_testcase)

// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
