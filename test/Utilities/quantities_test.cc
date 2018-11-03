/**
 * @file   test/Utilities/quantities_test.cc
 * @brief  Unit test for `quantities.h` header
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 30, 2018
 * @see    icaruscode/Utilities/quantities.h
 *
 */

// Boost libraries
#define BOOST_TEST_MODULE ( quantities_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "icaruscode/Utilities/quantities/spacetime.h"
#include "icaruscode/Utilities/quantities/frequency.h"

// C/C++ standard libraries
#include <string>
#include <type_traits> // std::decay_t<>

// -----------------------------------------------------------------------------
void test_quantities_sign() {
  
  util::quantities::microseconds t { -4.0 };
  
  BOOST_CHECK_EQUAL(t, -4.0);
  static_assert(
    std::is_same<decltype(+t), util::quantities::microseconds>(),
    "Positive sign converts to a different type!"
    );
  BOOST_CHECK_EQUAL(+t, -4.0);
  static_assert(
    std::is_same<decltype(-t), util::quantities::microseconds>(),
    "Negative sign converts to a different type!"
    );
  BOOST_CHECK_EQUAL(-t, 4.0);
  static_assert(
    std::is_same<decltype(t.abs()), util::quantities::microseconds>(),
    "Negative sign converts to a different type!"
    );
  BOOST_CHECK_EQUAL(t.abs(), 4.0);
  
} // test_quantities_sign()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_quantities_conversions() {
  // 
  // conversions to other scales
  //
  util::quantities::seconds t_s { 7.0 };
  
  util::quantities::microseconds t_us(t_s);
  BOOST_CHECK_EQUAL(t_us, 7000000.0);
  
  t_us = t_s;
  BOOST_CHECK_EQUAL(t_us, 7000000.0);
  
  util::quantities::seconds t(t_us);
  BOOST_CHECK_EQUAL(t, 7.0);
  
} // test_quantities_conversions()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_quantities_comparisons() {
  // 
  // comparisons between quantities
  // 
  util::quantities::microseconds t_us { 7.0 };
  BOOST_CHECK_EQUAL(t_us, t_us);
  BOOST_CHECK( (t_us == t_us));
  BOOST_CHECK(!(t_us != t_us));
  BOOST_CHECK( (t_us >= t_us));
  BOOST_CHECK( (t_us <= t_us));
  BOOST_CHECK(!(t_us >  t_us));
  BOOST_CHECK(!(t_us <  t_us));
  
  util::quantities::nanoseconds  t_ns { 7.0 };
  BOOST_CHECK_NE(t_us, t_ns);
  BOOST_CHECK(!(t_us == t_ns));
  BOOST_CHECK( (t_us != t_ns));
  BOOST_CHECK( (t_us >= t_ns));
  BOOST_CHECK(!(t_us <= t_ns));
  BOOST_CHECK( (t_us >  t_ns));
  BOOST_CHECK(!(t_us <  t_ns));
  
  util::quantities::nanoseconds t2_ns { 7000.0 };
  BOOST_CHECK_EQUAL(t_us, t2_ns);
  BOOST_CHECK( (t_us == t2_ns));
  BOOST_CHECK(!(t_us != t2_ns));
  BOOST_CHECK( (t_us >= t2_ns));
  BOOST_CHECK( (t_us <= t2_ns));
  BOOST_CHECK(!(t_us >  t2_ns));
  BOOST_CHECK(!(t_us <  t2_ns));
  
  BOOST_CHECK_NE(t_ns, t2_ns);
  BOOST_CHECK(!(t_ns == t2_ns));
  BOOST_CHECK( (t_ns != t2_ns));
  BOOST_CHECK(!(t_ns >= t2_ns));
  BOOST_CHECK( (t_ns <= t2_ns));
  BOOST_CHECK(!(t_ns >  t2_ns));
  BOOST_CHECK( (t_ns <  t2_ns));
  
  // 
  // comparisons between quantity and number
  // 
  BOOST_CHECK_NE(t_us, 8.0);
  BOOST_CHECK(!(t_us == 8.0));
  BOOST_CHECK( (t_us != 8.0));
  BOOST_CHECK(!(t_us >= 8.0));
  BOOST_CHECK( (t_us <= 8.0));
  BOOST_CHECK(!(t_us >  8.0));
  BOOST_CHECK( (t_us <  8.0));
  
  BOOST_CHECK_EQUAL(t_us, 7.0);
  BOOST_CHECK( (t_us == 7.0));
  BOOST_CHECK(!(t_us != 7.0));
  BOOST_CHECK( (t_us >= 7.0));
  BOOST_CHECK( (t_us <= 7.0));
  BOOST_CHECK(!(t_us >  7.0));
  BOOST_CHECK(!(t_us <  7.0));
  
  BOOST_CHECK_NE(t_us, 6.0);
  BOOST_CHECK(!(t_us == 6.0));
  BOOST_CHECK( (t_us != 6.0));
  BOOST_CHECK( (t_us >= 6.0));
  BOOST_CHECK(!(t_us <= 6.0));
  BOOST_CHECK( (t_us >  6.0));
  BOOST_CHECK(!(t_us <  6.0));
  
  BOOST_CHECK(!(8.0 == t_us));
  BOOST_CHECK( (8.0 != t_us));
  BOOST_CHECK( (8.0 >= t_us));
  BOOST_CHECK(!(8.0 <= t_us));
  BOOST_CHECK( (8.0 >  t_us));
  BOOST_CHECK(!(8.0 <  t_us));
  
  BOOST_CHECK( (7.0 == t_us));
  BOOST_CHECK(!(7.0 != t_us));
  BOOST_CHECK( (7.0 >= t_us));
  BOOST_CHECK( (7.0 <= t_us));
  BOOST_CHECK(!(7.0 >  t_us));
  BOOST_CHECK(!(7.0 <  t_us));
  
  BOOST_CHECK(!(6.0 == t_us));
  BOOST_CHECK( (6.0 != t_us));
  BOOST_CHECK(!(6.0 >= t_us));
  BOOST_CHECK( (6.0 <= t_us));
  BOOST_CHECK(!(6.0 >  t_us));
  BOOST_CHECK( (6.0 <  t_us));
  
} // test_quantities_conversions()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_quantities_multiply_scalar() {
  //
  // multiplication and division by scalar
  //
  
  using namespace util::quantities::time_literals;
  
//   5_s * 6_s; // ERROR
//   5_s * 6_us; // ERROR
  
  util::quantities::seconds const t { 3.0 };
  auto const twice_t = 2.0 * t;
  static_assert(
    std::is_same
      <std::decay_t<decltype(twice_t)>, util::quantities::seconds>(),
    "Multiplication by a scalar converts to a different type!"
    );
  BOOST_CHECK_EQUAL(twice_t, 6.0);
  
  auto const t_twice = t * 2.0;
  static_assert(
    std::is_same
      <std::decay_t<decltype(t_twice)>, util::quantities::seconds>(),
    "Multiplication by a scalar converts to a different type!"
    );
  BOOST_CHECK_EQUAL(twice_t, 6.0);
  
  static_assert(
    std::is_same<decltype(twice_t / 2.0), util::quantities::seconds>(),
    "Division by a scalar converts to a different type!"
    );
  BOOST_CHECK_EQUAL(twice_t / 2.0, 3.0);
  
  static_assert(
    std::is_same<decltype(twice_t / t), double>(),
    "Division by a scalar is not the base type!"
    );
  BOOST_CHECK_EQUAL(twice_t / t, 2.0);
  
  static_assert(
    std::is_same<decltype(t / 300_us), double>(),
    "Division by a scalar is not the base type!"
    );
  BOOST_CHECK_EQUAL(t / 300_us, 10'000.0);
  
} // test_quantities_multiply_scalar()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_quantities_addition() {
  
  using namespace util::quantities::time_literals;
  
  //
  // sum and difference
  //
  
//  5_s + 700_ms; // ERROR!
//  5_s + 0.7; // ERROR!
  
  static_assert(
    std::is_same<std::decay_t<decltype(45_s + 5_s)>, util::quantities::seconds>(),
    "Addition converts to a different type!"
    );
  BOOST_CHECK_EQUAL(45_s + 5_s, 50_s);
  
  static_assert(
    std::is_same<decltype(5_s - 55_s), util::quantities::seconds>(),
    "Subtraction converts to a different type!"
    );
  BOOST_CHECK_EQUAL(5_s - 55_s, -50_s);
  
  
} // test_quantities_addition()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_quantities_increment() {
  
  using namespace util::quantities::time_literals;
  
  //
  // increment and decrement by a quantity
  //
  util::quantities::seconds t { 0.05 };
  
  t += 0.05_s;
  static_assert(
    std::is_same<decltype(t += 0.05_s), util::quantities::seconds&>(),
    "Increment converts to a different type!"
    );
  BOOST_CHECK_EQUAL(t, 0.1_s);
  
  t -= 0.05_s;
  static_assert(
    std::is_same<decltype(t -= 0.05_s), util::quantities::seconds&>(),
    "Decrement converts to a different type!"
    );
  BOOST_CHECK_EQUAL(t, 0.05_s);
  
  t += 50_ms;
  static_assert(
    std::is_same<decltype(t += 50_ms), util::quantities::seconds&>(),
    "Increment converts to a different type!"
    );
  BOOST_CHECK_EQUAL(t, 0.1_s);
  
  t -= 50_ms;
  static_assert(
    std::is_same<decltype(t -= 50_ms), util::quantities::seconds&>(),
    "Decrement converts to a different type!"
    );
  BOOST_CHECK_EQUAL(t, 0.05_s);
  
} // test_quantities_multiply_scalar()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_quantities_scale() {
  util::quantities::microseconds t { 11.0 };
  //
  // scaling
  //
  t *= 2.0;
  static_assert(
    std::is_same<decltype(t *= 2.0), util::quantities::microseconds&>(),
    "Scaling converts to a different type!"
    );
  BOOST_CHECK_EQUAL(t, 22.0);
  
  t /= 2.0;
  static_assert(
    std::is_same<decltype(t /= 2.0), util::quantities::microseconds&>(),
    "Scaling (division) converts to a different type!"
    );
  BOOST_CHECK_EQUAL(t, 11.0);
  
} // test_quantities_scale()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_quantities_literals() {
  
  using namespace util::quantities::time_literals;
  
  constexpr util::quantities::second t1 = 7_s;
  static_assert(t1 == 7.0, "Literal conversion failed.");
  
  constexpr util::quantities::microsecond t2 = 7_s;
  static_assert(t2 == 7000000.0, "Literal conversion failed.");
  
  util::quantities::microsecond t3;
  t3 = 7.0_s;
  BOOST_CHECK_EQUAL(t3, 7000000.0);
  BOOST_CHECK_EQUAL(t3, 7000000_us);
  
  static_assert(7000000_us == 7_s, "Literal conversion failed.");
  
} // test_quantities_literals()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_quantities() {
  
  // ---------------------------------------------------------------------------
  // default constructor
  // 
//  BOOST_TEST_CHECKPOINT("Default constructor");
  util::quantities::microseconds t1; // can't do much with this except assigning
  
  // ---------------------------------------------------------------------------
  // assignment
  // 
//  t1 = 4.0; // error!
  t1 = util::quantities::microseconds { 4.0 };
  BOOST_CHECK_EQUAL(std::to_string(t1.unit()), "us");
  BOOST_CHECK_EQUAL(std::to_string(t1), "4.000000 us");
  BOOST_CHECK_EQUAL(t1, 4.0);
  
  // ---------------------------------------------------------------------------
  // value constructor
  // 
  util::quantities::microseconds t2 { 7.0 };
  BOOST_CHECK_EQUAL(t2, 7.0);
  
  
  // ---------------------------------------------------------------------------
  
} // test_quantities()


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void test_constexpr_operations() {
  
  constexpr util::quantities::microseconds t1 { 10.0 };
  constexpr util::quantities::microseconds t2 { 20.0 };
  constexpr util::quantities::nanoseconds t_ns { 500.0 };
  constexpr util::quantities::nanoseconds t1_ns = t1; // convert
  
  static_assert(t1.value() == 10.0, "value()");
  static_assert(double(t1) == 10.0, "explicit conversion to plain number");
  static_assert(t1 == 10.0, "implicit conversion to plain number");
  static_assert(+t1 == 10.0, "unary +");
  static_assert(-t1 == -10.0, "unary -");
  static_assert(t1.abs() == 10.0, "abs()");
  
  static_assert( (t1 == t1   ), "comparison");
  static_assert(!(t1 == t2   ), "comparison");
  static_assert(!(t1 == t_ns ), "comparison");
  static_assert( (t1 == t1_ns), "comparison"); // rounding?
  
  static_assert(!(t1 != t1   ), "comparison");
  static_assert( (t1 != t2   ), "comparison");
  static_assert( (t1 != t_ns ), "comparison");
  static_assert(!(t1 != t1_ns), "comparison"); // rounding?
 
  static_assert( (t1 >= t1   ), "comparison");
  static_assert(!(t1 >= t2   ), "comparison");
  static_assert( (t1 >= t_ns ), "comparison");
  static_assert( (t1 >= t1_ns), "comparison"); // rounding?
  
  static_assert(!(t1 <  t1   ), "comparison");
  static_assert( (t1 <  t2   ), "comparison");
  static_assert(!(t1 <  t_ns ), "comparison");
  static_assert(!(t1 <  t1_ns), "comparison"); // rounding?
  
  static_assert( (t1 <= t1   ), "comparison");
  static_assert( (t1 <= t2   ), "comparison");
  static_assert(!(t1 <= t_ns ), "comparison");
  static_assert( (t1 <= t1_ns), "comparison"); // rounding?
  
  static_assert(!(t1 >  t1   ), "comparison");
  static_assert(!(t1 >  t2   ), "comparison");
  static_assert( (t1 >  t_ns ), "comparison");
  static_assert(!(t1 >  t1_ns), "comparison"); // rounding?
  
  static_assert(t1 * 2.0 == 20.0, "scaling");
  static_assert(2.0 * t1 == 20.0, "scaling");
  static_assert(t1 / 2.0 == 5.0, "scalirng");
  
  
  // ---------------------------------------------------------------------------
  
} // test_constexpr_operations()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(quantities_testcase) {
  
  test_quantities();
  test_quantities_sign();
  test_quantities_multiply_scalar();
  test_quantities_addition();
  test_quantities_increment();
  test_quantities_scale();
  test_quantities_conversions();
  test_quantities_comparisons();
  
  test_quantities_literals();
  
  test_constexpr_operations();
  
} // BOOST_AUTO_TEST_CASE(quantities_testcase)

// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
