/**
 * @file   NonRandomCounter_test.cc
 * @brief  Unit test for `NonRandomCounter` random engine.
 * @date   February 14, 2020
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    `icaruscode/Utilities/NonRandomCounter.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE NonRandomCounterTest
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// ICARUS libraries
#include "icaruscode/Utilities/NonRandomCounter.h"


//------------------------------------------------------------------------------
void Test() {
  
  constexpr std::size_t N = 1024U;
  constexpr long seed = 3U;
  
  util::NonRandomCounter engine { seed };
  
  std::cout << "engine { " << seed << " }: ";
  engine.showStatus();
  
  BOOST_CHECK_EQUAL(engine.name(), "NonRandomCounter");
  
  // test that the numbers are growing
  std::array<double, N> numbers;
  engine.flatArray(numbers.size(), numbers.data());
  std::cout << "\nengine.flatArray(" << numbers.size() << ", "
    << ((void*) numbers.data()) << ") => ";
  engine.showStatus();
  
  for (std::size_t i = 1U; i < N; ++i)
    BOOST_CHECK_LT(numbers[i - 1U], numbers[i]);
  
  // test that the numbers are reproducible
  engine.setSeed(seed, 0);
  std::cout << "\nengine.setSeed(" << seed << ") => ";
  engine.showStatus();
  
  for (double const expected: numbers)
    BOOST_CHECK_EQUAL(engine.flat(), expected); // same math, same result
  std::cout << "\nengine.flat() x " << numbers.size() << " => ";
  engine.showStatus();
  
  // test the other way to set a seed
  engine.setSeeds(&seed, 1);
  std::cout << "\nengine.setSeeds(" << ((void*) &seed) << ", 1) => ";
  engine.showStatus();
  
  for (double const expected: numbers)
    BOOST_CHECK_EQUAL(engine.flat(), expected); // same math, same result
  std::cout << "\nengine.flat() x " << numbers.size() << " => ";
  engine.showStatus();
  
} // Test()


//------------------------------------------------------------------------------
//---  The tests
//---
BOOST_AUTO_TEST_CASE( TestCase ) {
  
  Test();
  
} // BOOST_AUTO_TEST_CASE( TestCase )

