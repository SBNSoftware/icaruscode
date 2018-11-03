/**
 * @file   test/Utilities/datasize_test.cc
 * @brief  Unit test for `icaruscode/Utilities/quantities/datasize.h` header.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 2, 2018
 * @see    `icaruscode/Utilities/quantities.h`
 *
 */

// Boost libraries
#define BOOST_TEST_MODULE ( datasize_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "icaruscode/Utilities/quantities/datasize.h"

// C/C++ standard libraries
#include <string>
#include <iostream>
#include <type_traits> // std::decay_t<>


// -----------------------------------------------------------------------------
void test_datasize_literals() {
  
  using namespace util::quantities::datasize_literals;
  
  constexpr auto s_B = 256_B;
  static_assert(std::is_same<decltype(s_B), util::quantities::byte const>());
  static_assert(s_B == 256);
  static_assert(s_B == 256_B);
  std::cout << "Tested " << s_B << std::endl;
  
  constexpr auto s_kiB = 4_kiB;
  static_assert
    (std::is_same<decltype(s_kiB), util::quantities::kibibyte const>());
  static_assert(s_kiB == 4_kiB);
  static_assert(s_kiB == 4);
  static_assert(s_kiB == 4096_B);
  std::cout << "Tested " << s_kiB << std::endl;
  
  constexpr auto s_MiB = 4_MiB;
  static_assert
    (std::is_same<decltype(s_MiB), util::quantities::mebibyte const>());
  static_assert(s_MiB == 4_MiB);
  static_assert(s_MiB == 4);
  static_assert(s_MiB == 4096_kiB);
  std::cout << "Tested " << s_MiB << std::endl;
  
  constexpr auto s_GiB = 4_GiB;
  static_assert
    (std::is_same<decltype(s_GiB), util::quantities::gibibyte const>());
  static_assert(s_GiB == 4_GiB);
  static_assert(s_GiB == 4);
  static_assert(s_GiB == 4096_MiB);
  std::cout << "Tested " << s_GiB << std::endl;
  
  constexpr auto s_TiB = 4_TiB;
  static_assert
    (std::is_same<decltype(s_TiB), util::quantities::tebibyte const>());
  static_assert(s_TiB == 4_TiB);
  static_assert(s_TiB == 4);
  static_assert(s_TiB == 4096_GiB);
  std::cout << "Tested " << s_TiB << std::endl;
  
  constexpr auto s_PiB = 4_PiB;
  static_assert
    (std::is_same<decltype(s_PiB), util::quantities::pebibyte const>());
  static_assert(s_PiB == 4_PiB);
  static_assert(s_PiB == 4);
  static_assert(s_PiB == 4096_TiB);
  std::cout << "Tested " << s_PiB << std::endl;
  
  constexpr auto s_EiB = 4_EiB;
  static_assert
    (std::is_same<decltype(s_EiB), util::quantities::exbibyte const>());
  static_assert(s_EiB == 4_EiB);
  static_assert(s_EiB == 4);
  static_assert(s_EiB == 4096_PiB);
  std::cout << "Tested " << s_EiB << std::endl;
  
} // test_datasize_literals()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(datasize_testcase) {
  
  test_datasize_literals();
  
} // BOOST_AUTO_TEST_CASE(datasize_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
