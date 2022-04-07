/**
 * @file KeyValueParser_test.cc
 * @brief Unit test `icarus::details::KeyValueParser`
 * @date May 18, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see icaruscode/PMT/Algorithms/KeyValueParser.h
 * 
 */

// ICARUS libraries
#include "icaruscode/PMT/Algorithms/KeyValueParser.h"
#include "icaruscode/Decode/DecoderTools/details/KeyValuesData.h"

// Boost libraries
#define BOOST_TEST_MODULE ( KeyValueParser_test )
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <cstdint> // std::size_t


// -----------------------------------------------------------------------------
// --- implementation detail tests
// -----------------------------------------------------------------------------
void SPR_test() {
  
  std::string const info { R"(

# SPR test input file

Description: "
This is a test for the key-values parser with default settings.
It is expected to be used to describe the Single Photoelectron Response.
"

Contact: Gianluca Petrillo (petrillo@slac.stanford.edu)

Gain: 9.7e6  # from amplitude 4 mV
Tick: '2 ns'
Samples: 0.0 1.0 2.5 \
         4.5 3.0 2.5
Samples:+1.8 1.6 1.2 0.8 0.8 0.7 0.7 0.6
)" };
  
  icarus::details::KeyValueParser const parser;
  
  std::istringstream in { info };
  icarus::KeyValuesData const data { parser(in) };
  
  std::cout << "Input:\n" << std::string(80, '-') << "\n"
    << info
    << "\n" << std::string(80, '-')
    << "\nParsed:\n" << data
    << std::endl;
  
  BOOST_TEST(data.hasItem("Description"));
  if (data.hasItem("Description")) {
    icarus::KeyValuesData::Item const& item = data.getItem("Description");
    BOOST_TEST(item.nValues() == 1);
    BOOST_TEST(item.values[0] == "\n"
      "This is a test for the key-values parser with default settings.\n"
      "It is expected to be used to describe the Single Photoelectron Response.\n"
      );
  }
  
  BOOST_TEST(data.hasItem("Contact"));
  if (data.hasItem("Contact")) {
    icarus::KeyValuesData::Item const& item = data.getItem("Contact");
    BOOST_TEST(item.nValues() == 3);
    BOOST_TEST(item.values[0] == "Gianluca");
    BOOST_TEST(item.values[1] == "Petrillo");
    BOOST_TEST(item.values[2] == "(petrillo@slac.stanford.edu)");
  }
  
  BOOST_TEST(data.hasItem("Gain"));
  if (data.hasItem("Gain")) {
    icarus::KeyValuesData::Item const& item = data.getItem("Gain");
    BOOST_TEST(item.nValues() == 1);
    BOOST_TEST
      (item.getNumber<double>(0) == 9.7e6, boost::test_tools::tolerance(0.001));
  }
  
  BOOST_TEST(data.hasItem("Tick"));
  if (data.hasItem("Tick")) {
    icarus::KeyValuesData::Item const& item = data.getItem("Tick");
    BOOST_TEST(item.nValues() == 1);
    BOOST_TEST(item.values[0] == "2 ns");
  }
  
  BOOST_TEST(data.hasItem("Samples"));
  if (data.hasItem("Samples")) {
    icarus::KeyValuesData::Item const& item = data.getItem("Samples");
    std::array const samples
      { 0.0, 1.0, 2.5, 4.5, 3.0, 2.5, 1.8, 1.6, 1.2, 0.8, 0.8, 0.7, 0.7, 0.6 };
    std::vector<double> values;
    for (std::size_t i = 0; i < item.nValues(); ++i)
      values.push_back(item.getNumber<double>(i));
    BOOST_CHECK_EQUAL_COLLECTIONS
      (values.cbegin(), values.cend(), samples.cbegin(), samples.cend());
  }
  
} // SPR_test()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(parse_testcase) {
  
  SPR_test();
  
} // BOOST_AUTO_TEST_CASE(parse_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
