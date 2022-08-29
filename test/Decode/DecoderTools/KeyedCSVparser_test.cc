/**
 * @file   test/Decode/DecoderTools/KeyedCSVparser_test.cc
 * @brief  Unit test for `KeyedCSVparser.h` header.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 4, 2022
 * @see    `icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h`
 *
 */

// ICARUS libraries
#include "icaruscode/Decode/DecoderTools/details/KeyedCSVparser.h"

// Boost libraries
#define BOOST_TEST_MODULE ( KeyedCSVparser_test )
#include <boost/test/unit_test.hpp>

// C/C++ standard library
#include <string>
#include <vector>
#include <cstdint> // std::uint32_t


// -----------------------------------------------------------------------------
// --- KeyedCSVparser tests
// -----------------------------------------------------------------------------

void KeyedCSVparser_documentation_test() {
  
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * using namespace std::string_view_literals;
   * icarus::details::KeyedCSVparser parser;
   * parser.addPatterns({
   *       { "Trigger Type", 1U } // expect one value (even if contains letters)
   *     , { "TriggerWindows", 1U } // expect one value (even if contains letters)
   *     , { "TPChitTimes", icarus::details::KeyedCSVparser::FixedSize }
   *          // the first value is an integer, count of how many other values
   *   });
   * 
   * icarus::KeyValuesData data = parser(
   *   "TriggerType, S5, Triggers, TriggerWindows, 0C0B,"
   *   " TPChits, 12, 130, 0, 0, TPChitTimes, 3, -1.1, -0.3, 0.1, PMThits, 8"sv
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  
  using namespace std::string_literals;
  icarus::details::KeyedCSVparser parser;
  parser.addPatterns({
        { "TriggerType", 1U } // expect one value (even if contains letters)
      , { "TriggerWindows", 1U } // expect one value (even if contains letters)
      , { "TPChitTimes", icarus::details::KeyedCSVparser::FixedSize }
           // the first value is an integer, count of how many other values
    });
  
  icarus::KeyValuesData data = parser(
    "TriggerType, S5, Triggers, TriggerWindows, 0C0B,"
    " TPChits, 12, 130, 0, 0, TPChitTimes, 3, -1.1, -0.3, 0.1, PMThits, 8"s
    );

  // deep test:
  std::vector<int> const expectedTPChits { 12, 130, 0, 0 };
  std::vector<float> const expectedTPCtimes { -1.1, -0.3, 0.1 };
  
  icarus::KeyValuesData::Item const* item = nullptr;
  
  std::cout << data << std::endl;
  
  BOOST_TEST(data.hasItem("TriggerType"));
  item = data.findItem("TriggerType");
  BOOST_TEST(item);
  BOOST_TEST(&(data.getItem("TriggerType")) == item);
  BOOST_TEST(item->nValues() == 1);
  BOOST_TEST(item->value() == "S5");
  
  BOOST_TEST(data.hasItem("Triggers"));
  item = data.findItem("Triggers");
  BOOST_TEST(item);
  BOOST_TEST(&(data.getItem("Triggers")) == item);
  BOOST_TEST(item->values().empty());
  BOOST_TEST(item->getVector<int>().empty());
  
  BOOST_TEST(data.hasItem("TriggerWindows"));
  item = data.findItem("TriggerWindows");
  BOOST_TEST(item);
  BOOST_TEST(&(data.getItem("TriggerWindows")) == item);
  BOOST_TEST(item->nValues() == 1);
  BOOST_TEST(item->value(0) == "0C0B");
  BOOST_TEST(item->getNumber<std::uint32_t>(0, 16) == 0x0C0B);
  
  BOOST_TEST(data.hasItem("TPChits"));
  item = data.findItem("TPChits");
  BOOST_TEST(item);
  BOOST_TEST(&(data.getItem("TPChits")) == item);
  BOOST_TEST(item->nValues() == 4);
  BOOST_TEST(item->value(0) == "12");
  BOOST_TEST(item->value(1) == "130");
  BOOST_TEST(item->value(2) == "0");
  BOOST_TEST(item->value(3) == "0");
  BOOST_TEST(item->getVector<int>() == expectedTPChits);
  
  BOOST_TEST(data.hasItem("TPChitTimes"));
  item = data.findItem("TPChitTimes");
  BOOST_TEST(item);
  BOOST_TEST(&(data.getItem("TPChitTimes")) == item);
  BOOST_TEST(item->nValues() == 4);
  BOOST_TEST(item->value(0) == "3");
  BOOST_TEST(item->value(1) == "-1.1");
  BOOST_TEST(item->value(2) == "-0.3");
  BOOST_TEST(item->value(3) == "0.1");
  BOOST_TEST(item->getSizedVector<float>() == expectedTPCtimes);
  
  BOOST_TEST(!data.hasItem("CRThits"));
  item = data.findItem("CRThits");
  BOOST_TEST(!item);
  
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * std::string triggerType = data.getItem("TriggerType").value(0);
   * std::vector<int> triggers = data.getItem("Triggers").getVector<int>();
   * std::uint32_t triggerWindowBits
   *  = data.getItem("TriggerWindows").getNumber<std::uint32_t>(0, 16); // base 16
   * std::vector<int> TPChits = data.getItem("TPChits").getVector<int>();
   * std::vector<float> TPCtimes
   *  = data.getItem("TPChitTimes").getSizedVector<float>();
   * std::vector<int> CRThits;
   * if (auto const* item = data.findItem("CRThits"))
   *   CRThits = item->getVector<int>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  
  std::string triggerType = data.getItem("TriggerType").value(0);
  BOOST_TEST(triggerType == "S5");
  
  std::vector<int> triggers = data.getItem("Triggers").getVector<int>();
  BOOST_TEST(triggers.empty());
  
  std::uint32_t triggerWindowBits
   = data.getItem("TriggerWindows").getNumber<std::uint32_t>(0, 16); // base 16
  BOOST_TEST(triggerWindowBits == 0x0C0B);
  
  std::vector<int> TPChits = data.getItem("TPChits").getVector<int>();
  BOOST_TEST(TPChits == expectedTPChits);
  
  std::vector<float> TPCtimes
   = data.getItem("TPChitTimes").getSizedVector<float>();
  BOOST_TEST(TPCtimes == expectedTPCtimes);
  
  std::vector<int> CRThits;
  BOOST_TEST(!data.hasItem("CRThits"));
  BOOST_TEST(!data.findItem("CRThits"));
  if (auto const* item = data.findItem("CRThits"))
    CRThits = item->getVector<int>();
  BOOST_TEST(CRThits.empty());
  
} // KeyedCSVparser_documentation_test()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(KeyedCSVparser_documentation_testcase) {
  
  KeyedCSVparser_documentation_test();
  
} // BOOST_AUTO_TEST_CASE(KeyedCSVparser_documentation_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
