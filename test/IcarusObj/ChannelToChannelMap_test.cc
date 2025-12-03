/**
 * @file   test/IcarusObj/ChannelToChannelMap_test.cc
 * @brief  Unit test for `ChannelToChannelMap.h` object.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 4, 2025
 * @see    `icaruscode/IcarusObj/ChannelToChannelMap.h`
 */

// ICARUS libraries
#include "icaruscode/IcarusObj/ChannelToChannelMap.h"
#include "larcorealg/CoreUtils/counter.h"

// Boost libraries
#define BOOST_TEST_MODULE ( ChannelToChannelMap_test )
#include <boost/test/unit_test.hpp>

// C/C++ standard library
#include <iterator> // std::distance()
#include <numeric> // std::iota()
#include <utility> // std::move()
#include <vector>
#include <iostream>


// -----------------------------------------------------------------------------
template <typename Key = int, typename Value = Key>
icarus::ChannelToChannelMap<Key, Value> createAdderMap() {
  using ChannelMap_t = icarus::ChannelToChannelMap<Key, Value>;
  
  ChannelMap_t map;
  
  typename ChannelMap_t::KeyChannel_t PMTchannel = 0;
  for (auto const channel
    : util::counter<typename ChannelMap_t::KeyChannel_t>(0xA000, 0xA018)
  ) {
    typename ChannelMap_t::MappedChannels_t toChannels(15);
    std::iota(
      begin(toChannels), end(toChannels),
      typename ChannelMap_t::MappedChannel_t{ PMTchannel }
      );
    PMTchannel += toChannels.size();
    map.addChannel(channel, std::move(toChannels));
  }
  
  return map;
} // createAdderMap()


// -----------------------------------------------------------------------------
// --- ChannelToChannelMap tests
// -----------------------------------------------------------------------------
void ChannelToChannelMap_empty_test() {
  
  icarus::ChannelToChannelMap<> const map;
  
  BOOST_CHECK(map.empty());
  BOOST_TEST(map.size() == 0);
  BOOST_TEST(std::distance(map.begin(),   map.end()  ) == 0);
  BOOST_TEST(std::distance(map.cbegin(),  map.cend() ) == 0);
  BOOST_TEST(std::distance(map.rbegin(),  map.rend() ) == 0);
  BOOST_TEST(std::distance(map.crbegin(), map.crend()) == 0);
  
  
  BOOST_CHECK(std::empty(map));
  BOOST_TEST(std::size(map) == 0);
  
} // ChannelToChannelMap_empty_test()


template <typename Key, typename Value>
void ChannelToChannelMap_test() {
  
  using ChannelMap_t = icarus::ChannelToChannelMap<Key, Value>;
  
  ChannelMap_t const map = createAdderMap<Key, Value>();
  
  BOOST_TEST(map.size() == 24);
  BOOST_TEST(std::distance(map.begin(),   map.end()  ) == 24);
  BOOST_TEST(std::distance(map.cbegin(),  map.cend() ) == 24);
  BOOST_TEST(std::distance(map.rbegin(),  map.rend() ) == 24);
  BOOST_TEST(std::distance(map.crbegin(), map.crend()) == 24);
  BOOST_TEST(map.lowestChannel() == 0xA000);
  BOOST_TEST(map.highestChannel() == 0xA017);
  
  for (auto const channel
    : util::counter<typename ChannelMap_t::KeyChannel_t>(0xA000, 0xA020)
  ) {
    BOOST_TEST_INFO_SCOPE("Channel: " << channel);
    
    if (channel < 0xA018) {
      typename ChannelMap_t::MappedChannels_t expected(15);
      std::iota(begin(expected), end(expected), (channel - 0xA000) * 15);
      BOOST_CHECK(map.contains(channel));
      BOOST_TEST(map.at(channel) == expected);
    }
    else {
      BOOST_CHECK(!map.contains(channel));
      BOOST_CHECK_THROW(map.at(0), std::out_of_range);
      BOOST_CHECK(map[channel].empty());
    }
  } // for channel
  
  
  BOOST_TEST(std::size(map) == 24);
  
} // ChannelToChannelMap_value_test()


void ChannelToChannelMap_documentation_CompareByChannel_test() {
  
  /* The promise:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const it = std::lower_bound
   *   (begin(map), end(map), channel, icarus::ChannelToChannelMap<>::CompareByChannel{});
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  
  icarus::ChannelToChannelMap const map = createAdderMap();
  
  int const channel = 0xA002;
  auto const it = std::lower_bound
    (begin(map), end(map), channel, icarus::ChannelToChannelMap<>::CompareByChannel{});
  
  BOOST_CHECK((it != map.end()));
  BOOST_TEST(it->first == channel);
  BOOST_TEST(it->second.size() == 15);
  
} // ChannelToChannelMap_documentation_CompareByChannel_test()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(ChannelToChannelMap_testcase) {
  
  ChannelToChannelMap_empty_test();
  BOOST_TEST_CONTEXT("Test with signed int keys") {
    ChannelToChannelMap_test<signed int, signed int>();
  }
  BOOST_TEST_CONTEXT("Test with unsigned int keys") {
    ChannelToChannelMap_test<unsigned int, unsigned int>();
  }
  
} // BOOST_AUTO_TEST_CASE(ChannelToChannelMap_testcase)


BOOST_AUTO_TEST_CASE(ChannelToChannelMap_documentation_testcase) {
  
  ChannelToChannelMap_documentation_CompareByChannel_test();
  
} // BOOST_AUTO_TEST_CASE(ChannelToChannelMap_documentation_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
