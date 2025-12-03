/**
 * @file   icaruscode/IcarusObj/ChannelToChannelMap.h
 * @brief  Object storing a channel map.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/IcarusObj/ChannelToChannelMap.tcc
 * 
 */

#ifndef ICARUSCODE_ICARUSOBJ_CHANNELTOCHANNELMAP_H
#define ICARUSCODE_ICARUSOBJ_CHANNELTOCHANNELMAP_H

// C/C++ standard libraries
#include <stdexcept> // std::out_of_range (convenience only)
#include <utility> // std::pair
#include <vector>


// -----------------------------------------------------------------------------
namespace icarus
  { template <typename Key, typename Value> class ChannelToChannelMap; }
/**
 * @brief Data structure connecting a channel ID to one or more other channel IDs.
 * @tparam Key (default: `int`) type of the channel ID used as key of the map
 * @tparam Value (default: `Key`) type of the channel ID used as mapped values
 * 
 * This data object is deliberately generic so that it can be used for different
 * contexts. For that reason it does not use a particular channel ID type.
 * 
 * The object stores its information sorted, and it provides O(log(N)) access.
 * 
 * The only way to populate this object is via `addChannel()` function calls,
 * one per channel.
 */
template <typename Key = int, typename Value = Key>
class icarus::ChannelToChannelMap {
  
    public:
  using KeyChannel_t = Key; ///< Type used for key channel ID.
  using MappedChannel_t = Value; ///< Type used for mapped channel IDs.
  
  ///< Type used for the collection of mapped channel IDs.
  using MappedChannels_t = std::vector<MappedChannel_t>;
  
  
  // the usual C++ container crap (but all constant):
  using key_type               = KeyChannel_t;
  using mapped_type            = MappedChannels_t;
  using key_compare            = std::less<KeyChannel_t>;
  
  using value_type             = std::pair<KeyChannel_t, MappedChannels_t>;
  using size_type              = typename std::vector<value_type>::size_type;
  using difference_type        = typename std::vector<value_type>::difference_type;
  using reference              = typename std::vector<value_type>::const_reference;
  using const_reference        = typename std::vector<value_type>::const_reference;
  using pointer                = typename std::vector<value_type>::const_pointer;
  using const_pointer          = typename std::vector<value_type>::const_pointer;
  using iterator               = typename std::vector<value_type>::const_iterator;
  using const_iterator         = typename std::vector<value_type>::const_iterator;
  using reverse_iterator       = typename std::vector<value_type>::const_reverse_iterator;
  using const_reverse_iterator = typename std::vector<value_type>::const_reverse_iterator;
  using allocator_type         = typename std::vector<value_type>::allocator_type;
  
  
  /**
   * @brief Returns the list of channels associated to `channel`.
   * @param channel the channel to get information about
   * @return list of channels, empty if `channel` is not available
   */
  MappedChannels_t const& operator[] (KeyChannel_t channel) const noexcept;
  
  /**
   * @brief Returns the list of channels associated to `channel`.
   * @param channel the channel to get information about
   * @return list of channels
   * @throw std::out_of_range if `channel` is not available
   */
  MappedChannels_t const& at(KeyChannel_t channel) const;
  
  /// Returns whether there is information on the specified channel.
  bool contains(KeyChannel_t channel) const noexcept;
  
  /// Returns the lowest channel in the map. Undefined behaviour on empty map.
  KeyChannel_t lowestChannel() const;
  
  /// Returns the highest channel in the map. Undefined behaviour on empty map.
  KeyChannel_t highestChannel() const;
  
  /**
   * @brief Adds or replaces an element to the map.
   * @param channel the single channel to map
   * @param toChannels all the channels to map `channel` to
   * @return this object (for chaining `addChannel()` calls)
   * 
   * The map channels are sorted at insertion time.
   */
  ChannelToChannelMap& addChannel
    (KeyChannel_t channel, std::vector<MappedChannel_t> toChannels);
  
  // --- BEGIN ---  Standard container interface  ------------------------------
  /**
   * @name Standard container interface
   * 
   * Well, `std::map`-style standard at least.
   */
  /// @{
  
  const_iterator begin() const noexcept { return fMap.cbegin(); }
  const_iterator cbegin() const noexcept { return fMap.cbegin(); }
  const_reverse_iterator rbegin() const noexcept { return fMap.crbegin(); }
  const_reverse_iterator crbegin() const noexcept { return fMap.crbegin(); }
  
  const_iterator end() const noexcept { return fMap.cend(); }
  const_iterator cend() const noexcept { return fMap.cend(); }
  const_reverse_iterator rend() const noexcept { return fMap.crend(); }
  const_reverse_iterator crend() const noexcept { return fMap.crend(); }
  
  bool empty() const noexcept { return fMap.empty(); }
  const_pointer data() const noexcept { return fMap.data(); }
  
  /// Returns the number of channels we have information about.
  std::size_t size() const noexcept { return fMap.size(); }
  
  /// @}
  // ---- END ----  Standard container interface  ------------------------------
  
  
  /**
   * @brief Convenience pair comparison.
   *
   * Public because sooner or later may be useful.
   * Use like in:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const it = std::lower_bound
   *   (begin(map), end(map), channel, icarus::ChannelToChannelMap<>::CompareByChannel{});
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  struct CompareByChannel {
    static auto channelOf(value_type const& v) noexcept { return v.first; }
    static auto channelOf(key_type k) noexcept { return k; }
    
    template <typename T, typename U>
    bool operator() (T&& a, U&& b) const noexcept
      { return channelOf(a) < channelOf(b); }
  }; // CompareByChannel
  
  
    private:
  
  /// Used to return non-existing channels.
  static MappedChannels_t const NoChannels;
  
  
  std::vector<value_type> fMap; ///< List of (channel, associated channels)
  
  /// Returns a pointer to the pair with `channel`, `nullptr` if not present.
  const_pointer find(KeyChannel_t channel) const noexcept;
  
}; // icarus::ChannelToChannelMap



namespace icarus {
  
  /// @name Convenience free functions for standard containers.
  /// @{
  
  template <typename Key, typename Value>
  auto begin  (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.begin();   }
  template <typename Key, typename Value>
  auto cbegin (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.cbegin();  }
  template <typename Key, typename Value>
  auto rbegin (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.rbegin();  }
  template <typename Key, typename Value>
  auto crbegin(icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.crbegin(); }
  template <typename Key, typename Value>
  auto end    (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.end();     }
  template <typename Key, typename Value>
  auto cend   (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.cend();    }
  template <typename Key, typename Value>
  auto rend   (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.rend();    }
  template <typename Key, typename Value>
  auto crend  (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.crend();   }
  
  template <typename Key, typename Value>
  auto size (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.size();  }
  template <typename Key, typename Value>
  auto empty(icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.empty(); }
  template <typename Key, typename Value>
  auto data (icarus::ChannelToChannelMap<Key, Value> const& map) noexcept { return map.data();  }
  
  /// @}
  
  
} // namespace icarus


// -----------------------------------------------------------------------------
// ---  Inline implementation
// -----------------------------------------------------------------------------

// template implementation
#include "icaruscode/IcarusObj/ChannelToChannelMap.tcc"


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_ICARUSOBJ_CHANNELTOCHANNELMAP_H
