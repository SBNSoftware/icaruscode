/**
 * @file icaruscode/Decode/ChannelMapping/PositionFinder.h
 * @brief Simple utility for finding multiple keys with a single code line.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 */

#ifndef ICARUSCODE_DECODE_CHANNELMAPPING_POSITIONFINDER_H
#define ICARUSCODE_DECODE_CHANNELMAPPING_POSITIONFINDER_H

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// C++ standard libraries
#include <array>
#include <limits> // std::numeric_limits<>
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace icarus::ns::util { template <typename Coll> class PositionFinder; }
/**
 * @brief Returns the index of strings in a set collection.
 * @tparam Coll type of collection of keys being parsed
 * 
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * PositionFinder const namePositions{ Names };
 * auto const [ L1, L5, L3 ] = namePositions("one", "five", "three");
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will set `L1`, `L5` and `L3` to `0`, `4` and `2` respectively.
 * Indices of names that are not present are returned as `NotAvailable`.
 */
template <typename Coll>
class icarus::ns::util::PositionFinder {
  
  Coll const& fKeys; ///< Reference to an external pool of keys.
  
  /// Inserts the position of `first` and `others...` starting at `*position`.
  template<typename String, typename... OtherStrings>
  void fillPositions(
    std::size_t* position, String const& first, OtherStrings const&... others
    ) const;
  
    public:
  
  /// Special value used as "index" for targets that are not in the pool.
  static constexpr std::size_t NotAvailable
     = std::numeric_limits<std::size_t>::max();
  
  /// Constructor: uses `keys` as pool to find elements from.
  PositionFinder(Coll const& keys): fKeys{ keys } {}
  
  /// Returns the position in the pool of `target` or `NotAvailable`.
  template <typename String>
  std::size_t find(String const& target) const;
  
  /// Returns in an array the position of all targets specified in `selection`.
  template <typename... Strings>
  std::array<std::size_t, sizeof...(Strings)> operator()
    (Strings const&... selection) const;
  
}; // icarus::ns::util::PositionFinder<>


// -----------------------------------------------------------------------------
// ---  icarus::ns::util::PositionFinder implementation
// -----------------------------------------------------------------------------
template<typename Coll>
template<typename String, typename... OtherStrings>
void icarus::ns::util::PositionFinder<Coll>::fillPositions(
  std::size_t* position, String const& first, OtherStrings const&... others
) const {
  *position = find(first);
  if constexpr(sizeof...(OtherStrings) > 0)
    fillPositions(++position, others...);
} // icarus::ns::util::PositionFinder<Coll>::fillPositions()


// -----------------------------------------------------------------------------
template <typename Coll>
template <typename String>
std::size_t icarus::ns::util::PositionFinder<Coll>::find
  (String const& target) const
{
  for (auto const& [ index, key ]: ::util::enumerate(fKeys))
    if (key == target) return index;
  return NotAvailable;
} // icarus::ns::util::PositionFinder<Coll>::find()


// -----------------------------------------------------------------------------
template <typename Coll>
template <typename... Strings>
std::array<std::size_t, sizeof...(Strings)>
icarus::ns::util::PositionFinder<Coll>::operator()
  (Strings const&... selection) const
{
  std::array<std::size_t, sizeof...(Strings)> positions;
  fillPositions(positions.data(), selection...);
  return positions;
} // icarus::ns::util::PositionFinder<Coll>::operator()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_DECODE_CHANNELMAPPING_POSITIONFINDER_H
