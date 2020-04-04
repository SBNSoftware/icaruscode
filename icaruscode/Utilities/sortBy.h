/**
 * @file   icaruscode/Utilities/sortBy.h
 * @brief  Provides `sortBy()` class of utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 3, 2020
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_SORTBY_H
#define ICARUSCODE_UTILITIES_SORTBY_H


// C/C++ standard libraries
#include <algorithm> // std::transform(), std::sort()
#include <functional> // std::less<>
#include <vector>
#include <iterator> // std::back_inserter()
#include <utility> // std::pair>?
#include <type_traits> // std::is_base_of_v


// -----------------------------------------------------------------------------
namespace util {
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns a vectors to pointers to `coll` elements, sorted by `key`.
   * @tparam BIter type of begin iterator to objects to be sorted
   * @tparam EIter type of end iterator to objects to be sorted
   * @tparam Key type of functor extracting the key from an element
   * @tparam Sorter (default: `std::less`) type of functor comparing two keys
   * @param coll collection of objects to be sorted
   * @param key functor extracting the key from an element
   * @param sorter (default: `std::less{}`) functor comparing two keys
   * @return a vector of pointers to `coll` elements, sorted by `key`
   * @see `sortCollBy()`
   * 
   * A vector of pointers to all elements of `coll` is returned.
   * The pointers are constant only if `Coll` is a constant type.
   * The order of the pointed elements is driven by `sorter` applied to the
   * key of each element (`key(item)`).
   * 
   * @note As an exception, if the elements of `Coll` are already C pointers,
   *       the returned collection is a copy of those pointers rather than
   *       pointers to them.
   */
  template <
    typename BIter, typename EIter,
    typename Key, typename Sorter = std::less<void>
    >
  auto sortBy(BIter begin, EIter end, Key key, Sorter sorter = {});

  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns a vectors to pointers to `coll` elements, sorted by `key`.
   * @tparam Coll type of collection of objects to be sorted
   * @tparam Key type of functor extracting the key from an element
   * @tparam Sorter (default: `std::less`) type of functor comparing two keys
   * @param coll collection of objects to be sorted
   * @param key functor extracting the key from an element
   * @param sorter (default: `std::less{}`) functor comparing two keys
   * @return a vector of pointers to `coll` elements, sorted by `key`
   * @see `sortBy()`
   * 
   * A vector of pointers to all elements of `coll` is returned.
   * The pointers are constant only if `Coll` is a constant type.
   * The order of the pointed elements is driven by `sorter` applied to the
   * key of each element (`key(item)`).
   * 
   * @note As an exception, if the elements of `Coll` are already C pointers,
   *       the returned collection is a copy of those pointers rather than
   *       pointers to them.
   */
  template <typename Coll, typename Key, typename Sorter = std::less<void>>
  auto sortCollBy(Coll& coll, Key key, Sorter sorter = {});
  
  // ---------------------------------------------------------------------------
  
} // namespace util



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
namespace util::details {
  
  template <typename Iter>
  constexpr bool is_random_access_iterator_v = std::is_base_of_v<
    std::random_access_iterator_tag,
    typename std::iterator_traits<Iter>::iterator_category
    >;
  
} // namespace util::details

// -----------------------------------------------------------------------------
template <
  typename BIter, typename EIter,
  typename Key, typename Sorter /* = std::less<void> */
  >
auto util::sortBy
  (BIter begin, EIter end, Key key, Sorter sorter /* = {} */)
{
  
  /*
   * 0. establish whether we are dealing with pointers or not
   * 1. create a collection of pairs { key, pointer to element }
   * 2. sort that collection on the first element (key)
   * 3. create a collection of the pointer to element (second element of the
   *    pairs from the collection just sorted), and return it
   */
  using value_type = typename BIter::value_type;
  
  //
  // 0. establish whether we are dealing with pointers or not
  //
  static constexpr bool isPointer = std::is_pointer_v<value_type>;
  
  using pointer_type = std::conditional_t
    <isPointer, value_type, typename std::iterator_traits<BIter>::pointer>;
  
  using Key_t = std::decay_t<decltype(key(*begin))>;
  using SortingPair_t = std::pair<Key_t, pointer_type>;
  
  //
  // 1. create a collection of pairs { key, pointer to element }
  //
  auto getPointer = [](auto& elem)
    { if constexpr(isPointer) return elem; else return &elem; };
  auto makePair = [&key, getPointer](auto&& item)
    { return SortingPair_t(key(item), getPointer(item)); };
  
  std::vector<SortingPair_t> sortingColl;
  
  // reserve size, but only if to discover the size is fast
  if constexpr(details::is_random_access_iterator_v<BIter>)
    sortingColl.reserve(std::distance(begin, end));
  std::transform(begin, end, back_inserter(sortingColl), makePair);
  
  //
  // 2. sort that collection on the first element (key)
  //
  auto pair_sorter
    = [&sorter](SortingPair_t const& a, SortingPair_t const& b)
      { return sorter(a.first, b.first); }
    ;
  std::sort(sortingColl.begin(), sortingColl.end(), pair_sorter);
  
  //
  // 3. create a collection of the pointer to element (second element of the
  //    pairs from the collection just sorted), and return it
  //
  std::vector<pointer_type> sortedColl;
  sortedColl.reserve(sortingColl.size());
  std::transform(
    sortingColl.begin(), sortingColl.end(), back_inserter(sortedColl),
    [](SortingPair_t const& pair){ return pair.second; }
    );
  
  return sortedColl;
  
} // util::sortBy()


//------------------------------------------------------------------------------
template <typename Coll, typename Key, typename Sorter /* = std::less<void> */>
auto util::sortCollBy(Coll& coll, Key key, Sorter sorter /* = {} */) {
  using std::begin, std::end;
  return sortBy(begin(coll), end(coll), key, sorter);
} // util::sortCollBy()


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_SORTBY_H
