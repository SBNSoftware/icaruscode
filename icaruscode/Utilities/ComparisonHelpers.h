/**
 * @file    icaruscode/Utilities/ComparisonHelpers.h
 * @brief   Small utilities to facilitate sorting and lookup.
 * @author  Gianluca Petrillo (SLAC, petrillo@slac.stanford.edu)
 * @date    March 5, 2023
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_COMPARISONHELPERS_H
#define ICARUSCODE_UTILITIES_COMPARISONHELPERS_H


// C/C++ libraries
#include <utility> // std::move(), std::get()


// -----------------------------------------------------------------------------
namespace util {
  template <typename Key, std::size_t N, typename Comp = std::less<>>
  class ElementComp;
}

/**
 * @brief Functor applying a specified function on a certain tuple element.
 * @tparam Key type of the key comparison
 * @tparam N the number of the element to apply the functor on
 * @tparam Comp the type of the functor
 * 
 * The `Key` type is the type of the element `N` of the tuples being compared.
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::vector<std::pair<int, std::string>> names;
 * 
 * std::sort(names.begin(), names.end()); // sorts by int, then std::string
 * 
 * int code = 7;
 * if (!std::binary_search(names.begin(), names.end(), code, ElementComp<int, 0>{}))
 *   std::cerr << "code=" << code << " not available" << std::endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * @note This comparison object may fail to compile when used with tuple
 *       objects that are convertible to the `Key` type.
 * 
 */
template <typename Key, std::size_t N, typename Comp /* = std::less<> */>
class util::ElementComp {
  
  Comp comp; ///< Local copy of the comparer object.
  
    public:
  
  using Key_t = Key; ///< Type of the value being compared.
  
  using Compare_t = Comp; ///< Type of comparison functor.
  
  /// Index of the tuple being compared.
  static constexpr std::size_t Index = N;
  
  /// Constructor: can specify the comparison between keys.
  ElementComp(Compare_t comp = {}): comp(std::move(comp)) {}
  
  // --- BEGIN -- Comparison operators -----------------------------------------
  /// @name Comparison operators.
  /// @{
  template <typename Tuple>
  constexpr bool operator()(Key_t const& a, Tuple const& b) const noexcept
    { return cmp(a, b, comp); }
  
  template <typename Tuple>
  constexpr bool operator()(Tuple const& a, Key_t const& b) const noexcept
    { return cmp(a, b, comp); }
  
  template <typename Tuple>
  constexpr bool operator()(Key_t const& a, Key_t const& b) const noexcept
    { return cmp(a, b, comp); }
  
  /// @}
  // --- END ---- Comparison operators -----------------------------------------
  
  
  // --- BEGIN -- Static comparisons -------------------------------------------
  /// @name Static comparisons and utilities.
  /// @{
  
  /// Returns the key of a tuple.
  template <typename Tuple>
  static auto getKey(Tuple const& tuple) -> decltype(auto)
    { using std::get; return get<Index>(tuple); }
  
  template <typename Tuple, typename Cmp = Comp>
  static constexpr bool cmp
    (Key_t const& a, Tuple const& b, Cmp comp = {}) noexcept
    { return comp(a, getKey(b)); }
  
  template <typename Tuple, typename Cmp = Comp>
  static constexpr bool cmp
    (Tuple const& a, Key_t const& b, Cmp comp = {}) noexcept
    { return comp(getKey(a), b); }
  
  template <typename Tuple, typename Cmp = Comp>
  static constexpr bool cmp
    (Key_t const& a, Key_t const& b, Cmp comp = {}) noexcept
    { return comp(a, b); }
  
  /// @}
  // --- END ---- Static comparisons -------------------------------------------
  
}; // util::ElementComp<>


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_COMPARISONHELPERS_H
