/**
 * @file    icaruscode/Utilities/ArtAssociationCaches.h
 * @brief   Caching utilities for _art_ associations.
 * @author  Gianluca Petrillo (SLAC, petrillo@slac.stanford.edu)
 * @date    March 5, 2023
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_ARTASSOCIATIONCACHES_H
#define ICARUSCODE_UTILITIES_ARTASSOCIATIONCACHES_H


// LArSoft and ICARUS libraries
#include "icaruscode/Utilities/ComparisonHelpers.h"
#include "larcorealg/CoreUtils/get_elements.h" // util::get_const_elements()

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/ProductID.h"

// C/C++ libraries
#include <algorithm> // std::inplace_merge(), std::binary_search(), ...
#include <vector>
#include <utility> // std::move(), ...


// -----------------------------------------------------------------------------
namespace util {
  template <typename Left, typename Right> class OneToOneAssociationCache;
} // namespace util
/**
 * @brief Cache with fast lookup for associations to _art_ pointers.
 * @tparam Left the type for the query
 * @tparam Right the associated type to discover
 * 
 * This object is functionally similar to `art::FindOneP` in that it presents
 * the information from an _art_ association.
 * The main difference is that `art::FindOneP` requires to know at construction
 * all the pointers we are interested in, and then it accesses them by the index
 * in the interesting pointers list. `OneToOneAssociationCache` instead reads in
 * all the available associations and keeps them available for lookup based on
 * the full _art_ pointer.
 * `art::FindOneP` should be preferred whenever it is suitable, since it is
 * faster. On the other end, `OneToOneAssociationCache` is more flexible in that
 * queries don't need to rely on the knowledge of the index of the pointer at
 * construction time, but it's slower (lookup takes logarithmic rather than
 * constant time).
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * util::OneToOneAssociationCache const cache
 *   { event.getProduct<art::Assns<recob::Track, anab::T0>>(T0tag) };
 * 
 * std::cout << selectedTracks.size() << " tracks selected:";
 * for (art::Ptr<recob::Track> const& trackPtr: selectedTracks) {
 *   art::Ptr<anab::T0> const& t0ptr = cache(trackPtr);
 *   std::cout << "\n  track ID=" << trackPtr->ID();
 *   if (t0ptr) std::cout << " pinned to time " << t0ptr->Time() << " ns";
 *   else std::cout << " time is not known";
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Note that, assuming `selectedTracks` is an object of type
 * `std::vector<art::Ptr<recob::Track>>`, this specific example could have been
 * better written with `art::FindOneP cache{ selectedTracks, event, T0tag };`
 * (or just `art::FindOne`), but the lookup should then track the index of the
 * track in `selectedTracks` itself.
 * 
 */
template <typename Left, typename Right>
class util::OneToOneAssociationCache {
  
    public:
  using Assns_t = art::Assns<Left, Right>;
  using Cache_t = OneToOneAssociationCache<Left, Right>;
  
  /**
   * @brief Default constructor: an empty cache.
   * 
   * The only supported ways to fill a default-constructed cache are assignment
   * and a call to `mergeCache()` from an existing cache object.
   */
  OneToOneAssociationCache() = default;
  
  /// Constructor: reads the associations from `assns`.
  OneToOneAssociationCache(Assns_t const& assns);
  
  /**
   * @brief Constructor from a sequence.
   * @tparam BIter type of begin iterator
   * @tparam EIter type of end iterator
   * @tparam begin begin forward iterator of the input sequence
   * @tparam end end forward iterator of the input sequence
   * 
   * This constructor copies the associations from the `begin`/`end` sequence.
   * The elements of the sequence shoudl be pairs of _art_ pointers, the first
   * with the left object and the second with the right object associated to it.
   * 
   * Requirements
   * -------------
   * 
   * The value pointed by the iterators should be such that applying to it
   * `get<0>()` and `get<1>()` returns an _art_ pointer to the left and right
   * type respectively (or something that can be converted into it).
   * 
   * Sequences of _art_ pointer pairs and `art::Assns` have iterators compatible
   * with these requirements.
   */
  template <typename BIter, typename EIter>
  OneToOneAssociationCache(BIter begin, EIter end);
  
  
  /// Copies the entries of the `other` cache into this one.
  void mergeCache(Cache_t const& other);
  
  
  // --- BEGIN -- Query interface ----------------------------------------------
  /// Returns whether the cache has no element in.
  bool empty() const;
  
  /// Returns whether a pointer with the given product `id` is in cache.
  bool hasProduct(art::ProductID const& id) const;

  /// Returns whether a pointer from the same data product as `ptr` is in cache.
  bool hasProduct(art::Ptr<Left> const& ptr) const;

  /// Returns the "right" object associated to the "left" `ptr`.
  art::Ptr<Right> lookup(art::Ptr<Left> const& ptr) const;
  
  /// Returns the "right" object associated to the "left" `ptr`.
  art::Ptr<Right> operator() (art::Ptr<Left> const& ptr) const;
  
  // --- END ---- Query interface ----------------------------------------------
  
    private:
  
  /// The current content of the cache.
  std::vector<typename Assns_t::assn_t> fCache;
  std::vector<art::ProductID> fProducts; ///< Product IDs featured in the cache.
  
  void prepareCache();
  
  /// Updates the internal record of product IDs.
  void updateProducts();
  
}; // class util::OneToOneAssociationCache


namespace util {
  // deduction guide: pick the types from the ones of the association
  template <typename Left, typename Right>
  OneToOneAssociationCache(art::Assns<Left, Right> const& assns)
    -> OneToOneAssociationCache<Left, Right>;
}


// -----------------------------------------------------------------------------
// --- Template implementation
// -----------------------------------------------------------------------------
// --- util::OneToOneAssociationCache
// -----------------------------------------------------------------------------
template <typename Left, typename Right>
util::OneToOneAssociationCache<Left, Right>::OneToOneAssociationCache
  (Assns_t const& assns)
  : OneToOneAssociationCache{ assns.begin(), assns.end() }
  {}


// -----------------------------------------------------------------------------
template <typename Left, typename Right>
template <typename BIter, typename EIter>
util::OneToOneAssociationCache<Left, Right>::OneToOneAssociationCache
  (BIter begin, EIter end)
{
  using std::get;
  for (auto it = begin; it != end; ++it)
    fCache.emplace_back(get<0>(*it), get<1>(*it));
  prepareCache();
} // util::OneToOneAssociationCache<>::OneToOneAssociationCache(Iter)


// -----------------------------------------------------------------------------
template <typename Left, typename Right>
void util::OneToOneAssociationCache<Left, Right>::mergeCache
  (Cache_t const& other)
{
  // associations
  fCache.reserve(fCache.size() + other.fCache.size());
  auto const itNew = fCache.end();
  fCache.insert(itNew, other.fCache.begin(), other.fCache.end());
  std::inplace_merge(fCache.begin(), itNew, fCache.end());
  // product IDs
  std::vector<art::ProductID> newProducts;
  std::set_union(
    fProducts.cbegin(), fProducts.cend(),
    other.fProducts.begin(), other.fProducts.end(),
    back_inserter(newProducts)
    );
  fProducts = std::move(newProducts);
} // util::OneToOneAssociationCache<>::mergeCache()


// -----------------------------------------------------------------------------
template <typename Left, typename Right>
bool util::OneToOneAssociationCache<Left, Right>::empty() const
  { return fCache.empty(); }


// -----------------------------------------------------------------------------
template <typename Left, typename Right>
bool util::OneToOneAssociationCache<Left, Right>::hasProduct
  (art::ProductID const& id) const
{
  return std::binary_search(fProducts.begin(), fProducts.end(), id);
}


template <typename Left, typename Right>
bool util::OneToOneAssociationCache<Left, Right>::hasProduct
  (art::Ptr<Left> const& ptr) const
  { return hasProduct(ptr.id()); }


// -----------------------------------------------------------------------------
template <typename Left, typename Right>
auto util::OneToOneAssociationCache<Left, Right>::lookup
  (art::Ptr<Left> const& ptr) const -> art::Ptr<Right>
{
  auto const cend = fCache.end();
  auto const it = std::lower_bound(
    fCache.begin(), cend, ptr,
    util::ElementComp<art::Ptr<Left>, 0>() // compare with Ptr from element 0
    );
  if ((it == cend) || (it->first != ptr)) return {};
  return it->second;
} // util::OneToOneAssociationCache<>::lookup()


template <typename Left, typename Right>
auto util::OneToOneAssociationCache<Left, Right>::operator()
  (art::Ptr<Left> const& ptr) const -> art::Ptr<Right>
  { return lookup(ptr); }


// -----------------------------------------------------------------------------
template <typename Left, typename Right>
void util::OneToOneAssociationCache<Left, Right>::prepareCache() {
  std::sort(fCache.begin(), fCache.end());
  updateProducts();
}


// -----------------------------------------------------------------------------
template <typename Left, typename Right>
void util::OneToOneAssociationCache<Left, Right>::updateProducts() {
  fProducts.clear();
  if (fCache.empty()) return;
  if (fCache.front().first == fCache.back().first) // worth a special case
    fProducts.push_back(fCache.front().first.id());
  else {
    auto it = fCache.cbegin();
    auto const cend = fCache.cend();
    fProducts.push_back(it->first.id());
    while (++it != cend) {
      if (it->first.id() != fProducts.back())
        fProducts.push_back(it->first.id());
    }
  }
} // util::OneToOneAssociationCache<>::updateProducts()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_UTILITIES_ARTASSOCIATIONCACHES_H
