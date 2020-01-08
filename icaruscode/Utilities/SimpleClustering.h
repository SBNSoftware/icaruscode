/**
 * @file   icaruscode/Utilities/SimpleClustering.h
 * @brief  Algorithms to cluster objects according to specified criteria.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 7, 2020
 *
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_SIMPLECLUSTERING_H
#define ICARUSCODE_UTILITIES_SIMPLECLUSTERING_H

// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::transform(), std::sort()
#include <utility> // std::pair, std::move(), std::declval()
#include <type_traits> // std::decay_t
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace util {

  /**
   * @brief Performs a simple clustering.
   * @tparam Coll type of collection of objects to cluster
   * @tparam KeyOp type of operation extracting the relevant key for clustering
   * @tparam CmpOp type of operation determining if object belongs to a cluster
   * @tparam RefOp type of operation extracting the object to store in cluster
   * @tparam KeySortOp type of operation sorting the clustering keys
   * @param coll collection of objects to cluster
   * @param keyFunc operation extracting the relevant key for clustering
   * @param sameGroup operation determining if an object belongs to a cluster
   * @param objRef operation extracting the object to store in the cluster
   * @param keySort operation sorting the clustering keys
   * @return a STL vector of clusters of object "references"
   *
   * The algorithm clusters objects whose key is compatible.
   * The key must be sortable (`keySort`). Each cluster is assigned a key value,
   * and all the unclustered objects whose key (`keyFunc`) is compatible
   * (`sameGroup`) with that cluster key are added to that cluster.
   * The objects from `coll` are considered in order of key value, with the
   * order being defined by `keySort`.
   * Each cluster contains a derivative of the original object (`objRef`), which
   * may be for example a pointer to the original object, a copy of it, or
   * in fact anything `objRef` returns.
   *
   * The return value is in the form of a `std::vector` in which each element
   * represents a cluster. Each of these clusters is a `std::vector` of
   * the derivative objects.
   *
   *
   * Requirements
   * -------------
   *
   * * `Coll`: an iterable collection of objects of type `Object_t`;
   * * `KeyOp`: a functor returning the key of an object, like
   *     `Key_t keyFunc(Object_t const&)`;
   * * `KeySortOp`: a functor like `bool keySort(Key_t a, Key_t b)`
   *     returning strict ordering between key values `a` and `b`;
   * * `CmpOp`: a functor like `bool keySort(Key_t a, Key_t b)` returning if
   *     an object with key `b` should belong to a cluster with key `a`;
   * * `RefOp`: a functor returning the object to store in the cluster starting
   *     from an object in the original collection: `ObjRef_t objRef(Object_t)`;
   *     note that the result should probably *not* be a C++ reference since
   *     they don't go along well with containers.
   *
   */
  template <
    typename Coll,
    typename KeyOp, typename CmpOp, typename RefOp, typename KeySortOp
    >
  auto clusterBy(
    Coll const& objs,
    KeyOp keyFunc, CmpOp sameGroup, RefOp objRef, KeySortOp keySort
    );

  /// A version of `clusterBy()` storing a copy of each object as "reference".
  template <typename Coll, typename KeyOp, typename CmpOp, typename KeySortOp>
  auto clusterBy
    (Coll const& objs, KeyOp keyFunc, CmpOp sameGroup, KeySortOp keySort);

  // ---------------------------------------------------------------------------

} // namespace util



// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
namespace util::details {

  template <std::size_t I, typename KeySort>
  struct TupleElementOp {
    KeySort sorter;

    TupleElementOp(KeySort keySort): sorter(keySort) {}

    template <typename A, typename B>
    auto operator() (A&& a, B&& b) const
      { return sorter(std::forward<A>(a), std::forward<B>(b)); }

  }; // TupleElementOp<>

  template <std::size_t I, typename KeySort>
  TupleElementOp<I, KeySort> makeTupleElementOp(KeySort keySort)
    { return { keySort }; }

} // namespace util::details


// -----------------------------------------------------------------------------
template <
  typename Coll,
  typename KeyOp, typename CmpOp, typename RefOp, typename KeySortOp
  >
auto util::clusterBy(
  Coll const& objs,
  KeyOp keyFunc, CmpOp sameGroup, RefOp objRef, KeySortOp keySort
) {

  /*
   * 1. create a list of object key (`keyFunc`) and object reference (`objRef`)
   * 2. sort that list by key (`keySort`)
   * 3. cluster object references with nearby keys (`sameGroup`)
   *
   * A new cluster is started every time for a new object `sameGroup` return
   * `false` when applied on the key of the first object in the current cluster
   * and the key of the new object.
   *
   */

  //
  // some definitions
  //
  using Object_t = typename Coll::value_type;
  using ObjectRef_t = std::decay_t<decltype(objRef(std::declval<Object_t>()))>;
  using Cluster_t = std::vector<ObjectRef_t>;
  using Clusters_t = std::vector<Cluster_t>;

  auto makeKeyAndObjRef = [keyFunc, objRef](auto&& obj)
    { return std::make_pair(keyFunc(obj), objRef(obj)); };

  using KeyAndObjRef_t = decltype(makeKeyAndObjRef(*(util::begin(objs))));

  //
  // working data: "map" of objects by key
  //
  std::vector<KeyAndObjRef_t> KeysAndObjRefs;
  KeysAndObjRefs.reserve(objs.size());
  std::transform(
    util::begin(objs), util::end(objs), std::back_inserter(KeysAndObjRefs),
    makeKeyAndObjRef
    );

  // sort the map by key
  std::sort(
    util::begin(KeysAndObjRefs), util::end(KeysAndObjRefs),
    details::makeTupleElementOp<0U>(keySort)
    );

  //
  // cluster the objects in the map by key proximity
  //
  Clusters_t clusters;
  if (KeysAndObjRefs.empty()) return clusters;

  auto iKeyAndObjRef = util::begin(KeysAndObjRefs);
  auto const end = util::end(KeysAndObjRefs);

  auto startCluster = [](auto& keyAndObjRef)
    {
      return
        std::make_pair(keyAndObjRef.first, Cluster_t{ keyAndObjRef.second });
    };
  auto addToCluster = [](auto& cluster, auto& keyAndObjRef)
    { cluster.second.push_back(std::move(keyAndObjRef.second)); };

  auto currentCluster = startCluster(*iKeyAndObjRef);
  while(++iKeyAndObjRef != end) {

    if (sameGroup(iKeyAndObjRef->first, currentCluster.first)) {
      addToCluster(currentCluster, *iKeyAndObjRef);
    }
    else {
      clusters.push_back(std::move(currentCluster.second));
      currentCluster = startCluster(*iKeyAndObjRef);
    }

  } // while
  clusters.push_back(std::move(currentCluster.second));

  return clusters;

} // util::clusterBy(Coll, KeyOp, CmpOp, RefOp, KeySortOp)


// -----------------------------------------------------------------------------
template <typename Coll, typename KeyOp, typename CmpOp, typename KeySortOp>
auto util::clusterBy
  (Coll const& objs, KeyOp keyFunc, CmpOp sameGroup, KeySortOp keySort)
{
  return clusterBy(objs,
    std::move(keyFunc), std::move(sameGroup),
    [](auto const& obj){ return obj; },
    std::move(keySort)
    );
} // util::clusterBy(Coll, KeyOp, CmpOp, KeySortOp)


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_UTILITIES_SIMPLECLUSTERING_H
