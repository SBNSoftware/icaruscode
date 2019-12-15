/**
 * @file   icaruscode/Utilities/DataProductPointerMap.h
 * @brief  Utilities to map data pointer elements to their `art::Ptr`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 13, 2019
 *
 * This library is header only.
 */

#ifndef ICARUSCODE_UTILITIES_DATAPRODUCTPOINTERMAP_H
#define ICARUSCODE_UTILITIES_DATAPRODUCTPOINTERMAP_H


// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Ptr.h"

// C/C++ standard libraries
#include <map>
#include <vector>
#include <type_traits> // std::is_same_v


namespace util {

  // ---------------------------------------------------------------------------
  namespace details {

    // this is just to stress that the fact this is a STL map may change.
    template <typename T>
    using DataProductPointerMap_t = std::map<T const*, art::Ptr<T>>;

  } // namespace details

  // ---------------------------------------------------------------------------
  /// Type of data in a _art_ handle to vector data product.
  template <typename Handle>
  using ArtHandleData_t = typename Handle::element_type::value_type;

  /// Type of map for data product pointers to _art_ pointers.
  template <typename T>
  using DataProductPointerMap_t = details::DataProductPointerMap_t<T>;


  // ---------------------------------------------------------------------------
  /**
   * @brief Creates a map from address of data product element to _art_ pointer
   *        to it.
   * @tparam Handle type of handle to data product (e.g. `art::ValidHandle`)
   * @param handle _art_ handle to the data product
   * @param event the _art_ event the data product belongs to
   * @return a map from data product element pointer to _art_ pointer
   *
   * Returns a map from the address of any of the objects in the data product
   * in `event` pointed by `handle` and the _art_ pointer to those objects.
   * Provided that the type in the data product is `T`, the returned map is an
   * object guaranteed to have:
   *  * standard mobility (copy and move constructors and assignment operators)
   *  * an unckecked access operator `operator[] (T const*) const` returning
   *    some form of `art::Ptr<T>` (`art::Ptr<T>`, `art::Ptr<T> const&`...);
   *  * a checked access operator `at(T const*) const` similar to `operator[]`
   *    but throwing an exception derived from `std::out_of_range` if the
   *    element is not found;
   *  * an `empty()`, a `size()` and a `clear()` method.
   *
   * Example: with `event` a `art::Event` object and `waveformHandle` an handle
   * to a `std::vector<raw::OpDetWaveform>` data product in `event`:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * DataProductPointerMap_t<raw::OpDetWaveform> const& opDetWavePtrs
   *   = util::mapDataProductPointers(event, waveformHandle);
   *
   * for (raw::OpDetWaveform const& waveform: *waveformHandle) {
   *
   *   art::Ptr<raw::OpDetWaveform> const& ptr = opDetWavePtrs.at(&waveform);
   *
   *   // ...
   *
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * in the loop the `art::Ptr` of each `waveform` is tracked without carrying
   * around the index of the waveform in the original data product.
   *
   * There are alternative approaches to this: one can rely on the contiguous
   * memory model of `std::vector` and carry around the address of the first
   * waveform, `&(waveformHandle->front())`, and a properly initialized
   * `art::PtrMaker` which can then invoked as
   * `ptrMaker(&waveform - waveform0ptr)`; or work directly with `art::Ptr`.
   * This approach is more factorised than the first alternative (knowledge
   * of the event is not required any more during the iteration) and arguably
   * faster than the second, where _art_ pointer dereferencing always has
   * a small overhead, but it uses more memory.
   *
   * @note The returned object is currently a C++ STL data container.
   *       Optimizations are possible using a different data structure
   *       (for example, a vector-based container with index based on the
   *       difference described in the example above) and they may in future
   *       be implemented. Because of this reason, it is recommended that the
   *       produced map is stored in variables declared with an `auto const&`
   *       type or as a `DataProductPointerMap_t` instance, rather than with
   *       the currently underlying data type.
   */
  template <typename Handle>
  DataProductPointerMap_t<ArtHandleData_t<Handle>> mapDataProductPointers
    (art::Event const& event, Handle const& handle);

  // ---------------------------------------------------------------------------

} // namespace util


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
//------------------------------------------------------------------------------
template <typename Handle>
auto util::mapDataProductPointers(art::Event const& event, Handle const& handle)
  -> DataProductPointerMap_t<ArtHandleData_t<Handle>>
{
  using Data_t = ArtHandleData_t<Handle>;
  using Map_t = DataProductPointerMap_t<Data_t>;

  static_assert(
    std::is_same_v<std::vector<Data_t>, typename Handle::element_type>,
    "mapDataProductPointers() requires handles to STL vectors of data"
    );

  Map_t map;
  art::PtrMaker<Data_t> makePtr { event, handle.id() };
  for (auto const& [ i, item ]: util::enumerate(*handle))
    map[&item] = makePtr(i);
  return map;
} // util::mapDataProductPointers()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_UTILITIES_DATAPRODUCTPOINTERMAP_H
