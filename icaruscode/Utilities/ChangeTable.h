/**
 * @file   icaruscode/Utilities/ChangeTable.h
 * @brief  Class and configuration to change one value into another.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 18, 2026
 */

#ifndef ICARUSCODE_UTILITIES_CHANGETABLE_H
#define ICARUSCODE_UTILITIES_CHANGETABLE_H

// framework libraries
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <cstddef> // std::size_t
#include <map>
#include <utility> // std::pair<>, std::move()
#include <vector>


namespace util::fhicl {
  template <typename FromType, typename ToType = FromType> class ChangeConfig;
  template
    <typename FromType, typename ToType = FromType, std::size_t N = -1ULL>
  class ChangeSequence;
}

/**
  * @brief Sequence of change tables for FHiCL configuration.
  * @tparam FromType type of the value to be changed
  * @tparam ToType (default: `FromType`) type of the new value
  * @tparam N (default: no limit) number of elements in the sequence
  *
  * Use this in a FHiCL configuration specification like:
  * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
  * struct Config {
  *
  *   using Name = fhicl::Name;
  *   using Comment = fhicl::Comment;
  *
  *   util::fhicl::ChangeSequence<> Changes {
  *     Name{ "Changes" },
  *     Comment{ "List of change instructions" }
  *   };
  *
  * };
  * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  *
  * This object also exposes a limited subset of `fhicl::Sequence` interface.
  *
  */
template
  <typename FromType, typename ToType /* = FromType */, std::size_t N /* =-1 */>
class util::fhicl::ChangeSequence {
    public:

  using From_t = FromType; ///< Type of the value to be changed.
  using To_t = ToType;     ///< Type of the new value.

  /// Specification of one change.
  using ChangeSpec_t = std::pair<From_t, To_t>;
  using ChangeMap_t = std::map<From_t, To_t>;

  /// Specifications of all changes.
  using ChangeSpecs_t = std::vector<ChangeSpec_t>;

  /// Constructor: assigns the name (mandatory) and comment to the sequence.
  ChangeSequence(::fhicl::Name&& name, ::fhicl::Comment&& comment = {});

  /// Constructor: assigns the name (mandatory) and comment to the sequence,
  /// and a conditional use clause.
  ChangeSequence(
    ::fhicl::Name&& name, ::fhicl::Comment&& comment,
    std::function<bool()> maybeUse
    );

  /// Returns the list of configured change specifications.
  ChangeSpecs_t operator() () const;

  /// Returns the configured changes as a map (old value) &rarr; (new value)
  ChangeMap_t map() const;

  // --- BEGIN ---  Pass-through interface  ------------------------------------
  /// @name Pass-through `fhicl::Sequence` interface
  /// @{

  //
  // sequence interface
  //
  bool empty() const noexcept { return fSeq.empty(); }

  std::size_t size() const noexcept { return fSeq.size(); }

  ChangeSpecs_t operator() (std::size_t index) const { return fSeq(index); }

  //
  // parameter interface
  //
  std::string const& key() const { return fSeq.key(); }
  std::string const& name() const { return fSeq.name(); }
  std::string const& comment() const { return fSeq.comment(); }

  bool is_conditional() const { return fSeq.is_conditional(); }
  bool should_use() const { return fSeq.should_use(); }

  //
  // voodoo
  //
  void prepare_elements_for_validation(std::size_t const n)
    { fSeq.prepare_elements_for_validation(n); }

  // void walk_elements(ParameterWalker<tt::const_flavor::require_non_const>& pw)
  //   { fSeq.walk_elements(pw); }
  void walk_elements
    (::fhicl::detail::ParameterWalker<tt::const_flavor::require_const>& pw) const
    { fSeq.walk_elements(pw); }

  /// @}
  // ---  END  ---  Pass-through interface  ------------------------------------


    protected:

  using ChangeConfig_t = ChangeConfig<From_t, To_t>;
  using ChangeTable_t = ::fhicl::TableAs<ChangeSpec_t, ChangeConfig_t>;

  /// Internal representation of the FHiCL sequence.
  ::fhicl::Sequence<ChangeTable_t, N> fSeq;

}; // util::fhicl::ChangeSequence<>


/**
  * @brief Configuration object for a single change specification.
  * @tparam FromType type of the value to be changed
  * @tparam ToType (default: `FromType`) type of the new value
  *
  * This configuration object includes two keys:
  *  * `From` (mandatory value of `FromType` type): value to be changed
  *  * `To` (mandatory value of `ToType` type): new value
  *
  * This class can be enclosed in `fhicl::Table` and similar:
  * ```cpp
  * using ChangeConfig_t = util::fhicl::ChangeConfig<int>;
  *
  * fhicl::Table<ChangeConfig> SingleChange{ fhicl::Name{ "Single" } };
  * fhicl::Sequence<fhicl::Table<ChangeConfig>> ManyChanges
  *   { fhicl::Name{ "Changes" } };
  * fhicl::TableAs<ChangeConfig_t::Pair_t, ChangeConfig> SingleAutoChange
  *  { fhicl::Name{ "Single" } };
  * ```
  */
template <typename FromType, typename ToType /* = FromType */>
struct util::fhicl::ChangeConfig {

  using From_t = FromType;
  using To_t = ToType;
  using Pair_t = std::pair<From_t, To_t>;

  ::fhicl::Atom<From_t> From
    { ::fhicl::Name{ "From" }, ::fhicl::Comment{ "value to change" } };
  ::fhicl::Atom<To_t> To
    { ::fhicl::Name{ "To"   }, ::fhicl::Comment{ "new value" }       };

}; // util::fhicl::ChangeConfig


// -----------------------------------------------------------------------------
namespace util::fhicl {

  /// Converts a `ChangeConfig` configuration into a from/to value pair.
  template <typename FromType, typename ToType>
  std::pair<FromType, ToType> convert
    (ChangeConfig<FromType, ToType> const& config)
    { return { config.From(), config.To() }; }

} // namespace util::fhicl


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
template <typename FromType, typename ToType /* = FromType */, std::size_t N>
util::fhicl::ChangeSequence<FromType, ToType, N>::ChangeSequence
  (::fhicl::Name&& name, ::fhicl::Comment&& comment)
  : fSeq{ std::move(name), std::move(comment), ChangeSpecs_t{} }
  {}

// -----------------------------------------------------------------------------
template <typename FromType, typename ToType /* = FromType */, std::size_t N>
util::fhicl::ChangeSequence<FromType, ToType, N>::ChangeSequence(
  ::fhicl::Name&& name, ::fhicl::Comment&& comment,
  std::function<bool()> maybeUse
)
  : fSeq{ std::move(name), std::move(comment), ChangeSpecs_t{}, maybeUse }
{}


// -----------------------------------------------------------------------------
template <typename FromType, typename ToType /* = FromType */, std::size_t N>
auto util::fhicl::ChangeSequence<FromType, ToType, N>::operator() () const
  -> ChangeSpecs_t
{
  return fSeq();
} // util::fhicl::ChangeSequence<>::operator() ()


// -----------------------------------------------------------------------------
template <typename FromType, typename ToType /* = FromType */, std::size_t N>
auto util::fhicl::ChangeSequence<FromType, ToType, N>::map() const
  -> ChangeMap_t
{
  ChangeMap_t map;
  for (ChangeSpec_t const& spec: fSeq()) map.insert(spec); // first wins
  return map;
}


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_CHANGETABLE_H
