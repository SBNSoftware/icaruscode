/**
 * @file   icaruscode/PMT/Trigger/Utilities/TriggerGateDataIteration.h
 * @brief  Utilities to iterate through trigger gata data objects.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 19, 2025
 * 
 * This library is header-only.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEDATAITERATION_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEDATAITERATION_H

// SBN/ICARUS libraries
#include "sbnobj/ICARUS/PMT/Trigger/Data/TriggerGateData.h" // for convenience

// C/C++ libraries
#include <cstdint> // std::ptrdiff_t
#include <iterator> // std::input_iterator_tag


//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  template <typename TriggerGateDataType> class triggergatedata_iterator;
  
  template <typename TriggerGateDataType>
  auto intervalsOver(
    TriggerGateDataType const& gateData,
    typename TriggerGateDataType::OpeningCount_t openingLevel
    );
  
} // namespace icarus::trigger


// -----------------------------------------------------------------------------
/**
 * @brief Iterator to a `icarus::trigger::TriggerGateData`.
 * @tparam TriggerGateDataType type of `icarus::trigger::TriggerGateData`
 *                             to iterate through
 * 
 * This is an active iterator class which applies a discrimination algorithm
 * on demand to find the next interval open to a specified level.
 * The level is specified on construction.
 * 
 * This is an _input iterator_: dereferenciation returns a temporary object
 * with start and stop ticks. It can be dereferenced multiple times (bonus
 * compared to the standard input iterator).
 * 
 * End iterators are of this same type, and they all compare and look the same
 * (even if they are the result of the end of iteration of different trigger
 * gate data objects).
 * The end iterator is `icarus::trigger::triggergatedata_iterator::end`.
 * 
 * @note This is a "constant" iterator which will not allow changes to the
 *       trigger gate data being iterated through.
 */
template <typename TriggerGateDataType>
class icarus::trigger::triggergatedata_iterator {
  
    public:
  
  /// Type of trigger gate data (`icarus::trigger::TriggerGateData`) iterated.
  using TriggerGateData_t = TriggerGateDataType;
  
  /// Type of this iterator.
  using iterator_t = triggergatedata_iterator<TriggerGateData_t>;
  
  /// Type of the tick point.
  using Tick_t = typename TriggerGateData_t::ClockTick_t;
  
  /// Type of gate opening level (e.g. `unsigned int`).
  using OpeningCount_t = typename TriggerGateData_t::OpeningCount_t;
  
  
  /// @name Iterator trait types.
  /// @{
  
  using difference_type = std::ptrdiff_t;
  
  /// Type returned by this iterator.
  struct value_type {
    Tick_t start = TriggerGateData_t::MinTick;
    Tick_t stop = TriggerGateData_t::MinTick;
    
    constexpr bool empty() const { return stop == start; }
    
    constexpr operator bool() const { return !empty(); }
    constexpr bool operator! () const { return empty(); }
    constexpr auto duration() const { return stop - start; }
  }; // value_type
  
  using reference = value_type; // not a C++ reference type
  using pointer = value_type const*;
  using iterator_category = std::input_iterator_tag;
  
  /// @}
  
  /// Nominal value of an end iterator.
  static constexpr iterator_t end{};
  
  
  /// Default constructor: invalid iterator.
  triggergatedata_iterator() = default;
  
  /// Constructor: seeks the first interval above `openingLevel`.
  triggergatedata_iterator
    (TriggerGateData_t const& gate, OpeningCount_t openingLevel);
  
  /// Returns the value the iterator is on.
  constexpr value_type operator*() const { return fNextValue; }
  
  /// Returns a pointer to the value for dereferenciation.
  constexpr pointer operator->() const { return &fNextValue; }
  
  /// Moves to the next interval, or to the end.
  /// @return the value of the iterator after the increment.
  iterator_t& operator++();
  
  /// Moves to the next interval, or to the end (postfix).
  /// @return the value of the iterator before the increment.
  iterator_t operator++(int);
  
  
  /// Returns whether `other` will produce the same sequence as this.
  constexpr bool operator== (iterator_t const& other) const;
  
  /// Returns whether `other` will not produce the same sequence as this.
  constexpr bool operator!= (iterator_t const& other) const;
  
    private:
  
  TriggerGateData_t const* fGate = nullptr; ///< Trigger gate being iterated.
  OpeningCount_t fOpeningLevel = 0; ///< The opening level being sought.
  
  value_type fNextValue; ///< The next value to be reported.
  
  /// Moves to the next value.
  void next();
  
  /// Sets this iterator to be a "end-iterator".
  void makeEnd();
  
}; // icarus::trigger::triggergatedata_iterator


//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
//---  details implementation
//------------------------------------------------------------------------------
namespace icarus::trigger::details {
  
  /// Class returning `begin()`/`end()` iterators on demand (range-for loop).
  template <typename BIter, typename EIter = BIter>
  struct IteratorBox {
    BIter b;
    EIter e;
    constexpr BIter begin() const { return b; }
    constexpr EIter end() const { return e; }
  }; // IteratorBox
  
} // namespace icarus::trigger::details


//------------------------------------------------------------------------------
//---  icarus::trigger::triggergatedata_iterator
//------------------------------------------------------------------------------
template <typename TriggerGateDataType>
icarus::trigger::triggergatedata_iterator<TriggerGateDataType>::triggergatedata_iterator
  (TriggerGateData_t const& gate, OpeningCount_t openingLevel)
  : fGate{ &gate }, fOpeningLevel{ openingLevel }
{
  next();
}


//------------------------------------------------------------------------------
template <typename TriggerGateDataType>
auto icarus::trigger::triggergatedata_iterator<TriggerGateDataType>::operator++()
  -> iterator_t&
{
  next();
  return *this;
}


//------------------------------------------------------------------------------
template <typename TriggerGateDataType>
auto icarus::trigger::triggergatedata_iterator<TriggerGateDataType>::operator++(int)
  -> iterator_t
{
  iterator_t old = *this;
  next();
  return old;
}


//------------------------------------------------------------------------------
template <typename TriggerGateDataType>
constexpr bool icarus::trigger::triggergatedata_iterator<TriggerGateDataType>::operator==
  (iterator_t const& other) const
{
  return (other.fGate == fGate) && (other.fOpeningLevel == fOpeningLevel)
    && (fNextValue.start == other.fNextValue.start);
}


//------------------------------------------------------------------------------
template <typename TriggerGateDataType>
constexpr bool icarus::trigger::triggergatedata_iterator<TriggerGateDataType>::operator!=
  (iterator_t const& other) const
{
  return (other.fGate != fGate) || (other.fOpeningLevel != fOpeningLevel)
    || (fNextValue.start != other.fNextValue.start);
}


//------------------------------------------------------------------------------
template <typename TriggerGateDataType>
void icarus::trigger::triggergatedata_iterator<TriggerGateDataType>::next() {
  value_type nextValue; // offer strong guarantee against exceptions
  if (fGate) {
    nextValue.start = fGate->findOpen(fOpeningLevel, fNextValue.stop);
    nextValue.stop = fGate->findClose(fOpeningLevel, nextValue.start);
  }
  fNextValue = nextValue;
  if (fNextValue.start == TriggerGateData_t::MaxTick) makeEnd();
}


//------------------------------------------------------------------------------
template <typename TriggerGateDataType>
void icarus::trigger::triggergatedata_iterator<TriggerGateDataType>::makeEnd() {
  *this = end;
}


//------------------------------------------------------------------------------
/**
 * @brief Function for iterating through `icarus::trigger::TriggerGateData`.
 * @tparam TriggerGateDataType type of trigger gate data to iterate through
 * @param gate the trigger gate data to iterate through
 * @param openingLevel opening level for the determination of intervals
 * @return an object suitable for range-for loops over
 * @see `icarus::trigger::triggergatedata_iterator`
 * 
 * This function allows the iteration on all intervals above the specified
 * `openingLevel`. For example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * void printInterval7(icarus::trigger::TriggerGateData<int, int> const& gate) {
 *   for (auto const& interval: intervalsOver(gate, 7))
 *     std::cout << interval.start << " - " << interval.stop << "; ";
 *   std::cout << std::endl;
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will print, one by one, all intervals with opening `7`.
 * The same value can be also immediately unpacked:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * void printInterval7(icarus::trigger::TriggerGateData<int, int> const& gate) {
 *   for (auto [ start, stop ]: intervalsOver(gate, 7))
 *     std::cout << start << " - " << stop << "; ";
 *   std::cout << std::endl;
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * The types of `start` and `stop` are tick values (`int` in the example above).
 */
template <typename TriggerGateDataType>
auto icarus::trigger::intervalsOver(
  TriggerGateDataType const& gateData,
  typename TriggerGateDataType::OpeningCount_t openingLevel
) {
  using iterator_t = triggergatedata_iterator<TriggerGateDataType>;
  return details::IteratorBox<iterator_t>
    { { gateData, openingLevel }, iterator_t::end };
} // icarus::trigger::intervalsOver()



//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_TRIGGERGATEDATAITERATION_H
