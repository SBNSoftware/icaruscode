/**
 * @file   icaruscode/PMT/Trigger/Data/ReadoutTriggerGate.h
 * @brief  A trigger gate data object associated to one or more channels.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 16, 2019
 * 
 * This is a header-only library.
 */

#ifndef ICARUSCODE_PMT_TRIGGER_DATA_READOUTTRIGGERGATE_H
#define ICARUSCODE_PMT_TRIGGER_DATA_READOUTTRIGGERGATE_H


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Data/TriggerGateData.h"

// LArSoft libraries

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <vector>
#include <initializer_list>
#include <string>
#include <stdexcept> // std::runtime_error
#include <utility> // std::move()
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------

//
// declarations
//
namespace icarus::trigger {
  
  /// All `ReadoutTriggerGate` template instances derive from this tag.
  struct ReadoutTriggerGateTag {};
  
  template <typename Tick, typename TickInterval, typename ChannelIDType>
  class ReadoutTriggerGate;
  
  template <typename Tick, typename TickInterval, typename ChannelIDType>
  std::ostream& operator<< (
    std::ostream&, ReadoutTriggerGate<Tick, TickInterval, ChannelIDType> const&
    );
  
  
  // --- BEGIN Exceptions ------------------------------------------------------
  /// @name Exceptions associated to `icarus::trigger::ReadoutTriggerGate` class
  /// @{
  
  /// Base class for all exceptions from `icarus::trigger::ReadoutTriggerGate`.
  struct ReadoutTriggerGateError: std::runtime_error
    { using std::runtime_error::runtime_error; };
  
  /// No channel associated to the readout gate.
  struct NoChannelError: ReadoutTriggerGateError
    { using ReadoutTriggerGateError::ReadoutTriggerGateError; };
  
  /// More than one channel associated to the readout gate.
  struct MoreThanOneChannelError: ReadoutTriggerGateError {
    MoreThanOneChannelError(std::string msg, std::size_t nChannels)
      : ReadoutTriggerGateError(std::move(msg)), nChannels(nChannels) {}
    using ReadoutTriggerGateError::ReadoutTriggerGateError;
    
    std::size_t nChannels; ///< Number of channels.
  }; // struct MoreThanOneChannelError
  
  // --- END Exceptions --------------------------------------------------------
  
  /// Type traits: `Gate` type derives from a `ReadoutTriggerGate` instance.
  template <typename Gate>
  struct isReadoutTriggerGate;
  
  /// Flag: `true` if `Gate` type derives from a `ReadoutTriggerGate` instance.
  template <typename Gate>
  constexpr bool isReadoutTriggerGate_v = isReadoutTriggerGate<Gate>::value;
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Logical multi-level gate associated to one or more readout channels.
 * @tparam Tick type used to count the ticks
 * @tparam TickInterval type used to quantify tick difference
 * @tparam ChannelIDType type of channel ID
 * 
 * This object is a trigger gate associated with one or more readout channels.
 * The channels are expressed as channel identifiers.
 * 
 * @note This object should be parametrized with optical ticks
 *       (`detinfo::timescales::optical_tick`). But currently the quantities
 *       (`util::quantity` derived objects) are not well suited to be serialized
 *       by ROOT, because
 *       1. they are defined in `lardataalg` rather than `lardataobj`
 *       2. writing all their serialization is daunting (see how ROOT dealt with
 *           GenVector vectors for an example of how to do it)
 *       So we chicken out here and use simple data types instead.
 */
template <typename Tick, typename TickInterval, typename ChannelIDType>
class icarus::trigger::ReadoutTriggerGate
  : public icarus::trigger::TriggerGateData<Tick, TickInterval>
  , public ReadoutTriggerGateTag
{
  /// Type of the base class.
  using Base_t = icarus::trigger::TriggerGateData<Tick, TickInterval>;
  
  /// Type of this class.
  using This_t
    = icarus::trigger::ReadoutTriggerGate<Tick, TickInterval, ChannelIDType>;
  
    public:
  
  using ClockTick_t = Tick; ///< Tick point.
  using ClockTicks_t = TickInterval; ///< Tick interval.
  using ChannelID_t = ChannelIDType; ///< Type of stored channel ID.
  
  /// Type for gate data access.
  using GateData_t = Base_t;
  
  /// Type of list of associated channels.
  using ChannelList_t = std::vector<ChannelID_t>;
  
  
  /// Constructor: a closed gate with no associated channels
  /// @see `ReadoutTriggerGate(std::initializer_list<ChannelID_t>)`,
  ///      `addChannel()`, `addChannels()`
  ReadoutTriggerGate() = default;
  
  ReadoutTriggerGate(ReadoutTriggerGate const&) = default;
  ReadoutTriggerGate(ReadoutTriggerGate&&) = default;
  ReadoutTriggerGate& operator= (ReadoutTriggerGate const&) = default;
  ReadoutTriggerGate& operator= (ReadoutTriggerGate&&) = default;
  
  /// Constructor: a closed gate associated to the specified `channels`.
  ReadoutTriggerGate(std::initializer_list<ChannelID_t> channels);
  
  
  //@{
  /// Copies/steals all the levels from the specified data.
  ReadoutTriggerGate& operator= (GateData_t const& data)
    { GateData_t::operator=(data); return *this; }
  ReadoutTriggerGate& operator= (GateData_t&& data)
    { GateData_t::operator=(std::move(data)); return *this; }
  //@}


  // --- BEGIN Gate query ------------------------------------------------------
  /// @name Gate query
  /// @{

  /// Access to the underlying gate level data.
  GateData_t const& gateLevels() const { return *this; }

  /// Access to the underlying gate level data (mutable).
  GateData_t& gateLevels() { return *this; }

  /// @}
  // --- END Gate query --------------------------------------------------------



  // --- BEGIN Access to channel information -----------------------------------
  /**
   * @brief Returns whether there is any channel id associated to the gate data.
   * @return whether channel ids are associated to the gate data
   * @see `hasChannel()`
   */
  bool hasChannels() const { return !fChannels.empty(); }
  
  /**
   * @brief Returns whether exactly one channel id associated to the gate data.
   * @return whether exactly one channel id associated to the gate data
   * @see `channel()`
   * 
   * If this methods returns `true`, `channel()` can safely be used.
   */
  bool hasChannel() const { return nChannels() == 1U; }
  
  /// Returns the number of associated channels.
  std::size_t nChannels() const { return fChannels.size(); }
  
  /**
   * @brief Returns the channels associated to the gate data.
   * @return an iterable object with all channel IDs associated to the gate data
   * @see `hasChannels()`, `channel()`
   */
  decltype(auto) channels() const;
  
  /**
   * @brief Returns the channel associated to the gate data.
   * @return the channel associated to the gate data
   * @throw MoreThanOneChannelError if more than one associated channel
   * @throw NoChannelError if no channel is associated
   * @see `hasChannel()`
   */
  ChannelID_t channel() const;
  
  
  /// Associates the specified channel to this readout gate.
  This_t& addChannel(ChannelID_t const channel);
  
  /// Associates the specified channels to this readout gate.
  This_t& addChannels(std::initializer_list<ChannelID_t> channels)
    { associateChannels(channels); return *this; }
  
  // --- END Access to channel information -------------------------------------
  
  
  
  // --- BEGIN Combination operations ------------------------------------------
  /// @name Combination operations
  /// @{
  
  /**
   * @brief Combines with a gate, keeping the minimum opening among the two.
   * @param other gate to combine to
   * @return this object
   * @see `Max()`
   * 
   * Multi-level equivalent of an _and_ logical operation.
   */
  ReadoutTriggerGate& Min(ReadoutTriggerGate const& other);

  /**
   * @brief Combines with a gate, keeping the maximum opening among the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Sum()`
   * 
   * Multi-level equivalent of an _or_ logical operation.
   */
  ReadoutTriggerGate& Max(ReadoutTriggerGate const& other);

  /**
   * @brief Combines with a gate, keeping the sum of openings of the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Max()`
   */
  ReadoutTriggerGate& Sum(ReadoutTriggerGate const& other);

  /**
   * @brief Combines with a gate, keeping the product of openings of the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Max()`, `Sum()`
   */
  ReadoutTriggerGate& Mul(ReadoutTriggerGate const& other);

  /**
   * @brief Returns a gate with the minimum opening between the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the minimum opening among `a` and `b`
   * @see `Max()`
   * 
   * Multi-level equivalent of an _and_ logical operation.
   */
  static ReadoutTriggerGate Min
    (ReadoutTriggerGate const& a, ReadoutTriggerGate const& b);

  /**
   * @brief Returns a gate with the maximum opening between the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the maximum opening among `a` and `b`
   * @see `Min()`, `Sum()`
   * 
   * Multi-level equivalent of an _or_ logical operation.
   */
  static ReadoutTriggerGate Max
    (ReadoutTriggerGate const& a, ReadoutTriggerGate const& b);

  /**
   * @brief Returns a gate with opening sum of the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the total opening of `a` and `b`
   * @see `Max()`
   */
  static ReadoutTriggerGate Sum
    (ReadoutTriggerGate const& a, ReadoutTriggerGate const& b);

  /**
   * @brief Returns a gate with opening product of the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the product of openings of `a` and `b`
   * @see `Max()`
   */
  static ReadoutTriggerGate Mul
    (ReadoutTriggerGate const& a, ReadoutTriggerGate const& b);


  /**
   * @brief Returns a gate combination of the openings of two other gates.
   * @tparam Op binary operation: `OpeningCount_t` (x2) to `OpeningCount_t`
   * @param op symmetric binary combination operation
   * @param a first gate
   * @param b second gate
   * @return gate with opening combination of `a` and `b`
   * 
   * For this algorithm to work, the operation needs to be symmetric, i.e.
   * `op(c1, c2) == op(c2, c1)` for every valid combinations of counts
   * `c1` and `c2`.
   * 
   */
  template <typename Op>
  static ReadoutTriggerGate SymmetricCombination(
    Op&& op, ReadoutTriggerGate const& a, ReadoutTriggerGate const& b,
    ClockTicks_t aDelay = ClockTicks_t{ 0 },
    ClockTicks_t bDelay = ClockTicks_t{ 0 }
    );
  
  /// @}
  // --- END Combination operations --------------------------------------------
  
  
    protected:
  // we allow some manipulation by the derived classes
  using GateEvolution_t = typename Base_t::GateEvolution_t;
  
  /// Protected constructor: set the data directly.
  ReadoutTriggerGate(GateEvolution_t&& gateLevel, ChannelList_t&& channels)
    : Base_t(std::move(gateLevel)), fChannels(std::move(channels))
    { normalizeChannels(); }
  
  /// Protected constructor: set the data directly.
  template <typename BIter, typename EIter>
  ReadoutTriggerGate(GateEvolution_t&& gateLevel, BIter b, EIter e)
    : Base_t(std::move(gateLevel)), fChannels(b, e)
    { normalizeChannels(); }
  
  
  /// Removes duplicate channel IDs and sorts the remaining ones.
  ChannelList_t& normalizeChannels() { return normalizeChannels(fChannels); }
  
  //@{
  /// Associates this data with the channels from the specified list.
  void associateChannels
    (std::initializer_list<ChannelID_t> const& moreChannels);
  void associateChannels(ChannelList_t const& moreChannels);
  //@}
  
  /// Associates this data with the channels from the `other` gate.
  void associateChannelsFromGate(ReadoutTriggerGate const& other)
    { associateChannels(other.channels()); }
  
  template <typename BIter, typename EIter>
  static ChannelList_t& mergeSortedChannelsInto
    (ChannelList_t& channels, BIter b, EIter e);


  /// Removes duplicate channel IDs.
  static ChannelList_t& normalizeSortedChannels(ChannelList_t& channels);

  //@{
  /// Removes duplicate channel IDs and sorts the remaining ones.
  static ChannelList_t& normalizeChannels(ChannelList_t& channels);
  static ChannelList_t normalizeChannels(ChannelList_t&& channels);
  //@}
  
  /// Adds channels from iterators `b` to `e` into `channels` (returned).
  template <typename BIter, typename EIter>
  static ChannelList_t& mergeChannelsInto
    (ChannelList_t& channels, BIter b, EIter e);
  
  
  /// Returns a merged list of channels from `a` and `b`.
  static ChannelList_t mergeChannels
    (ChannelList_t const& a, ChannelList_t const& b);
  
    private:
  
  /// List of readout channels associated to this data.
  ChannelList_t fChannels;
  
  
}; // class icarus::trigger::ReadoutTriggerGate


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------

#include "icaruscode/PMT/Trigger/Data/ReadoutTriggerGate.tcc"

//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_DATA_READOUTTRIGGERGATE_H
