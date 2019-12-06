/**
 * @file   icaruscode/PMT/Trigger/Data/TriggerGateData.h
 * @brief  A logical multilevel gate for triggering.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 1, 2019
 * @see    `icaruscode/PMT/Trigger/Data/TriggerGateData.tcc`
 * 
 * This is a header only library.
 * 
 * @bug This class is too much for a data product. Merging utilities,
 *      and maybe more things, should be moved out of it.
 * 
 */

#ifndef ICARUSCODE_PMT_TRIGGER_DATA_TRIGGERGATEDATA_H
#define ICARUSCODE_PMT_TRIGGER_DATA_TRIGGERGATEDATA_H

// C/C++ standard libraries
#include <vector>
#include <iosfwd> // std::ostream
#include <optional>
#include <limits>
#include <utility> // std::pair, std::move()
#include <type_traits> // std::make_signed_t


// --- BEGIN -- Preliminary declarations and definitions -----------------------
namespace icarus::trigger {
  
  //----------------------------------------------------------------------------
  
  namespace details {
    
    /// Type of events that can happen at a certain tick of a trigger gate.
    enum class TriggerGateEventType {
      Unknown, ///< Unknown state. This can't be good.
      Set,     ///< Set the level, overriding the previous level.
      Shift    ///< Relative shift, _adding_ an offset to the previous level.
    }; // TriggerGateEventType
    
    template
      <typename ClockTick = unsigned int, typename OpeningCount = unsigned int>
    struct TriggerGateStatus {
      
      static constexpr ClockTick MinTick
        = std::numeric_limits<ClockTick>::min();
      
      /// The event which yielded this status.
      TriggerGateEventType event { TriggerGateEventType::Unknown };
      
      /// When the status starts being valid.
      ClockTick tick { std::numeric_limits<ClockTick>::lowest() };
      
      OpeningCount opening { 0 }; ///< The total opening of the gate.
      
      // default constructor is needed by ROOT I/O
      constexpr TriggerGateStatus() = default;
      constexpr TriggerGateStatus
        (TriggerGateEventType event, ClockTick tick, OpeningCount opening)
        : event(event), tick(tick), opening(opening)
        {}
      
    }; // struct TriggerGateStatus
    
    
  } // namespace details
  
  //
  // declarations
  //
  template <typename Tick, typename TickInterval>
  class TriggerGateData;
  
  
  template <typename TK, typename TI>
  std::ostream& operator<< (std::ostream&, TriggerGateData<TK, TI> const&);
  
  template <typename TK, typename TI>
  std::ostream& operator<<
    (std::ostream&, typename TriggerGateData<TK, TI>::Status const&);
  
  
} // namespace icarus::trigger
// --- END -- Preliminary declarations and definitions -------------------------


//------------------------------------------------------------------------------
/**
 * @brief Logical multi-level gate.
 * @tparam Tick type used to count the ticks
 * @tparam TickInterval type used to quantify tick difference
 * 
 * A `TriggerGate` object tracks a logical multi-level gate, that is a gate
 * whose level can be not only `0` or `1`, but an arbitrary integral number.
 * The gate level is defined as a positive number. The general convention is
 * that a level of `0` means the gate is closed, while any other level
 * (typically `1`) means it's open. In a multi-level gate, one can choose to
 * interpret the gate as open only if it reaches a chosen level higher than `1`.
 * The time domain of the gate is measured in ticks, which in this context are
 * @ref DetectorClocksIntroClocks "optical detector clock ticks" defined on
 * the @ref DetectorClocksElectronicsTime "electronics time scale".
 * That means that tick #0 is set
 * @ref DetectorClocksElectronicsStartTime "at the origin of that time scale".
 * The gate is unlimited, i.e. it can accommodate changes on any time tick
 * from #0 on (up to the machine limitations, `MaxTick`). For practical
 * purposes, `lastTick()` is instrumental to find the tick of the last gate
 * change.
 * 
 * Internally, the gate is represented by "events". The first event, at tick
 * #0, sets the level of the gate to the lowest possible (`0`, _closed_).
 * Each event may change the level at a given tick and all the following ones
 * until `MaxTick`. An actual gate is therefore defined by two events, an
 * opening event (that is an increase in level) and a closing one (a shift of
 * the same amount in the opposite direction).
 * 
 */
template <typename Tick, typename TickInterval>
class icarus::trigger::TriggerGateData {
  
    public:
  
  // --- BEGIN -- Data type definitions ----------------------------------------
  using triggergatedata_t = TriggerGateData<Tick, TickInterval>; ///< This type.
  
  /// Type of a point in time, measured in ticks.
  using ClockTick_t = Tick;
  
  /// Type of a time interval, measured in ticks.
  using ClockTicks_t = TickInterval;
  
  /// Type of count of number of open channels.
  using OpeningCount_t = unsigned int;
  
  /// Type representing a variation of open channels.
  using OpeningDiff_t = std::make_signed_t<OpeningCount_t>;
  
  // --- END -- Data type definitions ------------------------------------------
  
  
  /// An unbearably small tick number.
  static constexpr ClockTick_t MinTick
    = std::numeric_limits<ClockTick_t>::min();
  
  /// An unbearably large tick number.
  static constexpr ClockTick_t MaxTick
    = std::numeric_limits<ClockTick_t>::max();
  
  
  /// Constructor: a closed gate for the channel in `waveform`.
  TriggerGateData(): fGateLevel({ NewGateStatus }) {}
  
  
  // --- BEGIN Query -----------------------------------------------------------
  /// @name Query
  /// @{
  
  /// Returns the number of ticks this gate covers.
  ClockTick_t lastTick() const { return fGateLevel.back().tick; }
  
  /// Returns the opening count of the gate at the specified `tick`.
  OpeningCount_t openingCount(ClockTick_t tick) const;
  
  /// Returns whether this gate never opened.
  bool alwaysClosed() const { return findOpenStatus() == fGateLevel.end(); }
  
  /// Returns whether the gate is open at all at the specified `tick`.
  bool isOpen(ClockTick_t tick) const { return openingCount(tick) > 0U; }
  
  /**
   * @brief Returns the tick at which the gate opened.
   * @param minOpening minimum count to consider the gate open (default: `1`)
   * @param start first tick to check _(by default, the first available one)_
   * @param end if getting to this tick, give up _(by default, to the end)_
   * @return the tick at which the gate opened, or `end` if never
   * 
   * If the count is already above (or equal to) `minOpening` at the specified
   * `start` tick, the current level is ignored and only the next change which
   * keeps the gate open is reported.
   */
  ClockTick_t findOpen(
    OpeningCount_t minOpening = 1U,
    ClockTick_t start = MinTick, ClockTick_t end = MaxTick
    ) const;
  
  /**
   * @brief Returns the tick at which the gate has the maximum opening.
   * @param start first tick to check _(by default, the first available one)_
   * @param end if getting to this tick, give up _(by default, to the end)_
   * @return the first tick with the maximum opening, or `end` if never
   */
  ClockTick_t findMaxOpen
    (ClockTick_t start = MinTick, ClockTick_t end = MaxTick) const;
  
  /**
   * @brief Returns the range of trigger opening values in the specified range.
   * @param start the first tick to consider
   * @param end the first tick to exclude
   * @return a pair, { lower, upper } opening values
   * 
   * This method returns the lowest opening value seen in the tick range, and
   * the lowest opening value above that one _not seen_ in the same range.
   * Therefore, a gate always closed (opening value `0`) in the full range will
   * return `{ 0, 1 }`, a gate always opened (opening value `1`) will return
   * `{ 1, 2 }`, a gate which opened or closed within the tick range will return
   * `{ 0, 2 }`. In case the tick range is empty, the returned value contains
   * the same value for both extremes.
   * 
   * This function does not describe where the two extreme counts were found.
   */
  std::pair<OpeningCount_t, OpeningCount_t> openingRange
    (ClockTick_t start, ClockTick_t end) const;
  
  // --- END Query -------------------------------------------------------------
  
  
  // --- BEGIN Gate opening and closing operations -----------------------------
  /// @name Gate opening and closing operations
  /// @{
  /// Changes the opening to match `openingCount` at the specified time.
  void setOpeningAt(ClockTick_t tick, OpeningCount_t openingCount);
  
  /// Open this gate at the specified time (increase the opening by `count`).
  void openAt(ClockTick_t tick, OpeningDiff_t count)
    { openBetween(tick, MaxTick, count); }
  /// Open this gate at the specified time.
  void openAt(ClockTick_t tick) { openAt(tick, 1); }
  
  /// Open this gate at specified `start` tick, and close it at `end` tick.
  void openBetween(ClockTick_t start, ClockTick_t end, OpeningDiff_t count = 1);
  
  /// Open this gate at the specified time, and close it `length` ticks later.
  void openFor(ClockTick_t tick, ClockTicks_t length, OpeningDiff_t count = 1)
    { openBetween(tick, tick + length, count); }

  /// Close this gate at the specified time (decrease the opening by `count`).
  void closeAt(ClockTick_t tick, OpeningDiff_t count)
    { return openAt(tick, -count); }
  /// Close this gate at the specified time.
  void closeAt(ClockTick_t tick) { closeAt(tick, 1); }
  
  /// Completely close this gate at the specified time.
  void closeAllAt(ClockTick_t tick) { setOpeningAt(tick, 0U); }
  
  
  /// @}
  // --- END Gate opening and closing operations -------------------------------
  
  
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
  triggergatedata_t& Min(triggergatedata_t const& other);

  /**
   * @brief Combines with a gate, keeping the maximum opening among the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Sum()`
   * 
   * Multi-level equivalent of an _or_ logical operation.
   */
  triggergatedata_t& Max(triggergatedata_t const& other);

  /**
   * @brief Combines with a gate, keeping the sum of openings of the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Max()`
   */
  triggergatedata_t& Sum(triggergatedata_t const& other);

  /**
   * @brief Combines with a gate, keeping the product of openings of the two.
   * @param other gate to combine to
   * @return this object
   * @see `Min()`, `Max()`, `Sum()`
   * 
   * Multi-level equivalent of an extended _and_.
   */
  triggergatedata_t& Mul(triggergatedata_t const& other);

  /**
   * @brief Returns a gate with the minimum opening between the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the minimum opening among `a` and `b`
   * @see `Max()`
   * 
   * Multi-level equivalent of an _and_ logical operation.
   */
  static triggergatedata_t Min
    (triggergatedata_t const& a, triggergatedata_t const& b);

  /**
   * @brief Returns a gate with the maximum opening between the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the maximum opening among `a` and `b`
   * @see `Min()`, `Sum()`
   * 
   * Multi-level equivalent of an _or_ logical operation.
   */
  static triggergatedata_t Max
    (triggergatedata_t const& a, triggergatedata_t const& b);

  /**
   * @brief Returns a gate with opening sum of the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the total opening of `a` and `b`
   * @see `Max()`
   */
  static triggergatedata_t Sum
    (triggergatedata_t const& a, triggergatedata_t const& b);

  /**
   * @brief Returns a gate with opening product of the specified two.
   * @param a first gate
   * @param b second gate
   * @return gate with at every tick the product of openings of `a` and `b`
   * @see `Max()`, `Sum()`
   */
  static triggergatedata_t Mul
    (triggergatedata_t const& a, triggergatedata_t const& b);


  /**
   * @brief Returns a gate combination of the openings of two other gates.
   * @tparam Op binary operation: `OpeningCount_t` (x2) to `OpeningCount_t`
   * @param op symmetric binary combination operation
   * @param a first gate
   * @param b second gate
   * @param aDelay ticks of delay to be added to the first gate
   * @param bDelay ticks of delay to be added to the second gate
   * @return gate with opening combination of `a` and `b`
   * 
   * For this algorithm to work, the operation needs to be symmetric, i.e.
   * `op(c1, c2) == op(c2, c1)` for every valid combinations of counts
   * `c1` and `c2`.
   * 
   */
  template <typename Op>
  static triggergatedata_t SymmetricCombination(
    Op&& op, triggergatedata_t const& a, triggergatedata_t const& b,
    ClockTicks_t aDelay = ClockTicks_t{},
    ClockTicks_t bDelay = ClockTicks_t{}
    );
  
  /// @}
  // --- END Combination operations --------------------------------------------
  
    protected:
  
  // --- BEGIN Gate data types -------------------------------------------------
  using EventType = details::TriggerGateEventType;
  using Status = details::TriggerGateStatus<ClockTick_t, OpeningCount_t>;
  
  /// Type to describe the time evolution of the gate.
  using GateEvolution_t = std::vector<Status>;
  
  /// Protected constructor: set the data directly.
  TriggerGateData(GateEvolution_t&& gateLevel)
    : fGateLevel(std::move(gateLevel)) {}
  
  
  friend std::ostream& operator<< <ClockTick_t, ClockTicks_t>
    (std::ostream&, triggergatedata_t const&);
  friend std::ostream& operator<< <ClockTick_t, ClockTicks_t>
    (std::ostream&, Status const&);
  
  
    private:
  
  // --- BEGIN Private gate data types -----------------------------------------
  using status_iterator = typename GateEvolution_t::iterator;
  using status_const_iterator = typename GateEvolution_t::const_iterator;
  
  /// Comparison by tick number.
  struct CompareTick {
    static constexpr ClockTick_t tickOf(ClockTick_t tick) { return tick; }
    static constexpr ClockTick_t tickOf(Status const& status)
      { return status.tick; }
    template <typename A, typename B>
    bool operator() (A&& a, B&& b) const { return tickOf(a) < tickOf(b); }
  }; // CompareTick
  
  // --- END Private gate data types -------------------------------------------
  
  
  // This can't be `constexpr` because it's of a nested class type and we can't
  // call its constructor(s) during class definition. It sucks.
  /// A new gate starts with this status: opening level set to 0.
  static Status const NewGateStatus;
  
  
  GateEvolution_t fGateLevel; ///< Evolution of the gate in time.
  
  
  /// Returns a const-iterator to the status current at `tick`, or no value.
  std::optional<status_const_iterator> findLastStatusFor
    (ClockTick_t tick) const;
  
  /// Returns an iterator to the status current at `tick`, or no value.
  std::optional<status_iterator> findLastStatusFor(ClockTick_t tick);

  /// Returns an iterator to the status current at `tick`.
  /// @throw cet::exception if tick can't be handled
  status_iterator findLastStatusForTickOrThrow(ClockTick_t tick);
  
  /// Returns a const-iterator to the status current at `tick`.
  /// @throw cet::exception if tick can't be handled
  status_const_iterator findLastStatusForTickOrThrow(ClockTick_t tick) const;
  
  /// Returns an iterator to the first open status (@see `findOpen()`).
  status_const_iterator findOpenStatus(
    OpeningCount_t minOpening = 1U,
    ClockTick_t start = MinTick, ClockTick_t end = MaxTick
    ) const;

  /// Returns an iterator to the maximum open status (@see `findMaxOpen()`).
  status_const_iterator findMaxOpenStatus
    (ClockTick_t start = MinTick, ClockTick_t end = MaxTick) const;

  /// Maintenance operation: removes unconsequential stati.
  void compact();
  
}; // class icarus::trigger::TriggerGateData<>


//------------------------------------------------------------------------------
//--- Template implementation
//------------------------------------------------------------------------------

#include "icaruscode/PMT/Trigger/Data/TriggerGateData.tcc"

//------------------------------------------------------------------------------



#endif // ICARUSCODE_PMT_TRIGGER_DATA_TRIGGERGATEDATA_H

