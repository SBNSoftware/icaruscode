/**
 * @file   icaruscode/PMT/Algorithms/PMTReadoutWindowMaker.h
 * @brief  Algorithm turning readout seeds into a list of readout windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 9, 2022
 * 
 * This library is header-only.
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_PMTREADOUTWINDOWMAKER_H
#define ICARUSCODE_PMT_ALGORITHMS_PMTREADOUTWINDOWMAKER_H

// ICARUS libraries
#include "icarusalg/Utilities/mfLoggingClass.h"

// C/C++ standard library
#include <ostream>
#include <vector>
#include <utility> // std::pair, std::move()
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  template <typename Time, typename Interval = typename Time::interval_t>
  class PMTReadoutWindowMaker;
}
/**
 * @brief Algorithm turning readout seeds into a list of readout windows.
 * @tparam Time type of the time point of the primitives
 * @tparam Interval type of the interval between time points
 * 
 * Given a definition of readout window in the configuration, this algorithm
 * builds a list of readout windows around "seeds".
 * A seed is a time when the trigger of the readout window happens, and each
 * readout window will include a pre-trigger buffer (`preBuffer()`) before
 * that time and a post-trigger buffer (`postBuffer()`) after it.
 * 
 * This class is designed to work with LArSoft timescale data types like
 * `detinfo::timescales::electronics_time`, but it can also be used with simple
 * data types (like `double`), provided that the `Interval` template parameter
 * is correctly set.
 */
template <typename Time, typename Interval /* = typename Time::interval_t */>
class icarus::opdet::PMTReadoutWindowMaker
  : public icarus::ns::util::mfLoggingClass
{
  
    public:
  
  using TimePoint_t = Time; ///< Type of time used.
  using TimeInterval_t = Interval; ///< Type of time interval used.
  
  /// A single contiguous readout window.
  struct Window_t: std::pair<TimePoint_t, TimePoint_t> {
    using Base_t = std::pair<TimePoint_t, TimePoint_t>;
    
    using Base_t::Base_t;
    
    constexpr TimePoint_t start() const noexcept { return Base_t::first; }
    constexpr TimePoint_t stop() const noexcept { return Base_t::second; }
    constexpr auto width() const noexcept { return stop() - start(); }
    
    constexpr bool empty() const noexcept { return start() >= stop(); }
    
    constexpr bool includes(TimePoint_t t) const noexcept
      { return (t >= start()) && (t < stop()); }
    
    /// Sets the start of this interval to `start`.
    void setStart(TimePoint_t start) { Base_t::first = start; }
    
    /// Sets the end of this interval to `stop`.
    void setStop(TimePoint_t stop) { Base_t::second = stop; }
    
    
    friend std::ostream& operator<< (std::ostream& out, Window_t const& w)
      {
        return
          out << w.start() << " -- " << w.stop() << " (" << w.width() << ")";
      }
    
  }; // Window_t
  
  
  /**
   * @brief Constructor: specifies all the algorithm parameters.
   * @param readoutWindow the total time of the readout window
   * @param preTriggerWindow the time of the window preceding a readout seed
   * @param logCategory name of the stream to send console messages to
   */
  PMTReadoutWindowMaker(
    TimeInterval_t readoutWindow, TimeInterval_t preTriggerWindow,
    std::string logCategory = "PMTReadoutWindowMaker"
    );
  
  // --- BEGIN -- Configuration query ------------------------------------------
  /// @name Configuration query
  /// @{
  
  /// Returns the total time of a readout buffer.
  TimeInterval_t buffer() const { return preBuffer() + postBuffer(); }
  
  /// Returns the time of the readout buffer before its seed.
  TimeInterval_t preBuffer() const { return fPreBuffer; }
  
  /// Returns the time of the readout buffer after its seed.
  TimeInterval_t postBuffer() const { return fPostBuffer; }
  
  /// @}
  // --- END ---- Configuration query ------------------------------------------
  
  
  // --- BEGIN -- Algorithm operations -----------------------------------------
  /// @name Algorithm operations
  /// @{
  
  //@{
  /**
   * @brief Returns a list of readout windows from the specified seeds.
   * @param seeds the list of seeds to be considered
   * @return a list of readout windows
   * 
   * Each seed generates a readout window with the size specified by the
   * algorithm configuration (`buffer()`) for each of the time seeds specified
   * in the list.
   *
   * Overlapping windows are merged.
   * 
   * The seeds are expected to be sorted in increasing time.
   */
  std::vector<Window_t> makeWindows
    (std::vector<TimePoint_t> const& seeds) const;
  
  std::vector<Window_t> operator()
    (std::vector<TimePoint_t> const& seeds) const { return makeWindows(seeds); }
  
  //@}
  
  /// @}
  // --- END ---- Algorithm operations -----------------------------------------
  
    private:
  
  /// Extension of readout buffer before trigger.
  TimeInterval_t const fPreBuffer;
  
  /// Extension of readout buffer after trigger.
  TimeInterval_t const fPostBuffer;
  
}; // icarus::opdet::PMTReadoutWindowMaker



// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename Time, typename Interval>
icarus::opdet::PMTReadoutWindowMaker<Time, Interval>::PMTReadoutWindowMaker(
  TimeInterval_t readoutWindow, TimeInterval_t preTriggerWindow,
  std::string logCategory /* = "PMTReadoutWindowMaker" */
)
  : icarus::ns::util::mfLoggingClass{ std::move(logCategory) }
  , fPreBuffer{ preTriggerWindow }
  , fPostBuffer{ readoutWindow - fPreBuffer }
  {}


// -----------------------------------------------------------------------------
template <typename Time, typename Interval>
auto icarus::opdet::PMTReadoutWindowMaker<Time, Interval>::makeWindows
  (std::vector<TimePoint_t> const& seeds) const
  -> std::vector<Window_t>
{
  std::vector<Window_t> windows;
  if (seeds.empty()) return windows;
  
  auto itSeed = seeds.begin();
  auto const send = seeds.end();
  
  windows.emplace_back(*itSeed - fPreBuffer, *itSeed + fPostBuffer);
  Window_t* lastWindow = &(windows.back());
  
  mfLogTrace() << "Processing " << seeds.size() << " readout seeds."
    << "\nFirst window from t=" << *itSeed << ": " << *lastWindow;
  
  while (++itSeed != send) {
    
    TimePoint_t const time = *itSeed;
    
    TimePoint_t const start = time - fPreBuffer;
    TimePoint_t const stop = time + fPostBuffer;
    
    mfLogTrace() << "Seed at: " << time << " (starts at " << start << ")";
    if (start <= lastWindow->stop()) { // merge
      assert(lastWindow->start() <= start);
      lastWindow->setStop(stop);
      
      mfLogTrace()
        << "  starts within last window: extended, now " << *lastWindow;
      
    }
    else { // new
      windows.emplace_back(start, stop);
      lastWindow = &(windows.back());
      mfLogTrace() << "  new window: " << *lastWindow;
    }
    
  } // while
  
  return windows;
  
} // icarus::opdet::PMTReadoutWindowMaker<>::makeWindows()


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_ALGORITHMS_PMTREADOUTWINDOWMAKER_H
