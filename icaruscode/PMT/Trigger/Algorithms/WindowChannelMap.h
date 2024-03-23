/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WindowChannelMap.h
 * @brief  Data structure enclosing information for trigger sliding windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/WindowChannelMap.cxx
 * 
 * Linkage to the implementation file is required when using `dump()` methods.
 */


#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWCHANNELMAP_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWCHANNELMAP_H


// LArSoft libraries
#include <cstdint>  // uint16_t in OpDetWaveform.h
#include "lardataobj/RawData/OpDetWaveform.h" // raw::Channel_t
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::CryostatID

// C/C++ standard libraries
#include <iosfwd> // std::ostream
#include <vector>
#include <string>
#include <utility> // std::move()
#include <limits> // std::numeric_limits<>
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
namespace icarus::trigger {
  class WindowChannelMap;
  std::ostream& operator<< (std::ostream& , WindowChannelMap const&);
} // namespace icarus::trigger

/**
 * @brief Information about composition and topology of trigger sliding windows.
 * 
 * This class collects optical detector windows defined as a collection of
 * channels, supposedly contiguous in space.
 * A rudimentary topology is described with each window having upstream,
 * downstream windows (with respect to the nominal beam direction) and a
 * single opposite window.
 * 
 * The class exposes a minimal collection interface and is therefore iterable.
 * Windows are identified by an index that matches their position in the
 * collection.
 * 
 * This is just a data structure, and algorithms building the composition and
 * topology of the windows are not part of it.
 */
class icarus::trigger::WindowChannelMap {
  
    public:
  
  using WindowIndex_t = std::size_t; ///< Type of window index
  
  /// Special index denoting an invalid window.
  static const WindowIndex_t InvalidWindowIndex
    = std::numeric_limits<WindowIndex_t>::max();
  
  /// Geometric location and composition of the window.
  struct WindowComposition_t {
    
    geo::Point_t center; ///< Center of the window.
    
    geo::CryostatID cryoid; ///< Which cryostat the channels are in.
    
    /// Optical detector channels covered by this window.
    std::vector<raw::Channel_t> channels;
    
    /// Returns whether the window is in a single, known cryostat.
    bool hasCryostat() const { return cryoid.isValid; }
    
  }; // struct WindowComposition_t
  
  /// Information of the identity and neighbourhood of a window.
  struct WindowTopology_t {
    
    WindowIndex_t index; ///< Index of the window this information is about.
    
    /// Index of the window opposite to this one.
    WindowIndex_t opposite { InvalidWindowIndex };
    
    /// Index of the window upstream of this one.
    WindowIndex_t upstream { InvalidWindowIndex };
    
    /// Index of the window downstream of this one.
    WindowIndex_t downstream { InvalidWindowIndex };
    
    /// Returns whether the main window has another upstream of it.
    bool hasUpstreamWindow() const { return isValidWindow(upstream); }
    
    /// Returns whether the main window has another downstream of it.
    bool hasDownstreamWindow() const { return isValidWindow(downstream); }
    
    /// Returns whether the main window has another downstream of it.
    bool hasOppositeWindow() const { return isValidWindow(opposite); }
    
  }; // struct WindowTopology_t
  
  /// Information of a single window.
  struct WindowInfo_t {
    
    WindowComposition_t composition;
    
    WindowTopology_t topology;
    
    /// Prints the information content (single line).
    void dump(std::ostream& out, std::string const& indent = "") const;
    
  }; // struct WindowInfo_t
  
  
  /// Construction: moves the `windows` information into the map.
  explicit WindowChannelMap(std::vector<WindowInfo_t>&& windows)
    : fWindows(std::move(windows)) {}
  
  // --- BEGIN Access to window information ------------------------------------
  /// @name Access to window information.
  /// @{
  
  /// Number of sliding windows.
  std::size_t nWindows() const { return fWindows.size(); }
  
  /// Returns whether a window with the specified index is present.
  bool hasWindow(WindowIndex_t index) const { return index < nWindows(); }
  
  /// Returns the information for the window with specified `index` (unchecked).
  WindowInfo_t const& info(WindowIndex_t index) const { return fWindows[index]; }
  
  /// @}
  // --- END Access to window information --------------------------------------

  // @{
  /// Prints the content of the full mapping.
  void dump(
    std::ostream& out, std::string const& indent, std::string const& firstIndent
    ) const;
  void dump(std::ostream& out, std::string const& indent = "") const
    { dump(out, indent, indent); }
  // @}
  

  
  /// @name Iterable standard interface
  /// @{
  using value_type = WindowInfo_t;
  using reference_type = WindowInfo_t const&;
  using pointer_type = WindowInfo_t const*;
  
  std::size_t size() const { return nWindows(); }
  bool empty() const { return fWindows.empty(); }
  auto begin() const { return fWindows.begin(); }
  auto end() const { return fWindows.end(); }
  auto cbegin() const { return begin(); }
  auto cend() const { return end(); }
  /// @}
  
  
  /// Returns whether the specified index is not the invalid window index.
  static bool isValidWindow(WindowIndex_t index)
    { return index != InvalidWindowIndex; }
  
  
    private:
  std::vector<WindowInfo_t> fWindows; /// Information for each window.
  
}; // icarus::trigger::WindowChannelMap


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWCHANNELMAP_H
