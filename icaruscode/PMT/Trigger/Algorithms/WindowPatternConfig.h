/**
 * @file   icaruscode/PMT/Trigger/Algorithms/WindowPatternConfig.h
 * @brief  FHiCL configuration structure for `icarus::trigger::WindowPattern`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 21, 2021
 * @see    icaruscode/PMT/Trigger/Algorithms/WindowPattern.cxx
 */

#ifndef ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWPATTERNCONFIG_H
#define ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWPATTERNCONFIG_H


// C/C++ standard libraries
#include <vector>
#include <string>
#include <limits>


// -----------------------------------------------------------------------------
namespace icarus::trigger::ns::fhicl {
  
  struct WindowPatternConfig;
  
  /**
   * @brief Function for conversion `WindowPatternConfig` -> `WindowPattern`.
   *
   * The conversion is implicitly applied by `fhicl::TableAs<>`.
   */
  icarus::trigger::WindowPattern convert(WindowPatternConfig const& config);
  
  //----------------------------------------------------------------------------
  /// Configuration element for a trigger window.
  using WindowPatternTable = 
    ::fhicl::TableAs<icarus::trigger::WindowPattern, WindowPatternConfig>;
  
  /// Configuration element for any number of trigger windows.
  using WindowPatternSequence = ::fhicl::Sequence<
    ::fhicl::TableAs<icarus::trigger::WindowPattern, WindowPatternConfig>
    >;

} // namespace icarus::trigger


//------------------------------------------------------------------------------
/**
 * @brief Configuration for a trigger window (`icarus::trigger::WindowPattern`).
 * @see `icarus::trigger::WindowPattern` and `convert()`
 * 
 * This FHiCL configuration class defines a trigger window, filling a
 * `icarus::trigger::WindowPattern` object (or a sequence).
 * 
 * The configuration describes the minimum requirement for the main window and
 * each of the neighboring windows, and whether their existence is required.
 * All parameters except for the requirement on the main window are optional,
 * and if omitted (and assigned the default values) will be always satisfied.
 * 
 * Example of usage:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * #include "icaruscode/PMT/Trigger/Algorithms/WindowPatternConfig.h"
 * #include "icaruscode/PMT/Trigger/Algorithms/WindowPattern.h"
 * 
 * class MyModule: public art::EDAnalyzer {
 *   
 *   icarus::trigger::WindowPattern const fOnePattern;
 *   icarus::trigger::WindowPatterns_t const fManyPatterns;
 *   
 *     public:
 *   struct Config {
 *     
 *     icarus::trigger::ns::fhicl::WindowPatternTable OneWindow {
 *       fhicl::Name("OneWindow"),
 *       fhicl::Comment("a single window configuration")
 *       };
 *     
 *     icarus::trigger::ns::fhicl::WindowPatternSequence ManyWindows {
 *       fhicl::Name("ManyWindows"),
 *       fhicl::Comment("many trigger window configurations [default: none]"),
 *       icarus::trigger::WindowPatterns_t{}
 *       };
 *     
 *   }; // struct Config
 *   using Parameters = art::EDAnalyzer::Table<Config>;
 *   
 *   MyModule(Parameters const& config);
 *   
 *   // ...
 * }; // class MyModule
 * 
 * 
 * MyModule::MyModule(Parameters const& config)
 *   : art::EDAnalyzer(config)
 *   , fOnePattern(config().OneWindow())
 *   , fManyPatterns(config().ManyWindows())
 * {
 *   // ...
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will support a configuration block like:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * physics.analyzers.mymodule: {
 *   module_type: MyModule
 *   OneWindow: {
 *     inMainWindow:     1
 *     inOppositeWindow: 1
 *   }
 *   ManyWindows: [
 *     {
 *       # M3D1
 *       inMainWindow:       3
 *       inDownstreamWindow: 1
 *     },
 *     {
 *       # M3D1req
 *       inMainWindow:               3
 *       inDownstreamWindow:         1
 *       requireDownstreamWindow: true
 *     },
 *     {
 *       # M5D2
 *       inMainWindow:       5
 *       inDownstreamWindow: 2
 *     },
 *     {
 *       # M5D2req
 *       inMainWindow:               5
 *       inDownstreamWindow:         2
 *       requireDownstreamWindow: true
 *     }
 *   ]
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
struct icarus::trigger::ns::fhicl::WindowPatternConfig {
  
  using Name = ::fhicl::Name;
  using Comment = ::fhicl::Comment;
  template <typename T> using Atom = ::fhicl::Atom<T>;
  
  Atom<unsigned int> inMainWindow {
    Name("inMainWindow"),
    Comment("minimum fired primitives in the main sliding window")
    };
  Atom<unsigned int> inUpstreamWindow {
    Name("inUpstreamWindow"),
    Comment(
      "minimum fired primitives in the sliding window upstream of main one"
      ),
    0U // default
    };
  Atom<unsigned int> inDownstreamWindow {
    Name("inDownstreamWindow"),
    Comment(
      "minimum fired primitives in the sliding window downstream of main one"
      ),
    0U // default
    };
  Atom<unsigned int> inOppositeWindow {
    Name("inOppositeWindow"),
    Comment(
      "minimum fired primitives in the sliding window opposite of main one"
      ),
    0U // default
    };
  
  Atom<bool> requireUpstreamWindow {
    Name("requireUpstreamWindow"),
    Comment("an upstream window must be present (no border main window)"),
    false
    };
  Atom<bool> requireDownstreamWindow {
    Name("requireDownstreamWindow"),
    Comment("a downstream window must be present (no border main window)"),
    false
    };
  
}; // icarus::trigger::ns::fhicl::WindowPatternConfig


// -----------------------------------------------------------------------------
// --- Inline implementation
// -----------------------------------------------------------------------------
namespace icarus::trigger::ns::fhicl {
  
  icarus::trigger::WindowPattern convert(WindowPatternConfig const& config) {
    return {
      config.inMainWindow(),           // minInMainWindow
      config.inUpstreamWindow(),       // minInUpstreamWindow
      config.inDownstreamWindow(),     // minInDownstreamWindow
      config.inOppositeWindow(),       // minInOppositeWindow
      config.requireUpstreamWindow(),  // requireUpstreamWindow
      config.requireDownstreamWindow() // requireDownstreamWindow
      };
  } // icarus::trigger::ns::fhicl::convert()
  
} // namespace icarus::trigger::ns::fhicl


// -----------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_ALGORITHMS_WINDOWPATTERNCONFIG_H
