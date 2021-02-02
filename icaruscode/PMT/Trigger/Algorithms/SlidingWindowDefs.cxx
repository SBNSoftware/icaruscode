/**
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.cxx
 * @brief  Definition for PMT sliding windows.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 26, 2021
 * @file   icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.h
 * 
 */


// library header
#include "icaruscode/PMT/Trigger/Algorithms/SlidingWindowDefs.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"

// C/C++ standard libraries
#include <ostream>


// -----------------------------------------------------------------------------
void icarus::trigger::printTriggerWindowChannels
  (std::ostream& out, TriggerWindowChannels_t const& window)
{
  /*
   * Format:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * std::cout << "Window ";
   * printTriggerWindowChannels(std::cout, window);
   * std::cout << " defined." << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * may yield:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Window 6 channels ( 1, 2, 3, 4, 5, 6 ) defined.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  out << window.size() << " channels";
  if (window.empty()) return;
  auto iChannel = window.begin();
  auto const cend = window.end();
  out << " ( " << *iChannel;
  while (++iChannel != cend) out << ", " << *iChannel;
  out << " )";
  
} // icarus::trigger::printTriggerWindowChannels()


// -----------------------------------------------------------------------------
void icarus::trigger::printTriggerWindowDefs
  (std::ostream& out, TriggerWindowDefs_t const& windows)
{
  /*
   * Format:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * std::cout << "Blah blah blah >";
   * printTriggerWindowDefs(std::cout, windows);
   * std::cout << "And blah blah." << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * may yield:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Blah blah blah >3 windows:
   *  #0: 6 channels ( 1, 2, 3, 4, 5, 6 )
   *  #1: 6 channels ( 4, 5, 6, 7, 8, 9 )
   *  #2: 6 channels ( 7, 8, 9, 10, 11, 12 )
   * And blah blah.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  out << windows.size() << " windows:";
  for (auto const& [ iWindow, window ]: util::enumerate(windows)) {
    out << "\n #" << iWindow << ": ";
    printTriggerWindowChannels(out, window);
  } // for
  
} // icarus::trigger::printTriggerWindowDefs()


// -----------------------------------------------------------------------------
