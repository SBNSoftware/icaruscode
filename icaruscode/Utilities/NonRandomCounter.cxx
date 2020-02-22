/**
 * @file   icaruscode/Utilities/NonRandomCounter.cxx
 * @brief  Non-random number engine for profiling purposes (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 28, 2020
 * @see    `icaruscode/Utilities/NonRandomCounter.h`
 *
 */

// library header
#include "icaruscode/Utilities/NonRandomCounter.h"

// C/C++ standard library
#include <fstream>
#include <iostream> // std::cout


// -----------------------------------------------------------------------------
void util::NonRandomCounter::saveStatus
  (const char filename[] /* = "NonRandomCounter.conf" */) const
{
  std::ofstream f { filename };
  if (!f) return; // shouldn't we complain? HepJamesRandom does not
  
  f << count;
  
} // util::NonRandomCounter::saveStatus()


// -----------------------------------------------------------------------------
void util::NonRandomCounter::restoreStatus
  (const char filename[] /* = "NonRandomCounter.conf" */)
{
  std::ifstream f { filename };
  if (!f) return; // shouldn't we complain? HepJamesRandom does not
  
  f >> count;
  
} // util::NonRandomCounter::restoreStatus()


// -----------------------------------------------------------------------------
void util::NonRandomCounter::showStatus() const {
  std::cout << "Counter: " << count << std::endl;
}


// -----------------------------------------------------------------------------
