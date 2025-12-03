/**
 * @file   icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.cxx
 * @brief  Filtering circuit response interface.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 25, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.h
 * 
 */

#undef NDEBUG // FIXME

// library header
#include "icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.h"

// ICARUS and framework libraries
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// C/C++ standard libraries
#include <ostream>


// -----------------------------------------------------------------------------
// --- icarus::trigger::FilteringCircuit
// -----------------------------------------------------------------------------
auto icarus::trigger::FilteringCircuit::makeFrequencies
  (Data_t period, std::size_t n) -> Frequencies_t
{
  std::size_t const K = n / 2 + 1;
  Frequencies_t s(K);
  
  auto const baseOmega = static_cast<Data_t>(2 * util::pi() / (n * period));
  
  for (std::size_t k = 0; k < K; ++k) s[k].imag(k * baseOmega);
  
  return s;
} // icarus::trigger::FilteringCircuit::makeFrequencies()


// -----------------------------------------------------------------------------
void icarus::trigger::FilteringCircuit::doDumpConfig
  (std::ostream& out, details::Indenter& nextLine) const
{
  // no not end the last line.
  out << nextLine << "FilteringCircuit: no configuration needed.";
}


// -----------------------------------------------------------------------------
