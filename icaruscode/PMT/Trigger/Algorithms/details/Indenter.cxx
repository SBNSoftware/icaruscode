/**
 * @file   icaruscode/PMT/Trigger/Algorithms/details/Indenter.cxx
 * @brief  Simple helper to track the indentation level.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 12, 2025
 * @see    icaruscode/PMT/Trigger/Algorithms/Indenter.h
 * 
 * Not much in this implementation file.
 * Mostly it was separated not to have to load STL stream library headers
 * (the header only forward-defines `std::ostream`).
 */

// library header
#include "icaruscode/PMT/Trigger/Algorithms/details/Indenter.h"

// C/C++ standard libraries
#include <ostream>


// -----------------------------------------------------------------------------
std::ostream& icarus::trigger::details::Indenter::insertNewLine
  (std::ostream& out)
{
  return out << nextLine();
}


// -----------------------------------------------------------------------------
auto icarus::trigger::details::Indenter::nested
  (std::string addIndent, std::string addFirst /* = "" */) const -> Indenter
{
  // member-by-member initialization
  Indenter nested; // count is reset
  nested.fIndent = fIndent + std::move(addIndent);
  nested.fFirstIndent = std::move(addFirst);
  return nested;
}


// -----------------------------------------------------------------------------
std::ostream& icarus::trigger::details::operator<<
  (std::ostream& out, Indenter& indenter)
{
  return indenter.insertNewLine(out);
}


// -----------------------------------------------------------------------------
