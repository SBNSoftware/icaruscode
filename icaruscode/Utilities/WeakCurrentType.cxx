/**
 * @file   icaruscode/Utilities/WeakCurrentType.cxx
 * @brief  A C++ type to describe the type of weak current (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 8, 2019
 * @see    `icaruscode/Utilities/WeakCurrentType.h`
 */


// library header
#include "icaruscode/Utilities/WeakCurrentType.h"

// C/C++ standard library
#include <algorithm> // std::transform()
#include <iterator> // std::back_inserter()
#include <stdexcept> // std::logic_error
#include <cctype>


//------------------------------------------------------------------------------
namespace {
  
  // ---------------------------------------------------------------------------
  bool anyOf(std::string const& value, std::initializer_list<std::string> keys)
    { return std::find(keys.begin(), keys.end(), value) != keys.end(); }
  
  
  // ---------------------------------------------------------------------------
  
} // local namespace


//------------------------------------------------------------------------------
std::string icarus::WeakCurrentType::name() const {
  
  switch (fType) {
    case CC:  return "charged";
    case NC:  return "neutral";
    case any: return "any";
  } // switch
  throw std::logic_error("icarus::WeakCurrentType::name()");
} // icarus::WeakCurrentType::name()


//------------------------------------------------------------------------------
std::string icarus::WeakCurrentType::shortName() const {
  
  switch (fType) {
    case CC:  return "CC";
    case NC:  return "NC";
    case any: return "any";
  } // switch
  throw std::logic_error("icarus::WeakCurrentType::shortName()");
  
} // icarus::WeakCurrentType::shortName()


//------------------------------------------------------------------------------
std::string icarus::WeakCurrentType::to_upper(std::string const& s) {
  
  std::string S;
  S.reserve(s.size());
  
  auto char_toupper = [](unsigned char c){ return std::toupper(c); };
  std::transform(s.begin(), s.end(), std::back_inserter(S), char_toupper);
  
  return S;
} // icarus::WeakCurrentType::to_upper()


//------------------------------------------------------------------------------
auto icarus::WeakCurrentType::parse(std::string const& spec) -> CurrentType {
  
  std::string const SPEC = to_upper(spec);
  
  if (anyOf(SPEC, { "CHARGED", "CC" })) return CC;
  if (anyOf(SPEC, { "NEUTRAL", "NC" })) return NC;
  if (anyOf(SPEC, { "", "ANY" }))       return any;
  
  throw cet::exception("WeakCurrentType")
    << "Invalid weak current specification: '" << spec << "'\n";
  
} // icarus::WeakCurrentType::parse()


//------------------------------------------------------------------------------
