/**
 * @file   icaruscode/Utilities/ArtDataProductSelectors.cxx
 * @brief  Helpers to pass to `art::Event::getMany()` and similar (implement.).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 30, 2022
 * @see    icaruscode/Utilities/ArtDataProductSelectors.h
 */

// library header
#include "icaruscode/Utilities/ArtDataProductSelectors.h"

// framework libraries
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Utilities/InputTag.h"

// C/C++ standard libraries
#include <regex>
#include <string>
#include <utility> // std::move()


// -----------------------------------------------------------------------------
util::RegexDataProductSelector::RegexDataProductSelector
  (std::vector<ProductRegex> patterns)
  : fPatterns(std::move(patterns))
{}


// -----------------------------------------------------------------------------
auto util::RegexDataProductSelector::makePattern
  (art::InputTag const& spec) -> ProductRegex
{
  
  auto const regexIfNotEmpty = [](std::string const& ptn)
    { return ptn.empty()? std::nullopt: std::optional{ std::regex{ ptn }}; };
  
  return {
      regexIfNotEmpty(spec.process())   // process
    , regexIfNotEmpty(spec.label())     // module
    , regexIfNotEmpty(spec.instance())  // instance
    };
  
} // util::RegexDataProductSelector::makePattern()


// -----------------------------------------------------------------------------
bool util::RegexDataProductSelector::doMatch
  (art::BranchDescription const& brDescr) const
{
  art::InputTag const& tag = brDescr.inputTag();
  for (ProductRegex const& pattern: fPatterns)
    if (matchPattern(tag, pattern)) return true;
  return false;
} // util::RegexDataProductSelector::doMatch()


// -----------------------------------------------------------------------------
std::string util::RegexDataProductSelector::doPrint
  (std::string const& indent) const
{
  return indent
    + "any of " + std::to_string(fPatterns.size()) + " patterns";
} // util::RegexDataProductSelector::doPrint()


// -----------------------------------------------------------------------------
bool util::RegexDataProductSelector::matchPattern
  (art::InputTag const& tag, ProductRegex const& ptn) const
{
  if (ptn.process  && !std::regex_match(tag.process(),  *(ptn.process) ))
    return false;
  if (ptn.module   && !std::regex_match(tag.label(),    *(ptn.module)  ))
    return false;
  if (ptn.instance && !std::regex_match(tag.instance(), *(ptn.instance)))
    return false;
  return true;
} // util::RegexDataProductSelector::matchPattern()


// -----------------------------------------------------------------------------

