/**
 * @file   icaruscode/Utilities/ArtDataProductSelectors.h
 * @brief  Helpers to pass to `art::Event::getMany()` and similar.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 30, 2022
 */

#ifndef ICARUSCODE_UTILITIES_ARTDATAPRODUCTSELECTORS_H
#define ICARUSCODE_UTILITIES_ARTDATAPRODUCTSELECTORS_H

// framework libraries
#include "art/Framework/Principal/SelectorBase.h"

// C/C++ standard libraries
#include <regex>
#include <vector>
#include <optional>


// -----------------------------------------------------------------------------
namespace art { class InputTag; class BranchDescription; }

// -----------------------------------------------------------------------------
namespace util { class RegexDataProductSelector; }
/**
 * @brief Matches products by regex on process, module label and instance name.
 * 
 * A pattern is matched if all specified subpatterns (process, module, instance)
 * match. Any component of the pattern may be unspecified.
 * A product is selected if it matches any of the configured patterns.
 * 
 * Patterns are matched with `std::regex_match()`, i.e. the whole string must
 * match.
 * 
 * The patterns may be created from one or a sequence of input-tag-like objects
 * using `makePattern()` and `makePatterns()` static functions, respectively.
 */
class util::RegexDataProductSelector: public art::SelectorBase {
    public:
  
  /// A pattern on a input tag (empty matches everything).
  struct ProductRegex {
    std::optional<std::regex> process;
    std::optional<std::regex> module;
    std::optional<std::regex> instance;
  };
  
  /// Initializes the selector with all the supported patterns.
  RegexDataProductSelector(std::vector<ProductRegex> patterns);
  
  
  /**
   * @brief Parses a input tag to create a single data product pattern.
   * @param spec pattern specification
   * @return a pattern object
   * 
   * A specification is in the form of an `art::InputTag`, but each element
   * is a regex pattern rather than a standard name.
   * 
   * If the input tag is created from a string, in the usual form
   * `<processName>:<moduleLabel>:<instanceName>` (with the usual omissions),
   * each element can't contain a `:` character (no escaping is supported).
   * 
   * An empty `tag` matches everything.
   */
  static ProductRegex makePattern(art::InputTag const& spec);
  
  
  /**
   * @brief Parses a sequence of input tags to create data product patterns.
   * @tparam SpecColl type of collection of specifications
   * @param specs sequence of pattern specifications
   * @return a sequence of pattern objects, one per input specification
   * 
   * Each specification in `specs` is converted into a pattern via
   * `makePattern()`. The specification is explicitly converted into an input
   * tag.
   */
  template <typename SpecColl>
  static std::vector<ProductRegex> makePatterns(SpecColl const& specs);
  
    private:
  
  std::vector<ProductRegex> fPatterns; ///< All the selection patterns.
  
  
  /// Returns whether data product described by `brDescr` matches.
  virtual bool doMatch(art::BranchDescription const& brDescr) const override;
  
  /// Part of the message used when no data product matches.
  virtual std::string doPrint(std::string const& indent) const override;
  
  
  /// Returns whether the input `tag` matches the pattern `ptn`.
  bool matchPattern(art::InputTag const& tag, ProductRegex const& ptn) const;
  
}; // class util::RegexDataProductSelector



// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
template <typename SpecColl>
auto util::RegexDataProductSelector::makePatterns
  (SpecColl const& specs) -> std::vector<ProductRegex>
{
  using std::size;
  std::vector<ProductRegex> ptns;
  ptns.reserve(size(specs));
  for (auto const& spec: specs)
    ptns.push_back(makePattern(art::InputTag{ spec }));
  return ptns;
} // util::RegexDataProductSelector::makePatterns()


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_ARTDATAPRODUCTSELECTORS_H
