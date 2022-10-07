/**
 * @file   FilterOnArtPathOutcome_module.cc
 * @brief  An _art_ filter picking its response from a previous one.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 4, 2022
 *
 */

// framework libraries
#include "art/Framework/Core/SharedFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "canvas/Persistency/Common/HLTGlobalStatus.h"
#include "canvas/Persistency/Common/HLTPathStatus.h"
#include "canvas/Persistency/Common/HLTenums.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <algorithm> // std::find()
#include <vector>
#include <string>
#include <optional>


//------------------------------------------------------------------------------
/**
 * @brief Filter module using a response an existing trigger path response.
 *
 * This filter emits the same response as for a _art_ path (not just a filter)
 * that was run previously.
 * 
 * 
 * Input
 * ------
 * 
 * Given a configured `Path`, this module reads from the event the data product
 * named `TriggerResults` (`art::TriggerResults` type) from the process
 * specified in that configuration parameter.
 * 
 *
 * Configuration parameters
 * =========================
 *
 * * `Path` (string, _mandatory_): the trigger path specification in the form
 *   `process:path name`.
 *
 */
class FilterOnArtPathOutcome: public art::SharedFilter {
  
    public:

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<std::string> Path {
      Name{ "Path" },
      Comment{ "trigger path specification: \"<process name:path name>\"" }
      };
    
    fhicl::Atom<bool> ResponseWhenNotRun {
      Name{ "ResponseWhenNotRun" },
      Comment{ "filter response when the reference path was not run" },
      false
      };
    
    fhicl::Atom<bool> ResponseOnError {
      Name{ "ResponseOnError" },
      Comment{ "filter response when the reference path threw an exception" },
      false
      };
    
  }; // Config
  
  using Parameters = art::SharedFilter::Table<Config>;
  
  
  /// Constructor. No surprise here.
  FilterOnArtPathOutcome(Parameters const& params, const art::ProcessingFrame&);
  
  
  virtual bool filter(art::Event& event, const art::ProcessingFrame&) override;
  
  
    private:
  
  struct TriggerSpec {
    art::InputTag const tag; ///< Input tag of the trigger result data product.
    std::string const path; ///< Name of the path to be read.
  }; // TriggerSpec;
  
  
  TriggerSpec const fTriggerSpec; ///< Trigger configuration.
  bool const fResponseWhenNotRun; ///< What to do if path was not run at all.
  bool const fResponseOnError; ///< What to do if path ended with error.
  
  
  /// Reads the needed trigger data product, and throws if not found.
  art::TriggerResults const& readTriggerResults
    (art::Event const& event, art::InputTag const& tag) const;

  /// Returns the status of the path with the required `name`.
  /// @throw cet::exception if not found
  art::HLTPathStatus const& findPath
    (art::TriggerResults const& results, std::string const& name) const;

  /// Returns the list of path names stored in `results`.
  /// @throw art::Exception (code: `art::errors::Unknown`) on logic errors
  std::vector<std::string> pathNames(art::TriggerResults const& results) const;

  /// Translates a path status into the response from this filter.
  bool responseFromPath(art::HLTPathStatus const& path) const;
  
  /// Parses a trigger path specification.
  static TriggerSpec parseTriggerSpec(std::string const& spec);
  
}; // class FilterOnArtPathOutcome


//------------------------------------------------------------------------------
FilterOnArtPathOutcome::FilterOnArtPathOutcome
  (Parameters const& params, art::ProcessingFrame const&)
  : art::SharedFilter{ params }
  , fTriggerSpec       { parseTriggerSpec(params().Path()) }
  , fResponseWhenNotRun{ params().ResponseWhenNotRun() }
  , fResponseOnError   { params().ResponseOnError() }
{
  
  async<art::InEvent>();
  
} // FilterOnArtPathOutcome::FilterOnArtPathOutcome()


//------------------------------------------------------------------------------
bool FilterOnArtPathOutcome::filter
  (art::Event& event, const art::ProcessingFrame&)
{
  
  art::TriggerResults const& trgResults
    = readTriggerResults(event, fTriggerSpec.tag);
  
  art::HLTPathStatus const& path = findPath(trgResults, fTriggerSpec.path);
  
  return responseFromPath(path);
  
} // FilterOnArtPathOutcome::filter()


//------------------------------------------------------------------------------
art::TriggerResults const& FilterOnArtPathOutcome::readTriggerResults
  (art::Event const& event, art::InputTag const& tag) const
{

  try {
    return event.getProduct<art::TriggerResults>(tag);
  }
  catch (art::Exception const& e) {
    if (e.categoryCode() != art::errors::ProductNotFound) throw;
    
    std::optional<std::vector<art::InputTag>> available;
    try {
      std::vector<art::Handle<art::TriggerResults>> const& triggerHandles
        = event.getMany<art::TriggerResults>();
      available.emplace();
      for (art::Handle<art::TriggerResults> const& handle: triggerHandles) {
        if (handle.isValid() && handle.provenance())
          available->push_back(handle.provenance()->inputTag());
      } // for
    }
    catch (...) {} // we tried to be nice, and failed.
    
    art::Exception msg { e.categoryCode(), "", e };
    msg << "Trigger data product '" << tag.encode() << "' not found.";
    if (available) {
      msg << "\n" << available->size() << " trigger products available";
      if (available->empty()) msg << ".";
      else {
        msg << ":";
        for (art::InputTag const& tag: *available)
          msg << "\n - '" << tag.encode() << "'";
      } // if ... else
    } // if available
    throw msg << "\n";
  }
  
} // FilterOnArtPathOutcome::readTriggerResults()


//------------------------------------------------------------------------------
art::HLTPathStatus const& FilterOnArtPathOutcome::findPath
  (art::TriggerResults const& results, std::string const& name) const
{
  std::vector<std::string> const names = pathNames(results);
  
  //
  // find the proper one...
  //
  
  auto iPathName = names.end();
  do { // quick-exit block
    
    // try verbatim first:
    iPathName = std::find(names.begin(), names.end(), name);
    if (iPathName != names.end()) break;
    
    // art has the custom of prepending an index to the name ("<index>:<name>");
    // let's try to ignore it
    iPathName = std::find_if(names.begin(), names.end(), 
      [name](std::string const& s)
        {
          std::size_t const iSep = s.find(':');
          return ((iSep == std::string::npos)? s: s.substr(iSep+1)) == name;
        }
      );
  
  } while (false);
  
  if (iPathName == names.end()) {
    cet::exception e{ "FilterOnArtPathOutcome" };
    e << "The trigger path '" << fTriggerSpec.tag.process()
      << "' does not include '" << name << "'! The " << names.size()
      << " available paths are:";
    for (std::string const& pathName: names)
      e << "\n - '" << pathName << "'";
    throw e << "\n";
  }
  
  // in principle we could also extract the index from the name; we don't.
  std::size_t const pathIndex = std::distance(names.begin(), iPathName);
  
  return results.at(pathIndex);
  
} // FilterOnArtPathOutcome::findPath()


//------------------------------------------------------------------------------
std::vector<std::string> FilterOnArtPathOutcome::pathNames
  (art::TriggerResults const& results) const
{
  //
  // list of paths for this event;
  // copied from art/Framework/Core/EventSelector.cc of art 3.9.3, `dataFor()`:
  //
  fhicl::ParameterSet pset;
  if (!fhicl::ParameterSetRegistry::get(results.parameterSetID(), pset)) {
    // "This should never happen";
    // just in case, I leave the message and blame to art
    throw art::Exception(art::errors::Unknown)
      << "FilterOnArtPathOutcome::findPath cannot find the trigger names for\n"
      << "a process for which the configuration has requested that the\n"
      << "OutputModule use TriggerResults to select events from.  This should\n"
      << "be impossible, please send information to reproduce this problem to\n"
      << "the art developers at artists@fnal.gov.\n";
  }
  auto const names = pset.get<std::vector<std::string>>("trigger_paths", {});
  if (names.size() != results.size()) {
    throw art::Exception(art::errors::Unknown)
      << "FilterOnArtPathOutcome::findPath: path names vector and\n"
      << "TriggerResults are different sizes (" << names.size()
      << " vs. " << results.size() << ").  This should be impossible,\n"
      << "please send information to reproduce this problem to\n"
      << "the art developers.\n";
  }
  
  return names;
  
} // FilterOnArtPathOutcome::pathNames()


//------------------------------------------------------------------------------
bool FilterOnArtPathOutcome::responseFromPath
  (art::HLTPathStatus const& path) const
{
  
  if (!path.wasrun()) { // path was not run at all
    return fResponseWhenNotRun;
  }
  else if (path.error()) { // path terminated with an exception
    return fResponseOnError;
  }
  else return path.accept();
  
} // FilterOnArtPathOutcome::responseFromPath()


//------------------------------------------------------------------------------
auto FilterOnArtPathOutcome::parseTriggerSpec(std::string const& spec)
  -> TriggerSpec
{
  std::string processName, pathName;
  auto const iSep = spec.find(':');
  if (iSep == std::string::npos) {
    pathName = spec;
  }
  else {
    processName = spec.substr(0, iSep);
    pathName = spec.substr(iSep + 1);
  }
  
  return { art::InputTag{ "TriggerResults", "", processName }, pathName };
  
} // FilterOnArtPathOutcome::parseTriggerSpec()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(FilterOnArtPathOutcome)


//------------------------------------------------------------------------------

