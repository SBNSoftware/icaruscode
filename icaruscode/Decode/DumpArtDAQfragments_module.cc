/**
 * @file   DumpArtDAQfragments_module.cc
 * @brief  Dumps on console the content of `artdaq::Fragment` data products.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 10, 2021
 */

// SBN libraries
#include "icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.h" // dumpFragment()

// LArSoft libraries
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "canvas/Persistency/Provenance/ProductToken.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <algorithm> // std::transform()
#include <vector>
#include <string>


//------------------------------------------------------------------------------
namespace sbn { class DumpArtDAQfragments; }

/**
 * @brief Dumps on console the content of `artdaq::Fragment` collections.
 * 
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<artdaq::Fragment>`: data to be dumped; the dump is binary,
 *   and it does not attempt an interpretation of the content of the fragment
 *   (in future, it could at least attempt unpacking the container fragments).
 * 
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description DumpArtDAQfragments`.
 * 
 * * `FragmentTags` (list of data product input tags): the tags identifying the
 *     data products to dump.
 * * `OutputCategory` (string, default: `DumpArtDAQfragments`): name of the
 *     message facility output stream to dump the information into
 * 
 */
class sbn::DumpArtDAQfragments: public art::EDAnalyzer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> FragmentTags {
      Name("FragmentTags"),
      Comment("tag of fragment collection data products to be dumped")
      };

    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("name of the category used for the output"),
      "DumpArtDAQfragments"
      };

  }; // struct Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit DumpArtDAQfragments(Parameters const& config);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Does nothing, but it is mandatory.
  virtual void analyze(art::Event const& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  /// Input data tokens.
  std::vector<art::InputTag> const fInputTags;
  
  /// Input data tokens.
  std::vector<art::ProductToken<artdaq::Fragments>> const fInputTokens;
  
  /// Category used for message facility stream.
  std::string const fOutputCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  /// Get art tokens for specified input data products.
  std::vector<art::ProductToken<artdaq::Fragments>> getFragmentTokens
    (std::vector<art::InputTag> const& tags);
  
  /// Dumps a single fragment collection.
  void dumpFragments(
    art::Event const& event,
    art::InputTag const& inputTag,
    art::ProductToken<artdaq::Fragments> const& inputToken
    ) const;
  
}; // sbn::DumpArtDAQfragments


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- sbn::DumpArtDAQfragments
//------------------------------------------------------------------------------
sbn::DumpArtDAQfragments::DumpArtDAQfragments
  (Parameters const& config)
  : art::EDAnalyzer(config)
  // configuration
  , fInputTags     (config().FragmentTags())
  , fInputTokens   (getFragmentTokens(fInputTags))
  , fOutputCategory(config().OutputCategory())
{
  assert(fInputTags.size() == fInputTokens.size());
} // sbn::DumpArtDAQfragments::DumpArtDAQfragments()


//------------------------------------------------------------------------------
void sbn::DumpArtDAQfragments::analyze(art::Event const& event) {
  
  mf::LogVerbatim(fOutputCategory) << event.id() << ":";
  
  for (auto const& [ tag, token ]: util::zip(fInputTags, fInputTokens))
    dumpFragments(event, tag, token);

} // sbn::DumpArtDAQfragments::analyze()


//------------------------------------------------------------------------------
void sbn::DumpArtDAQfragments::dumpFragments(
  art::Event const& event,
  art::InputTag const& inputTag,
  art::ProductToken<artdaq::Fragments> const& inputToken
) const {
  
  art::Handle<artdaq::Fragments> fragmentHandle;
  //bool found = event.getByToken(inputToken, fragmentHandle);
  
  if ( !(fragmentHandle = event.getHandle<artdaq::Fragments>(inputToken)) ) {
    mf::LogVerbatim(fOutputCategory)
      << "No fragment collection found as '" << inputTag.encode() << "'.";
    return;
  }
  
  artdaq::Fragments const& fragments = *fragmentHandle;
  mf::LogVerbatim log { fOutputCategory };
  log << "Data product '" << inputTag.encode() << "' has " << fragments.size()
    << " fragments.";
  for (auto const& [ iFragment, fragment ]: util::enumerate(fragments)) {
    // TODO special case for artdaq::ContainerFragment?
    log << "\n[#" << iFragment << "] " << sbndaq::dumpFragment(fragment);
  } // for
  
} // sbn::DumpArtDAQfragments::dumpFragments()


//------------------------------------------------------------------------------
std::vector<art::ProductToken<artdaq::Fragments>>
sbn::DumpArtDAQfragments::getFragmentTokens
  (std::vector<art::InputTag> const& tags)
{
  auto getToken = [this](art::InputTag const& tag)
    { return consumes<artdaq::Fragments>(tag); };
  
  std::vector<art::ProductToken<artdaq::Fragments>> tokens;
  std::transform(begin(tags), end(tags), back_inserter(tokens), getToken);
  return tokens;
} // sbn::DumpArtDAQfragments::getFragmentTokens()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::DumpArtDAQfragments)


//------------------------------------------------------------------------------
