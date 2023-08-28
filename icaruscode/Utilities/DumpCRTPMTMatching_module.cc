/**
 * @file   icaruscode/Utilities/DumpCRTPMTMatching_module.cc
 * @brief  Dumps on screen the content of the CRT/PMT matchings.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 8, 2023
 */

// LArSoft and ICARUS libraries
#include "icaruscode/Utilities/IcarusObjectSelectors.h"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcorealg/CoreUtils/enumerate.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

// C//C++ standard libraries
#include <optional>
#include <map>
#include <vector>
#include <string>


namespace { class ProductTagIndex; }

namespace icarus::crt { class DumpCRTPMTMatching; }
/**
 * @brief Prints the content of all the tracks on screen
 *
 * This analyser prints the content of all the tracks into the
 * LogInfo/LogVerbatim stream.
 * The associated objects are printed only if they were produced with the
 * same input tag as the tracks.
 *
 * Configuration parameters
 * -------------------------
 *
 * - `MatchingTag` (input tag, _required_): tag of the data product with the
 *   matchings to dump.
 * - *PrintFlashAssociations* (boolean, default: `true`): prints some
 *   information from the associated flash.
 * - *PrintCRTHitAssociations* (boolean, default: `true`): prints some
 *   information from the associated CRT hits.
 * - `LogCategory` (string, default: `"DumpCRTPMTMatching"`): the category
 *   used for the output (useful for filtering).
 *
 */
class icarus::crt::DumpCRTPMTMatching: public art::EDAnalyzer {
  
    public:
  using CRTPMTMatching = sbn::crt::CRTPMTMatching;
  using MatchedCRT = sbn::crt::MatchedCRT;
  /// Configuration object
  struct Config {
    using Comment = fhicl::Comment;
    using Name = fhicl::Name;

    fhicl::Atom<art::InputTag> MatchingTag{
      Name("MatchingTag"),
      Comment("input tag for the CRT/PMT matching information to be dumped")
    };
    
    fhicl::Atom<bool> PrintFlashAssociations{
      Name("PrintFlashAssociations"),
      Comment("prints some information from the associated flash"),
      true
    };
    
    fhicl::Atom<bool> PrintCRTHitAssociations{
      Name("PrintCRTHitAssociations"),
      Comment("prints some information from the associated CRT hits"),
      true
    };
    
    fhicl::Atom<std::string> LogCategory{
      Name("LogCategory"),
      Comment("name of the category used for message facility output"),
      "DumpCRTPMTMatching"
    };

  }; // Config

  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Default constructor
  explicit DumpCRTPMTMatching(Parameters const& config);

  /// Does the printing
  void analyze(art::Event const& evt) override;

    private:
  art::InputTag const fMatchingTag; ///< Input tag for matching input data product.
  bool const fPrintFlashAssociations; ///< Whether to print associated flashes.
  bool const fPrintCRTHitAssociations; ///< Whether to print associated CRT hits.
  std::string const fLogCategory;   ///< category for `mf::LogInfo` output.

  /// Dumps information about the specified matching.
  void DumpMatching(
    unsigned int iMatching,
    sbn::crt::CRTPMTMatching const& matchInfo,
    art::Ptr<recob::OpFlash> const& flash,
    std::vector<art::Ptr<sbn::crt::CRTHit>> const& CRThits,
    ProductTagIndex const& tagIndex,
    std::string const& indent = ""
    ) const;

}; // class icarus::crt::DumpCRTPMTMatching


//------------------------------------------------------------------------------
//---  implementation
//------------------------------------------------------------------------------
namespace {
  
  class ProductTagIndex {
    
    static art::InputTag const EmptyTag;
    
    art::Event const& fEvent;
    std::unordered_map<art::ProductID, art::InputTag> fTags;
    
    art::InputTag extractTag(art::ProductID id) const
      {
        auto const descr = fEvent.getProductDescription(id);
        return descr? descr->inputTag(): art::InputTag{};
      }
    
      public:
    ProductTagIndex(art::Event const& event): fEvent(event) {}
    
    void add(art::ProductID id) { fTags.emplace(id, extractTag(id)); }
    template <typename T>
    void add(art::Ptr<T> const& ptr) { add(ptr.id()); }
    
    art::InputTag const& tagOf(art::ProductID id) const
      {
        if (auto it = fTags.find(id); it != fTags.end()) return it->second;
        return EmptyTag;
      }
    
    template <typename T>
    art::InputTag const& tagOf(art::Ptr<T> const& ptr) const
      { return tagOf(ptr.id()); }
    
    art::InputTag const& operator() (art::ProductID id) const
      { return tagOf(id); }
    
    template <typename T>
    art::InputTag const& operator() (art::Ptr<T> const& ptr) const
      { return tagOf(ptr); }
    
  }; // class ProductTagIndex
  
  art::InputTag const ProductTagIndex::EmptyTag;
  
} // local namespace


//------------------------------------------------------------------------------
icarus::crt::DumpCRTPMTMatching::DumpCRTPMTMatching(Parameters const& params)
  : art::EDAnalyzer(params)
  , fMatchingTag(params().MatchingTag())
  , fPrintFlashAssociations(params().PrintFlashAssociations())
  , fPrintCRTHitAssociations(params().PrintCRTHitAssociations())
  , fLogCategory(params().LogCategory())
{}


//------------------------------------------------------------------------------
void icarus::crt::DumpCRTPMTMatching::analyze(art::Event const& event) {

  //
  // collect all the input
  //
  auto const& matchHandle
    = event.getValidHandle<std::vector<sbn::crt::CRTPMTMatching>>(fMatchingTag);

  mf::LogInfo(fLogCategory)
    << "The event contains " << matchHandle->size() << " '"
    << fMatchingTag.encode() << "' CRT/PMT matches.";

  auto const flashAssns = fPrintFlashAssociations
    ? std::make_optional
      (art::FindOneP<recob::OpFlash>(matchHandle, event, fMatchingTag))
    : std::nullopt
    ;
  if (flashAssns && !flashAssns->isValid()) {
    throw art::Exception(art::errors::ProductNotFound)
      << "No flash associated with '" << fMatchingTag.encode() << "' matches.\n";
  }

  auto const hitAssns = fPrintCRTHitAssociations
    ? std::make_optional
      (art::FindManyP<sbn::crt::CRTHit>(matchHandle, event, fMatchingTag))
    : std::nullopt
    ;
  if (hitAssns && !hitAssns->isValid()) {
    throw art::Exception(art::errors::ProductNotFound)
      << "No CRT hits associated with '" << fMatchingTag.encode() << "' matches.\n";
  }
  
  ProductTagIndex tagIndex{ event };

  for (auto const& [ iMatch, match ]: util::enumerate(*matchHandle)) {
    
    if (flashAssns) tagIndex.add(flashAssns->at(iMatch));
    if (hitAssns) {
      for (art::Ptr<sbn::crt::CRTHit> const& hitPtr: hitAssns->at(iMatch))
        tagIndex.add(hitPtr);
    }
    DumpMatching(iMatch, match
      , flashAssns? flashAssns->at(iMatch): art::Ptr<recob::OpFlash>{}
      , hitAssns? hitAssns->at(iMatch): std::vector<art::Ptr<sbn::crt::CRTHit>>{}
      , tagIndex
      );
    
  }   // for matches
  
} // icarus::crt::DumpCRTPMTMatching::analyze()


//---------------------------------------------------------------------------
void icarus::crt::DumpCRTPMTMatching::DumpMatching(
  unsigned int iMatching,
  sbn::crt::CRTPMTMatching const& matchInfo,
  art::Ptr<recob::OpFlash> const& flash,
  std::vector<art::Ptr<sbn::crt::CRTHit>> const& CRThits,
  ProductTagIndex const& tagIndex,
  std::string const& indent /* = "" */
) const {
  
  mf::LogVerbatim log(fLogCategory);
  
  log << "[#" << iMatching << "] match:";
  if (flash) log << " 1 flash";
  if (!CRThits.empty()) {
    if (flash) log << " and";
    log << " " << CRThits.size() << " CRT hits";
  }
  
  log << " (" << (
    icarus::crt::MatchTypeSelector.hasOption(matchInfo.flashClassification)
    ? icarus::crt::MatchTypeSelector.get(matchInfo.flashClassification).name()
    : "<unknown [#" + std::to_string(static_cast<int>(matchInfo.flashClassification)) + "]>"
    ) << ")";
  
  auto const hasCount
    = [](unsigned int count){ return count != CRTPMTMatching::NoCount; };
  
  if (hasCount(matchInfo.nTopCRTHitsBefore) || hasCount(matchInfo.nTopCRTHitsAfter))
  {
    log << " with top CRT hits";
    if (hasCount(matchInfo.nTopCRTHitsBefore))
      log << " " << matchInfo.nTopCRTHitsBefore << " before";
    if (hasCount(matchInfo.nTopCRTHitsAfter)) {
      if (hasCount(matchInfo.nTopCRTHitsBefore)) log << " and";
      log << " " << matchInfo.nTopCRTHitsAfter << " after";
    }
    log << " the flash";
  }
  if (hasCount(matchInfo.nSideCRTHitsBefore) || hasCount(matchInfo.nSideCRTHitsAfter))
  {
    if (hasCount(matchInfo.nTopCRTHitsBefore) || hasCount(matchInfo.nTopCRTHitsAfter))
      log << ",";
    log << " with side CRT hits";
    if (hasCount(matchInfo.nSideCRTHitsBefore))
      log << " " << matchInfo.nSideCRTHitsBefore << " before";
    if (hasCount(matchInfo.nSideCRTHitsAfter)) {
      if (hasCount(matchInfo.nSideCRTHitsBefore)) log << " and";
      log << " " << matchInfo.nSideCRTHitsAfter << " after";
    }
    log << " the flash";
  }
  
  // flash ID
  log << "\n" << indent << "  flash";
  if (matchInfo.flashID == CRTPMTMatching::NoID) log << " <unknown ID>";
  else                                           log << " ID=" << matchInfo.flashID;
  
  // flash times
  if (matchInfo.flashTime != CRTPMTMatching::NoTime)
    log << " at " << matchInfo.flashTime << " us";
  if (matchInfo.flashGateTime != CRTPMTMatching::NoTime)
    log << ", " << (matchInfo.flashGateTime*1000) << " ns after beam gate";
  
  // optical hit times
  if (matchInfo.firstOpHitStartTime != CRTPMTMatching::NoTime
    || matchInfo.firstOpHitPeakTime != CRTPMTMatching::NoTime)
  {
    log << " (";
    if (matchInfo.firstOpHitStartTime != CRTPMTMatching::NoTime)
      log << "first hit start at " << matchInfo.firstOpHitStartTime;
    if (matchInfo.firstOpHitPeakTime != CRTPMTMatching::NoTime) {
      if (matchInfo.firstOpHitStartTime != CRTPMTMatching::NoTime)
        log << ", ";
      log << "peak at " << matchInfo.firstOpHitPeakTime << " us";
    }
    log << ")";
  }
  
  // gates and p.e.
  log << "," << (matchInfo.flashInGate? "": " not") << " in beam gate,"
    << (matchInfo.flashInBeam? "": " not") << " in spill";
  if (matchInfo.flashPE >= 0.0)
    log << ", with " << matchInfo.flashPE << " p.e.";
  
  // position and width
  log << ", at " << matchInfo.flashPosition << " cm";
  if (matchInfo.flashYWidth != CRTPMTMatching::NoWidth
    || matchInfo.flashZWidth != CRTPMTMatching::NoWidth)
  {
    log << " (width ";
    if (matchInfo.flashZWidth != CRTPMTMatching::NoWidth)
      log << "z: " << matchInfo.flashZWidth << " cm";
    if (matchInfo.flashYWidth != CRTPMTMatching::NoWidth) {
      if (matchInfo.flashZWidth != CRTPMTMatching::NoWidth)
        log << ", ";
      log << "y: " << matchInfo.flashYWidth << " cm";
    }
    log << ")";
  }
  
  // CRT hits
  for (auto const& [ iCRT, matchedHit ]: util::enumerate(matchInfo.matchedCRTHits)) {
    log << "\n" << indent << "  CRT hit #" << iCRT;
    if (matchedHit.time != MatchedCRT::NoTime)
      log << " at " << matchedHit.time << " us";
    if (matchedHit.PMTTimeDiff != MatchedCRT::NoTime)
      log << " (" << (matchedHit.PMTTimeDiff * 1000) << " ns from flash)";
    log << " at " << matchedHit.position << " cm";
    if (matchedHit.sys != MatchedCRT::NoLocation
      || matchedHit.region != MatchedCRT::NoLocation)
    {
      log << "; sys/region: " << matchedHit.sys << "/" << matchedHit.region;
    }
  } // for
  
  if (flash) {
    log << "\n" << indent << "  associated flash: "
      << tagIndex(flash).encode() << "#"
      << flash.key() << " at " << flash->Time() << " us and (";
    if (flash->hasXCenter()) log << flash->XCenter();
    else log << "n/a";
    log << ", " << flash->YCenter() << ", " << flash->ZCenter() << ") cm";
  }
  if (!CRThits.empty()) {
    log << "\n" << indent << "  associated CRT hits:";
    for (art::Ptr<sbn::crt::CRTHit> const& hitPtr: CRThits) {
      log
        << "  " << tagIndex(hitPtr).encode() << "#" << hitPtr.key()
        << " at T0=" << hitPtr->ts0() << " ns, T1="
        << hitPtr->ts1()/1000.0 << " us, (" << hitPtr->x_pos << ", "
        << hitPtr->y_pos << ", " << hitPtr->z_pos << ") cm, " << hitPtr->peshit
        << " p.e., on " << hitPtr->tagger << ";";
    }
  }
  
} // icarus::crt::DumpCRTPMTMatching::DumpMatching()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::crt::DumpCRTPMTMatching)

//------------------------------------------------------------------------------
