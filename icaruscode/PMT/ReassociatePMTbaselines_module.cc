/**
 * @file   ReassociatePMTbaselines_module.cc
 * @brief  Associates baseline objects to pruned optical waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 6, 2023
 */


// ICARUS libraries
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker
#include "icaruscode/PMT/Data/WaveformRMS.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icarusalg/Utilities/AssnsCrosser.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "larcorealg/CoreUtils/counter.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <utility> // std::move()
#include <memory> // std::unique_ptr
#include <algorithm> // std::sort(), std::lower_bound()
#include <vector>
#include <string>
#include <optional>
#include <functional> // std::less
#include <cassert>


//------------------------------------------------------------------------------
namespace icarus::trigger { class ReassociatePMTbaselines; }

/**
 * @brief Associates baseline and RMS objects to pruned optical waveforms.
 * 
 * The issue at hand is that baselines and their RMS are computed and assigned
 * based on the full information available on a PMT channel, and then associated
 * to all the waveforms of that channel. Then, those waveforms are dropped but
 * some of them are _copied_ to a new data product. The issue is that this copy
 * does not replicate the associations of the waveforms being copied.
 * So we are left with all the baseline objects, all their associations to the
 * original waveforms, no original waveform, and some other waveforms which are
 * copies of part of the original ones. We are left to figure out which baseline
 * belongs to which copied waveform.
 * 
 * In addition, baseline and RMS information do not contain a channel number.
 * 
 * This module applies heuristic algorithm based on the assumptions:
 *  * associations between baselines and the original waveforms are still
 *    present (even if the pointers to those waveforms can't be resolved);
 *  * metadata of the original waveforms is still available and associated to
 *    the original waveforms.
 *  * RMS objects are in the same order as the baseline objects; it is assumed
 *    that baseline of index _i_ has RMS of index _i_.
 * 
 * So the plan is to:
 * 1. Identify the original waveform metadata matching the available waveform
 *    copies (based on channel and timestamp);
 * 2. Associate the metadata of the original waveforms with the baseline
 *    by hopping through their waveform association;
 * 3. Not to repeat the association traversal for RMS but rather use the result
 *    from the baseline one, assuming RMS to mirror baselines.
 * 4. Finally, create a new set of baseline and RMS objects associated
 *    to the available waveform copies, in the same order (similar protocol).
 * 
 * 
 * Output data products
 * =====================
 * 
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>`,
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformRMS>`: associations between
 *   each available waveform and its baseline and RMS.
 * * `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>` (if
 *   `RecreateMetaAssns` is set): associations between each available waveform
 *   and its metadata. Normally produced by `CopyBeamTimePMTwaveforms` module,
 *   which was buggy until <a href="https://github.com/SBNSoftware/icaruscode/pull/617"><tt>icaruscode</tt> PR #617</a>) was merged.
 * 
 * 
 * Input data products
 * =====================
 * 
 * * `std::vector<raw::OpDetWaveform>` (`WaveformTag`): the available waveforms
 *   which need to be associated to baselines and RMS.
 * * `std::vector<sbn::OpDetWaveformMeta>` (`OriginalWaveformTag`): the
 *   metadata of the _original_ waveforms (which may have been dropped).
 * * `art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>
 *   (`OriginalWaveformTag`): the associations of the metadata above with the
 *   original waveforms.
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>` (`BaselineTag`):
 *   associations of the baselines to the original waveforms
 * * `std::vector<icarus::WaveformRMS>` (optional, `RMSTag`): the RMS of
 *   each of the original waveforms (does not need to be available).
 * 
 * 
 * Service requirements
 * ---------------------
 * 
 * The following services are _required_:
 * * _art_ message facility
 * * an implementation of `DetectorClocksService` for the determination of the
 *   times for waveform metadata.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description ReassociatePMTbaselines`.
 * 
 * * `WaveformTag` (input tag): the data product containing all surviving
 *   optical detector waveforms
 * * `BaselineTag` (input tag): the data product containing one
 *   baseline per waveform in data product (from `OpticalWaveforms`)
 * * `RMSTag` (input tag, optional): the data product containing one baseline
 *   per waveform in data product (from `OpticalWaveforms`)
 * * `OriginalWaveformTag` (input tag): the data product creating the metadata
 *   of the original waveforms, and possibly the associations to them.
 * * `RecreateMetaAssns` (flag, default: `false`): if set, associations between
 *   the waveforms in `WaveformTag` and the metadata of the corresponding ones
 *   in the original waveform data product are (re)created.
 * * `LogCategory` (string, default: `"ReassociatePMTbaselines"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * 
 * 
 */
class icarus::trigger::ReassociatePMTbaselines: public art::SharedProducer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    
    fhicl::Atom<art::InputTag> WaveformTag{
      Name{ "WaveformTag" },
      Comment{ "input tag of digitized optical waveform data product" }
      };
    
    fhicl::Atom<art::InputTag> BaselineTag{
      Name{ "BaselineTag" },
      Comment{ "input tag for waveform baselines associations to waveforms" }
      };
    
    fhicl::OptionalAtom<art::InputTag> RMSTag{
      Name{ "RMSTag" },
      Comment{ "input tag for waveform RMS associations to waveforms" }
      };
    
    fhicl::Atom<art::InputTag> OriginalWaveformTag{
      Name{ "OriginalWaveformTag" },
      Comment{ "input tag for original waveform baselines and metadata" }
      };
    
    fhicl::Atom<bool> RecreateMetaAssns{
      Name{ "RecreateMetaAssns" },
      Comment{ "recreate the associations between metadata and waveforms" },
      false
      };
    
    fhicl::Atom<std::string> LogCategory{
      Name{ "LogCategory" },
      Comment{ "tag of the module output to console via message facility" },
      "ReassociatePMTbaselines"
      };
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  explicit ReassociatePMTbaselines
    (Parameters const& config, art::ProcessingFrame const& frame);
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  /// Creates the data products.
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  /// Sorting for `sbn::OpDetWaveformMeta` objects.
  /// Compares channel, then start time, then stop time.
  struct OpDetWaveformMetaSorter {
    bool operator()
      (sbn::OpDetWaveformMeta const& a, sbn::OpDetWaveformMeta const& b)
      const noexcept;
  };
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fWaveformTag; ///< Optical waveform input tag.
  
  art::InputTag const fBaselineTag; ///< Waveform baseline input tag.
  art::InputTag const fRMSTag; ///< Waveform RMS input tag.
  
  art::InputTag const fOriginalWaveformTag; ///< Original waveform input tag.
  
  bool const fRecreateMetaAssns; ///< Whether to recreate metadata associations.
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  // --- END Service variables -------------------------------------------------
  
  
  
}; // icarus::trigger::ReassociatePMTbaselines



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Keeps a map of the _art_ pointers of some objects.
   * @tparam T type of the tracked objects; required to be comparable
   * @tparam Comp (default: `std::less<T>`) type of functor comparison
   * 
   * The object looks up a _art_ pointer for an object which compares equal
   * to a requested one.
   * 
   * To create an `PtrFinder` object, use `makePtrFinder()` helper function.
   * 
   * Requirements
   * -------------
   * 
   * Object of type `Comp` must support:
   * * `bool operator() (const T& a, const T& b)` returning whether `a` is
   *   strictly before `b` (for whatever definition). This is a strict order
   *   comparison.
   * 
   */
  template <typename T, typename Comp = std::less<T>>
  class PtrFinder {
    
    using Data_t = T;
    
    using Ptr_t = art::Ptr<Data_t>;
    
    using Storage_t = std::vector<Ptr_t>;
    
    struct Comparer {
      
      Comp cmp;
      
      bool operator() (Ptr_t const& a, Ptr_t const& b) const
        { return cmp(*a, *b); }
      bool operator() (Ptr_t const& a, Data_t const& b) const
        { return cmp(*a, b); }
      bool operator() (Data_t const& a, Ptr_t const& b) const
        { return cmp(a, *b); }
      bool operator() (Data_t const& a, Data_t const& b) const
        { return cmp(a, b); }
      
    }; // Comparer
    
    
    Storage_t fMap; ///< Cached map of objects.
    
    Comparer fComp; ///< Comparison operator
    
    /// Hidden constructor: create an object with `makeIndexMap()`.
    PtrFinder(Storage_t sortedPtrs, Comp comp)
      : fMap{ std::move(sortedPtrs) }, fComp{ std::move(comp) } {}
    
      public:
    
    /// Null pointer returned on lookup failure.
    static art::Ptr<T> const NullPtr;
    
    /// Returns the _art_ pointer of the object equal to `obj`.
    art::Ptr<T> const& operator[] (T const& obj) const
      { return findObject(obj); }
    
    /// Returns the _art_ pointer of the object equal to `obj`.
    art::Ptr<T> const& findObject(T const& obj) const;
    
    
    // friends
    template <typename Handle, typename OtherComp>
    friend PtrFinder<typename Handle::element_type::value_type, OtherComp>
      makePtrFinder(Handle const&, OtherComp);
    
  }; // class PtrFinder
  
  
  // ---------------------------------------------------------------------------
  /// Creates a `IndexMap` object with the specified collection and comparison.
  template
    <typename Handle, typename Comp = std::less<typename Handle::value_type>>
  PtrFinder<typename Handle::element_type::value_type, Comp> makePtrFinder
    (Handle const& handle, Comp comp = Comp{})
  {
    using Data_t = typename Handle::element_type::value_type;
    
    using PtrFinder_t = PtrFinder<Data_t, Comp>;
    
    typename PtrFinder_t::Storage_t sortedPtrs;
    art::fill_ptr_vector(sortedPtrs, handle);
    
    std::sort(
      begin(sortedPtrs), end(sortedPtrs), typename PtrFinder_t::Comparer{ comp }
      );
    
    return PtrFinder_t{ std::move(sortedPtrs), std::move(comp) };
  } // makePtrFinder()
  
  
  // ---------------------------------------------------------------------------
  template <typename T>
  std::unique_ptr<T> moveToUniquePtr(T& obj)
    { return std::make_unique<T>(std::move(obj)); }
  
} // local namespace

namespace sbn {
  
  std::ostream& operator<<
    (std::ostream& out, sbn::OpDetWaveformMeta const& meta)
  {
    out << "Ch=" << meta.ChannelNumber() << " Time=[ " << meta.startTime
      << " ; " << meta.endTime << " ] (" << meta.nSamples << " samples)";
    // should also print flags... I don't care here
    return out;
  }
  
} // namespace sbn


//------------------------------------------------------------------------------
//--- PtrFinder
//------------------------------------------------------------------------------
namespace {
  
  template <typename T, typename Comp /* = std::less<T> */>
  art::Ptr<T> const PtrFinder<T, Comp>::NullPtr;
  
  
  template <typename T, typename Comp /* = std::less<T> */>
  art::Ptr<T> const& PtrFinder<T, Comp>::findObject(T const& obj) const {
    
    auto const it = std::lower_bound(fMap.begin(), fMap.end(), obj, fComp);
    // lower_bound() points to an object `*it` no smaller than `obj`,
    // so `comp(*it, obj)` is `false`: `*it` is greater or equal to `obj`;
    // if `comp(obj, *it)` is `false` too, `obj` is _equal_ to `*it` and we have
    // a match; otherwise, we still don't have a match.
    return ((it == fMap.end()) || fComp(obj, *it))? NullPtr: *it;
    
  } // PtrFinder<>::findObject()
  
  
} // local namespace

//------------------------------------------------------------------------------
//--- icarus::trigger::ReassociatePMTbaselines
//------------------------------------------------------------------------------
icarus::trigger::ReassociatePMTbaselines::ReassociatePMTbaselines
  (Parameters const& config, art::ProcessingFrame const& frame)
  : art::SharedProducer{ config }
  // configuration
  , fWaveformTag        { config().WaveformTag()                   }
  , fBaselineTag        { config().BaselineTag()                   }
  , fRMSTag             { config().RMSTag().value_or(fBaselineTag) }
  , fOriginalWaveformTag{ config().OriginalWaveformTag()           }
  , fRecreateMetaAssns  { config().RecreateMetaAssns()             }
  , fLogCategory        { config().LogCategory()                   }
  // service cache
{
  async<art::InEvent>();
  
  //
  // configuration post-processing
  //
  {
    mf::LogInfo log(fLogCategory);
    log << "Configuration:"
      "\n * match waveforms '" << fWaveformTag.encode()
        << "' to waveform metadata '" << fOriginalWaveformTag.encode()
        << "' via baselines '" << fBaselineTag.encode()
        << "' and RMS '" << fRMSTag.encode()
      ;
    if (fRecreateMetaAssns) {
      log
        << "\n * recreate associations between matched waveforms and metadata";
    }
  } // nameless block
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  consumes<std::vector<sbn::OpDetWaveformMeta>>(fOriginalWaveformTag);
  consumes<art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>>
    (fOriginalWaveformTag);
  consumes<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>
    (fBaselineTag);
  consumes<art::Assns<raw::OpDetWaveform, icarus::WaveformRMS>>(fRMSTag);
  // consumes<std::vector<icarus::WaveformBaseline>>(fBaselineTag);

  //
  // declaration of output
  //
  
  produces<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>();
  produces<art::Assns<raw::OpDetWaveform, icarus::WaveformRMS>>();
  if (fRecreateMetaAssns)
    produces<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>();
  
} // icarus::trigger::ReassociatePMTbaselines::ReassociatePMTbaselines()


//------------------------------------------------------------------------------
void icarus::trigger::ReassociatePMTbaselines::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
 /* 
  * So the plan is to:
  * 1. associate the metadata of the original waveforms with the baseline and RMS
  *    by hopping through their waveform association;
  * 2. identify the original waveform metadata matching the available waveform
  *    copies (based on channel and timestamp);
  * 3. finally, create a new set of baseline and RMS objects associated
  *    to the available waveform copies, in the same order (similar protocol).
  */
  
  detinfo::DetectorTimings const detTimings{
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  
  //
  // prepare input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fWaveformTag);
  
  auto const& waveformRMShandle
    = event.getValidHandle<std::vector<icarus::WaveformRMS>>(fRMSTag);
  
  auto const& origMetaHandle
    = event.getValidHandle<std::vector<sbn::OpDetWaveformMeta>>
      (fOriginalWaveformTag)
    ;
  
  // associate original metadata with (original) baselines:
  auto const metaToBaseline = makeAssnsCrosser(event,
      icarus::ns::util::startFrom<sbn::OpDetWaveformMeta>{}
    , icarus::ns::util::hopTo    <raw::OpDetWaveform>{ fOriginalWaveformTag }
    , icarus::ns::util::hopTo    <icarus::WaveformBaseline>{ fBaselineTag }
    );
  
  // find which original metadata we need
  art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline> waveToBaseline;
  art::Assns<raw::OpDetWaveform, icarus::WaveformRMS> waveToRMS;
  art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> waveToMeta;
  
  PtrFinder const metadataPointers
    = makePtrFinder(origMetaHandle, OpDetWaveformMetaSorter{});
  
  sbn::OpDetWaveformMetaMaker const makeMetadata{ detTimings };
 
  for (auto const iWave: util::counter(waveformHandle->size())) {
    
    art::Ptr<raw::OpDetWaveform> wavePtr{ waveformHandle, iWave };
    raw::OpDetWaveform const& waveform{ *wavePtr };
    
    // recreate the metadata for this waveform and match it to the known one
    sbn::OpDetWaveformMeta const meta{ makeMetadata(waveform) };
    art::Ptr<sbn::OpDetWaveformMeta> origMetaPtr
      = metadataPointers.findObject(meta);
    if (!origMetaPtr) {
      // this warning suggests either a bug in the algorithm, a violation of
      // the assumptions or some wrong configuration.
      mf::LogProblem{ fLogCategory }
        << "Waveform '" << fWaveformTag.encode() << "' #" << wavePtr.key()
        << " (" << meta << ") couldn't be matched to any original waveform.";
      continue;
    }
    
    art::Ptr<icarus::WaveformBaseline> baselinePtr
      = metaToBaseline.assPtr(origMetaPtr);
    if (!baselinePtr) {
      // this warning suggests either a violation of the assumptions or wrong
      // configuration.
      mf::LogProblem{ fLogCategory }
        << "Waveform '" << fWaveformTag.encode() << "' #" << wavePtr.key()
        << " (" << meta << ") couldn't be matched to any baseline.";
      continue;
    }
    
    art::Ptr<icarus::WaveformRMS> RMSptr
      { waveformRMShandle, baselinePtr.key() };
    
    waveToBaseline.addSingle(wavePtr, std::move(baselinePtr));
    waveToRMS.addSingle(wavePtr, std::move(RMSptr));
    if (fRecreateMetaAssns)
      waveToMeta.addSingle(wavePtr, std::move(origMetaPtr));
    
  } // for
  
  //
  // send the data products to the event
  //
  event.put(moveToUniquePtr(waveToBaseline));
  event.put(moveToUniquePtr(waveToRMS));
  if (fRecreateMetaAssns) event.put(moveToUniquePtr(waveToMeta));
  
} // icarus::trigger::ReassociatePMTbaselines::produce()


//------------------------------------------------------------------------------
bool
icarus::trigger::ReassociatePMTbaselines::OpDetWaveformMetaSorter::operator()
  (sbn::OpDetWaveformMeta const& a, sbn::OpDetWaveformMeta const& b)
  const noexcept
{
  if (a.channel != b.channel) return a.channel < b.channel;
  if (a.startTime != b.startTime) return a.startTime < b.startTime;
  return a.endTime < b.endTime;
} // icarus::trigger::ReassociatePMTbaselines::OpDetWaveformMetaSorter::operator()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::ReassociatePMTbaselines)


//------------------------------------------------------------------------------
