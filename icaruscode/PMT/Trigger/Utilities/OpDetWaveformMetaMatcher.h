/**
 * @file   icaruscode/PMT/Trigger/Utilities/OpDetWaveformMetaMatcher.h
 * @brief  Utilities for matching `raw::OpDetWaveform` and their
 *         `sbn::OpDetWaveformMeta`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 10, 2022
 */

#ifndef ICARUSCODE_PMT_TRIGGER_UTILITIES_OPDETWAVEFORMMETAMATCHER_H
#define ICARUSCODE_PMT_TRIGGER_UTILITIES_OPDETWAVEFORMMETAMATCHER_H


// ICARUS libraries
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h" // sbn::OpDetWaveformMeta
#include "sbnobj/ICARUS/PMT/Trigger/Data/OpticalTriggerGate.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Persistency/Provenance/ProductID.h"

// C/C++ standard libraries
#include <algorithm> // std::lower_bound()
#include <unordered_map>
// #include <utility> // std::move()


// -----------------------------------------------------------------------------
namespace icarus::trigger {
  template <typename Event> class OpDetWaveformMetaMatcher;
}
/**
 * @brief Object to facilitate the discovery of the `raw::OpDetWaveform` a
 *        `sbn::OpDetWaveformMeta` objects comes from.
 * @tparam Event type of the framework event to read data products from
 *
 * This algorithm will look in the associated _art_ event for the
 * `raw::OpDetWaveform` associated to any specified `sbn::OpDetWaveformMeta`.
 * 
 * The `sbn::OpDetWaveformMeta` object must be specified by _art_ pointer.
 * The algorithm assumes that the module that created the
 * `sbn::OpDetWaveformMeta` object also created an association between it and
 * the original `raw::OpDetWaveform`, and will read that association.
 *
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * icarus::trigger::OpDetWaveformMetaMatcher waveformMetaMatcher{ event };
 * art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform> waveAssns;
 * for (auto const [ gatePtr, metaPtr ]: metaAssns)
 *       waveAssns.addSingle(gatePtr, waveformMetaMatcher(metaPtr));
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will fill `waveAssns` with one gate/waveform association for each
 * gate/metadata association in `metaAssns`.
 * Note that, in this simple example, if no waveform pointer is matched to a
 * metadata item, a null pointer to waveform is stored in the association.
 * 
 */
template <typename Event>
class icarus::trigger::OpDetWaveformMetaMatcher {
  
  using Event_t = Event; ///< Type of the framework event to read data from.
  
  /// Type of associations used in this algorithm,
  using WaveformMetaAssns_t
    = art::Assns<sbn::OpDetWaveformMeta, raw::OpDetWaveform>;
  
  /// Compares a sbn::OpDetWaveformMeta` and an element of `WaveformMetaAssns_t`
  struct CmpFirst {
    
    using Assn_t = typename WaveformMetaAssns_t::assn_t; // a pair
    using MetaPtr_t = typename Assn_t::first_type; // an art pointer
    
    template <typename A, typename B>
    bool operator() (A const& a, B const& b) const { return cmp(a, b); }
    
    static bool cmp(Assn_t const& a, Assn_t const& b)
      { return cmp(a.first, b.first); }
    
    static bool cmp(Assn_t const& a, MetaPtr_t const& b)
      { return cmp(a.first, b); }
    
    static bool cmp(MetaPtr_t const& a, Assn_t const& b)
      { return cmp(a, b.first); }
    
    static bool cmp(MetaPtr_t const& a, MetaPtr_t const& b)
      { return a.key() < b.key(); }
    
  }; // struct CmpFirst
  
  
  Event_t const& fEvent; ///< Event to read associations from.
  
  /// All associations discovered so far.
  std::unordered_map<art::ProductID, WaveformMetaAssns_t const*> fAssns;
  
  
  /// Loads the association including `pid`.
  /// @return pointer to the associations containing `pid`, `nullptr` if n/a.
  WaveformMetaAssns_t const* loadAssociations(art::ProductID const& pid);
  
  /// Returns the waveform associated to `meta` if in `assns` or a null pointer.
  art::Ptr<raw::OpDetWaveform> findAssociatedWaveform(
    art::Ptr<sbn::OpDetWaveformMeta> const& meta,
    WaveformMetaAssns_t const& assns
    ) const;
  
    public:
  
  /// Constructor: associates to the specified art event.
  OpDetWaveformMetaMatcher(Event const& event);
  
  // @{
  /// Returns the waveform associated to `meta`, or a null pointer if not found.
  art::Ptr<raw::OpDetWaveform> fetchAssociatedWaveform
    (art::Ptr<sbn::OpDetWaveformMeta> const& meta);
  art::Ptr<raw::OpDetWaveform> operator() // alias
    (art::Ptr<sbn::OpDetWaveformMeta> const& meta)
    { return fetchAssociatedWaveform(meta); }
  // @}
  
  
}; // class icarus::trigger::OpDetWaveformMetaMatcher



//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
template <typename Event>
icarus::trigger::OpDetWaveformMetaMatcher<Event>::OpDetWaveformMetaMatcher
  (Event const& event)
  : fEvent(event)
  {}


//------------------------------------------------------------------------------
template <typename Event>
art::Ptr<raw::OpDetWaveform>
icarus::trigger::OpDetWaveformMetaMatcher<Event>::fetchAssociatedWaveform
  (art::Ptr<sbn::OpDetWaveformMeta> const& meta)
{
  
  WaveformMetaAssns_t const* assns = loadAssociations(meta.id());
  return assns
    ? findAssociatedWaveform(meta, *assns): art::Ptr<raw::OpDetWaveform>{};
  
} // icarus::trigger::OpDetWaveformMetaMatcher<>::fetchAssociatedWaveform()


//------------------------------------------------------------------------------
template <typename Event>
auto icarus::trigger::OpDetWaveformMetaMatcher<Event>::loadAssociations
  (art::ProductID const& pid) -> WaveformMetaAssns_t const*
{
  
  // if already there, the better
  if (auto const itAssns = fAssns.find(pid); itAssns != fAssns.end())
    return itAssns->second;
  
  // let's start optimistic: we won't find it
  auto& assnCache = (fAssns[pid] = nullptr);
  
  // find the product label of pid
  cet::exempt_ptr<art::BranchDescription const> metaDescr
    = fEvent.getProductDescription(pid);
  if (!metaDescr) return nullptr;
  
  // shall we check that the product type is the right one? no, we shan't
  art::InputTag const tag = metaDescr->inputTag();
  
  auto const assnsHandle = fEvent.template getHandle<WaveformMetaAssns_t>(tag);
  if (!assnsHandle.isValid()) return nullptr;
  
  return assnCache = assnsHandle.product();
  
} // icarus::trigger::OpDetWaveformMetaMatcher<>::loadAssociations()


//------------------------------------------------------------------------------
template <typename Event>
art::Ptr<raw::OpDetWaveform>
icarus::trigger::OpDetWaveformMetaMatcher<Event>::findAssociatedWaveform(
  art::Ptr<sbn::OpDetWaveformMeta> const& meta,
  WaveformMetaAssns_t const& assns
) const {
  
  auto const key = meta.key();
  
  // try our luck: 1-1 association?
  if ((assns.size() > key) && (assns.at(key).first == meta))
    return assns[key].second;
  
  // nope, go binary search
  if (auto itAssn = std::lower_bound(assns.begin(), assns.end(), meta, CmpFirst{});
      itAssn != assns.end()
  ) {
    if (itAssn->first == meta) return itAssn->second;
  }
  
  // still nope, go linear search
  for (auto const& assn: assns) if (assn.first == meta) return assn.second;
  
  // still nope, not found!
  return {};
  
} // icarus::trigger::OpDetWaveformMetaMatcher::fetchAssociatedWaveform()


//------------------------------------------------------------------------------

#endif // ICARUSCODE_PMT_TRIGGER_UTILITIES_OPDETWAVEFORMMETAMATCHER_H
