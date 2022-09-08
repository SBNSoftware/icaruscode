/**
 * @file  icaruscode/Analysis/trigger/details/TriggerResponseManager.cxx
 * @brief Helper managing the trigger response part of a TTree (implementation).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date November 24, 2021
 * @see icaruscode/Analysis/trigger/details/TriggerResponseManager.h
 *
 */

// library header
#include "icaruscode/Analysis/trigger/details/TriggerResponseManager.h"

// LArSoft libraries
#include "larcorealg/CoreUtils/counter.h"

// framework libraries
#include "art/Framework/Principal/Handle.h"
// #include "canvas/Persistency/Common/FindOneP.h"
// #include "canvas/Persistency/Common/Assns.h"
// #include "canvas/Persistency/Common/Ptr.h"

// ROOT libraries
#include "TTree.h"
#include "TBranch.h"

// C/C++ libraries
#include <vector>
#include <string>
#include <cassert>


// -----------------------------------------------------------------------------
// ---  sbn::details::TriggerResponseManager::Extractors
// -----------------------------------------------------------------------------
sbn::details::TriggerResponseManager::Extractors::Extractors
  (art::Event const& event, std::vector<TriggerInfoBranch_t>& branchInfoList)
  : fBranchInfo(&branchInfoList)
{
  
  for (TriggerInfoBranch_t const& branchInfo: *fBranchInfo)
    fInputData.emplace_back(buildInputData(event, branchInfo));
  
} // sbn::details::TriggerResponseManager::Extractors::Extractors()


// -----------------------------------------------------------------------------
void sbn::details::TriggerResponseManager::Extractors::fetch(std::size_t iEntry)
{
  assert(fBranchInfo);
  
  std::vector<TriggerInfoBranch_t> const& branchInfo = *fBranchInfo;
  assert(branchInfo.size() == fInputData.size());
  
  for (std::size_t iBranch: util::counter(branchInfo.size()))
    fetchBranch(branchInfo[iBranch], fInputData[iBranch], iEntry);
  
} // sbn::details::TriggerResponseManager::Extractors::fetch()


// -----------------------------------------------------------------------------
void sbn::details::TriggerResponseManager::Extractors::consumesInputData
  (art::ConsumesCollector& collector, TriggerInfoBranch_t const& branchInfo)
{
  
  collector.consumes<std::vector<raw::Trigger>>(branchInfo.triggerTag);
  
} // sbn::details::TriggerResponseManager::Extractors::consumesInputData()


// -----------------------------------------------------------------------------
auto sbn::details::TriggerResponseManager::Extractors::buildInputData
  (art::Event const& event, TriggerInfoBranch_t const& branchInfo) const
  -> TriggerInputData_t
{
  
  return TriggerInputData_t{
    event.getValidHandle<std::vector<raw::Trigger>>(branchInfo.triggerTag)
      .product()  // triggers
    };
  
} // sbn::details::TriggerResponseManager::Extractors::buildInputData()


// -----------------------------------------------------------------------------
void sbn::details::TriggerResponseManager::Extractors::fetchBranch(
  TriggerInfoBranch_t const& info,
  TriggerInputData_t& inputData,
  std::size_t iEntry
) {
  assert(inputData.triggers);
  
  raw::Trigger const& trigger = inputData.triggers->at(iEntry);
  TriggerInfo_t& data = *(info.data);
  
  data.fired = (trigger.TriggerBits() != 0);
  data.triggerTime = trigger.TriggerTime();
  data.gateTime = trigger.BeamGateTime();
  
} // sbn::details::TriggerResponseManager::Extractors::fetchBranch()


// -----------------------------------------------------------------------------
// ---  sbn::details::TriggerResponseManager
// -----------------------------------------------------------------------------
std::string const&
sbn::details::TriggerResponseManager::TriggerInfo_t::TriggerResponseBranchStructure()
{
  // using a static std::string constant data member caused segmentation fault
  // at the closure of the program; so we get just a tiny bit less static.
  static std::string const specs { "time/D:gateStart/D:fired/O" };
  return specs;
}


// -----------------------------------------------------------------------------
sbn::details::TriggerResponseManager::TriggerResponseManager(
  std::vector<TriggerInputSpec_t> const& triggerSpecs,
  art::ConsumesCollector& collector,
  TTree& tree
)
  : fBranchInfo{ buildTriggerResponseBranches(tree, triggerSpecs) }
{
  
  declareConsumables(collector);
  
} // sbn::details::TriggerResponseManager::TriggerResponseManager()


// -----------------------------------------------------------------------------
auto sbn::details::TriggerResponseManager::extractorsFor
  (art::Event const& event) -> Extractors
{
  return Extractors{ event, fBranchInfo };
} // sbn::details::TriggerResponseManager::extractorsFor()


// -----------------------------------------------------------------------------
void sbn::details::TriggerResponseManager::declareConsumables
  (art::ConsumesCollector& collector) const
{
  for (TriggerInfoBranch_t const& branchInfo: fBranchInfo)
    Extractors::consumesInputData(collector, branchInfo);
} // sbn::details::TriggerResponseManager::declareConsumables()


// -----------------------------------------------------------------------------
auto sbn::details::TriggerResponseManager::buildTriggerResponseBranches
  (TTree& tree, std::vector<TriggerInputSpec_t> const& triggerSpecs) const
  -> std::vector<TriggerInfoBranch_t>
{
  std::vector<TriggerInfoBranch_t> branchInfoList;
  for (TriggerInputSpec_t const& spec: triggerSpecs)
    branchInfoList.push_back(buildTriggerResponseBranch(tree, spec));
  return branchInfoList;
} // sbn::details::TriggerResponseManager::buildTriggerResponseBranches()


// -----------------------------------------------------------------------------
auto sbn::details::TriggerResponseManager::buildTriggerResponseBranch
  (TTree& tree, TriggerInputSpec_t const& spec) const -> TriggerInfoBranch_t
{
  TriggerInfoBranch_t branchInfo {
      spec.name     // name
    , spec.inputTag // triggerTag
    // the rest is default-constructed
    }; // TriggerInfoBranch_t
  
  branchInfo.branch = tree.Branch(
    spec.name.c_str(), branchInfo.data.get(),
    TriggerInfo_t::TriggerResponseBranchStructure().c_str()
    );
  
  return branchInfo;
} // sbn::details::TriggerResponseManager::buildTriggerResponseBranches()


// -----------------------------------------------------------------------------
