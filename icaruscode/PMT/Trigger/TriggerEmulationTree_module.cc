/**
 * @file    TriggerEmulationTree_module.cc
 * @authors Gianluca Petrillo (SLAC, petrillo@slac.stanford.edu)
 * @date    February 18, 2022
 * 
 * Borrowed from `TriggerEmulationTree` module used for trigger efficiency
 * analysis in data.
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Utilities/TriggerResponseManager.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"

// LArSoft libraries
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger
#include "lardataobj/Simulation/BeamGateInfo.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TTree.h"

// C/C++ libraries
#include <vector>
#include <string>



// -----------------------------------------------------------------------------
namespace sbn { class TriggerEmulationTree; }
/**
 * @brief Fills a ROOT tree with trigger emulation results.
 * 
 * 
 * Trigger information
 * ====================
 * 
 * There are two distinct groups of information about triggers in the tree.
 * The first one is the hardware trigger, that is the trigger that was actually
 * used to decide to record the event. 
 * The other group of information is from trigger logic "emulation" on the
 * collected data (optical detector waveforms): there may be multiple branches
 * in this group, one for each emulated trigger logic, threshold and gate.
 * 
 * 
 * Simulated trigger information
 * ------------------------------
 * 
 * A branch is added for each configured trigger logic.
 * This module actually ignores the details of the logic yielding the results:
 * a trigger result is a group of data products under the same tag.
 * For example, a tag `M1` (presumably, a very loose trigger logic requiring
 * just one PMT pair above threshold anywhere in the detector, with no hint of
 * which that threshold is) will produce a branch with name `"M1"` including
 * information from the `M1` data products (so far, only a collection of
 * `raw::Trigger` is read).
 * Note that a trigger object is expected to be present _regardless whether
 * the trigger fired or not_. It is assumed that the trigger fired if at least
 * one of the `raw::Trigger::TriggerBits()` is set, not fired otherwise.
 * 
 * Each branch mirrors a data record, `TriggerInfo_t`, reporting the response
 * for one simulated trigger logic.
 * The structure and interpretation of the trigger response branch is described
 * in ROOT TTree code by the string in
 * `TriggerInfo_t::TriggerResponseBranchStructure`, and its meaning is:
 * * `fired` (boolean): whether this trigger has triggered at all during the
 *   considered trigger gate.
 * * `time` (double): time of the trigger, taken directly from
 *   `raw::Trigger::TriggerTime()`; it is expected to be measured in
 *   microseconds and relative to the nominal time reference of the event,
 *   which for data is the hardware trigger time and for simulation is the
 *   "hardware" beam gate time (however that is defined).
 * * `gate` (double): start time of the gate where the trigger logic was
 *   simulated. It is in the same time scale as the trigger time, and
 *   directly taken from `raw::Trigger::BeamGateTime()`.
 * 
 */
class sbn::TriggerEmulationTree: public art::EDAnalyzer {
  
  using TriggerInputSpec_t
    = sbn::details::TriggerResponseManager::TriggerInputSpec_t;
  
    public:
  
  /// Data record the trigger response branch is based on.
  using TriggerInfo_t = sbn::details::TriggerResponseManager::TriggerInfo_t;
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    /// Information on a single trigger source for the tree.
    struct TriggerSpecConfig {
      
      fhicl::Atom<std::string> Name {
        fhicl::Name("Name"),
        Comment("name of the trigger (e.g. `\"M5O3\"`), used for branch name")
        };
      fhicl::Atom<art::InputTag> TriggerTag {
        fhicl::Name("TriggerTag"),
        Comment("tag of the input trigger info data product")
        };
      
    }; // TriggerSpecConfig
    using TriggerSpecConfigTable
      = fhicl::TableAs<TriggerInputSpec_t, TriggerSpecConfig>;
    
    fhicl::Atom<art::InputTag> BeamGateProducer {
      Name("BeamGateProducer"),
      Comment("tag of beam gate information")
      // mandatory
      };
    
    fhicl::Atom<art::InputTag> TriggerProducer {
      Name("TriggerProducer"),
      Comment("tag of hardware trigger information")
      // mandatory
      };
    
    fhicl::Sequence<TriggerSpecConfigTable> EmulatedTriggers {
      Name("EmulatedTriggers"),
      Comment("the emulated triggers to include in the tree")
    };
    
    fhicl::Atom<std::string> TreeName {
      Name("TreeName"),
      Comment("name of the output ROOT tree"),
      "TriggerEmu" // default
    };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("label for output messages of this module instance"),
      "TriggerEmulationTree" // default
      };
    
  }; // Config
  
  using Parameters = art::EDAnalyzer::Table<Config>;
  
  explicit TriggerEmulationTree(Parameters const& config);

  void analyze(art::Event const& event) override;
  
  void endJob() override;

    private:
  
  // --- BEGIN -- Tree data structures -----------------------------------------
  
  /// Data structure for the beam gate data in the tree.
  struct BeamInfo_t {
    static constexpr char branchDef[]
      = "beamGateStart/l:beamGateDuration/F:beamGateType/i";
    
    std::uint64_t beamGateSimStart = 0;
    float beamGateDuration = -1.0;
    unsigned int beamGateType = 999;
  }; // BeamInfo_t
  
  /// Data structure for the global trigger data in the tree.
  struct GlobalTriggerInfo_t {
    static constexpr char branchDef[]
      = "beamType/i:triggerTime/l:beamGateTime/l:triggerID/i:gateID/i";
    
    unsigned int beamType = 0U;
    std::uint64_t triggerTime = 0ULL;
    std::uint64_t beamGateTime = 0ULL;
    unsigned int triggerID = 0U;
    unsigned int gateID = 0U;
  }; // GlobalTriggerInfo_t
  
  // --- END ---- Tree data structures -----------------------------------------
  
  // --- BEGIN -- Configuration parameters -------------------------------------

  art::InputTag const fBeamGateProducer;
  art::InputTag const fTriggerProducer;
  std::string const fLogCategory;
  
  // --- END ---- Configuration parameters -------------------------------------

  // --- BEGIN -- Tree buffers -------------------------------------------------
  
  struct {
    unsigned int run;
    unsigned int subRun;
    unsigned int event;
  } fEventID; // anonymous
  
  BeamInfo_t fBeamInfo;
  GlobalTriggerInfo_t fGlobalTriggerInfo; ///< Data for global trigger in tree.
  
  // --- END ---- Tree buffers -------------------------------------------------
  
  TTree *fStoreTree = nullptr;

  unsigned int fTotalProcessed = 0U;
  
  
  // `convert()` needs to be a free function
  friend TriggerInputSpec_t convert(Config::TriggerSpecConfig const& config);
  
  
  // --- BEGIN -- Trigger response branches ------------------------------------
  
  ///< Manages filling of trigger result branches.
  sbn::details::TriggerResponseManager fTriggerResponses;
  
  // --- END ---- Trigger response branches ------------------------------------
  
}; // sbn::TriggerEmulationTree


// -----------------------------------------------------------------------------
namespace sbn {
  
  TriggerEmulationTree::TriggerInputSpec_t convert
    (TriggerEmulationTree::Config::TriggerSpecConfig const& config)
  {
    return {
        config.Name()       // name
      , config.TriggerTag() // inputTag
      };
  } // convert(sbn::TriggerSpecConfig)
  
} // namespace sbn


// -----------------------------------------------------------------------------
sbn::TriggerEmulationTree::TriggerEmulationTree(Parameters const& config)
  : art::EDAnalyzer{ config }
  // configuration
  , fBeamGateProducer { config().BeamGateProducer() }
  , fTriggerProducer  { config().TriggerProducer() }
  , fLogCategory      { config().LogCategory() }
  // algorithms
  , fStoreTree {
      art::ServiceHandle<art::TFileService>()->make<TTree>
        (config().TreeName().c_str(), "Trigger emulation results")
      }
  , fTriggerResponses
      { config().EmulatedTriggers(), consumesCollector(), *fStoreTree }
{
  
  //
  // declaration of input
  //
  
  consumes<std::vector<sim::BeamGateInfo>>(fBeamGateProducer);

  //
  // tree population
  //
  fStoreTree->Branch("run", &fEventID.run);
  fStoreTree->Branch("subrun", &fEventID.subRun);
  fStoreTree->Branch("event", &fEventID.event);
  fStoreTree->Branch("beamInfo", &fBeamInfo, BeamInfo_t::branchDef);
  fStoreTree->Branch
    ("globalTrigger", &fGlobalTriggerInfo, GlobalTriggerInfo_t::branchDef);
  
} // sbn::TriggerEmulationTree::TriggerEmulationTree()


void sbn::TriggerEmulationTree::analyze(art::Event const& event)
{
  // Implementation of required member function here.
  fEventID = { event.run(), event.subRun(), event.event() };

  //
  // beam gate
  //
  fBeamInfo = {};
  
  auto const& beamgate
    = event.getProduct<std::vector<sim::BeamGateInfo>>(fBeamGateProducer);
  if(beamgate.empty())
    mf::LogWarning(fLogCategory) << "No Beam Gate Information!";
  if(beamgate.size() > 1U)
    mf::LogWarning(fLogCategory) << "Event has multiple beam gate info labels! (maybe this changes later to be standard)";
  fBeamInfo.beamGateSimStart = beamgate.front().Start();
  fBeamInfo.beamGateDuration = beamgate.front().Width();
  fBeamInfo.beamGateType = beamgate.front().BeamType();

  //
  // hardware trigger
  //
  auto const& triggerinfo
    = event.getProduct<sbn::ExtraTriggerInfo>(fTriggerProducer);
  fGlobalTriggerInfo.beamType
    = static_cast<unsigned int>(triggerinfo.sourceType);
  fGlobalTriggerInfo.triggerTime = triggerinfo.triggerTimestamp;
  fGlobalTriggerInfo.beamGateTime = triggerinfo.beamGateTimestamp;
  fGlobalTriggerInfo.triggerID = triggerinfo.triggerID;
  fGlobalTriggerInfo.gateID = triggerinfo.gateID;
  
  // get an extractor bound to this event
  sbn::details::TriggerResponseManager::Extractors triggerResponseExtractors
    = fTriggerResponses.extractorsFor(event);
  
  triggerResponseExtractors.fetch(0); // only the first trigger in each list

  ++fTotalProcessed;
  fStoreTree->Fill();

} // sbn::TriggerEmulationTree::analyze()


// -----------------------------------------------------------------------------
void sbn::TriggerEmulationTree::endJob() {
  
  mf::LogInfo(fLogCategory) << "Processed " << fTotalProcessed << " events.";
  
} // sbn::TriggerEmulationTree::endJob()


// -----------------------------------------------------------------------------

DEFINE_ART_MODULE(sbn::TriggerEmulationTree)
