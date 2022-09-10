/**
 * @file  icaruscode/Analysis/trigger/details/TriggerResponseManager.h
 * @brief Helper managing the trigger response part of a TTree.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date November 24, 2021
 *
 *
 */

#ifndef ICARUSCODE_ANALYSIS_TRIGGER_DETAILS_TRIGGERRESPONSEMANAGER_H
#define ICARUSCODE_ANALYSIS_TRIGGER_DETAILS_TRIGGERRESPONSEMANAGER_H


// LArSoft libraries
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"
#include "lardataobj/RawData/TriggerData.h" // raw::Trigger

// framework libraries
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

// C/C++ libraries
#include <vector>
#include <string>
#include <memory> // std::make_unique<>
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
// forward declarations
class TTree;
class TBranch;

// -----------------------------------------------------------------------------
namespace sbn::details { class TriggerResponseManager; }
/**
 * @brief Manages extraction of trigger results and filling of their branches.
 * 
 * This class defines and manages the structure of the branches pertaining the
 * trigger simulation in a ROOT tree.
 * 
 * It supports multiple trigger algorithms in a single ROOT tree entry
 * (which may be defined for example as a physics event or as a particle within
 * a physics event).
 * 
 * 
 * Configuration
 * --------------
 * 
 * The manager is configured with a list of trigger specifications
 * (`TriggerInputSpec_t`).
 * One simple branch (with multiple leaves) is created for each specifications.
 * 
 * 
 * Usage pattern
 * --------------
 * 
 * This manager class encloses all the steps needed to fill a tree with trigger
 * response data. Users will have to take very few steps to integrate its
 * functionality in a module.
 * 
 * 1. A single instance of this class is created associated to a ROOT tree,
 *    configured with a set of trigger specifications. A class data member is
 *    recommended for this pattern.
 * 2. On each _art_ event, an extractor object is created (`extractFor()`)
 *    which will manage the extraction of trigger data from that event.
 * 3. On each entry in the tree, `Extractor::fetch()` is called for that entry,
 *    which will fill all the relevant branch data from the available
 *    information. This single steps makes all trigger data ready for filling.
 * 4. When all additional data is ready, the tree can be filled with the entry
 *    data (`TTree::Fill()`).
 * 
 * The configuration of this class is performed via a custom configuration data
 * structure. The module may read it directly from its own (FHiCL) configuration
 * or fill it in any other way.
 * 
 * 
 * @note This class owns the memory associated to the tree branches.
 *       As such, it needs to exist as long as the buffers are needed, and those
 *       buffers must not change memory location.
 *       For this reason, the class is not copyable, but it is moveable (moving
 *       will preserve the address of the buffers).
 *
 */
class sbn::details::TriggerResponseManager
  : private lar::UncopiableClass // this class is moveable but not copyable
{
  
  struct TriggerInfoBranch_t;
  
    public:
  
  /**
   * @brief Information about a single trigger logic (hardware or emulated).
   * 
   * This data structure is the base for the tree branch.
   * Each instance of the data structure represents a single trigger logic
   * response.
   * 
   * Default constructor represents a trigger that did not fire at all.
   */
  struct TriggerInfo_t {
    
    /// ROOT TTree specification for this data structure (leaf list).
    static std::string const& TriggerResponseBranchStructure();
    
    /// Mnemonic value for absence of trigger time information
    static constexpr double NotTriggeredTime = -999999.0;
    
    /// Time of the trigger
    /// (@ref DetectorClocksElectronicsTime "electronics time scale").
    double triggerTime = NotTriggeredTime;

    /// Time of the opening of the gate for the trigger evaluation.
    double gateTime = NotTriggeredTime;
    
    bool fired = false; ///< Whether this trigger fired.
    
  }; // TriggerInfo_t
  
  
  /// Configuration specifications for the emulation of a trigger logic.
  struct TriggerInputSpec_t {
    std::string name;
    art::InputTag inputTag;
  }; // TriggerInputSpec_t

  /// Trigger information extractors tied to an event.
  class Extractors {
    
    /// Data pertaining a single branch.
    struct TriggerInputData_t {
      std::vector<raw::Trigger> const* triggers; ///< Trigger results.
    }; // struct TriggerInputData_t
    
    /// Pointer to the information for each supported branch.
    std::vector<TriggerInfoBranch_t> const* fBranchInfo { nullptr };
    
    /// Data for each branch in this event (index parallel to `fBranchInfo`).
    std::vector<TriggerInputData_t> fInputData;
    
    
    /// Returns all data necessary to `branchInfo` extracted from `event`.
    TriggerInputData_t buildInputData
      (art::Event const& event, TriggerInfoBranch_t const& branchInfo) const;
    
    /// Fills the specified trigger logic branch with information from its
    /// trigger results using the `inputData` entry with index `iEntry`.
    void fetchBranch(
      TriggerInfoBranch_t const& info,
      TriggerInputData_t& inputData,
      std::size_t iEntry
      );
    
      public:
    
    /// Reads all data from `events` needed for branches in `branchInfoList`.
    Extractors(
      art::Event const& event,
      std::vector<TriggerInfoBranch_t>& branchInfoList
      );

    /// Fills all branch data with information from trigger entries with index
    /// `iEntry`.
    void fetch(std::size_t iEntry);
    
    
    /// Declares all data products that need to be read for `branchInfo`.
    static void consumesInputData(
      art::ConsumesCollector& collector, TriggerInfoBranch_t const& branchInfo
      );
    
  }; // class Extractors
  
  
  /// Initializes `tree` to accommodate the specified trigger information.
  TriggerResponseManager(
    std::vector<TriggerInputSpec_t> const& triggerSpecs,
    art::ConsumesCollector& collector,
    TTree& tree
    );
  
  
  /// Returns an object to extract trigger information from `event`.
  Extractors extractorsFor(art::Event const& event);
  
  
    private:
  
  /// Data for a single trigger logic output branch.
  struct TriggerInfoBranch_t {
    std::string name;
    art::InputTag triggerTag;
    std::unique_ptr<TriggerInfo_t> data = std::make_unique<TriggerInfo_t>();
    TBranch* branch = nullptr;
  }; // TriggerInfoBranch_t
  
  
  /// Data structures for the tree.
  std::vector<TriggerInfoBranch_t> fBranchInfo;
  
  
  /// Sets up the tree branches and returns the branch information structures.
  std::vector<TriggerInfoBranch_t> buildTriggerResponseBranches
    (TTree& tree, std::vector<TriggerInputSpec_t> const& triggerSpecs) const;
  
  /// Sets up a tree branch and returns its branch information structure.
  TriggerInfoBranch_t buildTriggerResponseBranch
    (TTree& tree, TriggerInputSpec_t const& spec) const;

  /// Declares all the data products we are going to read.
  void declareConsumables(art::ConsumesCollector& collector) const;
  
}; // class sbn::details::TriggerResponseManager

// -----------------------------------------------------------------------------


#endif // ICARUSCODE_ANALYSIS_TRIGGER_DETAILS_TRIGGERRESPONSEMANAGER_H
