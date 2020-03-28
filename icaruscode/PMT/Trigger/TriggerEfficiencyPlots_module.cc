/**
 * @file   TriggerEfficiencyPlots_module.cc
 * @brief  Plots of trigger efficiency.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 9, 2020
 */


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Algorithms/TriggerTypes.h" // ADCCounts_t
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // FillTriggerGates()
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icaruscode/PMT/Trigger/Utilities/ROOTutils.h" // util::ROOT
#include "icaruscode/PMT/Trigger/Data/OpticalTriggerGate.h"
#include "icaruscode/PMT/Trigger/Data/MultiChannelOpticalTriggerGate.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/TensorIndices.h" // util::MatrixIndices
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/MCDumpers/MCDumperUtils.h" // sim::TruthInteractionTypeName
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds, ...
#include "lardataalg/Utilities/quantities/energy.h" // megaelectronvolt, ...
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/get_elements.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#if 0
// #include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::::static_assert_on<>
#include "lardataobj/RawData/OpDetWaveform.h"
#endif // 0

// nutools libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "cetlib_except/exception.h"
#if 0
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"
#endif // 0

// ROOT libraries
#include "TEfficiency.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

// C/C++ standard libraries
#include <ostream>
#include <atomic>
#include <map>
#include <vector>
#include <string>
#include <functional> // std::function<>, std::reference_wrapper<>
#include <memory> // std::unique_ptr
#include <utility> // std::pair<>, std::move()
#include <cstddef> // std::size_t


//------------------------------------------------------------------------------
using MeV = util::quantities::megaelectronvolt;
using GeV = util::quantities::gigaelectronvolt;

using namespace util::quantities::time_literals; // ""_ns ...


//------------------------------------------------------------------------------
/// Information about the event.
struct EventInfo_t {
  
  /// Constructor. As if nobody noticed.
  EventInfo_t() { fInteractions.fill(0U); }
  
  // --- BEGIN Query interface -----------------------------------------------
  /// @name Query interface
  /// @{

  /// Returns the number of weak charged current interactions in the event.
  unsigned int nWeakChargedCurrentInteractions() const
    { return fInteractions[itWCC]; }

  /// Returns the number of weak neutral current interactions in the event.
  unsigned int nWeakNeutralCurrentInteractions() const
    { return fInteractions[itWNC]; }

  /// Returns the number of weak current interactions in the event.
  unsigned int nWeakCurrentInteractions() const
    {
      return
        nWeakChargedCurrentInteractions() + nWeakNeutralCurrentInteractions();
    }

  /// Returns whether the event is generated as a neutrino CC interaction.
  bool isWeakChargedCurrent() const
    { return nWeakChargedCurrentInteractions() > 0U; }

  /// Returns whether the event is generated as a neutrino NC interaction.
  bool isWeakNeutralCurrent() const
    { return nWeakNeutralCurrentInteractions() > 0U; }

  /// Returns whether the event is generated as a neutrino interaction.
  bool isNeutrino() const
    { return isWeakChargedCurrent() || isWeakNeutralCurrent(); }
  
  /// Returns which neutrino flavor is present in an event
  bool isNu_mu() const { return nu_mu; }
  bool isNu_e() const { return nu_e; }

  /// Returns the neutrino PDG code
  int NeutrinoPDG() const { return fNeutrinoPDG; }

  /// Returns the total energy deposited in the detector during the event [GeV]
  GeV DepositedEnergy() const { return fEnergyDepTotal; }
  
  /// Returns the total energy deposited in the detector during beam [GeV]
  GeV DepositedEnergyInSpill() const { return fEnergyDepSpill; }

  /// Returns the energy deposited in the active volume during the event [GeV]
  GeV DepositedEnergyInActiveVolume() const { return fEnergyDepActive; }
  
  /// Returns the energy deposited in the active volume during the beam [GeV]
  GeV DepositedEnergyInSpillInActiveVolume() const
    { return fEnergyDepSpillActive; }

  // Returns the neutrino energy [GeV]
  GeV NeutrinoEnergy() const { return fNeutrinoEnergy; }

  // Returns the lepton energy [GeV]
  GeV LeptonEnergy() const { return fLeptonEnergy; }

  // Returns the interaction type
  int InteractionType() const { return fInteractionType; }

  /// Returns whether this type of event has a known vertex.
  bool hasVertex() const { return !fVertices.empty(); }
  
  /// Returns whether there is an interaction within the active volume.
  bool isInActiveVolume() const { return fInActiveVolume; }
  
  /// @}
  // --- END Query interface -------------------------------------------------


  // --- BEGIN Set interface -------------------------------------------------
  /// @name Set interface
  /// @{

  /// Marks this event as including _n_ more weak charged current interactions.
  void AddWeakChargedCurrentInteractions(unsigned int n = 1U)
    { fInteractions[itWCC] += n; }

  /// Marks this event as including _n_ more weak neutral current interactions.
  void AddWeakNeutralCurrentInteractions(unsigned int n = 1U)
    { fInteractions[itWNC] += n; }

  /// Marks the neutrino flavor
  void SetNu_mu(bool numu) { nu_mu = numu; }
  void SetNu_e(bool nue) { nu_e = nue; }

  /// Marks this event's neutrino type
  void SetNeutrinoPDG(int NU) { fNeutrinoPDG = NU; }

  /// Sets the total deposited energy of the event [GeV]
  void SetDepositedEnergy(GeV e) { fEnergyDepTotal = e; }

  /// Sets the energy of the event deposited during beam gate [GeV]
  void SetDepositedEnergyInSpill(GeV e) { fEnergyDepSpill = e; }

  /// Sets the total deposited energy of the event in active volume [GeV]
  void SetDepositedEnergyInActiveVolume(GeV e) { fEnergyDepActive = e; }

  /// Sets the energy of the event deposited during beam gate in active volume
  /// [GeV]
  void SetDepositedEnergyInActiveVolumeInSpill(GeV e)
    { fEnergyDepSpillActive = e; }

  // Sets the neutrino energy.
  void SetNeutrinoEnergy(GeV eNu) { fNeutrinoEnergy = eNu; }

  // Sets the lepton energy.
  void SetLeptonEnergy(GeV eL) { fLeptonEnergy = eL; }

  // Sets the interaction type
  void SetInteractionType(int type) { fInteractionType = type; }

  /// Set whether the event has relevant activity in the active volume.
  void SetInActiveVolume(bool active = true) { fInActiveVolume = active; }
  
  /// Adds a point to the list of interaction vertices in the event.
  void AddVertex(geo::Point_t const& vertex) { fVertices.push_back(vertex); }
  
  std::vector<geo::Point_t> fVertices; ///< Position of all vertices.

  /// @}
  // --- END Set interface ---------------------------------------------------

  /// Prints the content of the object into a stream.
  void dump(std::ostream& out) const
    {
      out << "Event contains:";
      if (isNeutrino()) {
        if (nWeakChargedCurrentInteractions())
          out << " " << nWeakChargedCurrentInteractions() << " CC";
        if (nWeakNeutralCurrentInteractions())
          out << " " << nWeakNeutralCurrentInteractions() << " NC";
        if (isNu_mu()) out << " nu_mu";
        if (isNu_e()) out << " nu_e";
        out << "\nThe first neutrino has E=" << NeutrinoEnergy()
          << " and becomes a lepton with E=" << LeptonEnergy()
          << " with a " << sim::TruthInteractionTypeName(InteractionType())
          << " interaction"
          ;
      }
      else {
        out << " no neutrino interaction";
      }
      out << "\nTotal deposited energy: " << DepositedEnergy()
        << ", of which in spill " << DepositedEnergyInSpill()
        << ", in active volume " << DepositedEnergyInActiveVolume()
        << ", in active volume and in spill "
          << DepositedEnergyInSpillInActiveVolume();
      if (fVertices.empty()) {
        out << "\nNo interaction vertex found.";
      }
      else {
        auto iVertex = fVertices.begin();
        auto const vend = fVertices.end();
        out
          << "\n" << fVertices.size() << " interaction vertices: " << *iVertex;
        while (++iVertex != vend) out << "; " << *iVertex;
        out << ".";
      }
      out << "\nThe event is" << (isInActiveVolume()? "": " NOT")
        << " marked as in the active volume of the detector.";
      out << "\n";
      out.flush();
    } // dump()

    private:

  // --- BEGIN interaction type constants ------------------------------------
  
  static constexpr std::size_t itWCC { 0U }; ///< Charged weak current.
  static constexpr std::size_t itWNC { 1U }; ///< Neutral weak current.
  static constexpr std::size_t NInteractionTypes { 2U };

  // --- END interaction type constants --------------------------------------

  std::array<unsigned int, NInteractionTypes> fInteractions;

  GeV fEnergyDepTotal       { 0.0 }; ///< Total deposited energy.
  GeV fEnergyDepSpill       { 0.0 }; ///< Energy deposited in spill.
  GeV fEnergyDepActive      { 0.0 }; ///< Energy deposited in active volume.
  /// Energy deposited in active volume in spill.
  GeV fEnergyDepSpillActive { 0.0 };

  int fNeutrinoPDG { 0 };
  int fInteractionType { 0 };

  GeV fNeutrinoEnergy { 0.0 };
  GeV fLeptonEnergy { 0.0 };
  //GeV fNucleonEnergy { 0.0 };
  
  bool nu_mu { false };
  bool nu_e { false };
  
  /// Whether the event has activity inside the active volume.
  bool fInActiveVolume { false };
  
  //std::vector<geo::Point_t> fVertices; ///< Position of all vertices.
  
}; // struct EventInfo_t

std::ostream& operator<< (std::ostream& out, EventInfo_t const& info)
  { info.dump(out); return out; }


//------------------------------------------------------------------------------
class PlotCategory {
  
    public:
  
  /// Type of test function.
  using QualifyFunc_t = std::function<bool(EventInfo_t const&)>;
  
  
  /// Constructor from category name and test function.
  PlotCategory(
    std::string name, std::string descr = {},
    QualifyFunc_t&& test = [](EventInfo_t const&){ return true; }
    )
    : fName(std::move(name)), fDescr(std::move(descr)), fTest(std::move(test))
    {}
  
  /// Returns the name of the category.
  std::string const& name() const { return fName; }
  
  /// Returns the description of the category.
  std::string const& description() const { return fDescr; }
  
  /// Returns whether the event belong to this category.
  bool test(EventInfo_t const& info) const { return fTest(info); }
  
  operator std::string() const { return name(); }
  bool operator() (EventInfo_t const& info) const { return test(info); }
  
    private:
  
  std::string fName;
  std::string fDescr;
  QualifyFunc_t fTest;
  
}; // class PlotCategory
using PlotCategories_t = std::vector<PlotCategory>;

/**
 * @brief List of event categories.
 * 
 * category name  | condition
 * -------------- | ------------------------------------------------------------
 * `All`          | any event
 *  ---Nu_mu      |
 *  ---Nu_e       |
 * `NuCC`         | at least one generated charged current neutrino interaction
 *  ---Nu_mu      |
 *  ---Nu_e       |
 * `NuNC`         | at least one generated neutral current neutrino interaction
 *  ---Nu_mu      |
 *  ---Nu_e       |
 * 
 * 
 */
PlotCategories_t const PlotCategories {

  PlotCategory{
    "All"
    },

  PlotCategory{
	"All nu_mu", "nu_mu",
    [](EventInfo_t const& info){ return info.isNu_mu(); }
	},

  PlotCategory{
	"All nu_e", "nu_e",
    [](EventInfo_t const& info){ return info.isNu_e(); }
	},

  PlotCategory{
    "NuCC", "CC",
    [](EventInfo_t const& info){ return info.isWeakChargedCurrent(); }
    },

  PlotCategory{
    "NuCC_mu", "CC_mu",
    [](EventInfo_t const& info){ return (info.isWeakChargedCurrent() & info.isNu_mu()); }
    },

  PlotCategory{
    "NuCC_e", "CC_e",
    [](EventInfo_t const& info){ return (info.isWeakChargedCurrent() & info.isNu_e()); }
    },

  PlotCategory{
    "NuNC", "NC",
    [](EventInfo_t const& info){ return info.isWeakNeutralCurrent(); }
    },

  PlotCategory{
    "NuNC_mu", "NC_mu",
    [](EventInfo_t const& info){ return (info.isWeakNeutralCurrent() & info.isNu_mu()); }
    },

  PlotCategory{
    "NuNC_e", "NC_e",
    [](EventInfo_t const& info){ return (info.isWeakNeutralCurrent() & info.isNu_e()); }
    }

}; // PlotCategories[]



// --- BEGIN -- ROOT tree helpers ----------------------------------------------
struct TreeHolder {

  TreeHolder() = default;
  TreeHolder(TTree& tree): fTree(&tree) {}

  TTree& tree() { return *fTree; }
  TTree const& tree() const { return *fTree; }

    private:
  TTree* fTree = nullptr;

}; // struct TreeHolder


/**
 * @brief Class managing the serialization of event ID in a simple ROOT tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assignEvent()`, the branch addresses are assigned the values from the
 * event ID. The tree is not `Fill()`-ed.
 *
 * The tree structure is: `Run/i:SubRun/i:Event/i`, with a single branch per
 * element.
 */
struct EventIDTree: public TreeHolder {

  /// Creates the required branches and assigns addresses to them.
  EventIDTree(TTree& tree);

  /// Fills the information of the specified event.
  void assignID(art::EventID const& id);
  void assignEvent(art::Event const& event) { assignID(event.id()); }

  UInt_t fRunNo;
  UInt_t fSubRunNo;
  UInt_t fEventNo;

}; // struct EventIDTree


/**
 * @brief Class managing the serialization of plot information in a simple ROOT
 *        tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assign()`, the branch addresses are assigned the values from the
 * arguments. The tree is not `Fill()`-ed.
 *
 * The tree structure is:
 * `InPlots/O`,
 * with a single branch per element.
 *
 * Branches:
 *  * `InPlots` (bool): the event was _not_ filtered out before plotting
 *    (it may still belong to no category and eventually appear in no plot)
 *
 */
struct PlotInfoTree: public TreeHolder {

  /// Creates the required branches and assigns addresses to them.
  PlotInfoTree(TTree& tree);

  /**
   * @brief Fills the information of the specified event.
   * @param inPlots whether this event is plotted (as opposed to filtered out)
   */
  void assign(bool inPlots);

  Bool_t fInPlots;

}; // struct PlotInfoTree


/**
 * @brief Class managing the serialization of event information in a simple ROOT
 *        tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assignEvent()`, the branch addresses are assigned the values from the
 * event information. The tree is not `Fill()`-ed.
 *
 * The tree structure is:
 * `CC/i:NC/i:IntType:I/NuE:D/OutLeptE:D/TotE/D:SpillE/D:InActive/O`,
 * with a single branch per element.
 *
 * Branches:
 *  * `CC` (unsigned integer): number of neutrino CC interactions in the event
 *  * `NC` (unsigned integer): number of neutrino NC interactions in the event
 *  * `IntType` (integer): code of interaction type (see `simb::int_type_`)
 *  * `NuE` (double): energy of the generated initial state neutrino [GeV]
 *  * `OutLeptE` (double): energy of the generated final state lepton [GeV]
 *  * `TotE` (double): total deposited energy in the event [GeV]
 *  * `SpillE` (double): total deposited energy during the beam gate [GeV]
 *  * `InActive` (bool): whether an interaction happened in active volume'
 *      this requires an interaction vertex (e.g. cosmic rays are out)
 *
 */
struct EventInfoTree: public TreeHolder {

  /// Creates the required branches and assigns addresses to them.
  EventInfoTree(TTree& tree);

  /**
   * @brief Fills the information of the specified event.
   * @param info event information to fill the tree with
   * @param inPlots whether this event is plotted (as opposed to filtered out)
   */
  void assignEvent(EventInfo_t const& info);

  UInt_t fCC;
  UInt_t fNC;
  Int_t fIntType;
  Double_t fNuE;
  Double_t fOutLeptE;
  Double_t fTotE;
  Double_t fSpillE;
  Double_t fActiveE;
  Double_t fSpillActiveE;
  
  Bool_t fInActive;
  
}; // struct EventInfoTree


/**
 * @brief Class managing the serialization of trigger responses in a simple ROOT
 *        tree.
 *
 * The tree is supplied by the caller.
 * This object will create the proper branches into the tree and assign
 * addresses to them. Then it will assume they will stay assigned.
 *
 * On `assignResponse()`, the proper branch address is assigned the specified
 * trigger response (`true` or `false`).
 *
 * The branch structure is: a `RespTxxRxx/O` branch for each threshold (numeric)
 * and requirement (also numeric).
 * with a single branch per element.
 *
 */
struct ResponseTree: public TreeHolder {

  /// Constructor: accommodates that many thresholds and requirements.
  template <typename Thresholds, typename Requirements>
  ResponseTree
    (TTree& tree, Thresholds const& thresholds, Requirements const& minReqs);

  /// Assigns the response for the specified trigger.
  void assignResponse(std::size_t iThr, std::size_t iReq, bool resp);

  // `std::vector<bool>` is too special for us. Let's pack it to go.
  util::MatrixIndices indices;
  std::unique_ptr<bool[]> RespTxxRxx;

}; // struct ResponseTree


// --- END -- ROOT tree helpers ------------------------------------------------


//------------------------------------------------------------------------------
namespace icarus::trigger { class TriggerEfficiencyPlots; }
/**
 * @brief Produces plots about trigger simulation and trigger efficiency.
 * 
 * This module produces sets of plots based on trigger primitives given in
 * input.
 * 
 * The following documentation mainly deals with the most standard configuration
 * and operations in ICARUS. This module and the ones upstream of it have quite
 * some knobs that can be manipulated to test unorthodox configurations.
 * 
 * 
 * Overview of trigger simulation steps
 * =====================================
 * 
 * This simulation only trigger primitives derived from the optical detector
 * via hardware (V1730B boards).
 * The trigger simulation branches from the standard simulation of optical
 * detector waveforms (`icarus::SimPMTIcarus` module).
 * From there, multiple steps are needed.
 * 
 * 1. Produce single-PMT-channel discriminated waveforms: the discrimination
 *    threshold can be chosen ("ADC threshold" or just "threshold" in the
 *    following), and the result is one binary discriminated waveform that
 *    has value `0` when the waveform is under threshold and `1` when it's
 *    above threshold, with the same time discretization as the PMT waveform.
 *    Each discriminated waveform spans the whole event readout time and
 *    therefore it may merge information from multiple readout waveforms
 *    (`raw::OpDetWaveform`), but all from the same optical detector channel.
 *    This step can be performed with the module
 *    `icarus::trigger::DiscriminatePMTwaveforms`.
 * 2. Combine the discriminated waveforms in pairs, in the same fashion as the
 *    V1730B readout board does to produce LVDS output. The output of this step
 *    is roughly half as many discriminated outputs as there were from the
 *    previous step (some channels are not paired). This output will be called
 *    _trigger primitives_ because it is what the trigger logics is based on.
 *    We say that a trigger primitive is "on" when its level is `1`.
 *    This step can be performed with the module `icarus::trigger::LVDSgates`.
 * 3. Simulate the trigger logic based on the trigger primitives _(see below)_.
 *    This usually includes the requirement of coincidence with the beam gate.
 * 
 * Trigger logic may be complex, being implemented in a FPGA.
 * Many options are available, including:
 * 
 * * coincidence of a minimum number of trigger primitives: the event is trigger
 *   if at least _N_ trigger primitives are on at the same time;
 * * sliding window: primitives are grouped depending on their location, and
 *   there is a requirement on the minimum number of primitives in "on" level
 *   within one or more of these groups. For example, any sliding window of 30
 *   channels (corresponding to 16 trigger primitives in standard ICARUS PMT
 *   readout installation) has to contain at least 10 trigger primitives "on"
 *   at the same time; or there must be two sliding windows with that condition;
 *   etc. Sliding windows can be obtained from further processing of the LVDS
 *   trigger primitives, for example with the module
 *   `icarus::trigger::SlidingWindowTrigger`.
 * 
 * This module (`icarus::trigger::TriggerEfficiencyPlots`) is designed for
 * implementing a trigger logic based on coincidence of trigger primitives,
 * and to plot trigger efficiency and related information based on a requirement
 * of a minimum number of trigger primitives, in that enabling the simulation
 * in the fashion of the first option above. More details on the steps performed
 * by this module are documented @ref TriggerEfficiencyPlots_Algorithm "later".
 * 
 * For the implementation of sliding window, different logic needs to be
 * implemented, possibly keeping track of the location of each primitive, and
 * different plots are needed as well. This module may provide a template for
 * such implementation, but it hasn't been designed with enough flexibility to
 * accommodate that directly.
 * 
 * 
 * @anchor TriggerEfficiencyPlots_Data
 * Data objects for discriminated waveforms
 * -----------------------------------------
 *
 * A discriminated waveform is the information whether the level of a waveform
 * is beyond threshold, as function of time.
 * A discriminated waveform may be binary, i.e. with only levels `0` and `1`
 * based on a single threshold, or with multiple levels.
 * Also the numerical _addition_ of two binary discriminated waveforms
 * yields a multi-level waveform (in fact, three levels -- `0`, `1` and `2`).
 * 
 * We represent this data in the forms of "events": an event is a change of
 * level happening at a certain time. The class holding this information,
 * `icarus::trigger::TriggerGateData`, covers the whole time, starting with a
 * ground level `0`. The next event will be a "opening" that increases the
 * level, usually to `1`. Other changing events may follow, and typically the
 * last one will bring the level back to `0`.
 * 
 * This information is joined by a list of _channel numbers_ in order to
 * represent a discriminated waveform e.g. from the optical detector.
 * There may be one or more channels associated to a discriminated waveform,
 * but for us there should always be at least one.
 * The discriminated waveforms from PMT readout (`raw::OpDetWaveform`) are
 * associated to a single channel, while LVDS trigger primitives are most often
 * associated to two channels (some channels are not paired and will have only
 * one channel). A sliding window will have as many channels as the PMT it
 * covers. The global trigger should have _all_ channels, while in ICARUS each
 * of the two discriminated wabeforms from a cryostat should have half the
 * channels.
 * This information is represented in the class
 * `icarus::trigger::ReadoutTriggerGate`, which inherits from
 * `icarus::trigger::TriggerGateData`.
 * This class is generic and can hold any representation for the time of the
 * level changing events, for the levels, and for the identifiers of the
 * channels. ICARUS specifies a data type for each of these quantities, and
 * the resulting `icarus::trigger::ReadoutTriggerGate` class instance is called
 * `icarus::trigger::OpticalTriggerGateData_t`.
 * 
 * @note The class `icarus::trigger::OpticalTriggerGateData_t` is the one that
 *       gets written to disk in the _art_ ROOT files.
 *       That is not a class by itself, but rather an alias of
 *       `icarus::trigger::ReadoutTriggerGate`, and many coding tools will call
 *       it in the latter way.
 * 
 * The class `icarus::trigger::OpticalTriggerGate` is currently the most
 * commonly used in the code. It adds to the information of
 * `icarus::trigger::ReadoutTriggerGate`, from which it derives, a list of
 * optical waveforms (`raw::OpDetWaveform`) it originates from.
 * 
 * Finally, the classes `icarus::trigger::SingleChannelOpticalTriggerGate` and
 * `icarus::trigger::MultiChannelOpticalTriggerGate` do not add any information
 * to the `icarus::trigger::OpticalTriggerGate` they derive from, but they
 * have an interface explicitly tuned for discriminated waveforms from a single
 * channel or from multiple channels, respectively (for example, the former
 * provides a `channel()` method returning a single channel identifier, while
 * the latter provides a `channels()` method returning a list of channels).
 * 
 * These three higher level classes, `icarus::trigger::OpticalTriggerGate` and
 * derivatives, _can't be directly saved_ in _art_ ROOT files.
 * There are utilities available in
 * `icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h` that can convert them
 * into a collection of `icarus::trigger::OpticalTriggerGateData_t` objects
 * plus a collection of _art_ associations (for writing), and the other way
 * around (for reading). The module `icarus::trigger::LVDSgates` uses both sides
 * and can be used as an illustration of the functionality.
 * 
 * A module is provided, called `icarus::trigger::DumpTriggerGateData`, which
 * dumps on screen or on text file the information from a collection of
 * discriminated waveforms in a (sort-of) human readable format.
 * An example configuration for this module is provided in `icaruscode`, called
 * `dump_triggergatedata_icarus.fcl`.
 * 
 * The main functions to manipulate the trigger gates are defined in the very
 * base class, `icarus::trigger::TriggerGateData`: these allow e.g. to find
 * events and query the level at a given time.
 * Another important set of features is also present in
 * `icarus::trigger::TriggerGateData` and replicated in the higher levels:
 * the combination of trigger gates by sum (`OR`), multiplication (_AND_),
 * minimum and maximum value.
 * 
 * @note Combining a multi-level gate via `Min()` with a binary gate results
 *       into a binary gate which is logic AND of the two gates.
 * 
 * 
 * On the terminology
 * -------------------
 * 
 * In this documentation, the "discriminated waveforms" are sometimes called
 * "trigger gates" and sometimes "trigger primitives".
 * Although these concepts are, strictly, different, usually the difference
 * does not affect the meaning and they end up being exchanged carelessly.
 * 
 * 
 * @anchor TriggerEfficiencyPlots_Algorithm
 * Trigger logic algorithm
 * ========================
 * 
 * This section describes the trigger logic algorithm used in
 * `icarus::trigger::TriggerEfficiencyPlots` and its assumptions.
 * 
 * The algorithm treats all the trigger primitives equally, whether they
 * originate from one or from two channels (or 10 or 30), and wherever their
 * channels are in the detector.
 * The trigger primitives are combined in a multi-level gate by adding them,
 * so that the level of the resulting gate describes at any time matches how
 * many trigger primitives are on at that time.
 * 
 * This multi-level gate is set in coincidence with the beam gate by multiplying
 * the multi-level and the beam gates.
 * The beam gate opens at a time configured in `DetectorClocks` service provider
 * (`detinfo::DetectorClocks::BeamGateTime()`) and has a duration configured
 * in this module (`BeamGateDuration`).
 * 
 * At this point, the trigger gate is a multi-level gate suppressed everywhere
 * except than during the beam gate.
 * The algorithm handles multiple trigger definitions, or requirements.
 * Each requirement is simply how many trigger primitives must be open at the
 * same time for the trigger to fire. The values of these requirements are
 * set in the configuration (`MinimumPrimitives`).
 * To determine whether a trigger with a given requirement, i.e. with a required
 * minimum number of trigger primitives open at the same time, has fired, the
 * gate is scanned to find _the first tick_ where the level of the gate reaches
 * or passes this minimum required number. If such tick exists, the trigger is
 * considered to have fired, and at that time.
 * 
 * As a consequence, there are for each event as many different trigger
 * responses as how many different requirements are configured in
 * `MinimumPrimitives`, _times_ how many ADC thresholds are provided in input,
 * configured in `Thresholds`.
 * 
 * While there _is_ a parameter describing the time resolution of the trigger
 * (`TriggerTimeResolution`), this is currently only used for aesthetic purposes
 * to choose the binning of some plots: the resolution is _not_ superimposed
 * to the gates (yet).
 * 
 * This algorithm is currently implemented in
 * `detinfo::trigger::TriggerEfficiencyPlots::plotResponses()`.
 * 
 * 
 * Event categories
 * =================
 * 
 * Each event is assigned to a set of categories, and its information
 * contributes to the plots of all those categories.
 * The categories are defined in the `PlotCategories` list (hard-coded).
 * 
 * 
 * Module usage
 * =============
 * 
 * 
 * Input data products
 * --------------------
 * 
 * This module uses the following data products as input:
 * 
 * * trigger primitives:
 *     * `std::vector<icarus::trigger::OpticalTriggerGateData_t>` (labels out of
 *       `TriggerGatesTag` and `Thresholds`): full sets of discriminated
 *       waveforms, each waveform possibly covering multiple optical channels,
 *       and their associations to optical waveforms. One set per threshold;
 *     * associations with `raw::OpDetWaveform` (currently superfluous);
 * * event characterization:
 *     * `std::vector<simb::MCTruth>`: generator information (from
 *       `GeneratorTags`; multiple generators are allowed at the same time);
 *       it is used mostly to categorize the type of event (background,
 *       weak charged current, electron neutrino, etc.);
 *     * `std::vector<simb::MCParticle>`: particles propagating in the detector
 *       (from `PropagatedParticles`); currently not used;
 *     * `std::vector<sim::SimEnergyDeposit>`: energy deposited in the active
 *       liquid argon volume (from `EnergyDepositTags`); it is used to
 *       quantify the energy available to be detected in the event.
 * 
 * 
 * Output plots
 * -------------
 * 
 * The module produces a standard set of plots for each configured ADC threshold
 * and for each event category.
 * The plots are saved via _art_ service `TFileService` in a ROOT file and they
 * are organized in nested ROOT directories under the module label directory
 * which is assigned by _art_:
 * 
 * * `Thr###` (outer level) describes the ADC threshold on the discriminated
 *   waveforms: the threshold is `###` from the baseline;
 * * `\<Category\>` (inner level) describes the category of events included in
 *   the plots (e.g. `All`, `NuCC`; see `PlotCategories`).
 * 
 * Each of the inner ROOT directories contains a full set of plots, whose name
 * is the standard plot name followed by its event category and threshold
 * (e.g. `Eff_NuCC_Thr15` for the trigger efficiency plot on neutrino charged
 * current events with 15 ADC counts as PMT threshold).
 * 
 * Each set of plots, defined in
 * `icarus::trigger::TriggerEfficiencyPlots::initializePlotSet()`
 * (and, within, `initializeEventPlots()`
 * and `initializeEfficiencyPerTriggerPlots()`), is contributed only by the
 * events in the set category.
 * 
 * There are different "types" of plots. Some
 * @ref TriggerEfficiencyPlots_SelectionPlots "do not depend on triggering at all",
 * like the deposited energy distribution. Others
 * @ref TriggerEfficiencyPlots_MultiTriggerPlots "cross different trigger definitions",
 * like the trigger efficiency as function of trigger requirement. Others still
 * @ref TriggerEfficiencyPlots_SingleTriggerPlots "assume a single trigger definition":
 * this is the case of trigger efficiency plots versus energy. Finally, there are
 * @ref TriggerEfficiencyPlots_SingleTriggerResponsePlots "plots that depend on a specific trigger definition and outcome":
 * this is the case of all the plots including only triggering or non-triggering
 * events.
 * 
 * A list of plots follows for each plot type.
 * All the plots are always relative to a specific optical detector channel
 * threshold (ADC) and a broad event category.
 * 
 * 
 * ### Plots independent of the triggers (selection plots)
 * 
 * @anchor TriggerEfficiencyPlots_SelectionPlots
 * 
 * These plots are stored directly in a threshold/category folder:
 * 
 * * `EnergyInSpill`: total energy deposited in the detector during the time the
 *   beam gate is open. It is proportional to the amount of scintillation light
 *   in the event;
 * * `EnergyInSpillActive`: energy deposited in the active volume of the
 *   detector during the time the beam gate is open; the active volume is
 *   defined as the union of all TPC drift voulmes;
 * * plots specific to neutrino interactions (if the event is not a neutrino
 *   interaction, it will not contribute to them); if not otherwise specified,
 *   only the first neutrino interaction in the event is considered:
 *     * `InteractionType`: code of the interaction type, as in
 *       `sim::TruthInteractionTypeName`;
 *     * `NeutrinoEnergy`: generated energy of the interacting neutrino;
 *     * `LeptonEnergy`: generated energy of the lepton out of the _first_
 *       neutrino interaction;
 *     * `InteractionVertexYZ`: coordinates of the location of all interactions
 *       in the event, in world coordinates, projected on the anode plane.
 * 
 * @note In fact, these plots usually do not even depend on the ADC threshold
 *       of the optical detector channels. Nevertheless, they are stored in the
 *       folders under specific thresholds, and they are duplicate.
 * 
 * 
 * ### Plots including different trigger requirements
 * 
 * @anchor TriggerEfficiencyPlots_MultiTriggerPlots
 * 
 * These plots collect information from scenarios with different trigger
 * requirements, but still with the same underlying optical detector channel
 * ADC threshold.
 * They are also stored in the threshold/category folder:
 * 
 * * `Eff`: trigger efficiency defined as number of triggered events over the
 *   total number of events, as function of the minimum number of trigger
 *   primitives (`MinimumPrimitives`) to define a firing trigger; uncertainties
 *   are managed by `TEfficiency`.
 * * `TriggerTick`: distribution of the time of the earliest trigger for the
 *   event, as function of the minimum number of trigger primitives (as in
 *   `Eff`). It may happen that the event is such that there is e.g. a
 *   20-primitive flash, then subsiding, and then another 30-primitive flash.
 *   In such a case, in the trigger requirement "&geq; 15 primitives" such event
 *   will show at the time of the 20-primitive flash, while in the trigger
 *   requirement "&geq; 25 primitives" it will show at the time of the 
 *   30-primitive flash. Each event appears at most once for each trigger
 *   requirement, and it may not appear at all if does not fire a trigger.
 * * `NPrimitives`: the maximum number of primitives "on" at any time.
 * 
 * 
 * ### Plots depending on a specific trigger definition
 * 
 * @anchor TriggerEfficiencyPlots_SingleTriggerPlots
 * 
 * The following plots depend on the full definition of the trigger, including
 * PMT thresholds _and_ other requirements.
 * They are hosted each in a subfolder of the threshold/category folder, with
 * a name encoding the requirement: `ReqXX` for trigger with minimum required
 * primitives `XX`.
 * 
 * * `EffVsEnergyInSpill`: trigger efficiency as function of the total energy
 *   deposited during the beam gate;
 * * `EffVsEnergyInSpillActive`: trigger efficiency as function of the energy
 *   deposited in the TPC active volumes during the beam gate;
 * * `EffVsNeutrinoEnergy`, `EffVsLeptonEnergy`: trigger efficiency as function
 *   of the true energy of the first interacting neutrino, and of the outgoing
 *   lepton in the final state of the interaction, respectively;
 * * `TriggerTick`: time of the earliest trigger. Only triggering events
 *   contribute.
 * 
 * The parameters are defined in the same way as in the
 * @ref TriggerEfficiencyPlots_SelectionPlots "selection plots", unless stated
 * otherwise.
 *
 *
 * ### Plots depending on a specific trigger definition and response
 *
 * @anchor TriggerEfficiencyPlots_SingleTriggerResponsePlots
 *
 * These plots depend on the trigger definition, as the ones in the previous
 * type, and on the actual trigger response.
 * They are hosted each in a subfolder of the threshold/category/requirement
 * folder, with a name encoding the response: `triggering` for triggering
 * events, `nontriggering` for the others.
 * 
 * Their event pool is filtered to include only the events in the current
 * category which also pass, or do not pass, the trigger requirement.
 * 
 * The plots are:
 * * `EnergyInSpill`, `EnergyInSpillActive`, `InteractionType`,
 *   `NeutrinoEnergy`,`LeptonEnergy`, `InteractionVertexYZ`:
 *   these are defined in the same way as the
 *   @ref TriggerEfficiencyPlots_SelectionPlots "selection plots"
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description TriggerEfficiencyPlots`.
 * 
 * * `TriggerGatesTag` (string, mandatory): name of the module
 *     instance which produced the trigger primitives to be used as input;
 *     it must not include any instance name, as the instance names will be
 *     automatically added from `Thresholds` parameter.
 *     The typical trigger primitives used as input may be LVDS discriminated
 *     output (e.g. from `icarus::trigger::LVDSgates` module) or combinations
 *     of them (e.g. from `icarus::trigger::SlidingWindowTrigger` module).
 * * `Thresholds` (list of integers, mandatory): list of the discrimination
 *     thresholds to consider, in ADC counts. A data product containing a
 *     digital signal is read for each one of the thresholds, and the tag of the
 *     data product is expected to be the module label `TriggerGatesTag` with as
 *     instance name the value of the threshold (e.g. for a threshold of 60 ADC
 *     counts the data product tag might be `LVDSgates:60`).
 * * `GeneratorTags` (list of input tags, default: `[ generator ]`): a list of
 *     data products containing the particles generated by event generators;
 * * `DetectorParticleTag` (input tag, default: `largeant`): data product
 *     containing the list of particles going through the argon in the detector;
 * * `EnergyDepositTags`
 *     (list of input tags, default: `[ "largeant:TPCActive" ]`): a list of
 *     data products with energy depositions;
 * * `BeamGateDuration` (time, _mandatory_): the duration of the beam
 *     gate; _the time requires the unit to be explicitly specified_: use
 *     `"1.6 us"` for BNB, `9.5 us` for NuMI (also available as
 *     `BNB_settings.spill_duration` and `NuMI_settings.spill_duration` in
 *     `trigger_icarus.fcl`);
 * * `MinimumPrimitives` (list of integers, _mandatory_): a list of alternative
 *     requirements for the definition of a trigger; each value is the number
 *     of trigger primitives needed to be "on" at the same time for the trigger
 *     to fire;
 * * `TriggerTimeResolution` (time, default: `8 ns`): time resolution for the
 *     trigger primitives;
 * * `EventTreeName` (optional string): if specified, a simple ROOT tree is
 *     created with the information from each event (see `EventInfoTree` for
 *     its structure);
 * * `EventDetailsLogCategory` (optional string): if specified, information for
 *     each single event is output into the specified stream category; if
 *     the string is specified empty, the default module stream is used, as
 *     determined by `LogCategory` parameter;
 * * `LogCategory` (string, default `TriggerEfficiencyPlots`): name of category
 *     used to stream messages from this module into message facility.
 * 
 * An example job configuration is provided as `maketriggerplots_icarus.fcl`.
 * 
 * 
 * Technical description of the module
 * ====================================
 * 
 * This module reads
 * @ref TriggerEfficiencyPlots_Data "trigger gate data products" from different
 * ADC thresholds, and for each threahold it combines the data into a trigger
 * response depending on a ADC threshold.
 * It then applies different requirement levels and for each one it fills plots.
 *
 * All the input sets (each with its own ADC treshold) are treated independently
 * from each other.
 *
 * We can separate the plots in two categories: the ones which depend on the
 * trigger response, and the ones which do not.
 * The latter include event-wide information like the total energy deposited in
 * the detector by one event. They also include plots that depend on the ADC
 * threshold, e.g. the largest number of trigger primitives available at any
 * time during the event. Therefore there are more precisely three plot
 * categories:
 *
 * 1. plots which do not depend on the trigger primitives at all;
 * 2. plots which depend on the ADC threshold of the trigger primitives,
 *    but not on the trigger response;
 * 3. plots which depend on the trigger response (and therefore also on the
 *    ADC threshold of the trigger primitives).
 *
 * In addition, each event is assigned to event categories. These categories
 * do not depend on trigger primitives. For example, a electron neutrino neutral
 * current interaction event would be assigned to the neutral current neutrino
 * interaction category, to the electron neutrino category, and so on.
 * A set of plot is produced at once for each of the event categories.
 *
 * This module is designed as follows:
 *
 * * each _art_ event is treated independently from the others
 *   (this is the norm; method: `analyze()`);
 * * information for the plots which do not depend on the trigger primitives
 *   are extracted once per event (method: `extractEventInfo()`);
 * * the event is assigned its categories (method: `selectedPlotCategories()`,
 *   with the categories efined in `PlotCategories` global variable);
 * * for each input primitive set (i.e. for each ADC threshold):
 *     * trigger logic is applied (method: `combineTriggerPrimitives()`):
 *         * all trigger primitives are combined in a single nulti-level gate
 *         * the level of the gate is the count of trigger primitives on which
 *           the trigger response is based
 *     * the beam gate is imposed in coincidence (method: `applyBeamGate()`)
 *     * control is passed to method `plotResponses()` for plots, together with
 *       the list of all the categories the event belongs:
 *     * for each trigger requirement, that is a minimum number of
 *       primitives, the response to that requirement is plotted (as
 *       efficiency) and the time of the resulting trigger is added into an
 *       histogram
 *     * response-independent plots are also filled
 *     * the event is plotted for all the event categories which apply
 *
 * Note that the plots depending on the response may require multiple fillings.
 * For example, the trigger efficiency plot is a single plot in the plot set
 * which requires a filling for every requirement level.
 * In a similar fashion, the plot of trigger time requires one filling for
 * each requirement level that is met.
 *
 *
 * @anchor TriggerEfficiencyPlots_OrganizationOfPlots
 * Organization of the plots
 * --------------------------
 *
 * Plots are written on disk via the standard _art_ service `TFileService`,
 * which puts them in a ROOT subdirectory named as the instance of this module
 * being run.
 *
 * As described above, there are two "dimensions" for the plot sets: there is
 * one plot set for each ADC threshold and for each event category, the total
 * number of plot sets being the product of the options in the two dimensions.
 *
 * The code is structured to work by threshold by threshold, because the trigger
 * response depends on the threshold but not by the event category: the response
 * is computed once per threshold, then all the plots related to that response
 * (including all the event categories) are filled.
 *
 * The structure in the file reflects this logic, and there are two levels of
 * ROOT directories inside the one assigned by _art_ to this module:
 * * the outer level pins down the ADC threshold of the plots inside it;
 *   the name of the directory follows the pattern `Thr###` (### being the ADC
 *   threshold), which is the "threshold tag";
 * * the inner level is the category of events included in the plots, and the
 *   name of the directories reflect the name of the category as defined in
 *   the corresponding `PlotCategory` object (`PlotCategory::name()`); this
 *   defines the "category tag".
 *
 * In each inner directory, a complete set of plots is contained.
 * The name of each plot is a base name for that plot (e.g. `Eff` for the
 * efficiency plot) with the event category tag and the threshold tag appended
 * (separated by an underscore character, `"_"`). The title of the plot is also
 * modified to include the description of the plot category (as defined in
 * `PlotCategory::description()`) and the ADC threshold.
 * Therefore, within the same module directory, all plot names are different.
 *
 *
 * Adding a plot
 * --------------
 *
 * When adding a plot, two actions are necessary:
 *
 * 1. initialize the new plot
 * 2. fill the new plot
 *
 * The initialization happens in
 * `icarus::trigger::TriggerEfficiencyPlots::initializePlotSet()` method.
 * A request must be issued to the
 * @ref TriggerEfficiencyPlots_PlotSandboxes "plot sandbox" to "make" the plot.
 * In general it can be expected that all the arguments `args` in a call
 * `plots.make<PlotType>(args)` are forwarded to the constructor of `PlotType`,
 * so that to make a new `PlotType` object one has to consult only the
 * constructor of that class. Be aware though that the library expects the first
 * argument to be the ROOT name of the new object and the second to be its ROOT
 * title.
 * It is possible to work around this requirement with additional coding.
 *
 * The method performing the plotting is `plotResponses()`.
 * The plots should be filled inside loops going through all plot sets.
 * The existing code already does that in two loops, one for the plots
 * depending on the trigger response requirement, and the other for the plots
 * _not_ depending on it.
 *
 *
 * @anchor TriggerEfficiencyPlots_PlotSandboxes
 * ### About plot sandboxes
 *
 * For the sake of this module, a plot sandbox is an object similar to a ROOT
 * `TDirectory`, which can mangle the objects it contains to change ("process")
 * their name and title, and that interacts properly with `art::TFileDirectory`
 * so make `TFileService` happy.
 * The processing of the name and title is used to add the category and
 * threshold tags to the plot ROOT name and the corresponding annotations to
 * the plot title, as described
 * @ref TriggerEfficiencyPlots_OrganizationOfPlots "above".
 *
 *
 * @anchor TriggerEfficiencyPlots_AddingCategory
 * Adding an event category
 * -------------------------
 *
 * Event categories are listed in `PlotCategories`: each one is described by
 * a simple object `PlotCategory` which contains a name (used in ROOT
 * directories, plot name tags and to univocally identify the category),
 * a description (used in plot titles) and a condition for the event to fulfill
 * in order to qualify for the category. The condition is a function object
 * that accepts an `EventInfo_t` object with the relevant event information,
 * and returns its decision (`true` if the event belongs to this category).
 *
 * To add a category, it is enough to append an entry at the end of the
 * `PlotCategories` list, taking the existing ones as example.
 *
 * The condition function adjudicates based on the information in `EventInfo_t`.
 * It is likely though that if a new category is needed, the information
 * currently available in `EventInfo_t` is not enough to decide the response.
 * In that case, `EventInfo_t` must be extended to hold the newly required
 * information. In addition, the method `extractEventInfo()` must be modified
 * to set the new information. The documentation of that method expands on this
 * topic.
 *
 *
 * Changing the trigger logic
 * ---------------------------
 *
 * The trigger logic is hard-coded in `analyze()`, but it effectively has a
 * follow up in the plotting code.
 * In fact, in `analyze()` the primitives are combined into a single object,
 * but it is in `plotResponses()` that this object is interpreted to decide
 * whether the trigger fired according to a requirement.
 *
 * To implement a radically different trigger logic, several components need
 * to be coordinated:
 *
 * 1. configuration of trigger logic parameters
 * 2. combination of trigger primitives;
 * 3. the data structure the combination is stored in;
 * 4. interpretation of that data structure in terms of trigger firing;
 * 5. configuration of the trigger logic from user parameters.
 *
 * Because of the design of this module, the first and the third elements are
 * in different places (causing the necessity of the shared data structure that
 * is the second element) and changing one typically imply changing the other.
 *
 * This module simulates a trigger counting the number of primitives that are
 * "on" and imposing a required minimum on that number for the trigger to fire.
 * This is how _this_ implementation includes those elements:
 *
 * 1. combination of trigger primitives: primitives are summed into a single
 *    multilevel gate (with beam gate on top);
 * 2. the data structure the combination is stored in: a
 *    @ref TriggerEfficiencyPlots_Data "trigger gate object";
 * 3. interpretation of that data structure in terms of trigger firing:
 *    the logic requires a minimum number of primitives being on at the same
 *    time; the interpretation uses the level of the combined trigger gate as
 *    the number of primitives on at every time, and decides that the first
 *    time this number meets the requirement, the trigger fires;
 * 4. configuration of the trigger logic from user parameters:
 *    it's a parameter of the module itself (`MinimumPrimitives`).
 *
 * An example of how a different trigger logic _could_ be implemented following
 * a similar model. Let's assume we want a trigger based on two sliding windows.
 * A window is a group of LVDS signals from neighboring channels. We can split
 * the whole optical detector channels into many windows, and we can make the
 * windows overlap: for example, a window includes the first 30 channels behind
 * one anode plane, the second window includes 30 channels starting from the
 * fifteenth, the third window includes 30 channels starting from the thirtieth,
 * and so on. In ICARUS, 30 channels are served in 16 LVDS discriminated
 * waveforms, i.e. 16 trigger primitives. We may require as a loose trigger
 * that any window contains at least 10 primitives on at the same time (that is
 * about 20 PMT channels beyond ADC threshold), and as a tighter trigger that
 * there is a window with at least 8 trigger primitives on, _and_ the
 * corresponding window in the opposite TPC has at least 6 primitives on
 * at the same time. Note that a splitting 7/7 won't do: one of the two windows
 * must have 8 or more primitives on.
 * The input primitives from the point of view of the module now become the
 * "sliding window" combination of the LVDS trigger primitives, as produced e.g.
 * by `icarus::trigger::SlidingWindowTrigger` module.
 * How this can be implemented following the pattern of this module:
 *
 * 1. combination of trigger primitives: we prepare two primitives:
 *    one describes the maximum level among all the windows, and the other
 *    describes the level of the window opposite to the maximum (the way to
 *    achieve this is quite complicate);
 * 2. the data structure the combination is stored in: two
 *    @ref TriggerEfficiencyPlots_Data "trigger gate objects";
 * 3. interpretation of that data structure in terms of trigger firing:
 *    for the loose trigger, we look for the first time the maximum window level
 *    reaches the designated threshold; for the tighter trigger, some
 *    additional logic is required, discriminating the maximum window level
 *    primitive to the higher threshold, the opposite window level on the lower
 *    threshold, and setting them in coincidence; then the time the resulting
 *    gate opens, if any, is the trigger time;
 * 4. configuration of trigger logic parameters: the two requirements (single
 *    window, double window with their minimum number of LVDS primitives)
 *    are read from the module FHiCL configuration; for computing efficiency,
 *    either in the constructor or in a `beginJob()` method we may also prepare
 *    the associations (pairing) between opposite windows or just channels;
 *
 * @todo add ability to discriminate a trigger gate object?
 *       discrimination on level _N_ (>=N -> 1, \<N -> 0) can be obtained with
 *       adding a _-N_ uniform level, flooring on level 0 and (if needed) maxing
 *       on level 1.
 *
 * 
 * Code style: quantities with units
 * ----------------------------------
 * 
 * To avoid issues with missing or wrong conversions, this code often uses
 * LArSoft quantities library. A variable with a `Quantity` type is represented
 * with a specific unit (e.g. microseconds) and can undergo only limited
 * manipulation. The allowed manipulation should guarantee that the unit of
 * the quantity is correctly preserved. For example, it is not possible to
 * add to a `microseconds` interval a pure number (e.g. `9.0`), but rather
 * it is only possible to add time quantities (e.g. another `microseconds`
 * variable, or a `nanoseconds` variable, or a literal value with unit, like
 * `9.0_ns`), and those quantities are properly converted, with the caveat that
 * rounding errors may happen that a more careful explicit handling could avoid.
 * Also, there are two types of variables that can feature the same unit,
 * intervals and points. A point can't be scaled (e.g. you can't "double" the
 * time the beam gate opens, while you can double the beam gate _duration_)
 * and it can't be added to another point (the difference between two points
 * is an interval).
 * 
 * To avoid mistakes in the conversion between different time scales, the
 * LArSoft utility library `detinfo::DetectorTimings` is used, which is a
 * wrapper of `DetectorClocks` service provider that makes the interface to
 * convert times easier and uniform. Here we use it to convert time points
 * from the simulation (which are expressed in nanoseconds and are starting
 * at trigger time) into optical detector readout ticks, and vice versa.
 * The values returned by this library have encoded in them which time scale
 * they belong to, and in which unit they are measured.
 * 
 */
class icarus::trigger::TriggerEfficiencyPlots: public art::EDAnalyzer {

  using microseconds = util::quantities::intervals::microseconds;
  using nanoseconds = util::quantities::intervals::nanoseconds;

    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Sequence<art::InputTag> GeneratorTags {
      Name("GeneratorTags"),
      Comment("labels of the event generators"),
      std::vector<art::InputTag>{ "generator" }
      };

    fhicl::Atom<art::InputTag> DetectorParticleTag {
      Name("DetectorParticleTag"),
      Comment("label of simulated particles through the detector"),
      "largeant" // tradition demands
      };

    fhicl::Sequence<art::InputTag> EnergyDepositTags {
      Name("EnergyDeposits"),
      Comment("label of energy deposition data product(s) in the detector"),
      std::vector<art::InputTag>{ "largeant:TPCActive" }
      };

    fhicl::Atom<std::string> TriggerGatesTag {
      Name("TriggerGatesTag"),
      Comment("label of the input trigger gate data product (no instance name)")
      };

    fhicl::Sequence<raw::ADC_Count_t> Thresholds {
      Name("Thresholds"),
      Comment("thresholds to consider [ADC counts]")
      };

    fhicl::Atom<microseconds> BeamGateDuration {
      Name("BeamGateDuration"),
      Comment("length of time interval when optical triggers are accepted")
      };

    fhicl::Sequence<unsigned int> MinimumPrimitives {
      Name("MinimumPrimitives"),
      Comment("minimum required number of trigger primitives for a trigger")
      };
    
    fhicl::Atom<nanoseconds> TriggerTimeResolution {
      Name("TriggerTimeResolution"),
      Comment("resolution of trigger in time"),
      8_ns
      };
    
    fhicl::Atom<bool> PlotOnlyActiveVolume {
      Name("PlotOnlyActiveVolume"),
      Comment
        ("only events within TPC active volume are plot (if that makes sense)"),
      true
      };
    
    fhicl::OptionalAtom<std::string> EventTreeName {
      Name("EventTreeName"),
      Comment("name of a ROOT tree where to store event-by-event information")
      };

    fhicl::OptionalAtom<std::string> EventDetailsLogCategory {
      Name("EventDetailsLogCategory"),
      Comment("name of the category used for event information output")
      };

    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category used for the output"),
      "SlidingWindowTrigger" // default
      };

  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;
  // --- END Configuration -----------------------------------------------------


  // --- BEGIN Constructors ----------------------------------------------------
  explicit TriggerEfficiencyPlots(Parameters const& config);

  // Plugins should not be copied or assigned.
  TriggerEfficiencyPlots(TriggerEfficiencyPlots const&) = delete;
  TriggerEfficiencyPlots(TriggerEfficiencyPlots&&) = delete;
  TriggerEfficiencyPlots& operator=(TriggerEfficiencyPlots const&) = delete;
  TriggerEfficiencyPlots& operator=(TriggerEfficiencyPlots&&) = delete;

  // --- END Constructors ------------------------------------------------------


  // --- BEGIN Framework hooks -------------------------------------------------

  /// Fills the plots. Also extracts the information to fill them with.
  virtual void analyze(art::Event const& event) override;
  
  /// Prints end-of-job summaries.
  virtual void endJob() override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  /// Type of gate data without channel information.
  using TriggerGateData_t
    = icarus::trigger::OpticalTriggerGateData_t::GateData_t;
  
  using PlotSandbox = icarus::trigger::PlotSandbox; ///< Import type.
  
  /// List of references to plot sandboxes.
  using PlotSandboxRefs_t
    = std::vector<std::reference_wrapper<PlotSandbox const>>;
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  std::vector<art::InputTag> fGeneratorTags; ///< Generator data product tags.
  art::InputTag fDetectorParticleTag; ///< Input simulated particles.
  /// Energy deposition data product tags.
  std::vector<art::InputTag> fEnergyDepositTags;
  
  /// ADC thresholds to read, and the input tag connected to their data.
  std::map<icarus::trigger::ADCCounts_t, art::InputTag> fADCthresholds;
  
  /// Duration of the gate during with global optical triggers are accepted.
  microseconds fBeamGateDuration;
  
  /// Minimum number of trigger primitives for a trigger to happen.
  std::vector<unsigned int> fMinimumPrimitives;
  
  nanoseconds fTriggerTimeResolution; ///< Trigger resolution in time.
  
  bool fPlotOnlyActiveVolume; ///< Plot only events in active volume.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;

  std::string fLogEventDetails; ///< Steam where to print event info.
  
  // --- END Configuration variables -------------------------------------------


  // --- BEGIN Service variables -----------------------------------------------

  geo::GeometryCore const& fGeom;
  detinfo::DetectorClocks const& fDetClocks;
  detinfo::DetectorTimings fDetTimings;

  /// ROOT directory where all the plots are written.
  art::TFileDirectory fOutputDir;

  // --- END Service variables -------------------------------------------------

  
  // --- BEGIN Internal variables ----------------------------------------------
  
  /// Gate representing the time we expect light from beam interactions.
  icarus::trigger::OpticalTriggerGate const fBeamGate;
  
  /// Beam gate start and stop tick in optical detector scale.
  std::pair
    <detinfo::timescales::optical_tick, detinfo::timescales::optical_tick>
    const fBeamGateOpt;

  /// Beam gate start and stop time in simulation scale.
  std::pair
    <detinfo::timescales::simulation_time, detinfo::timescales::simulation_time>
    const fBeamGateSim;

  /// All plots, one set per ADC threshold.
  std::vector<PlotSandbox> fThresholdPlots;

  std::unique_ptr<EventIDTree> fIDTree; ///< Handler of ROOT tree output.
  std::unique_ptr<PlotInfoTree> fPlotTree; ///< Handler of ROOT tree output.
  std::unique_ptr<EventInfoTree> fEventTree; ///< Handler of ROOT tree output.
  std::unique_ptr<ResponseTree> fResponseTree; ///< Handler of ROOT tree output.
  
  std::atomic<unsigned int> nEvents { 0U }; ///< Count of seen events.
  std::atomic<unsigned int> nPlottedEvents { 0U }; ///< Count of plotted events.

  // --- END Internal variables ------------------------------------------------


  /// Initializes all the plot sets, one per threshold.
  void initializePlots(PlotCategories_t const& categories);

  /// Initializes full set of plots for (ADC threshold + category) into `plots`.
  void initializePlotSet(PlotSandbox& plots) const;

  /// Initializes set of plots per complete trigger definition into `plots`.
  void initializeEfficiencyPerTriggerPlots(PlotSandbox& plots) const;
  
  /// Initializes a single, trigger-independent plot set into `plots`.
  void initializeEventPlots(PlotSandbox& plots) const;
  
  
  /// Returns whether an event with the specified information should be included
  /// in the plots at all (it's a filter).
  bool shouldPlotEvent(EventInfo_t const& eventInfo) const;


  /// Returns the names of the plot categories event qualifies for.
  std::vector<std::string> selectPlotCategories
    (EventInfo_t const& info, PlotCategories_t const& categories) const;

  /**
   * @brief Extracts from `event` the relevant information on the physics event.
   *
   * This returns a `EventInfo_t` object filled as completely as possible.
   * The result is used in the context of
   * @ref TriggerEfficiencyPlots_AddingCategory "assigning the `event` to categories".
   *
   * This method can read any data product in the event, and has access to the
   * module configuration. If information is needed from a data product that
   * is not read yet, the following actions are needed:
   *
   * 1. decide the input tag of the data product(s) to be read; this is usually
   *    delegated to the user via module configuration (see `Config` structure);
   * 2. the required input tags need to be stored into `TriggerEfficiencyPlots`
   *    object;
   * 3. _art_ needs to be informed of the intention of reading the data products
   *    via `consumes()` or equivalent calls, in `TriggerEfficiencyPlots`
   *    constructor;
   * 4. the data product can be read with the usual means in this method.
   *
   * For an example of this procedure, see how the configuration parameter
   * `DetectorParticleTag` is treated: read by `Config::DetectorParticleTag`,
   * stored in `fDetectorParticleTag`, declared as "consumed" in
   * `TriggerEfficiencyPlots` constructor. Likewise a known list of data
   * products can be used: see `GeneratorTags` parameter, read by the sequence
   * `Config::GeneratorTags`, stored in the class field `fGeneratorTags`,
   * declared as "consumed" in a loop in `TriggerEfficiencyPlots` constructor,
   * and read in `extractEventInfo()` to be used within it.
   *
   */
  EventInfo_t extractEventInfo(art::Event const& event) const;
  
  /// Returns the TPC `point` is within, `nullptr` if none.
  geo::TPCGeo const* pointInTPC(geo::Point_t const& point) const;
  
  /// Returns in which TPC `point` is within the _active volume_ of;
  /// `nullptr` if none.
  geo::TPCGeo const* pointInActiveTPC(geo::Point_t const& point) const;
  
  /// Computes the trigger response from primitives with the given `threshold`.
  std::vector<TriggerGateData_t> combineTriggerPrimitives(
    art::Event const& event,
    icarus::trigger::ADCCounts_t const threshold,
    art::InputTag const& dataTag
    ) const;

  /// Applies the beam gate coincidence to the specified trigger primitive.
  template <typename GateObject>
  GateObject applyBeamGate(GateObject gate) const;

  /// Adds all the `responses` (one per threshold) to the plots.
  void plotResponses(
    std::size_t iThr, icarus::trigger::ADCCounts_t const threshold,
    PlotSandboxRefs_t const& plots, EventInfo_t const& info,
    std::vector<TriggerGateData_t> const& primitiveCounts
    ) const;

  /// Reads a set of input gates from the `event`.
  std::vector<std::vector<icarus::trigger::MultiChannelOpticalTriggerGate>> ReadTriggerGates
    (art::Event const& event, art::InputTag const& dataTag) const;
  
  static std::string thrAndCatName
    (std::string const& boxName, std::string const& category)
    { return boxName + "_" + category; }
  static std::string thrAndCatName
    (PlotSandbox const& box, std::string const& category)
    { return thrAndCatName(box.name(), category); }
  
  /// Returns a gate that is `Max()` of all the specified `gates`.
  template <typename TrigGateColl>
  static auto computeMaxGate(TrigGateColl const& gates);

  
}; // icarus::trigger::TriggerEfficiencyPlots


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
//--- icarus::trigger::TriggerEfficiencyPlots
//------------------------------------------------------------------------------
icarus::trigger::TriggerEfficiencyPlots::TriggerEfficiencyPlots
  (Parameters const& config)
  : art::EDAnalyzer(config)
  // configuration
  , fGeneratorTags        (config().GeneratorTags())
  , fDetectorParticleTag  (config().DetectorParticleTag())
  , fEnergyDepositTags    (config().EnergyDepositTags())
  , fBeamGateDuration     (config().BeamGateDuration())
  , fMinimumPrimitives    (config().MinimumPrimitives())
  , fTriggerTimeResolution(config().TriggerTimeResolution())
  , fPlotOnlyActiveVolume (config().PlotOnlyActiveVolume())
  , fLogCategory          (config().LogCategory())
  // services
  , fGeom      (*lar::providerFrom<geo::Geometry>())
  , fDetClocks (*lar::providerFrom<detinfo::DetectorClocksService>())
  , fDetTimings(fDetClocks)
  , fOutputDir (*art::ServiceHandle<art::TFileService>())
  // cached
  , fBeamGate(icarus::trigger::BeamGateMaker{ fDetClocks }(fBeamGateDuration))
  , fBeamGateOpt(
      fDetTimings.toOpticalTick(fDetTimings.BeamGateTime()),
      fDetTimings.toOpticalTick(fDetTimings.BeamGateTime() + fBeamGateDuration)
    )
  , fBeamGateSim(
      fDetTimings.toSimulationTime(fBeamGateOpt.first),
      fDetTimings.toSimulationTime(fDetTimings.BeamGateTime())
        + fBeamGateDuration
    )
{
  //
  // more complex parameter parsing
  //
  if (config().EventDetailsLogCategory(fLogEventDetails)) {
    // the parameter is somehow set, so fLogEventDetails won't be empty;
    // but if the parameter value is specified empty, we set it to fLogCategory
    if (fLogEventDetails.empty()) fLogEventDetails = fLogCategory;
  } // if EventDetailsLogCategory is specified

  std::string const discrModuleLabel = config().TriggerGatesTag();
  for (raw::ADC_Count_t threshold: config().Thresholds()) {
    fADCthresholds[icarus::trigger::ADCCounts_t{threshold}]
      = art::InputTag{ discrModuleLabel, util::to_string(threshold) };
  }

  if (config().EventTreeName.hasValue()) {
    std::string treeName;
    config().EventTreeName(treeName);

    fIDTree = std::make_unique<EventIDTree>
      (*(fOutputDir.make<TTree>(treeName.c_str(), "Event information")));
    fEventTree = std::make_unique<EventInfoTree>(fIDTree->tree());
    fPlotTree = std::make_unique<PlotInfoTree>(fIDTree->tree());
    fResponseTree = std::make_unique<ResponseTree>(
      fIDTree->tree(),
      util::get_elements<0U>(fADCthresholds), fMinimumPrimitives
      );

  } // if make tree

  //
  // input data declaration
  //
  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // event information
  for (art::InputTag const& inputTag: fGeneratorTags)
    consumes<std::vector<simb::MCTruth>>(inputTag);
  for (art::InputTag const& inputTag: fEnergyDepositTags)
    consumes<std::vector<sim::SimEnergyDeposit>>(inputTag);
//   consumes<std::vector<simb::MCParticle>>(fDetectorParticleTag);
  
  // trigger primitives
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    consumes<std::vector<OpticalTriggerGateData_t>>(inputDataTag);
    consumes<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (inputDataTag);
  } // for

  //
  // initialization of plots and plot categories
  //
  initializePlots(::PlotCategories);

  {
    mf::LogInfo log(fLogCategory);
    log << "\nConfigured " << fADCthresholds.size() << " thresholds:";
    for (auto const& [ threshold, dataTag ]: fADCthresholds)
      log << "\n * " << threshold << " ADC (from '" << dataTag.encode() << "')";
    log << "\nBeam gate is " << fBeamGate << " (" << fBeamGateSim.first
      << " -- " << fBeamGateSim.second << ")";
  } // local block

} // icarus::trigger::TriggerEfficiencyPlots::TriggerEfficiencyPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::analyze(art::Event const& event) {

  /*
   * 1. find out the features of the event and the categories it belongs to
   * 2. for each threshold:
   *   1. read the trigger primitives
   *   2. apply the beam gate on the primitives
   *   3. (optional) combine the trigger primitives
   *   4. generate the trigger response
   *   5. pick the plots to be filled
   *   6. add the response to all the plots
   *
   */
  
  ++nEvents;
  
  //
  // 1. find out the features of the event and the categories it belongs to
  //
  EventInfo_t const eventInfo = extractEventInfo(event);
  
  bool const bPlot = shouldPlotEvent(eventInfo);
  if (bPlot) ++nPlottedEvents;
  
  if (fIDTree) fIDTree->assignEvent(event);
  if (fPlotTree) fPlotTree->assign(bPlot);
  if (fEventTree) fEventTree->assignEvent(eventInfo);
  
  std::vector<std::string> selectedPlotCategories
    = selectPlotCategories(eventInfo, ::PlotCategories);
  {
    mf::LogTrace log(fLogCategory);
    log
      << "Event " << event.id() << " falls in " << selectedPlotCategories.size()
      << " categories:"
      ;
    for (std::string const& name: selectedPlotCategories)
      log << " \"" << name << "\"";
    // print the information on the event
  } // local block
  if (!fLogEventDetails.empty()) {
    mf::LogTrace(fLogEventDetails)
      << "Event " << event.id() << ": " << eventInfo;
  }

  //
  // 2. for each threshold:
  //
  for (auto&& [ iThr, thrPair, thrPlots ]
    : util::enumerate(fADCthresholds, fThresholdPlots)
  ) {

    auto const& [ threshold, dataTag ] = thrPair;

    //
    // 1.1. read the trigger primitives
    // 1.2. (optional) combine the trigger primitives
    // 1.3. apply the beam gate on the primitives
    // 1.4. generate the trigger response
    //

    std::vector<TriggerGateData_t> primitiveCounts;
    std::vector<TriggerGateData_t> combinedTriggerPrimitives = combineTriggerPrimitives(event, threshold, dataTag);

    for (TriggerGateData_t combinedPrimitive : (combinedTriggerPrimitives)) {
      TriggerGateData_t const primitiveCount
        = applyBeamGate(combinedPrimitive);

      primitiveCounts.push_back(primitiveCount);
    }
    
    //
    // 1.5. pick the plots to be filled
    //
    PlotSandboxRefs_t selectedPlots;
    
    if (bPlot) {
      for (std::string const& name: selectedPlotCategories)
        selectedPlots.emplace_back(*(thrPlots.findSandbox(name)));
    }
    
    // 1.6. add the response to the appropriate plots
    plotResponses(iThr, threshold, selectedPlots, eventInfo, primitiveCounts);
    
  } // for thresholds


  //
  // store information in output tree if any
  //
  if (fIDTree) fIDTree->tree().Fill();

} // icarus::trigger::TriggerEfficiencyPlots::analyze()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::endJob() {
  
  mf::LogInfo(fLogCategory)
    << nPlottedEvents << "/" << nEvents << " events plotted."
    ;
  
} // icarus::trigger::TriggerEfficiencyPlots::endJob()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::initializePlots
  (PlotCategories_t const& categories)
{
  using namespace std::string_literals;
  
  for (icarus::trigger::ADCCounts_t const threshold
    : util::get_elements<0U>(fADCthresholds))
  {
    // create a plot sandbox inside `fOutputDir` with a name/prefix `Thr###`
    auto const thr = threshold.value();
    icarus::trigger::PlotSandbox thrPlots { fOutputDir,
      "Thr"s + util::to_string(thr), "(thr: "s + util::to_string(thr) + ")"s };
    
    // create a subbox for each plot category
    for (PlotCategory const& category: categories) {
      PlotSandbox& plots = thrPlots.addSubSandbox(
        category.name(),
        category.description()
        );
      
      initializePlotSet(plots);
    }
    fThresholdPlots.push_back(std::move(thrPlots));
  } // for thresholds
  
  mf::LogTrace log(fLogCategory);
  log << "Created " << fThresholdPlots.size() << " plot boxes:\n";
  for (auto const& box: fThresholdPlots) {
    box.dump(log, "  ");
  } // for
  
} // icarus::trigger::TriggerEfficiencyPlots::initializePlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::initializePlotSet
  (PlotSandbox& plots) const
{
  
  // a variable binning for the required number of trigger primitives
  auto [ minimumPrimBinning, minimumPrimBinningLabels ]
    = util::ROOT::makeVariableBinningAndLabels(fMinimumPrimitives);
  assert(minimumPrimBinning.size() == minimumPrimBinningLabels.size() + 1U);

  {
    mf::LogTrace log(fLogCategory);
    log << "TriggerEfficiencyPlots (plots '"
      << plots.name() << "') variable binning including the "
      << fMinimumPrimitives.size() << " points {";
    for (auto value: fMinimumPrimitives) log << " " << value;
    log << " } => " << minimumPrimBinningLabels.size() << " bins: ";
    for (auto const& [ value, label ]
      : util::zip<1U>(minimumPrimBinning, minimumPrimBinningLabels))
    {
      log << " " << value << " (\"" << label << "\") =>";
    } // for
    log << " " << minimumPrimBinning.back();
  } // debug output block

  //
  // Triggering efficiency vs. threshold.
  //
  detinfo::timescales::optical_time_ticks const triggerResolutionTicks
    { fDetTimings.toOpticalTicks(fTriggerTimeResolution) };
  auto* TrigTime = plots.make<TH2F>(
    "TriggerTick",
    "Trigger time tick"
      ";minimum requested number of trigger primitives"
      ";optical time tick [ /" + util::to_string(triggerResolutionTicks) + " ]",
    minimumPrimBinning.size() - 1U, minimumPrimBinning.data(),
//    fMinimumPrimitives.back(), 0, fMinimumPrimitives.back() + 1
    (fBeamGateOpt.second - fBeamGateOpt.first) / triggerResolutionTicks,
    fBeamGateOpt.first.value(), fBeamGateOpt.second.value()
    );
  
  util::ROOT::applyAxisLabels(TrigTime->GetXaxis(), minimumPrimBinningLabels);
  
  auto* Eff = plots.make<TEfficiency>(
    "Eff",
    "Efficiency of triggering"
      ";minimum requested number of trigger primitives"
      ";trigger efficiency",
    minimumPrimBinning.size() - 1U, minimumPrimBinning.data()
//    fMinimumPrimitives.back(), 0, fMinimumPrimitives.back() + 1
    );
  
  // people are said to have earned hell for things like this;
  // but TEfficiency really does not expose the interface to assign labels to
  // its axes, which supposedly could be done had we chosen to create it by
  // histograms instead of directly as recommended.
  util::ROOT::applyAxisLabels(
    const_cast<TH1*>(Eff->GetTotalHistogram())->GetXaxis(),
    minimumPrimBinningLabels
    );
  
  //
  // plots independent of the trigger primitive requirements
  //
  plots.make<TH1F>(
    "NPrimitives",
    "Number of trigger primitives (\"channels firing at once\")"
    ";maximum trigger primitives at the same time"
    ";events",
    192, 0.0, 192.0 // large number, zoom in presentations!
    );
  
  //
  // Selection-related plots
  //
  initializeEventPlots(plots);
  

  //
  // Plots per trigger setting, split in triggering and not triggering events;
  // the plot set is the same as the "global" one.
  //
  using SS_t = std::pair<std::string, std::string>;
  std::array<SS_t, 2U> const classes {
    SS_t{ "triggering", "triggering events" },
    SS_t{ "nontriggering", "non-triggering events" }
    };
  for (auto minCount: fMinimumPrimitives) {
    
    std::string const minCountStr = std::to_string(minCount);
    
    // this defines a specific trigger, with its thresholds and settings
    PlotSandbox& reqBox = plots.addSubSandbox
      ("Req" + minCountStr, minCountStr + " channels required");
    
    initializeEfficiencyPerTriggerPlots(reqBox);
    
    for (auto const& [ name, desc ]: classes) {
      
      PlotSandbox& box = reqBox.addSubSandbox(name, desc);
      
      initializeEventPlots(box);
      
    } // for triggering requirement
  } // for triggering classes
  
  //
  // No trigger related plots
  //
  plots.make<TH1F>(
    "NeutrinoEnergy_NoTrig",
    "True Neutrino Energy of Non-Triggered Event"
      ";neutrino energy [GeV]"
      ";events",
    120,0.0,6.0 // GeV
  );
  plots.make<TH1F>(
    "EnergyInSpill_NoTrig",
    "Energy eposited during the beam gate opening of Non-Triggered Event"
      ";deposited energy  [ GeV ]"
      ";events  [ / 50 MeV ]",
    120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
    );
  plots.make<TH1I>(
    "InteractionType_NoTrig",
    "Interaction type of Non-Triggered Event"
      ";Interaction Type"
      ";events  [ / 50 MeV ]",
    200, 999.5, 1199.5 
    );
  plots.make<TH1F>(
    "LeptonEnergy_NoTrig",
    "Energy of outgoing lepton of Non-Triggered Event"
      ";lepton generated energy  [ GeV ]"
      ";events  [ / 50 MeV ]",
    120, 0.0, 6.0 
    );
 
} // icarus::trigger::TriggerEfficiencyPlots::initializePlotSet()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::initializeEfficiencyPerTriggerPlots
  (PlotSandbox& plots) const
{
  
  detinfo::timescales::optical_time_ticks const triggerResolutionTicks
    { fDetTimings.toOpticalTicks(fTriggerTimeResolution) };
  
  //
  // Triggering efficiency vs. something else
  //
  plots.make<TEfficiency>(
    "EffVsEnergyInSpill",
    "Efficiency of triggering vs. energy deposited in spill"
      ";energy deposited in spill  [ GeV ]"
      ";trigger efficiency  [ / 50 GeV ]",
    120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
    );
  
  plots.make<TEfficiency>(
    "EffVsEnergyInSpillActive",
    "Efficiency of triggering vs. energy deposited in active volume"
      ";energy deposited in active volume in spill  [ GeV ]"
      ";trigger efficiency  [ / 50 GeV ]",
    120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
    );
  
  plots.make<TEfficiency>(
    "EffVsNeutrinoEnergy",
    "Efficiency of triggering vs. neutrino energy"
      ";neutrino true energy  [ GeV ]"
      ";trigger efficiency  [ / 50 GeV ]",
    120, 0.0, 6.0 // 6 GeV is not that much for NuMI, but we should be ok
    );
  
  plots.make<TEfficiency>(
    "EffVsLeptonEnergy",
    "Efficiency of triggering vs. outgoing lepton energy"
      ";final state lepton true energy  [ GeV ]"
      ";trigger efficiency  [ / 50 GeV ]",
    120, 0.0, 6.0
    );
  
  plots.make<TH1F>(
    "TriggerTick",
    "Trigger time tick"
      ";optical time tick [ /" + util::to_string(triggerResolutionTicks) + " ]",
    (fBeamGateOpt.second - fBeamGateOpt.first) / triggerResolutionTicks,
    fBeamGateOpt.first.value(), fBeamGateOpt.second.value()
    );
  
} // initializeEfficiencyPerTriggerPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::initializeEventPlots
  (PlotSandbox& plots) const
{
  
  // a variable binning for the required number of trigger primitives
  auto [ minimumPrimBinning, minimumPrimBinningLabels ]
    = util::ROOT::makeVariableBinningAndLabels(fMinimumPrimitives);
  assert(minimumPrimBinning.size() == minimumPrimBinningLabels.size() + 1U);

  {
    mf::LogTrace log(fLogCategory);
    log << "TriggerEfficiencyPlots (plots '"
      << plots.name() << "') variable binning including the "
      << fMinimumPrimitives.size() << " points {";
    for (auto value: fMinimumPrimitives) log << " " << value;
    log << " } => " << minimumPrimBinningLabels.size() << " bins: ";
    for (auto const& [ value, label ]
      : util::zip<1U>(minimumPrimBinning, minimumPrimBinningLabels))
    {
      log << " " << value << " (\"" << label << "\") =>";
    } // for
    log << " " << minimumPrimBinning.back();
  } // debug output block

  //
  // Selection-related plots
  //
  plots.make<TH1F>(
    "NeutrinoEnergy",
    "True Neutrino Energy"
      ";neutrino energy [GeV]"
      ";events",
    120, 0.0, 6.0 // GeV
  );
  plots.make<TH1F>(
    "EnergyInSpill",
    "Energy deposited during the beam gate opening"
      ";energy deposited in spill [ GeV ]"
      ";events  [ / 50 MeV ]",
    120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
    );
  plots.make<TH1F>(
    "EnergyInSpillActive",
    "Energy deposited during the beam gate opening in active volume"
      ";energy deposited in active volume in spill [ GeV ]"
      ";events  [ / 50 MeV ]",
    120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
    );
  plots.make<TH1I>(
    "InteractionType",
    "Interaction type"
      ";Interaction Type"
      ";events",
    200, 999.5, 1199.5
    );
  plots.make<TH1F>(
    "LeptonEnergy",
    "Energy of outgoing lepton"
      ";deposited energy  [ GeV ]"
      ";events  [ / 50 MeV ]",
    120, 0.0, 6.0
    );
  plots.make<TH2F>(
    "InteractionVertexYZ",
    "Vertex of triggered interaction"
      ";Z [ / 20 cm ]"
      ";Y [ / 5 cm ]",
    120, -1200., +1200.,
    100,  -250.,  +250.
    );
  
} // icarus::trigger::TriggerEfficiencyPlots::initializeEventPlots()


//------------------------------------------------------------------------------
bool icarus::trigger::TriggerEfficiencyPlots::shouldPlotEvent
  (EventInfo_t const& eventInfo) const
{
  if (fPlotOnlyActiveVolume
    && eventInfo.hasVertex() && !eventInfo.isInActiveVolume())
  {
    return false;
  }
  
  return true;
} // icarus::trigger::TriggerEfficiencyPlots::shouldPlotEvent()


//------------------------------------------------------------------------------
std::vector<std::string>
icarus::trigger::TriggerEfficiencyPlots::selectPlotCategories
  (EventInfo_t const& info, PlotCategories_t const& categories) const
{
  std::vector<std::string> selected;
  
  for (auto const& category: categories)
    if (category(info)) selected.push_back(category);
  
  return selected;
  
} // icarus::trigger::TriggerEfficiencyPlots::selectPlotCategories()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlots::extractEventInfo
  (art::Event const& event) const -> EventInfo_t
{
  
  EventInfo_t info;
  
  //
  // generator information
  //
  for (art::InputTag const& inputTag: fGeneratorTags) {
  
    auto const& truthRecords
      = *(event.getValidHandle<std::vector<simb::MCTruth>>(inputTag));
    
    for (simb::MCTruth const& truth: truthRecords) {
      
      if (truth.NeutrinoSet()) {
        //
        // interaction flavor (nu_mu, nu_e)
        // interaction type (CC, NC)
        //

        simb::MCParticle const& nu = truth.GetNeutrino().Nu();
        info.SetNeutrinoPDG(nu.PdgCode());
        info.SetInteractionType(truth.GetNeutrino().InteractionType());

        info.SetNeutrinoEnergy(GeV{ nu.E() });
        info.SetLeptonEnergy(GeV{ truth.GetNeutrino().Lepton().E() });
        //info.SetNucleonEnergy(truth.GetNeutrino().HitNuc().E());

        switch (nu.PdgCode()) {
          case 14:
          case -14:
            info.SetNu_mu(true);
            break;
          case 12:
          case -12:
            info.SetNu_e(true);
            break;
        }

        switch (truth.GetNeutrino().CCNC()) {
          case simb::kCC: info.AddWeakChargedCurrentInteractions(); break;
          case simb::kNC: info.AddWeakNeutralCurrentInteractions(); break;
          default:
            mf::LogWarning(fLogCategory)
              << "Event " << event.id() << " has unexpected NC/CC flag ("
              << truth.GetNeutrino().CCNC() << ")";
        } // switch   
        
        // we do not trust the vertex (`GvX()`) of the neutrino particle,
        // since GenieHelper does not translate the vertex
        // of some of the particles from GENIE to detector frame;
        // trajectory is always translated:
        geo::Point_t const vertex { nu.EndX(), nu.EndY(), nu.EndZ() };
        info.AddVertex(vertex);
        
        geo::TPCGeo const* tpc = pointInActiveTPC(vertex);
        if (tpc) info.SetInActiveVolume();
        
      } // if neutrino event
      
    } // for truth records
    
  } // for generators
  
  //
  // propagation in the detector
  //
  GeV totalEnergy { 0.0 }, inSpillEnergy { 0.0 };
  GeV activeEnergy { 0.0 }, inSpillActiveEnergy { 0.0 };
  
  for (art::InputTag const& edepTag: fEnergyDepositTags) {
    
    auto const& energyDeposits
      = *(event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(edepTag));
    mf::LogTrace(fLogCategory)
      << "Event " << event.id() << " has " << energyDeposits.size()
      << " energy deposits recorded in '" << edepTag.encode() << "'";
    
    for (sim::SimEnergyDeposit const& edep: energyDeposits) {
      
      MeV const e { edep.Energy() }; // assuming it's stored in MeV
      
      detinfo::timescales::simulation_time const t { edep.Time() };
      bool const inSpill
        = (t >= fBeamGateSim.first) && (t <= fBeamGateSim.second);
      
      totalEnergy += e;
      if (inSpill) inSpillEnergy += e;
      
      if (pointInActiveTPC(edep.MidPoint())) {
        activeEnergy += e;
        if (inSpill) inSpillActiveEnergy += e;
      }
      
    } // for all energy deposits in the data product
    
  } // for all energy deposit tags
  
  info.SetDepositedEnergy(totalEnergy);
  info.SetDepositedEnergyInSpill(inSpillEnergy);
  info.SetDepositedEnergyInActiveVolume(activeEnergy);
  info.SetDepositedEnergyInActiveVolumeInSpill(inSpillActiveEnergy);
  
  mf::LogTrace(fLogCategory) << "Event " << event.id() << ": " << info;
  
  return info;
} // icarus::trigger::TriggerEfficiencyPlots::extractEventInfo()


//------------------------------------------------------------------------------
geo::TPCGeo const* icarus::trigger::TriggerEfficiencyPlots::pointInTPC
  (geo::Point_t const& point) const
{
  return fGeom.PositionToTPCptr(point);
} // icarus::trigger::TriggerEfficiencyPlots::pointInTPC()


//------------------------------------------------------------------------------
geo::TPCGeo const* icarus::trigger::TriggerEfficiencyPlots::pointInActiveTPC
  (geo::Point_t const& point) const
{
  geo::TPCGeo const* tpc = pointInTPC(point);
  return
    (tpc && tpc->ActiveBoundingBox().ContainsPosition(point))? tpc: nullptr;
} // icarus::trigger::TriggerEfficiencyPlots::pointInActiveTPC()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlots::combineTriggerPrimitives(
  art::Event const& event,
  icarus::trigger::ADCCounts_t const threshold,
  art::InputTag const& dataTag
) const -> std::vector<TriggerGateData_t>{

  //
  // simple count
  //

  std::vector<std::vector<icarus::trigger::MultiChannelOpticalTriggerGate>> const& cryoGates
    = ReadTriggerGates(event, dataTag);

  std::vector<TriggerGateData_t> cryoCombinedGate;

  for (auto const& gates: (cryoGates)) {
    mf::LogTrace(fLogCategory)
      << "Simulating trigger response with ADC threshold " << threshold
      << " from '" << dataTag.encode() << "' (" << gates.size() << " primitives)";

    TriggerGateData_t combinedGate;

    auto itGate = gates.begin();
    auto const gend = gates.end();
    if (itGate == gend) { // this is unexpected...
      mf::LogWarning(fLogCategory)
        << "  No trigger primitive found for threshold " << threshold
        << " ('" << dataTag.encode() << "')";
      return {};
    } // if no gates

    combinedGate = *itGate;
    while (++itGate != gend) combinedGate.Sum(*itGate);

    cryoCombinedGate.push_back(combinedGate) ;
  }

  return cryoCombinedGate;
} // icarus::trigger::TriggerEfficiencyPlots::combineTriggerPrimitives()


//------------------------------------------------------------------------------
template <typename GateObject>
GateObject icarus::trigger::TriggerEfficiencyPlots::applyBeamGate
  (GateObject gate) const
{
  return std::move(gate.Mul(fBeamGate));
} // icarus::trigger::TriggerEfficiencyPlots::applyBeamGate()


//------------------------------------------------------------------------------
// out-of-order definition, needs to be before `plotResponses()`
template <typename TrigGateColl>
auto icarus::trigger::TriggerEfficiencyPlots::computeMaxGate
  (TrigGateColl const& gates)
{
  
  // if `gates` is empty return a default-constructed gate of the contained type
  if (empty(gates)) return decltype(*begin(gates)){};
  
  auto iGate = cbegin(gates);
  auto const gend = cend(gates);
  auto maxGate = *iGate;
  while (++iGate != gend) maxGate.Max(*iGate);
  
  return maxGate;
} // icarus::trigger::TriggerEfficiencyPlots::computeMaxGate()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlots::plotResponses(
  std::size_t iThr,
  icarus::trigger::ADCCounts_t const threshold,
  PlotSandboxRefs_t const& plotSets,
  EventInfo_t const& eventInfo,
  std::vector<TriggerGateData_t> const& primitiveCounts
) const {
  
  /*
   * This function plots according to the configured minimum number of trigger
   * primitives: for each requirement of minimum number of primitives, the
   * earliest time where that requirement is met is found, and that is
   * considered as the trigger time.
   * 
   * The following quantities are drawn per ADC threshold and per plot category:
   * 
   * * per minimum number of primitives:
   *    * trigger efficiency (i.e. whether there was _any time_ a number of
   *      primitives fulfilling the requirement)
   *    * trigger time in ticks (distribution as 2D histogram)
   * * maximum number of trigger primitives present at any time
   * * deposited energy during beam spill
   * 
   */
  using namespace std::string_literals;
  
  using ClockTick_t = TriggerGateData_t::ClockTick_t;
  using OpeningCount_t = TriggerGateData_t::OpeningCount_t;
  
  using PrimitiveCount_t = std::pair<ClockTick_t, OpeningCount_t>;
  
  // largest number of trigger primitives at any time
  assert(!primitiveCounts.empty());
  auto const primitiveCount = computeMaxGate(primitiveCounts);
  
  auto const maxPrimitiveTime { primitiveCount.findMaxOpen() };
  PrimitiveCount_t const maxPrimitives
    { maxPrimitiveTime, primitiveCount.openingCount(maxPrimitiveTime) };

  mf::LogTrace(fLogCategory)
    << "Max primitive count in " << threshold << ": "
    << maxPrimitives.second << " at tick " << maxPrimitives.first << " ("
    << fDetTimings.toElectronicsTime
      (detinfo::DetectorTimings::optical_tick{ maxPrimitives.first })
    << ")"
    ;
  
  /*
   * Fill all the histograms for all the minimum primitive requirements
   * (filling the information whether or not the trigger fired),
   * for all the qualifying categories.
   * Note that in this type of plots each event appears in all bins
   * (may be with "fired" or "not fired" on each bin)
   */
  PrimitiveCount_t lastMinCount { TriggerGateData_t::MinTick, 0 };
  bool fired = true; // the final trigger response (changes with requirement)
  
  class HistGetter { // helper, since this seems "popular"
    PlotSandbox const& plots;
    
      public:
    HistGetter(PlotSandbox const& plots): plots(plots) {}
    
    PlotSandbox const& box() const { return plots; }
    
    TH1& Hist(std::string const& name) const { return plots.demand<TH1>(name); }
    TH2& Hist2D(std::string const& name) const { return plots.demand<TH2>(name); }
    TEfficiency& Eff(std::string const& name) const
      { return plots.demand<TEfficiency>(name); }
    
  }; // class HistGetter
  
  
  for (auto [ iReq, minCount ]: util::enumerate(fMinimumPrimitives)) {
    
    if (fired && (lastMinCount.second < minCount)) {
      // if we haven't passed this minimum yet
      ClockTick_t const time = primitiveCount.findOpen(minCount);
      if (time == TriggerGateData_t::MaxTick) {
        mf::LogTrace(fLogCategory)
          << "Never got at " << minCount << " primitives or above.";
        fired = false;
      }
      else {
        lastMinCount = { time, primitiveCount.openingCount(time) };
        mf::LogTrace(fLogCategory)
          << "Reached " << minCount << " primitives or above ("
          << lastMinCount.second << ") at " << lastMinCount.first << ".";
      }
    } // if
    
    // at this point we know we have minCount or more trigger primitives,
    // and the time of this one is in lastMinCount.first (just in case)
    
    if (fResponseTree) fResponseTree->assignResponse(iThr, iReq, fired);
    
    std::string const minCountStr { "Req" + std::to_string(minCount) };
    
    // go through all the plot categories this event qualifies for
    for (icarus::trigger::PlotSandbox const& plotSet: plotSets) {
      
      HistGetter const get { plotSet };
      
      // simple efficiency
      get.Eff("Eff"s).Fill(fired, minCount);
      
      // trigger time (if any)
      if (fired) {
        get.Hist2D("TriggerTick"s).Fill(minCount, lastMinCount.first);
        if (minCount == 1) {
          for (auto point : eventInfo.fVertices) {
            get.Hist2D("InteractionVertexYZ"s).Fill(point.Z(), point.Y());
          }
        }
      }
      
      //
      // plotting specific to this trigger definition
      //
      
      HistGetter const getTrigEff { plotSet.demandSandbox(minCountStr) };
      
      // efficiency plots
      getTrigEff.Eff("EffVsEnergyInSpill").Fill
        (fired, double(eventInfo.DepositedEnergyInSpill()));
      getTrigEff.Eff("EffVsEnergyInSpillActive").Fill
        (fired, double(eventInfo.DepositedEnergyInSpillInActiveVolume()));
      getTrigEff.Eff("EffVsNeutrinoEnergy").Fill
        (fired, double(eventInfo.NeutrinoEnergy()));
      getTrigEff.Eff("EffVsLeptonEnergy").Fill
        (fired, double(eventInfo.LeptonEnergy()));
      if (fired) {
        getTrigEff.Hist("TriggerTick"s).Fill(lastMinCount.first);
      }

      
      //
      // plotting split for triggering/not triggering events
      //
      
      HistGetter const getTrig
        { getTrigEff.box().demandSandbox(fired? "triggering": "nontriggering") };
      
      getTrig.Hist("EnergyInSpill"s).Fill(double(eventInfo.DepositedEnergyInSpill()));
      getTrig.Hist("EnergyInSpillActive"s).Fill(double(eventInfo.DepositedEnergyInSpillInActiveVolume()));
      if (eventInfo.isNeutrino()) {
        getTrig.Hist("NeutrinoEnergy"s).Fill(double(eventInfo.NeutrinoEnergy()));
        getTrig.Hist("InteractionType"s).Fill(eventInfo.InteractionType());
        getTrig.Hist("LeptonEnergy"s).Fill(double(eventInfo.LeptonEnergy()));
      } // if neutrino event
      TH2& vertexHist = getTrig.Hist2D("InteractionVertexYZ"s);
      for (auto const& point: eventInfo.fVertices)
        vertexHist.Fill(point.Z(), point.Y());
      

      //
      // non triggered events
      //
      if (!fired && minCount == 1 ) { // I only am interested in events that aren't triggered when there is a low multiplicity requirement
        get.Hist("EnergyInSpill_NoTrig"s).Fill(double(eventInfo.DepositedEnergyInSpill()));
        get.Hist("NeutrinoEnergy_NoTrig"s).Fill(double(eventInfo.NeutrinoEnergy()));
        get.Hist("InteractionType_NoTrig"s).Fill(eventInfo.InteractionType());
        get.Hist("LeptonEnergy_NoTrig"s).Fill(double(eventInfo.LeptonEnergy()));
        //getHist("NucleonEnergy_NoTrig"s)->Fill(double(eventInfo.NucleonEnergy())); 
      }

    } // for all qualifying plot categories
    
  } // for all thresholds
  
  /*
   * Now fill the plots independent of the trigger response:
   * the same value is plotted in all plot sets.
   * 
   */
  for (PlotSandbox const& plotSet: plotSets) {
    
    HistGetter const get(plotSet);
    
    // number of primitives
    get.Hist("NPrimitives"s).Fill(maxPrimitives.second);
    
    // selection-related plots:
    get.Hist("EnergyInSpill"s).Fill(double(eventInfo.DepositedEnergyInSpill()));
    get.Hist("EnergyInSpillActive"s).Fill(double(eventInfo.DepositedEnergyInSpillInActiveVolume()));
    
    if (eventInfo.isNeutrino()) {
      get.Hist("NeutrinoEnergy"s).Fill(double(eventInfo.NeutrinoEnergy()));
      get.Hist("InteractionType"s).Fill(eventInfo.InteractionType());
      get.Hist("LeptonEnergy"s).Fill(double(eventInfo.LeptonEnergy()));
    } // if neutrino event
    TH2& vertexHist = get.Hist2D("InteractionVertexYZ"s);
    for (auto const& point: eventInfo.fVertices)
      vertexHist.Fill(point.Z(), point.Y());
    
  } // for 

} // icarus::trigger::TriggerEfficiencyPlots::plotResponses()


//------------------------------------------------------------------------------
std::vector<std::vector<icarus::trigger::MultiChannelOpticalTriggerGate>>
icarus::trigger::TriggerEfficiencyPlots::ReadTriggerGates
  (art::Event const& event, art::InputTag const& dataTag) const
{

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // currently the associations are a waste (we include it not to have
  // potentially broken `MultiChannelOpticalTriggerGate` objects, but at the
  // time this module was written that was not important)
  auto const& gates
    = *(event.getValidHandle<std::vector<OpticalTriggerGateData_t>>(dataTag));
  auto const& gateToWaveforms = *(
    event.getValidHandle
      <art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>(dataTag)
    );
  try {
    auto check = icarus::trigger::FillTriggerGates<icarus::trigger::MultiChannelOpticalTriggerGate>
    (gates, gateToWaveforms);
  }
  catch (cet::exception const& e) {
    throw cet::exception("TriggerEfficiencyPlots", "", e)
      << "Error encountered while reading data products from '"
      << dataTag.encode() << "'\n";
  }

  auto allGates = icarus::trigger::FillTriggerGates<icarus::trigger::MultiChannelOpticalTriggerGate>
  (gates, gateToWaveforms);

  std::vector<geo::CryostatID> fChannelCryostat;

  fChannelCryostat.reserve(fGeom.NOpChannels());
  for (auto const opChannel: util::counter(fGeom.NOpChannels())) {
    if (!fGeom.IsValidOpChannel(opChannel)) continue;
    fChannelCryostat[opChannel] = fGeom.OpDetGeoFromOpChannel(opChannel).ID();
  } // for all channels
  
  std::vector<std::vector<icarus::trigger::MultiChannelOpticalTriggerGate>> gatesPerCryostat{ fGeom.Ncryostats() };
  for (auto& gate: allGates) {
    assert(gate.hasChannels());
    gatesPerCryostat[fChannelCryostat[gate.channels().front()].Cryostat].push_back(std::move(gate));
  }
  return gatesPerCryostat;

} // icarus::trigger::TriggerEfficiencyPlots::ReadTriggerGates()


//------------------------------------------------------------------------------
//--- EventIDTree
//------------------------------------------------------------------------------
EventIDTree::EventIDTree(TTree& tree): TreeHolder(tree) {

  this->tree().Branch("RunNo", &fRunNo);
  this->tree().Branch("SubRunNo", &fSubRunNo);
  this->tree().Branch("EventNo", &fEventNo);

} // EventIDTree::EventIDTree()


//------------------------------------------------------------------------------
void EventIDTree::assignID(art::EventID const& id) {

  fRunNo = id.run();
  fSubRunNo = id.subRun();
  fEventNo = id.event();

} // EventIDTree::assignID()


//------------------------------------------------------------------------------
//--- PlotInfoTree
//------------------------------------------------------------------------------
PlotInfoTree::PlotInfoTree(TTree& tree): TreeHolder(tree) {

  this->tree().Branch("InPlots", &fInPlots);

} // PlotInfoTree::PlotInfoTree()


//------------------------------------------------------------------------------
void PlotInfoTree::assign(bool inPlots) {

  fInPlots = static_cast<Bool_t>(inPlots);

} // PlotInfoTree::assignEvent()


//------------------------------------------------------------------------------
//--- EventInfoTree
//------------------------------------------------------------------------------
EventInfoTree::EventInfoTree(TTree& tree): TreeHolder(tree) {

  this->tree().Branch("CC",       &fCC);
  this->tree().Branch("NC",       &fNC);
  this->tree().Branch("IntType",  &fIntType);
  this->tree().Branch("NuE",      &fNuE);
  this->tree().Branch("OutLeptE", &fOutLeptE);
  
  this->tree().Branch("TotE",         &fTotE);
  this->tree().Branch("SpillE",       &fSpillE);
  this->tree().Branch("ActiveE",      &fActiveE);
  this->tree().Branch("SpillActiveE", &fSpillActiveE);
  this->tree().Branch("InActive",     &fInActive);
  
} // EventInfoTree::EventInfoTree()


//------------------------------------------------------------------------------
void EventInfoTree::assignEvent(EventInfo_t const& info) {

  fCC       = info.nWeakChargedCurrentInteractions();
  fNC       = info.nWeakNeutralCurrentInteractions();
  fIntType  = info.InteractionType();
  fNuE      = static_cast<Double_t>(info.NeutrinoEnergy());
  fOutLeptE = static_cast<Double_t>(info.LeptonEnergy());

  fTotE         = static_cast<Double_t>(info.DepositedEnergy());
  fSpillE       = static_cast<Double_t>(info.DepositedEnergyInSpill());
  fActiveE      = static_cast<Double_t>(info.DepositedEnergyInActiveVolume());
  fSpillActiveE = static_cast<Double_t>(info.DepositedEnergyInSpillInActiveVolume());
  fInActive     = static_cast<Bool_t>(info.isInActiveVolume());
  
} // EventInfoTree::assignEvent()


//------------------------------------------------------------------------------
//--- ResponseTree
//------------------------------------------------------------------------------
template <typename Thresholds, typename Requirements>
ResponseTree::ResponseTree
  (TTree& tree, Thresholds const& thresholds, Requirements const& minReqs)
  : TreeHolder(tree)
  , indices(std::size(thresholds), std::size(minReqs))
  , RespTxxRxx{ std::make_unique<bool[]>(indices.size()) }
{

  for (auto [ iThr, threshold]: util::enumerate(thresholds)) {
    std::string const thrStr = util::to_string(raw::ADC_Count_t(threshold));

    for (auto [ iReq, req ]: util::enumerate(minReqs)) {

      std::string const branchName
        = "RespT" + thrStr + "R" + util::to_string(req);

      this->tree().Branch
        (branchName.c_str(), &(RespTxxRxx[indices(iThr, iReq)]));

    } // for all requirements

  } // for all thresholds

} // ResponseTree::ResponseTree()


//------------------------------------------------------------------------------
void ResponseTree::assignResponse
  (std::size_t iThr, std::size_t iReq, bool resp)
{
  RespTxxRxx[indices(iThr, iReq)] = resp;
} // ResponseTree::assignResponse()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::TriggerEfficiencyPlots)


//------------------------------------------------------------------------------
