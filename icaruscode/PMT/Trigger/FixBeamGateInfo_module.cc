/**
 * @file   FixBeamGateInfo_module.cc
 * @brief  Rewrites a collection of sim::BeamGateInfo.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 9, 2022
 */

// LArSoft libraries
#include "lardataalg/Utilities/intervals_fhicl.h" // microseconds from FHiCL
#include "lardataalg/Utilities/quantities/spacetime.h" // microseconds
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // simulation_time
#include "lardataalg/Utilities/MultipleChoiceSelection.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

// framework libraries
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <optional>
#include <memory> // std::make_unique()


//------------------------------------------------------------------------------
namespace icarus::trigger { class FixBeamGateInfo; }

/**
 * @brief Rewrites a set collection of beam gates into each event.
 * 
 * This module allows to perform fixed transformations to an existing beam
 * gate data product (`sim::BeamGateInfo`), producing a new set of gates.
 * 
 * 
 * Output data products
 * =====================
 *
 * * `std::vector<sim::BeamGateInfo>`: a collection of as many
 *    `sim::BeamGateInfo` as in the input.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse online description of the parameters is printed by running
 * `lar --print-description FixBeamGateInfo`.
 * 
 * * `BeamGateTag` (input tag, mandatory): data product with the beam gates to
 *     be "fixed"
 * * `Changes` (list of tables): a list of beam gate changes to be performed;
 *     the changes are performed in the order of this list, all the ones which
 *     apply (e.g. not stopping at the first that applies). Each item must
 *     specify at least an action to be taken (changing either start or width
 *     of the gate). Each element in the list is a table including:
 *     * `Select` (optional table): if specified, the following changes apply
 *       only to the gates that match the constraints as follows.
 *         * `Types` (list of gate types): if specified, the changes will be
 *           applied only to gates of type specified in the list; the gate type
 *           is specified by name (e.g. `"BNB"`, `"NuMI"`). If not specified,
 *           the changes are applied to all types of gates.
 *     * `Width` (optional table): a table with the operations affecting the
 *         duration of the gate. The end of the gate is modified accordingly,
 *         while the start point is not changed. The following options are
 *         exclusive:
 *         * `SetTo` (time): duration is set to this interval;
 *         * `Change` (time): this value is added to the duration (negative
 *             decreases it, but never below `0`).
 *     * `Start` (optional table): a table with the operations affecting the
 *         start of the gate. The duration of the gate is always preserved.
 *         The following options are exclusive:
 *         * `SetTo` (time): start of the beam gate is moved to this time point,
 *           set on the electronics time scale;
 *         * `Change` (time): this value is added to the start time (negative
 *           anticipates it).
 * * `KeepInstanceName` (flag, default: `false`): if set, the output beam gate
 *     data product will have the same instance name as the input one specified
 *     in `BeamGateTag`.
 * * `OutputInstanceName` (string, default: empty): the instance name of the
 *     output data product (exclusive with `KeepInstanceName`).
 * * `LogCategory` (string, default: `FixBeamGateInfo`): name of the output
 *     stream category for console messages (managed by MessageFacility
 *     library).
 *
 * Times must be specified as strings with their unit, e.g. `"0.6 us"` for 0.6
 * microseconds.
 */
class icarus::trigger::FixBeamGateInfo: public art::SharedProducer {
  
  using simulation_time = detinfo::timescales::simulation_time; // alias
  
  /// All directions to change a beam gate.
  struct BeamChangeRecipe {
    
    struct GateSelector_t {
      
      std::vector<sim::BeamType_t> beamTypes; ///< Match these beam types.
      
      bool empty() const { return beamTypes.empty(); }
      bool isValid() const { return true; }
    }; // GateSelector_t
    
    /// Set of instructions for a change.
    template <typename P, typename I = typename P::interval_t>
    struct ChangeRecipe_t {
      std::optional<P> setValue; ///< Value to set.
      std::optional<I> addValue; ///< Value to add.
      
      bool empty() const { return !setValue && !addValue; }
      bool valid() const { return !(setValue && addValue); }
    }; // ChangeRecipe_t
    
    GateSelector_t selectGates; ///< Which gates to apply this recipe on.
    
    /// Instructions on how to change the gate start.
    std::optional<ChangeRecipe_t<simulation_time>> const start;
    
    /// Instructions on how to change the gate width.
    std::optional<ChangeRecipe_t<simulation_time::interval_t>> const width;
    
    bool empty() const
      { return (!width || width->empty()) && (!start || start->empty()); }
    bool valid() const
      {
        // validity checks
        if (!selectGates.isValid()) return false;
        if (start && !start->valid()) return false;
        if (width && !width->valid()) return false;
        // minimum requirement checks
        if (!start && !width) return false;
        return true;
      }
  }; // struct BeamChangeRecipe
  
  
    public:
  
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    enum class BeamType_t { // we need to translate enum into a strong type
        kUnknown = sim::kUnknown
      , kBNB     = sim::kBNB
      , kNuMI    = sim::kNuMI
    };
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    using microsecond = util::quantities::points::microsecond;
    using microseconds = util::quantities::intervals::microseconds;
    
    /// Settings to change a gate.
    struct ChangeGate {
      
      template <typename P, typename I = typename P::interval_t>
      struct ChangeConfig {
        fhicl::OptionalAtom<P> SetTo {
          Name{ "SetTo" },
          Comment{ "set to the specified time, regardless of the previous one" }
          };
        fhicl::OptionalAtom<I> Add {
          Name{ "Add" },
          Comment{ "add the specified time to the existing one" }
          };
      }; // ChangeConfig<>
      
      /// Configuration to select a gate to be changed.
      struct SelectGateConfig {
        
        fhicl::Sequence<std::string> Types {
          Name("Types"),
          Comment(
            "gate type to apply the changes on: "
            + BeamTypeSelector.optionListString()
            ),
          std::vector<std::string>
            { BeamTypeSelector.get(BeamType_t::kUnknown).name() }
          };
        
        std::vector<sim::BeamType_t> getBeamTypes() const
          { return Config::getBeamTypes(Types); }
        
      }; // struct SelectGateConfig
      
      fhicl::Table<SelectGateConfig> Select {
        Name{ "Select" },
        Comment
          { "apply these settings only to gates satisfying these criteria" }
        };
      
      fhicl::OptionalTable<ChangeConfig<simulation_time>> Start {
        Name{ "Start" },
        Comment{ "changes to the start of the beam gate" }
        };
      
      fhicl::OptionalTable<ChangeConfig<microseconds>> Width {
        Name{ "Width" },
        Comment{ "changes to the width of the beam gate" }
        };
      
      
      /// Converts a configuration table into a `ChangeRecipe_t`.
      template <typename P, typename I>
      static std::optional<BeamChangeRecipe::ChangeRecipe_t<P, I>> convert
        (std::optional<ChangeConfig<P, I>> const& config);
      
      /// Converts a configuration table into a `GateSelector_t`.
      static BeamChangeRecipe::GateSelector_t convert
        (SelectGateConfig const& config);
      
    }; // struct ChangeGate
    
    
    fhicl::Atom<art::InputTag> BeamGateTag {
      Name{ "BeamGateTag" },
      Comment{ "beam gate data product to operate on" }
      // mandatory
      };

    fhicl::Sequence<fhicl::TableAs<BeamChangeRecipe, ChangeGate>> Changes {
      Name{ "Changes" },
      Comment{ "sets of changes to apply to the beam gate" }
      // mandatory
      };

    fhicl::Atom<bool> KeepInstanceName {
      Name{ "KeepInstanceName" },
      Comment
        { "output data product has the same instance name as in BeamGateTag" },
      false // default
      };

    fhicl::Atom<std::string> OutputInstanceName {
      Name{ "OutputInstanceName" },
      Comment{ "instance name for the output data product" },
      ""
      };

    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "name of the category used for the output" },
      "FixBeamGateInfo" // default
      };
    
    
    /// Selector for `Type` parameter.
    static util::MultipleChoiceSelection<BeamType_t> const BeamTypeSelector;
    
    /// Converts a FHiCL atom into a beam type.
    static sim::BeamType_t getBeamType(fhicl::Atom<std::string> const& type);
    
    /// Converts a FHiCL sequence into a vector of beam types.
    static std::vector<sim::BeamType_t> getBeamTypes
      (fhicl::Sequence<std::string> const& type);
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit FixBeamGateInfo
    (Parameters const& config, art::ProcessingFrame const&);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  
  virtual void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fBeamGateTag; ///< Input beam gate data product.
  
  /// Changes on beam gate.
  std::vector<BeamChangeRecipe> const fChanges;
  
  std::string const fInstanceName; ///< Instance name for the output product.
  
  /// Message facility stream category for output.
  std::string const fLogCategory;
  
  // --- END Configuration variables -------------------------------------------
  
  
  /// Returns a "fixed" beam gate based on the input `beamGate` one.
  sim::BeamGateInfo fixBeamGate(sim::BeamGateInfo const& beamGate) const;
  
  
  /// Returns whether `gate` passes the specified `selection`.
  static bool acceptGate(
    sim::BeamGateInfo const& gate,
    BeamChangeRecipe::GateSelector_t const& selection
    );
  
  /// Applies the changes in `recipe` on the `target` value.
  template <typename T, typename P, typename I>
  static T& applyRecipe(
    T& target,
    std::optional<BeamChangeRecipe::ChangeRecipe_t<P, I>> const& recipe
    );
  

  friend BeamChangeRecipe convert(Config::ChangeGate const& config);
  friend struct dumpRecipe;
  
}; // class icarus::trigger::FixBeamGateInfo


//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  struct dumpRecipe {
    icarus::trigger::FixBeamGateInfo::BeamChangeRecipe const& recipe;
    std::string indent, firstIndent;
    dumpRecipe(
      icarus::trigger::FixBeamGateInfo::BeamChangeRecipe const& recipe,
      std::string indent, std::string firstIndent
      )
      : recipe(recipe)
      , indent(std::move(indent)), firstIndent(std::move(firstIndent))
      {}
    dumpRecipe(
      icarus::trigger::FixBeamGateInfo::BeamChangeRecipe const& recipe,
      std::string const& indent = ""
      )
      : recipe(recipe), indent(indent), firstIndent(indent)
      {}
    
  }; // dumpRecipe
  
  std::ostream& operator<< (std::ostream& out, dumpRecipe const& dr);
  
} // icarus::trigger


//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  // configuration conversions automatically picked by fhicl::TableAs<>
  
  FixBeamGateInfo::BeamChangeRecipe convert
    (FixBeamGateInfo::Config::ChangeGate const& config)
  {
    using ChangeGate = FixBeamGateInfo::Config::ChangeGate;
    return FixBeamGateInfo::BeamChangeRecipe{
        ChangeGate::convert(config.Select())  // selectGates
      , ChangeGate::convert(config.Start())   // start
      , ChangeGate::convert(config.Width())   // width
      };
  } // convert(FixBeamGateInfo::Config::ChangeGate)
  
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Returns whether any of the elements of `coll` compares equal to `value`.
  template <typename Coll, typename T>
  bool contains(Coll const& coll, T const& value) {
    using std::begin, std::end;
    auto const b = begin(coll); auto const e = end(coll);
    return std::find(b, e, value) != e;
  } // contains()
  
} // local namespace
//------------------------------------------------------------------------------
namespace icarus::trigger {
  
  util::MultipleChoiceSelection<FixBeamGateInfo::Config::BeamType_t> const
  FixBeamGateInfo::Config::BeamTypeSelector
    {
        { BeamType_t::kUnknown, "unknown" }
      , { BeamType_t::kBNB,     "BNB" }
      , { BeamType_t::kNuMI,    "NuMI" }
    };
  
} // namespace icarus::trigger


//------------------------------------------------------------------------------
icarus::trigger::FixBeamGateInfo::FixBeamGateInfo
  (Parameters const& config, art::ProcessingFrame const&)
  : art::SharedProducer(config)
  // configuration
  , fBeamGateTag{ config().BeamGateTag() }
  , fChanges{ config().Changes() }
  , fInstanceName{
      config().KeepInstanceName()
      ? fBeamGateTag.instance(): config().OutputInstanceName()
    }
  , fLogCategory(config().LogCategory())
{
  using namespace util::quantities::time_literals;
  
  //
  // parameter validation
  //
  if (config().KeepInstanceName() && !config().OutputInstanceName().empty()) {
    throw art::Exception{ art::errors::Configuration }
      << "Can't set both '" << config().KeepInstanceName.name()
        << " (" << std::boolalpha << config().KeepInstanceName() << ")"
      << "' and '" << config().OutputInstanceName.name() << "' ('"
        << config().OutputInstanceName() << "') at the same time.\n";
  }
  
  for (auto const& [ iChange, change ]: util::enumerate(fChanges)) {
    if (change.valid()) continue;
    throw art::Exception{ art::errors::Configuration }
      << "Incompatible requests for '" << config().Changes.name()
      << "' item [" << iChange << "] " << dumpRecipe(change, "  ", "")
      << "\n";
  } // for changes
  
  //
  // input data declaration
  //
  consumes<std::vector<sim::BeamGateInfo>>(fBeamGateTag);
  
  //
  // output data declaration
  //
  produces<std::vector<sim::BeamGateInfo>>(fInstanceName);
  
  async<art::InEvent>();
  
  //
  // configuration report (short)
  //
  {
    mf::LogInfo log { fLogCategory };
    log << "New beam gates based on '" << fBeamGateTag.encode()
      << "' following " << fChanges.size() << " rules:";
    for (auto const& [iChange, change]: util::enumerate(fChanges)) {
      log << "\n[" << iChange << "] " << dumpRecipe(change, "  ", "");
    } // for all change sets
  }
  
} // icarus::trigger::FixBeamGateInfo::FixBeamGateInfo()


//------------------------------------------------------------------------------
void icarus::trigger::FixBeamGateInfo::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  //
  // read input data product
  //
  auto const& beamGates
    = event.getProduct<std::vector<sim::BeamGateInfo>>(fBeamGateTag);
  
  //
  // "fix"
  //
  std::vector<sim::BeamGateInfo> fixedBeamGates;
  fixedBeamGates.reserve(beamGates.size());
  for (sim::BeamGateInfo const& beamGate: beamGates)
    fixedBeamGates.push_back(fixBeamGate(beamGate));
  
  //
  // put into the event
  //
  event.put(
    std::make_unique<std::vector<sim::BeamGateInfo>>
      (std::move(fixedBeamGates)),
    fInstanceName
    );
  
} // icarus::trigger::FixBeamGateInfo::produce()


//------------------------------------------------------------------------------
sim::BeamGateInfo icarus::trigger::FixBeamGateInfo::fixBeamGate
  (sim::BeamGateInfo const& beamGate) const
{
  using nanosecond = util::quantities::points::nanosecond;
  using nanoseconds = util::quantities::intervals::nanoseconds;
  
  // beamGate is supposed to be stored in simulation_time units (nanoseconds)
  simulation_time start { nanosecond{ beamGate.Start() }};
  simulation_time::interval_t width { nanoseconds{ beamGate.Width() }};
  
  for (BeamChangeRecipe const& change: fChanges) {
    if (!acceptGate(beamGate, change.selectGates)) continue;
    applyRecipe(start, change.start);
    applyRecipe(width, change.width);
  } // for all changes
  
  return sim::BeamGateInfo{
    nanosecond{ start }.value(),
    nanoseconds{ width }.value(),
    beamGate.BeamType()
    };
  
} // icarus::trigger::FixBeamGateInfo::fixBeamGate()


//------------------------------------------------------------------------------
bool icarus::trigger::FixBeamGateInfo::acceptGate
  (sim::BeamGateInfo const& gate, BeamChangeRecipe::GateSelector_t const& selection)
{
  if (!selection.beamTypes.empty()) {
    if (!contains(selection.beamTypes, gate.BeamType())) return false;
  }
  return true;
} // icarus::trigger::FixBeamGateInfo::acceptGate()


//------------------------------------------------------------------------------
template <typename T, typename P, typename I>
T& icarus::trigger::FixBeamGateInfo::applyRecipe(
  T& target, std::optional<BeamChangeRecipe::ChangeRecipe_t<P, I>> const& recipe
) {
  if (recipe) {
    if (recipe->setValue) target = *recipe->setValue;
    if (recipe->addValue) target += *recipe->addValue;
  }
  return target;
} // icarus::trigger::FixBeamGateInfo::applyRecipe()


//------------------------------------------------------------------------------
sim::BeamType_t icarus::trigger::FixBeamGateInfo::Config::getBeamType
  (fhicl::Atom<std::string> const& type)
{
  try {
    return static_cast<sim::BeamType_t>
      (BeamTypeSelector.parse(type()).value());
  }
  catch (util::MultipleChoiceSelectionBase::UnknownOptionError const& e)
  {
    throw art::Exception(art::errors::Configuration)
      << "Invalid value for '" << type.name()
      << "' parameter: '" << e.label() << "'; valid options: "
      << BeamTypeSelector.optionListString() << ".\n";
  }
} // icarus::trigger::FixBeamGateInfo::Config::getBeamType()


//------------------------------------------------------------------------------
auto icarus::trigger::FixBeamGateInfo::Config::getBeamTypes
  (fhicl::Sequence<std::string> const& types) -> std::vector<sim::BeamType_t>
{
  std::vector<sim::BeamType_t> beamTypes;
  beamTypes.reserve(types.size());
  for (std::string const& type: types()) {
    try {
      beamTypes.push_back
        (static_cast<sim::BeamType_t>(BeamTypeSelector.parse(type).value()));
    }
    catch (util::MultipleChoiceSelectionBase::UnknownOptionError const& e)
    {
      throw art::Exception(art::errors::Configuration)
        << "Invalid value for '" << types.name()
        << "' parameter: '" << e.label() << "'; valid options: "
        << BeamTypeSelector.optionListString() << ".\n";
    }
  }
  return beamTypes;
} // icarus::trigger::FixBeamGateInfo::Config::getBeamTypes()


//------------------------------------------------------------------------------
template <typename P, typename I>
auto icarus::trigger::FixBeamGateInfo::Config::ChangeGate::convert
  (std::optional<Config::ChangeGate::ChangeConfig<P, I>> const& config)
  -> std::optional<BeamChangeRecipe::ChangeRecipe_t<P, I>>
{
  return config
    ? std::optional{BeamChangeRecipe::ChangeRecipe_t<P, I>{
      config->SetTo(), config->Add()
      }}
    : std::nullopt
    ;
} // icarus::trigger::FixBeamGateInfo::Config::ChangeGate::convert(ChangeConfig)


//------------------------------------------------------------------------------
auto icarus::trigger::FixBeamGateInfo::Config::ChangeGate::convert
  (Config::ChangeGate::SelectGateConfig const& config)
  -> BeamChangeRecipe::GateSelector_t
{
  return BeamChangeRecipe::GateSelector_t { config.getBeamTypes() };
} // icarus::trigger::FixBeamGateInfo::..::ChangeGate::convert(SelectGateConfig)


//------------------------------------------------------------------------------
std::ostream& icarus::trigger::operator<<
  (std::ostream& out, dumpRecipe const& dr)
{
  using namespace util::quantities::time_literals;
  
  auto printChangeRecipe = [&out, &indent=dr.indent]
    (std::string const& name, auto const& recipe)
    {
      if (!recipe || recipe->empty()) return;
      out << "\n" << indent << "- " << name << ":";
      if (recipe->setValue) out << " set to " << (*recipe->setValue);
      if (recipe->addValue) {
        if (*recipe->addValue < 0_ns) out << " remove " << (-*recipe->addValue);
        else                          out << " add " << (*recipe->addValue);
      } // if adding
    }; // printChangeRecipe()
  
  // selection
  if (dr.recipe.empty()) {
    out << "(no action)";
  }
  else {
    out << dr.firstIndent;
    if (dr.recipe.selectGates.empty()) out << "(always)";
    else {
      out << "(only on gates of type:";
      for (sim::BeamType_t const gate: dr.recipe.selectGates.beamTypes) {
        // to query BeamTypeSelector we need its own strong index type as key...
        out << " " << FixBeamGateInfo::Config::BeamTypeSelector
          .get(static_cast<FixBeamGateInfo::Config::BeamType_t>(gate)).name();
      }
      out << ")";
    }
    
    printChangeRecipe("start", dr.recipe.start);
    printChangeRecipe("width", dr.recipe.width);
  }
  
  return out;
} // icarus::trigger::operator<<(dumpRecipe)


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::trigger::FixBeamGateInfo)


//------------------------------------------------------------------------------
