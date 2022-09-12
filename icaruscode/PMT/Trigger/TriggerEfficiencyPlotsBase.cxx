/**
 * @file   icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.cxx
 * @brief  Base class for _art_modules plotting trigger efficiencies.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 30, 2020
 * @see    icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h
 */

// library header
#include "icaruscode/PMT/Trigger/TriggerEfficiencyPlotsBase.h"


// ICARUS libraries
#include "icaruscode/PMT/Trigger/Algorithms/BeamGateMaker.h"
#include "icaruscode/PMT/Trigger/Utilities/TriggerDataUtils.h" // FillTriggerGates()
#include "icaruscode/PMT/Trigger/Utilities/PlotSandbox.h"
#include "icaruscode/Utilities/DetectorClocksHelpers.h" // makeDetTimings()
#include "icarusalg/Utilities/BinningSpecs.h"
#include "icarusalg/Utilities/ROOTutils.h" // util::ROOT
#include "icarusalg/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // optical_time_ticks..
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/values.h" // util::const_values()
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::to_string()
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
// #include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::::static_assert_on<>

// nutools libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h" // also simb::kCC, ...
#include "nusimdata/SimulationBase/MCParticle.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ConsumesCollector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h" // mf namespace
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TEfficiency.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TTree.h"

// C/C++ standard libraries
#include <ostream>
#include <vector>
#include <array>
#include <string>
#include <memory> // std::make_unique()
#include <utility> // std::pair<>, std::move()
#include <cassert>


//------------------------------------------------------------------------------
namespace {
  
  /// Returns a sorted copy of the specified collection.
  template <typename Coll>
  Coll sortCollection(Coll cont) {
    std::sort(cont.begin(), cont.end());
    return cont;
  }
  
} // local namespace

//------------------------------------------------------------------------------
//---  icarus::trigger::details::TriggerPassCounters
//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::create
  (Threshold_t const& threshold, std::string const& patternName) -> IndexPair_t
{
  
  std::size_t thrIndex = thresholdIndex(threshold);
  if (thrIndex == NoIndex) thrIndex = registerThreshold(threshold);
  
  std::size_t patIndex = patternIndex(patternName);
  if (patIndex == NoIndex) patIndex = registerPattern(patternName);
  
  if (thrIndex >= nThresholds())
    fCounters.resize(thrIndex + 1U, std::vector<Counter_t>{ nPatterns() });
  
  if (patIndex >= nPatterns()) {
    for (auto& thrCounters: fCounters) thrCounters.resize(patIndex + 1U);
  }
  
  assert(hasThreshold(thrIndex));
  assert(hasPattern(patIndex));
  return { thrIndex, patIndex };
  
} // icarus::trigger::details::TriggerPassCounters::create()


//------------------------------------------------------------------------------
std::size_t icarus::trigger::details::TriggerPassCounters::nThresholds() const
  { return fCounters.size(); }

  
//------------------------------------------------------------------------------
std::size_t icarus::trigger::details::TriggerPassCounters::nPatterns() const
  { return fCounters.empty()? 0U: fCounters.front().size(); }

  
//------------------------------------------------------------------------------
std::size_t icarus::trigger::details::TriggerPassCounters::registerThreshold
  (Threshold_t const& threshold)
{
  assert(!hasThreshold(threshold));
  
  std::size_t const newIndex = nThresholds();
  fThresholdIndex[threshold] = newIndex;
  
  assert(hasThreshold(threshold));
  return newIndex;
} // icarus::trigger::details::TriggerPassCounters::registerThreshold()


//------------------------------------------------------------------------------
std::size_t icarus::trigger::details::TriggerPassCounters::registerPattern
  (std::string const& name)
{
  assert(!hasPattern(name));
  
  std::size_t const newIndex = nPatterns();
  fPatternIndex[name] = newIndex;
  
  assert(hasPattern(name));
  return newIndex;
} // icarus::trigger::details::TriggerPassCounters::registerPattern()


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (Threshold_t const& threshold, std::string const& patternName) const
  -> Counter_t const&
{
  IndexPair_t const index
    { thresholdIndex(threshold), patternIndex(patternName) };
  if (index.first == NoIndex)
    throw std::out_of_range{ threshold };
  if (index.second == NoIndex)
    throw std::out_of_range{ patternName };
  return counter(index);
} // icarus::trigger::details::TriggerPassCounters::counter(Threshold_t, string)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (std::size_t threshold, std::string const& patternName) const
  -> Counter_t const&
{
  IndexPair_t const index { threshold, patternIndex(patternName) };
  if (index.first == NoIndex)
    throw std::out_of_range{ std::to_string(threshold) };
  if (index.second == NoIndex)
    throw std::out_of_range{ patternName };
  return counter(index);
} // icarus::trigger::details::TriggerPassCounters::counter(size_t, string)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (Threshold_t const& threshold, std::size_t pattern) const
  -> Counter_t const&
{
  IndexPair_t const index
    { thresholdIndex(threshold), pattern };
  if (index.first == NoIndex)
    throw std::out_of_range{ threshold };
  if (index.second == NoIndex)
    throw std::out_of_range{ std::to_string(pattern) };
  return counter(index);
} // icarus::trigger::details::TriggerPassCounters::counter(Threshold_t, size_t)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (std::size_t threshold, std::size_t pattern) const
  -> Counter_t const&
{
  IndexPair_t const index { threshold, pattern };
  if (index.first == NoIndex)
    throw std::out_of_range{ std::to_string(threshold) };
  if (index.second == NoIndex)
    throw std::out_of_range{ std::to_string(pattern) };
  return counter(index);
} // icarus::trigger::details::TriggerPassCounters::counter(size_t, size_t)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (IndexPair_t indices) const -> Counter_t const&
{
  if (!hasThreshold(indices.first))
    throw std::out_of_range{ std::to_string(indices.first) };
  
  auto const& thrCounters = fCounters[indices.first];
  if (indices.second >= thrCounters.size())
    throw std::out_of_range{ std::to_string(indices.second) };
  
  return thrCounters[indices.second];
} // icarus::trigger::details::TriggerPassCounters::counter(IndexPair_t)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (Threshold_t const& threshold, std::string const& patternName) -> Counter_t&
{
  IndexPair_t const index
    { thresholdIndex(threshold), patternIndex(patternName) };
  if (index.first == NoIndex)
    throw std::out_of_range{ threshold };
  if (index.second == NoIndex)
    throw std::out_of_range{ patternName };
  return counter(index);
} // icarus::trigger::details::TriggerPassCounters::counter(Threshold_t, string)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (std::size_t threshold, std::string const& patternName) -> Counter_t&
{
  IndexPair_t const index { threshold, patternIndex(patternName) };
  if (index.first == NoIndex)
    throw std::out_of_range{ std::to_string(threshold) };
  if (index.second == NoIndex)
    throw std::out_of_range{ patternName };
  return counter(index);
} // icarus::trigger::details::TriggerPassCounters::counter(size_t, string)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (Threshold_t const& threshold, std::size_t pattern) -> Counter_t&
{
  IndexPair_t const index
    { thresholdIndex(threshold), pattern };
  if (index.first == NoIndex)
    throw std::out_of_range{ threshold };
  if (index.second == NoIndex)
    throw std::out_of_range{ std::to_string(pattern) };
  return counter(index);
} // icarus::trigger::details::TriggerPassCounters::counter(Threshold_t, size_t)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter
  (std::size_t threshold, std::size_t pattern) -> Counter_t&
{
  IndexPair_t const index { threshold, pattern };
  if (index.first == NoIndex)
    throw std::out_of_range{ std::to_string(threshold) };
  if (index.second == NoIndex)
    throw std::out_of_range{ std::to_string(pattern) };
  return counter(index);
} // icarus::trigger::details::TriggerPassCounters::counter(size_t, size_t)


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::counter(IndexPair_t indices)
  -> Counter_t&
{
  if (!hasThreshold(indices.first))
    throw std::out_of_range{ std::to_string(indices.first) };
  
  auto& thrCounters = fCounters[indices.first];
  if (indices.second >= thrCounters.size())
    throw std::out_of_range{ std::to_string(indices.second) };
  
  return thrCounters[indices.second];
} // icarus::trigger::details::TriggerPassCounters::counter(IndexPair_t)


//------------------------------------------------------------------------------
bool icarus::trigger::details::TriggerPassCounters::hasThreshold
  (Threshold_t const& threshold) const
  { return fThresholdIndex.find(threshold) != fThresholdIndex.end(); }


//------------------------------------------------------------------------------
bool icarus::trigger::details::TriggerPassCounters::hasThreshold
  (std::size_t thresholdIndex) const
  { return thresholdIndex < nThresholds(); }


//------------------------------------------------------------------------------
bool icarus::trigger::details::TriggerPassCounters::hasPattern
  (std::string const& patternName) const
  { return fPatternIndex.find(patternName) != fPatternIndex.end(); }


//------------------------------------------------------------------------------
bool icarus::trigger::details::TriggerPassCounters::hasPattern
  (std::size_t patternIndex) const
  { return patternIndex < nPatterns(); }


//------------------------------------------------------------------------------
std::size_t icarus::trigger::details::TriggerPassCounters::thresholdIndex
  (Threshold_t const& threshold) const
{
  auto const iIndex = fThresholdIndex.find(threshold);
  return (iIndex == fThresholdIndex.end())? NoIndex: iIndex->second;
} // icarus::trigger::details::TriggerPassCounters::thresholdIndex()


//------------------------------------------------------------------------------
std::size_t icarus::trigger::details::TriggerPassCounters::patternIndex
  (std::string const& patternName) const
{
  auto const iIndex = fPatternIndex.find(patternName);
  return (iIndex == fPatternIndex.end())? NoIndex: iIndex->second;
} // icarus::trigger::details::TriggerPassCounters::patternIndex()


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::threshold
  (std::size_t index) const -> Threshold_t const&
{
  // reverse lookup: slow...
  for (auto const& [ threshold, thrIndex ]: fThresholdIndex)
    if (thrIndex == index) return threshold;
  throw std::out_of_range{ std::to_string(index) };
} // icarus::trigger::details::TriggerPassCounters::threshold()


//------------------------------------------------------------------------------
auto icarus::trigger::details::TriggerPassCounters::patternName
  (std::size_t index) const -> std::string const&
{
  // reverse lookup: slow...
  for (auto const& [ name, patIndex ]: fPatternIndex)
    if (patIndex == index) return name;
  throw std::out_of_range{ std::to_string(index) };
} // icarus::trigger::details::TriggerPassCounters::patternName()


//------------------------------------------------------------------------------
void icarus::trigger::details::TriggerPassCounters::dump
  (std::ostream& out) const
{
  out << "Triggers for " << nThresholds() << " thresholds and " << nPatterns()
    << " patterns:";
  for (auto const iThr: util::counter(nThresholds())) {
    
    assert(hasThreshold(iThr));
    out << "\n  threshold " << threshold(iThr) << " [#" << iThr << "]:";
    unsigned int nonEmptyPatterns = 0U;
    for (auto const iPat: util::counter(nPatterns())) {
      assert(hasPattern(iPat));
      auto const& counts = counter(iThr, iPat);
      if (counts.empty()) continue;
      out << "\n    " << patternName(iPat) << " [#" << iPat << "]: "
        << counts.passed() << " / " << counts.total();
      ++nonEmptyPatterns;
    } // for patterns
    if (nonEmptyPatterns == 0) out << " no events";
    
  } // for threshold
  out << "\n";
} // icarus::trigger::details::TriggerPassCounters::dump()


//------------------------------------------------------------------------------
std::ostream& icarus::trigger::details::operator<<
  (std::ostream& out, TriggerPassCounters const& counters)
  { counters.dump(out); return out; }


//------------------------------------------------------------------------------
//--- icarus::trigger::details::PlotInfoTree
//------------------------------------------------------------------------------
icarus::trigger::details::PlotInfoTree::PlotInfoTree(TTree& tree)
  : TreeHolder(tree)
{

  this->tree().Branch("InPlots", &fInPlots);

} // icarus::trigger::details::PlotInfoTree::PlotInfoTree()


//------------------------------------------------------------------------------
void icarus::trigger::details::PlotInfoTree::assign(bool inPlots) {

  fInPlots = static_cast<Bool_t>(inPlots);

} // icarus::trigger::details::PlotInfoTree::assignEvent()


//------------------------------------------------------------------------------
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
icarus::trigger::TriggerEfficiencyPlotsBase::PlotCategories_t const
icarus::trigger::TriggerEfficiencyPlotsBase::DefaultPlotCategories {

  PlotCategory{
    "All"
    },

  PlotCategory{
    "All nu_mu", "nu_mu",
    [](EventInfo_t const& info){ return info.hasGenerated() && info.isNu_mu(); }
  },

  PlotCategory{
    "All nu_e", "nu_e",
    [](EventInfo_t const& info){ return info.hasGenerated() && info.isNu_e(); }
  },

  PlotCategory{
    "NuCC", "CC",
    [](EventInfo_t const& info)
      { return info.hasGenerated() && info.isWeakChargedCurrent(); }
    },

  PlotCategory{
    "NuCC_mu", "CC_mu",
    [](EventInfo_t const& info)
      { return info.hasGenerated() && info.isWeakChargedCurrent() && info.isNu_mu(); }
    },

  PlotCategory{
    "NuCC_e", "CC_e",
    [](EventInfo_t const& info)
      { return info.hasGenerated() && info.isWeakChargedCurrent() && info.isNu_e(); }
    },

  PlotCategory{
    "NuNC", "NC",
    [](EventInfo_t const& info)
      { return info.hasGenerated() && info.isWeakNeutralCurrent(); }
    },

  PlotCategory{
    "NuNC_mu", "NC_mu",
    [](EventInfo_t const& info)
      { return info.hasGenerated() && info.isWeakNeutralCurrent() && info.isNu_mu(); }
    },

  PlotCategory{
    "NuNC_e", "NC_e",
    [](EventInfo_t const& info)
      { return info.hasGenerated() && info.isWeakNeutralCurrent() && info.isNu_e(); }
    },

  PlotCategory{
    "NoActivity", "no energy deposited in active volume during beam gate",
    [](EventInfo_t const& info)
      {
        using namespace util::quantities::energy_literals;
        return info.hasDepEnergy()
          && (info.DepositedEnergyInSpillInActiveVolume() == 0.0_GeV);
      }
    }

}; // icarus::trigger::TriggerEfficiencyPlotsBase::DefaultPlotCategories[]


//------------------------------------------------------------------------------
icarus::trigger::TriggerEfficiencyPlotsBase::TriggerEfficiencyPlotsBase
  (Config const& config, art::ConsumesCollector& consumer)
  // configuration
  : fDetectorParticleTag  (config.DetectorParticleTag())
  , fBeamGateDuration     (config.BeamGateDuration())
  , fBeamGateStart        (config.BeamGateStart())
  , fPreSpillWindow       (config.PreSpillWindow())
  , fPreSpillStart
      (fBeamGateStart - config.PreSpillWindowGap() - fPreSpillWindow)
  , fTriggerTimeResolution(config.TriggerTimeResolution())
  , fPlotOnlyActiveVolume (config.PlotOnlyActiveVolume())
  , fOnlyPlotCategories   (sortCollection(config.OnlyPlotCategories()))
  , fLogCategory          (config.LogCategory())
  // services
  , fGeom      (*lar::providerFrom<geo::Geometry>())
  , fOutputDir (*art::ServiceHandle<art::TFileService>())
  // cached
  , fEventInfoExtractorMaker(
      config.GeneratorTags(),              // truthTags
      makeEdepTag(config.EnergyDepositTags, config.EnergyDepositSummaryTag),
                                           // edepTags
      fGeom,                               // geom
      nullptr,                             // detProps
      nullptr,                             // detTimings
      fLogCategory,                        // logCategory
      consumer                             // consumesCollector
    )
  , fChannelCryostat(makeChannelCryostatMap(fGeom))
{
  //
  // more complex parameter parsing
  //
  if (config.EventDetailsLogCategory(fLogEventDetails)) {
    // the parameter is somehow set, so fLogEventDetails won't be empty;
    // but if the parameter value is specified empty, we set it to fLogCategory
    if (fLogEventDetails.empty()) fLogEventDetails = fLogCategory;
  } // if EventDetailsLogCategory is specified

  std::string const discrModuleLabel = config.TriggerGatesTag();
  for (std::string const& threshold: config.Thresholds())
    fADCthresholds[threshold] = art::InputTag{ discrModuleLabel, threshold };
  

  if (config.EventTreeName.hasValue()) {
    
    std::string treeName;
    config.EventTreeName(treeName);

    fIDTree = std::make_unique<details::EventIDTree>
      (*(fOutputDir.make<TTree>(treeName.c_str(), "Event information")));
    fEventTree = std::make_unique<details::EventInfoTree>
      (fIDTree->tree(), useGen(), useEDep());
    fPlotTree = std::make_unique<details::PlotInfoTree>(fIDTree->tree());

  } // if make tree

  //
  // input data declaration
  //
  using icarus::trigger::OpticalTriggerGateData_t; // for convenience
  
  // trigger primitives
  for (art::InputTag const& inputDataTag: util::const_values(fADCthresholds)) {
    consumer.consumes<std::vector<OpticalTriggerGateData_t>>(inputDataTag);
    consumer.consumes<art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>
      (inputDataTag);
  } // for

  {
    mf::LogInfo log(fLogCategory);
    log << "\nConfigured " << fADCthresholds.size() << " thresholds (ADC):";
    for (auto const& [ thresholdTag, dataTag ]: fADCthresholds)
      log << "\n * " << thresholdTag << " (from '" << dataTag.encode() << "')";
    
  } // local block
  
  
  if (!useGen()) {
    mf::LogVerbatim(fLogCategory)
      << "Generation information will not be produced.";
  }
  if (!useEDep()) {
    mf::LogVerbatim(fLogCategory)
      << "Energy deposition information will not be produced.";
  }
  
} // icarus::trigger::TriggerEfficiencyPlots::TriggerEfficiencyPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::process
  (art::Event const& event)
{

  /*
   * 1. find out the features of the event and the categories it belongs to
   * 2. for each threshold:
   *   1. read the trigger primitives
   *   2. pick the plots to be filled
   *   3. combine the trigger primitives, apply the beam gate,
   *      generate the trigger response, add the response to all the plots
   *      (delegated)
   *
   */
  
  ++nEvents;
  
  //
  // 1. find out the features of the event and the categories it belongs to
  //
  
  auto const [ detTimings, beamGate, preSpillWindow ] = makeGatePack(&event);
  
  if (auto oldGate = fBeamGateChangeCheck(beamGate); oldGate) {
    mf::LogDebug(fLogCategory)
      << "Beam gate has changed from " << oldGate->asOptTickRange()
      << " to " << beamGate.asOptTickRange() << " (optical tick)!";
  }
  
  EventInfo_t const eventInfo = fEventInfoExtractorMaker
    (beamGate.asSimulationRange(), preSpillWindow.asSimulationRange())(event);

  
  bool const bPlot = shouldPlotEvent(eventInfo);
  if (bPlot) ++nPlottedEvents;
  
  if (fIDTree) fIDTree->assignID(event.id());
  if (fPlotTree) fPlotTree->assign(bPlot);
  if (fEventTree) fEventTree->assignEvent(eventInfo);
  
  std::vector<std::string> selectedPlotCategories
    = selectPlotCategories(eventInfo, fPlotCategories);
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
  // 2. for each PMT threshold:
  //
  auto const clockData
   = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
  for (auto&& [ iThr, thrPair, thrPlots ]
    : util::enumerate(fADCthresholds, fThresholdPlots)
  ) {

    auto const& [ thresholdTag, dataTag ] = thrPair;

    //
    // 2.1. read the trigger primitives
    //
    
    TriggerGatesPerCryostat_t const& cryoGates
      = splitByCryostat(readTriggerGates(event, dataTag));
    
    //
    // 2.2. pick the plots to be filled
    //
    PlotSandboxRefs_t selectedPlots;
    
    if (bPlot) {
      for (std::string const& name: selectedPlotCategories)
        selectedPlots.emplace_back(*(thrPlots.findSandbox(name)));
    }
    
    //
    // 2.3. combine the trigger primitives, apply the beam gate,
    //      generate the trigger response, add the response to all the plots
    //
    simulateAndPlot(
      iThr, // settings index
      cryoGates,
      eventInfo,
      clockData,
      selectedPlots
      );

  } // for thresholds

  //
  // store information in output tree if any
  //
  if (fIDTree) fIDTree->tree().Fill();

} // icarus::trigger::TriggerEfficiencyPlotsBase::process()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::printSummary() const {
  
  mf::LogInfo log(fLogCategory);
  log << nPlottedEvents << "/" << nEvents << " events plotted.";
  
  log << "\n" << fPassCounters;
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::printSummary()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::initializePlots
  (PlotCategories_t categories, std::vector<SettingsInfo_t> const& settings)
{
  using namespace std::string_literals;
  
  auto const beamGate = icarus::trigger::makeBeamGateStruct
    (icarus::ns::util::makeDetTimings(), fBeamGateDuration, fBeamGateStart);

  fBeamGateChangeCheck(beamGate);
  
  if (fOnlyPlotCategories.empty()) fPlotCategories = std::move(categories);
  else {
    auto const plotThisCategory = [this](std::string const& name)
      {
        return std::binary_search
          (fOnlyPlotCategories.begin(), fOnlyPlotCategories.end(), name);
      };
    
    fPlotCategories.clear();
    for (auto&& plotCategory: categories) {
      if (!plotThisCategory(plotCategory.name())) continue;
      fPlotCategories.push_back(std::move(plotCategory));
    } // for
  }
  
  
  {
    mf::LogTrace(fLogCategory)
      << "Beam gate:"
      << "\n - electronics time: " << beamGate.asGate()
      << "\n - simulation time: " << beamGate.asSimulationRange()
      << "\n - optical ticks: " << beamGate.asOptTickRange()
      ;
    
    mf::LogInfo log(fLogCategory);
    log << "\nConfigured " << fADCthresholds.size() << " thresholds (ADC):";
    for (auto const& [ thresholdTag, dataTag ]: fADCthresholds)
      log << "\n * " << thresholdTag << " (from '" << dataTag.encode() << "')";
    log << "\nBeam gate for plots is " << beamGate.asSimulationRange();
    
    log << "\nConfigured " << fPlotCategories.size() << " plot categories"
      << (fPlotCategories.empty()? '.': ':');
    for (auto const& plotCategory: fPlotCategories) {
      log << "\n ['" << plotCategory.name() << "'] "
        << plotCategory.description();
    } // for
    
  } // local block
  
  
  for (std::string const& thresholdTag: util::get_elements<0U>(fADCthresholds))
  {
    // create a plot sandbox inside `fOutputDir` with a name/prefix `Thr###`
    icarus::trigger::PlotSandbox thrPlots
      { fOutputDir, "Thr"s + thresholdTag, "(thr: "s + thresholdTag + ")"s };
    
    // create a subbox for each plot category
    for (PlotCategory const& category: fPlotCategories) {
      PlotSandbox& plots = thrPlots.addSubSandbox(
        category.name(),
        category.description()
        );
      
      initializePlotSet(plots, settings);
    } // for plot category
    fThresholdPlots.push_back(std::move(thrPlots));
  } // for thresholds
  
  mf::LogTrace log(fLogCategory);
  log << "Created " << fThresholdPlots.size() << " plot boxes:\n";
  for (auto const& box: fThresholdPlots) {
    box.dump(log, "  ");
  } // for
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::initializePlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::initializePlotSet
  (PlotSandbox& plots, std::vector<SettingsInfo_t> const& settings) const
{
  
  //
  // Selection-related plots
  //
  initializeEventPlots(plots);
  
  initializePMTplots(plots);
  
  //
  // Plots per trigger setting, split in triggering and not triggering events;
  // the plot set is the same as the "global" one.
  //
  using SS_t = std::pair<std::string, std::string>;
  std::array<SS_t, 2U> const classes {
    SS_t{ "triggering", "triggering events" },
    SS_t{ "nontriggering", "non-triggering events" }
    };
  for (auto const& settingsDesc: settings) {
    
    // this defines a specific trigger, with its thresholds and settings
    PlotSandbox& reqBox
      = plots.addSubSandbox(settingsDesc.tag, settingsDesc.description);
    
    initializeEfficiencyPerTriggerPlots(reqBox);
    
    for (auto const& [ name, desc ]: classes) {
      
      PlotSandbox& box = reqBox.addSubSandbox(name, desc);
      
      initializeEventPlots(box);
      
      initializePMTplots(box);
      
    } // for triggering requirement
  } // for triggering classes
  
 
} // icarus::trigger::TriggerEfficiencyPlotsBase::initializePlotSet()


//------------------------------------------------------------------------------
void
icarus::trigger::TriggerEfficiencyPlotsBase::initializeEfficiencyPerTriggerPlots
  (PlotSandbox& plots) const
{
  
  auto const [ detTimings, beamGate, preSpillWindow ] = makeGatePack();

  detinfo::timescales::optical_time_ticks const triggerResolutionTicks
    { detTimings.toOpticalTicks(fTriggerTimeResolution) };
  
  auto const PreSpillDuration = preSpillWindow.asSimulationRange().duration();
  
  //
  // Triggering efficiency vs. something else
  //
  if (useEDep()) {
    plots.make<TEfficiency>(
      "EffVsEnergyInSpill",
      "Efficiency of triggering vs. energy deposited in spill"
        ";energy deposited in spill  [ GeV ]"
        ";trigger efficiency  [ / 50 MeV ]",
      120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
      );
    
    plots.make<TEfficiency>(
      "EffVsEnergyInSpillActive",
      "Efficiency of triggering vs. energy deposited in active volume"
        ";energy deposited in active volume in spill  [ GeV ]"
        ";trigger efficiency  [ / 50 MeV ]",
      120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
      );
    
    plots.make<TEfficiency>(
      "EffVsEnergyInPreSpill",
      (
        "Efficiency of triggering vs. energy deposited in pre-spill ("
        + to_string(PreSpillDuration) + ")"
        ";energy deposited in pre-spill  [ GeV ]"
        ";trigger efficiency  [ / 100 MeV ]"
      ).c_str(),
      120, 0.0, 12.0
      );
    
    plots.make<TEfficiency>(
      "EffVsEnergyInPreSpillActive",
      (
        "Efficiency of triggering vs. energy deposited in active volume"
          " (pre-spill: " + to_string(PreSpillDuration) + ")"
          ";energy deposited in active volume in pre-spill  [ GeV ]"
          ";trigger efficiency  [ / 100 MeV ]"
      ).c_str(),
      120, 0.0, 12.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
      );
  } // if plots with deposited energy
  
  if (useGen()) {
    plots.make<TEfficiency>(
      "EffVsNeutrinoEnergy",
      "Efficiency of triggering vs. neutrino energy"
        ";neutrino true energy  [ GeV ]"
        ";trigger efficiency  [ / 50 MeV ]",
      120, 0.0, 6.0 // 6 GeV is not that much for NuMI, but we should be ok
      );
    
    plots.make<TEfficiency>(
      "EffVsLeptonEnergy",
      "Efficiency of triggering vs. outgoing lepton energy"
        ";final state lepton true energy  [ GeV ]"
        ";trigger efficiency  [ / 50 MeV ]",
      120, 0.0, 6.0
      );
  } // if plots with generated info
  
  auto const& beamGateOpt = beamGate.asOptTickRange();
  plots.make<TH1F>(
    "TriggerTick",
    "Trigger time tick"
      ";optical time tick [ /" + util::to_string(triggerResolutionTicks) + " ]",
    beamGateOpt.duration() / triggerResolutionTicks,
    beamGateOpt.start().value(), beamGateOpt.end().value()
    );

  // plots will be relative to the beam gate:
  constexpr util::quantities::intervals::microseconds beamPlotPadding { 4_us };
  icarus::ns::util::BinningSpecs const beamGateBinning = alignBinningTo(
    icarus::ns::util::BinningSpecs{
      (
        std::min(
          preSpillWindow.asElectronicsTimeRange().start(),
          beamGate.asElectronicsTimeRange().start()
        )
        - beamGate.asElectronicsTimeRange().start()
        - beamPlotPadding
      ).value(),
      (std::max(
        preSpillWindow.asElectronicsTimeRange().end()
          - beamGate.asElectronicsTimeRange().start(),
        beamGate.duration()
        ) + beamPlotPadding).value(),
      fTriggerTimeResolution.value()
      },
      0.0
    );
  
  plots.make<TH1F>(
    "OpeningTimes",
    "Times at which trigger logic was satisfied"
      ";trigger time (relative to beam gate opening)  [ us ]"
      ";opened trigger gates",
    beamGateBinning.nBins(), beamGateBinning.lower(), beamGateBinning.upper()
    );
  
  plots.make<TH1F>(
    "TriggerTime",
    "Time of the trigger"
      ";trigger time (relative to beam gate opening)  [ us ]"
      ";opened trigger gates",
    beamGateBinning.nBins(), beamGateBinning.lower(), beamGateBinning.upper()
    );

  //
  // plots independent of the trigger primitive requirements
  //

} // icarus::trigger::TriggerEfficiencyPlotsBase::initializeEfficiencyPerTriggerPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::initializeEventPlots
  (PlotSandbox& plots) const
{
  
  auto const [ detTimings, beamGate, preSpillWindow ] = makeGatePack();
  
  auto const BeamGateDuration = beamGate.asSimulationRange().duration();
  auto const PreSpillDuration = preSpillWindow.asSimulationRange().duration();
  
  //
  // Selection-related plots
  //
  if (useGen()) {
    plots.make<TH1F>(
      "NeutrinoEnergy",
      "True Neutrino Energy"
        ";neutrino energy [GeV]"
        ";events",
      120, 0.0, 6.0 // GeV
    );
  }
  
  if (useEDep()) {
    plots.make<TH1F>(
      "EnergyInSpill",
      "Energy deposited during the beam gate opening"
        ";energy deposited in spill [ GeV ]"
        ";events  [ / 50 MeV ]",
      120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
      );
    plots.make<TH1F>(
      "EnergyInPreSpill",
      (
        "Energy deposited during the pre-spill window ("
          + to_string(PreSpillDuration) + ")"
        ";energy deposited in pre-spill [ GeV ]"
        ";events  [ / 100 MeV ]"
      ).c_str(),
      120, 0.0, 12.0
      );
    plots.make<TH1F>(
      "EnergyInSpillActive",
      "Energy deposited during the beam gate opening in active volume"
        ";energy deposited in active volume in spill [ GeV ]"
        ";events  [ / 50 MeV ]",
      120, 0.0, 6.0 // 6 GeV should be enough for a MIP crossing 20 m of detector
      );
    plots.make<TH1F>(
      "EnergyInPreSpillActive",
      (
        "Energy deposited in active volume during the pre-spill window ("
          + to_string(PreSpillDuration) + ")"
        ";energy deposited in active volume in pre-spill [ GeV ]"
        ";events  [ / 100 MeV ]"
      ).c_str(),
      120, 0.0, 12.0
      );
    plots.make<TH2F>(
      "EnergyInPreSpillVsSpillActive",
      (
        "Energy deposited in active volume" 
        ";energy in spill window (" + to_string(BeamGateDuration) + ")  [ GeV ]"
        ";energy in pre-spill window (" + to_string(PreSpillDuration)
          + ")  [ GeV ]"
      ).c_str(),
      120, 0.0, 6.0, 120, 0.0, 12.0
      );
  } // if use energy deposition
  
  if (useGen()) {
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
        ";beam direction (z)  [ / 20 cm ]"
        ";vertical direction (y)  [ / 5 cm ]",
      120, -1200., +1200.,
      100,  -250.,  +250.
      );

    plots.make<TH2F>(
      "InteractionTypeNeutrinoEnergy",
      "Interaction Type vs Neutrino Energy"
      ";InteractionType"
      ";Neutrino Energy",
      200,999.5,1199.5,
      120, 0.0, 6.0
      );

  } // if generated information
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::initializeEventPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::initializePMTplots
  (PlotSandbox& plots) const
{
  
  unsigned int const nOpChannels = fGeom.NOpChannels();
  
  //
  // plots independent of the trigger primitive requirements
  //
  plots.make<TH1I>(
    "ActivePMT",
    "PMT channels contributing to the trigger"
    ";channel with opened trigger gate"
    ";events",
    nOpChannels, // large number, zoom in presentations!
    0.0, static_cast<double>(nOpChannels)
    );
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::initializePMTplots()


//------------------------------------------------------------------------------
bool icarus::trigger::TriggerEfficiencyPlotsBase::shouldPlotEvent
  (EventInfo_t const& eventInfo) const
{
  if (fPlotOnlyActiveVolume
    && eventInfo.hasVertex() && !eventInfo.isInActiveVolume())
  {
    return false;
  }
  
  return true;
} // icarus::trigger::TriggerEfficiencyPlotsBase::shouldPlotEvent()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::fillEventPlots
  (EventInfo_t const& eventInfo, PlotSandbox const& plots) const
{
  
  using namespace std::string_literals;
  
  HistGetter const getTrig { plots };
  
  if (useEDep()) {
    assert(eventInfo.hasDepEnergy());
    getTrig.Hist("EnergyInSpill"s).Fill(double(eventInfo.DepositedEnergyInSpill()));
    getTrig.Hist("EnergyInSpillActive"s).Fill(double(eventInfo.DepositedEnergyInSpillInActiveVolume()));
    getTrig.Hist("EnergyInPreSpill"s)
      .Fill(double(eventInfo.DepositedEnergyInPreSpill()));
    getTrig.Hist("EnergyInPreSpillActive"s)
      .Fill(double(eventInfo.DepositedEnergyInPreSpillInActiveVolume()));
    getTrig.Hist2D("EnergyInPreSpillVsSpillActive"s).Fill(
      double(eventInfo.DepositedEnergyInSpillInActiveVolume()),
      double(eventInfo.DepositedEnergyInPreSpillInActiveVolume())
      );
  }
  if (useGen()) {
    if (eventInfo.isNeutrino()) {
      assert(eventInfo.hasGenerated());
      getTrig.Hist("NeutrinoEnergy"s).Fill(double(eventInfo.NeutrinoEnergy()));
      getTrig.Hist("InteractionType"s).Fill(eventInfo.InteractionType());
      getTrig.Hist("LeptonEnergy"s).Fill(double(eventInfo.LeptonEnergy()));
      getTrig.Hist("InteractionTypeNeutrinoEnergy"s).Fill(double(eventInfo.InteractionType()), double(eventInfo.NeutrinoEnergy()));
    } // if neutrino event
    TH2& vertexHist = getTrig.Hist2D("InteractionVertexYZ"s);
    for (auto const& point: eventInfo.Vertices())
      vertexHist.Fill(point.Z(), point.Y());
  } // if use generated information
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::fillEventPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::fillPMTplots
  (PMTInfo_t const& PMTinfo, PlotSandbox const& plots) const
{
  
  using namespace std::string_literals;
  
  HistGetter const getTrig { plots };
  
  auto& activePMThist = getTrig.Hist("ActivePMT"s);
  for (raw::Channel_t const channel: PMTinfo.activeChannels())
    activePMThist.Fill(channel);
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::fillPMTplots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::fillEfficiencyPlots(
  EventInfo_t const& eventInfo,
  TriggerInfo_t const& triggerInfo,
  PlotSandbox const& plots
) const {
  
  using namespace std::string_literals;
  using OpeningInfo_t = icarus::trigger::details::TriggerInfo_t::OpeningInfo_t;

  auto const detTimings = icarus::ns::util::makeDetTimings();

  HistGetter const getTrigEff { plots };
  
  bool const fired = triggerInfo.fired();

  // efficiency plots
  if (useEDep()) {
    getTrigEff.Eff("EffVsEnergyInSpill"s).Fill
      (fired, double(eventInfo.DepositedEnergyInSpill()));
    getTrigEff.Eff("EffVsEnergyInPreSpill"s).Fill
      (fired, double(eventInfo.DepositedEnergyInPreSpill()));
    getTrigEff.Eff("EffVsEnergyInSpillActive"s).Fill
      (fired, double(eventInfo.DepositedEnergyInSpillInActiveVolume()));
    getTrigEff.Eff("EffVsEnergyInPreSpillActive"s).Fill
      (fired, double(eventInfo.DepositedEnergyInPreSpillInActiveVolume()));
  } // if use energy deposits
  if (useGen()) {
    if (eventInfo.isNeutrino()) {
      getTrigEff.Eff("EffVsNeutrinoEnergy"s).Fill
        (fired, double(eventInfo.NeutrinoEnergy()));
      getTrigEff.Eff("EffVsLeptonEnergy"s).Fill
        (fired, double(eventInfo.LeptonEnergy()));
    }
  } // if use generated information
  
  if (fired) {
    detinfo::timescales::electronics_time const beamGateTime
      = detTimings.BeamGateTime();
    
    getTrigEff.Hist("TriggerTick"s).Fill(triggerInfo.atTick().value());
    
    // converts the tick in the argument into electronics time:
    auto openingTime = [&detTimings](OpeningInfo_t const& info)
      { return detTimings.toElectronicsTime(info.tick); };

    getTrigEff.Hist("TriggerTime"s).Fill
      ((openingTime(triggerInfo.main()) - beamGateTime).value());

    std::vector<OpeningInfo_t> const& allTriggerOpenings = triggerInfo.all();

    for (OpeningInfo_t const& opening : allTriggerOpenings) {
      getTrigEff.Hist("OpeningTimes"s).Fill
        ((openingTime(opening) - beamGateTime).value());
    } // for all trigger openings
    
  } // if fired
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::fillEfficiencyPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::fillAllEfficiencyPlots(
  EventInfo_t const& eventInfo,
  PMTInfo_t const& PMTinfo,
  TriggerInfo_t const& triggerInfo,
  PlotSandbox const& plots
) const {
  
  fillEfficiencyPlots(eventInfo, triggerInfo, plots);
  
  // plotting split for triggering/not triggering events
  fillEventPlots(
    eventInfo,
    plots.demandSandbox(triggerInfo.fired()? "triggering": "nontriggering")
    );
  
  fillPMTplots(
    PMTinfo,
    plots.demandSandbox(triggerInfo.fired()? "triggering": "nontriggering")
    );
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::fillAllEfficiencyPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::deleteEmptyPlots()
{
  
  for (auto& thrPlots: fThresholdPlots) deleteEmptyPlots(thrPlots);
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::deleteEmptyPlots()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlotsBase::createCountersForPattern
  (std::string const& patternName) -> std::size_t
{
  
  std::size_t patternIndex = fPassCounters.NoIndex;
  std::size_t iThr [[maybe_unused]] = 0U;
  for (std::string const& thresholdTag: util::get_elements<0U>(fADCthresholds))
  {
    
    auto const indices = fPassCounters.create(thresholdTag, patternName);
    if (patternIndex == fPassCounters.NoIndex) patternIndex = indices.second;
    else assert(indices.second == patternIndex);
    
  } // for thresholds
  
  return patternIndex;
} // icarus::trigger::TriggerEfficiencyPlotsBase::createCountersForPattern()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::registerTriggerResult
  (std::size_t threshold, std::size_t pattern, bool fired)
  { fPassCounters(threshold, pattern).add(fired); }


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::registerTriggerResult
  (std::size_t threshold, std::size_t pattern, TriggerInfo_t const& triggerInfo)
  { registerTriggerResult(threshold, pattern, triggerInfo.fired()); }


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlotsBase::makeGatePack
  (art::Event const* event /* = nullptr */) const -> GatePack_t
{
  auto const detTimings = icarus::ns::util::makeDetTimings(event);
  return GatePack_t{
    detTimings,
    icarus::trigger::makeBeamGateStruct
      (detTimings, fBeamGateDuration, fBeamGateStart),
    icarus::trigger::makeBeamGateStruct
      (detTimings, fPreSpillStart, fPreSpillStart + fPreSpillWindow)
    };
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::makeGatePack()


//------------------------------------------------------------------------------
std::vector<std::string>
icarus::trigger::TriggerEfficiencyPlotsBase::selectPlotCategories
  (EventInfo_t const& info, PlotCategories_t const& categories) const
{
  std::vector<std::string> selected;
  
  for (auto const& category: categories)
    if (category(info)) selected.push_back(category);
  
  return selected;
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::selectPlotCategories()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlotsBase::readTriggerGates
  (art::Event const& event, art::InputTag const& dataTag) const
  -> TriggerGates_t
{

  using icarus::trigger::OpticalTriggerGateData_t; // for convenience

  // currently the associations are a waste of time memory...
  auto const& gates
    = event.getProduct<std::vector<OpticalTriggerGateData_t>>(dataTag);
  auto const& gateToWaveforms = event.getProduct
      <art::Assns<OpticalTriggerGateData_t, sbn::OpDetWaveformMeta>>(dataTag);
  
  try {
    return icarus::trigger::FillTriggerGates(gates, gateToWaveforms);
  }
  catch (cet::exception const& e) {
    throw cet::exception("TriggerEfficiencyPlots", "", e)
      << "Error encountered while reading data products from '"
      << dataTag.encode() << "'\n";
  }

} // icarus::trigger::TriggerEfficiencyPlotsBase::readTriggerGates()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlotsBase::splitByCryostat
  (TriggerGates_t&& gates) const -> TriggerGatesPerCryostat_t
{

  TriggerGatesPerCryostat_t gatesPerCryostat{ fGeom.Ncryostats() };
  
  for (auto& gate: gates) {
    gatesPerCryostat[fChannelCryostat.at(gate.channels().front()).Cryostat]
      .push_back(std::move(gate));
  } // for gates
  
  return gatesPerCryostat;

} // icarus::trigger::TriggerEfficiencyPlotsBase::splitByCryostat()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlotsBase::extractActiveChannels
  (TriggerGatesPerCryostat_t const& cryoGates) -> std::vector<ChannelID_t>
{

  //
  // get channels contributing to gates in a fired event
  //

  std::vector<ChannelID_t> channelList;
  for (auto const& gates: cryoGates) {
    for (auto const& gate: icarus::trigger::gatesIn(gates)) {

      if (gate.alwaysClosed()) continue;
      for (auto const channel: gate.channels()) {
        channelList.push_back(channel);
      }
    } // for gates
  } // for

  // remove duplicates
  std::sort(channelList.begin(), channelList.end());
  auto const firstDuplicate
    = std::unique(channelList.begin(), channelList.end());
  channelList.erase(firstDuplicate, channelList.end());
  return channelList;
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::extractActiveChannels()


//------------------------------------------------------------------------------
bool icarus::trigger::TriggerEfficiencyPlotsBase::deleteEmptyPlots
  (PlotSandbox& plots) const
{
  TDirectory* baseDir = plots.getDirectory();
  if (!baseDir) return true; // no content, nothing to do
  
  // our plots first
  unsigned int nEntries = 0U, nDirectories = 0U, nDeleted = 0U;
  for (TObject* obj: *(baseDir->GetList())) {
    
    // we do not deal with directories (except for subboxes below)
    if (dynamic_cast<TDirectory*>(obj)) {
      ++nDirectories;
      continue;
    }
    
    ++nEntries;
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (auto hist = dynamic_cast<TH1 const*>(obj)) {
      if (hist->GetEntries() > 0) continue;
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    else if (auto graph = dynamic_cast<TGraph const*>(obj)) {
      if (graph->GetN() > 0) continue;
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    else if (auto* eff = dynamic_cast<TEfficiency const*>(obj)) {
      auto const* hist = eff->GetTotalHistogram();
      if (hist && hist->GetEntries() > 0) continue;
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    else if (auto* tree = dynamic_cast<TTree const*>(obj)) {
      if (tree->GetEntries() > 0) continue;
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // add here more supported object types
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    else continue; // we don't delete unknown objects
    
    mf::LogTrace(fLogCategory)
      << "Deleting empty " << obj->IsA()->GetName() << "['" << obj->GetName()
      << "'] from " << plots.name();
    delete obj;
    ++nDeleted;
    
  } // for objects in directory
  
  // if we have found no more directories than the ones expected
  // from the subboxes and all the other entries have been deleted,
  // this box might be empty
  
  bool empty
    = (nDeleted == nEntries) && (nDirectories <= plots.nSubSandboxes());
  
  // we can't delete the sandboxes while iterating on them...
  std::vector<std::string> toBeDeleted;
  for (PlotSandbox& subbox: plots.subSandboxes()) {
    if (!deleteEmptyPlots(subbox)) continue;
    toBeDeleted.push_back(subbox.name());
    mf::LogTrace(fLogCategory)
      << "Scheduling empty " << plots.name() << "/" << toBeDeleted.back() << " for deletion";
  } // for subboxes
  if (toBeDeleted.size() != plots.nSubSandboxes()) empty = false;
  for (std::string const& subName: toBeDeleted) {
    if (!plots.deleteSubSandbox(subName)) continue;
    mf::LogTrace(fLogCategory)
      << "Deleted box " << plots.name() << "/" << subName;
  } // for
  
  return empty;
} // icarus::trigger::TriggerEfficiencyPlotsBase::deleteEmptyPlots()


//------------------------------------------------------------------------------
auto icarus::trigger::TriggerEfficiencyPlotsBase::makeChannelCryostatMap
  (geo::GeometryCore const& geom) -> std::vector<geo::CryostatID>
{
  
  auto const nOpChannels = geom.NOpChannels();
  
  std::vector<geo::CryostatID> channelCryostatMap(nOpChannels);
  
  for (auto const opChannel: util::counter(nOpChannels)) {
    if (!geom.IsValidOpChannel(opChannel)) continue;
    channelCryostatMap.at(opChannel)
      = geom.OpDetGeoFromOpChannel(opChannel).ID();
  } // for all channels
  
  return channelCryostatMap;
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::makeChannelCryostatMap()


//------------------------------------------------------------------------------
icarus::trigger::details::EventInfoExtractor::EDepTags_t
icarus::trigger::TriggerEfficiencyPlotsBase::makeEdepTag(
  fhicl::Sequence<art::InputTag> const& EnergyDepositTags,
  fhicl::OptionalAtom<art::InputTag> const& EnergyDepositSummaryTag
) {
  
  if (auto summaryTag = util::fhicl::getOptionalValue(EnergyDepositSummaryTag))
  {
    return {
      icarus::trigger::details::EventInfoExtractor::SimEnergyDepositSummaryInputTag
        { *summaryTag }
      };
  }
  else {
    
    return { EnergyDepositTags() };
    
  }
  
} // icarus::trigger::MakeTriggerSimulationTree::makeEdepTag()


//------------------------------------------------------------------------------
