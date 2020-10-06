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
#include "icaruscode/PMT/Trigger/Utilities/ROOTutils.h" // util::ROOT
#include "icarusalg/Utilities/DetectorClocksHelpers.h" // makeDetTimings()

// LArSoft libraries
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
  , fLogCategory          (config.LogCategory())
  // services
  , fGeom      (*lar::providerFrom<geo::Geometry>())
  , fOutputDir (*art::ServiceHandle<art::TFileService>())
  // cached
  , fEventInfoExtractorMaker(
      config.GeneratorTags(),              // truthTags
      config.EnergyDepositTags(),          // edepTags
      fGeom,                               // geom
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
  for (raw::ADC_Count_t threshold: config.Thresholds()) {
    fADCthresholds[icarus::trigger::ADCCounts_t{threshold}]
      = art::InputTag{ discrModuleLabel, util::to_string(threshold) };
  }

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
    log << "\nConfigured " << fADCthresholds.size() << " thresholds:";
    for (auto const& [ threshold, dataTag ]: fADCthresholds)
      log << "\n * " << threshold << " ADC (from '" << dataTag.encode() << "')";
    
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
    mf::LogWarning(fLogCategory)
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

    auto const& [ threshold, dataTag ] = thrPair;

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
  
  mf::LogInfo(fLogCategory)
    << nPlottedEvents << "/" << nEvents << " events plotted."
    ;
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::printSummary()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::initializePlots
  (PlotCategories_t categories, std::vector<SettingsInfo_t> const& settings)
{
  using namespace std::string_literals;
  
  auto const beamGate = icarus::trigger::makeBeamGateStruct
    (icarus::ns::util::makeDetTimings(), fBeamGateDuration, fBeamGateStart);

  fBeamGateChangeCheck(beamGate);
  
  {
    mf::LogTrace(fLogCategory)
      << "Beam gate:"
      << "\n - electronics time: " << beamGate.asGate()
      << "\n - simulation time: " << beamGate.asSimulationRange()
      << "\n - optical ticks: " << beamGate.asOptTickRange()
      ;
    
    mf::LogInfo log(fLogCategory);
    log << "\nConfigured " << fADCthresholds.size() << " thresholds:";
    for (auto const& [ threshold, dataTag ]: fADCthresholds)
      log << "\n * " << threshold << " ADC (from '" << dataTag.encode() << "')";
    log << "\nBeam gate for plots is " << beamGate.asSimulationRange();
    
  } // local block
  
  
  fPlotCategories = std::move(categories);
  
  for (icarus::trigger::ADCCounts_t const threshold
    : util::get_elements<0U>(fADCthresholds))
  {
    // create a plot sandbox inside `fOutputDir` with a name/prefix `Thr###`
    auto const thr = threshold.value();
    icarus::trigger::PlotSandbox thrPlots { fOutputDir,
      "Thr"s + util::to_string(thr), "(thr: "s + util::to_string(thr) + ")"s };
    
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
  } // if generated information
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::initializeEventPlots()


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
    assert(eventInfo.hasGenerated());
    if (eventInfo.isNeutrino()) {
      getTrig.Hist("NeutrinoEnergy"s).Fill(double(eventInfo.NeutrinoEnergy()));
      getTrig.Hist("InteractionType"s).Fill(eventInfo.InteractionType());
      getTrig.Hist("LeptonEnergy"s).Fill(double(eventInfo.LeptonEnergy()));
    } // if neutrino event
    TH2& vertexHist = getTrig.Hist2D("InteractionVertexYZ"s);
    for (auto const& point: eventInfo.Vertices())
      vertexHist.Fill(point.Z(), point.Y());
  } // if use generated information
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::fillEventPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::fillEfficiencyPlots(
  EventInfo_t const& eventInfo,
  TriggerInfo_t const& triggerInfo,
  PlotSandbox const& plots
) const {
  
  using namespace std::string_literals;
  
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
    getTrigEff.Hist("TriggerTick"s).Fill(triggerInfo.atTick().value());
  }
  
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::fillEfficiencyPlots()


//------------------------------------------------------------------------------
void icarus::trigger::TriggerEfficiencyPlotsBase::fillAllEfficiencyPlots(
  EventInfo_t const& eventInfo,
  TriggerInfo_t const& triggerInfo,
  PlotSandbox const& plots
) const {
  
  fillEfficiencyPlots(eventInfo, triggerInfo, plots);
  
  // plotting split for triggering/not triggering events
  fillEventPlots(
    eventInfo,
    plots.demandSandbox(triggerInfo.fired()? "triggering": "nontriggering")
    );
  
} // icarus::trigger::TriggerEfficiencyPlotsBase::fillAllEfficiencyPlots()


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
    = *(event.getValidHandle<std::vector<OpticalTriggerGateData_t>>(dataTag));
  auto const& gateToWaveforms = *(
    event.getValidHandle
      <art::Assns<OpticalTriggerGateData_t, raw::OpDetWaveform>>(dataTag)
    );
  
  try {
    return icarus::trigger::FillTriggerGates<InputTriggerGate_t>
      (gates, gateToWaveforms);
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
    assert(gate.hasChannels());
    gatesPerCryostat[fChannelCryostat.at(gate.channels().front()).Cryostat]
      .push_back(std::move(gate));
  } // for gates
  
  return gatesPerCryostat;

} // icarus::trigger::TriggerEfficiencyPlotsBase::splitByCryostat()


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
