/**
 * @file   PMTWaveformBaselinesFromReadoutConfiguration_module.cc
 * @brief  Extracts PMT channel baselines from PMT readout configuration.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 24, 2021
 */

// ICARUS libraries
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h"
#include "sbnobj/Common/PMT/Data/V1730Configuration.h"
#include "sbnobj/Common/PMT/Data/V1730channelConfiguration.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"
// #include "icaruscode/Utilities/DataProductPointerMap.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TH2F.h"

// C/C++ standard libraries
#include <vector>
#include <memory> // std::make_unique()
#include <string>
#include <tuple> // std::tie()
#include <utility> // std::move(), std::pair<>
#include <limits> // std::numeric_limits<>
#include <cassert>

//------------------------------------------------------------------------------
namespace icarus { class PMTWaveformBaselinesFromReadoutConfiguration; }

/**
 * @brief Extracts PMT baseline settings from PMT readout configuration.
 * 
 * This module produces a baseline data product for each optical detector
 * waveform. The content is the same for all the events in each run, and it is
 * read from the specified PMT configuration data product.
 * 
 * Each waveform is associated with the baseline configured for its channel.
 * If there is no information for a given channel and that channel is requested,
 * an exception is thrown.
 * 
 * This module is interchangeable with `PMTWaveformBaselines`, which instead
 * extracts the baseline dynamically from each waveform.
 * 
 * 
 * Output data products
 * =====================
 * 
 * * a data product of type `std::vector<icarus::WaveformBaseline>`, with one
 *   baseline per input waveform; the baselines are guaranteed to be in the same
 *   order as the waveforms in the input collection;
 * * an _art_ association, one baseline per input optical waveform; normally
 *   this is not needed as long as the original PMT waveform data product is
 *   available; the type of the association is
 *   `art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform>`.
 * 
 * 
 * Output plots
 * -------------
 * 
 * * `Baselines`: baseline distribution, per channel; one entry [ADC]
 *     per configured channel per run.
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<raw::OpDetWaveform>`: a single waveform for each recorded
 *      optical detector activity; the activity belongs to a single channel, but
 *      there may be multiple waveforms on the same channel.
 * 
 * @note This module needs only the information of the channel of each waveform;
 *       because _art_ is not capable of selective readout, though, the whole
 *       optical detector readout data is loaded. When this module is the first
 *       one in the job accessing the PMT waveforms, the data loading time is
 *       accounted to it; the time the module needs to perform its specific
 *       operations is tiny in comparison. Anyway, since it is likely that
 *       another module in the job will then use PMT waveform data, this loading
 *       time is not wasted.
 * 
 * 
 * Service requirements
 * ---------------------
 * 
 * * `TFileService` and `Geometry` if `PlotBaselines` is enabled
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description PMTWaveformBaselinesFromReadoutConfiguration`.
 * 
 * * `OpticalWaveforms` (input tag, mandatory): the data product containing all
 *   optical detector waveforms.
 * * `PMTconfigurationTag` (input tag, mandatory): the run-level data product
 *   containing the full PMT readout configuration (including baseline per
 *   channel).
 * * `OutputCategory` (string, default:
 *   `"PMTWaveformBaselinesFromReadoutConfiguration"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 * * `PlotBaselines` (flag, default: `true`): whether to produce distributions
 *   of the configured baselines (_not implemented yet_).
 * * `PrintBaselines` (flag, default: `true`): if set to `true`, on each run it
 *   will print the baseline of all the configured channels with a LArSoft
 *   channel ID.
 * 
 */
class icarus::PMTWaveformBaselinesFromReadoutConfiguration
  : public art::EDProducer
{
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> OpticalWaveforms {
      Name("OpticalWaveforms"),
      Comment("label of input digitized optical waveform data product")
      };
    
    fhicl::Atom<art::InputTag> PMTconfigurationTag {
      Name("PMTconfigurationTag"),
      Comment("label of PMT readout configuration run-level data product")
      };
    
    fhicl::Atom<bool> PlotBaselines {
      Name("PlotBaselines"),
      Comment("produce plots on the extracted baseline"),
      true
      };
    
    fhicl::Atom<bool> PrintBaselines {
      Name("PrintBaselines"),
      Comment("prints baseline values for each configured channel on each run"),
      false
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("tag of the module output to console via message facility"),
      "PMTWaveformBaselinesFromReadoutConfiguration"
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  explicit PMTWaveformBaselinesFromReadoutConfiguration
    (Parameters const& config);
  
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  /// Prepares the plots to be filled.
  virtual void beginJob() override;
  
  /// Reads the PMT readout configuration.
  virtual void beginRun(art::Run& run) override;
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  /// Type for baseline (same as V1730channelConfiguration::baseline)
  using Baseline_t = signed short int;
  
  /// Mnemonic value for channels without baseline information.
  static constexpr Baseline_t NoBaseline
    = std::numeric_limits<Baseline_t>::min();
  
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fOpDetWaveformTag; ///< Input optical waveform tag.
  art::InputTag const fPMTconfigurationTag; ///< Input PMT readout config tag.
  bool const fPlotBaselines; ///< Whether to produce plots.
  bool const fPrintBaselines; ///< Whether to print baselines on each run.
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  // --- BEGIN Service variables -----------------------------------------------
  
  // --- END Service variables -------------------------------------------------
  
  unsigned int fConfigured = 0U; ///< Number of channels in PMT configuration.
  
  /// PMT baselines configured in the current run, indexed by channel.
  std::vector<short signed int> fBaselines;
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  // --- END Algorithms --------------------------------------------------------
  
  TH2* fHBaselines = nullptr; ///< All baselines, per channel.
  
  
  /// Creates all the plots to be filled by the module.
  void setupPlots();
  
  
  /// Returns the configured baseline.
  /// @throw cet::exception (category: "PMTWaveformBaselinesFromReadoutConfiguration")
  ///        if `channel` is not present in the configuration
  Baseline_t getBaseline(raw::Channel_t channel) const;
  
  /// Returns the number of channels in configuration and a baseline map.
  std::pair<unsigned int, std::vector<Baseline_t>>
  extractBaselinesFromConfiguration
    (sbn::PMTconfiguration const& PMTconfig) const;
  
  /// Returns the number of channels currently configured with a baseline.
  unsigned int nChannelsWithBaseline() const;
  
  /// Prints the current baselines on maesage facility (INFO level).
  void printBaselines() const;
  
}; // icarus::PMTWaveformBaselinesFromReadoutConfiguration



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Extracts the median of the collection between the specified iterators.
  template <typename BIter, typename EIter>
  auto median(BIter begin, EIter end) {
    
    using value_type = typename BIter::value_type;
    
    std::vector<value_type> data{ begin, end };
    assert(!data.empty());
    
    auto const middle = data.begin() + data.size() / 2;
    std::nth_element(data.begin(), middle, data.end());
    
    return *middle;
    
  } // median()
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::PMTWaveformBaselinesFromReadoutConfiguration
//------------------------------------------------------------------------------
icarus::PMTWaveformBaselinesFromReadoutConfiguration::PMTWaveformBaselinesFromReadoutConfiguration
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fOpDetWaveformTag(config().OpticalWaveforms())
  , fPMTconfigurationTag(config().PMTconfigurationTag())
  , fPlotBaselines(config().PlotBaselines())
  , fPrintBaselines(config().PrintBaselines())
  , fLogCategory(config().OutputCategory())
{
  //
  // optional configuration parameters
  //
  
  //
  // configuration report (currently, more like a placeholder)
  //
  mf::LogInfo(fLogCategory)
    << "Using configured baselines from '" << fPMTconfigurationTag.encode()
    << "', waveform by waveform, on '" << fOpDetWaveformTag.encode() << "'.";
  
  //
  // declaration of input
  //
  consumes<sbn::PMTconfiguration, art::InRun>(fPMTconfigurationTag);
  consumes<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  
  //
  // declaration of output
  //
  produces<std::vector<icarus::WaveformBaseline>>();
  produces<art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform>>();
  
} // icarus::PMTWaveformBaselinesFromReadoutConfiguration::PMTWaveformBaselinesFromReadoutConfiguration()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromReadoutConfiguration::beginJob() {
  
  //
  // set up the plots, if needed
  //
  if (fPlotBaselines) setupPlots();
  
} // icarus::PMTWaveformBaselinesFromReadoutConfiguration::beginJob()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromReadoutConfiguration::beginRun
  (art::Run& run)
{
  auto const& PMTconfig
    = run.getProduct<sbn::PMTconfiguration>(fPMTconfigurationTag);
  
  std::vector<Baseline_t> newBaselines;
  std::tie(fConfigured, newBaselines)
    = extractBaselinesFromConfiguration(PMTconfig);
  
  bool const changed = (fBaselines != newBaselines);
  fBaselines = std::move(newBaselines);
  
  if (fHBaselines) {
    for (auto const& [ channel, baseline ]: util::enumerate(fBaselines)) {
      if (baseline != NoBaseline)
        fHBaselines->Fill(double(channel), double(baseline));
    } // for
  } // if plots
  
  if (fPrintBaselines && changed) printBaselines();
  
} // icarus::PMTWaveformBaselinesFromReadoutConfiguration::beginRun()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromReadoutConfiguration::produce
  (art::Event& event)
{
  
  //
  // fetch input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  auto const& waveforms = *waveformHandle;
  
  /*
   * this may be needed in a future where processing happens per channel
   * rather than per waveform:
  // map address of waveform to art pointer to that waveform
  auto const& opDetWavePtrs
    = util::mapDataProductPointers(event, waveformHandle);
  */
  
  
  //
  // compute all the baselines
  //
  std::vector<icarus::WaveformBaseline> baselines;
  baselines.reserve(waveforms.size());
  
  art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform> baselineToWaveforms;
  
  art::PtrMaker<icarus::WaveformBaseline> const makeBaselinePtr(event);
  
  for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
    assert(iWaveform == baselines.size());
    
    icarus::WaveformBaseline const baseline
      { static_cast<float>(getBaseline(waveform.ChannelNumber())) };
    
    baselines.push_back(baseline);
    
    baselineToWaveforms.addSingle(
      makeBaselinePtr(iWaveform),
      art::Ptr<raw::OpDetWaveform>(waveformHandle, iWaveform)
      );
    
  } // for waveforms
  
  
  //
  // output
  //
  event.put(
    std::make_unique<std::vector<icarus::WaveformBaseline>>
      (std::move(baselines))
    );
  event.put(
    std::make_unique<art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform>>
      (std::move(baselineToWaveforms))
    );
  
  
} // icarus::PMTWaveformBaselinesFromReadoutConfiguration::produce()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromReadoutConfiguration::setupPlots() {
  
  auto const& tfs = *(art::ServiceHandle<art::TFileService>());
  auto const& geom = *(lar::providerFrom<geo::Geometry const>());
  
  unsigned int const nChannels = geom.NOpChannels();
  
  fHBaselines = tfs.make<TH2F>(
    "Baselines",
    "PMT baseline;channel;baseline per channel [ / 8 ADC ]",
    nChannels, 0.0, double(nChannels),
    256, 13312.0, 15360.0
    );
  
} // icarus::PMTWaveformBaselinesFromReadoutConfiguration::setupPlots()


//------------------------------------------------------------------------------
auto icarus::PMTWaveformBaselinesFromReadoutConfiguration::getBaseline
  (raw::Channel_t channel) const -> Baseline_t
{
  
  auto const channelSlot = static_cast<std::size_t>(channel);
  if (channelSlot < fBaselines.size()) {
    Baseline_t const baseline = fBaselines[channelSlot];
    if (baseline != NoBaseline) return baseline;
  }
  throw cet::exception("PMTWaveformBaselinesFromReadoutConfiguration")
    << "No baseline configured for channel " << channel << ".\n";
  
} // icarus::PMTWaveformBaselinesFromReadoutConfiguration::getBaseline()



//------------------------------------------------------------------------------
auto
icarus::PMTWaveformBaselinesFromReadoutConfiguration::extractBaselinesFromConfiguration
  (sbn::PMTconfiguration const& PMTconfig) const
  -> std::pair<unsigned int, std::vector<Baseline_t>>
{
  constexpr auto NoChannelID = sbn::V1730channelConfiguration::NoChannelID;
  
  unsigned int nConfigured = 0U;
  std::vector<Baseline_t> baselines;
  
  for (sbn::V1730Configuration const& boardConfig: PMTconfig.boards) {
    for (sbn::V1730channelConfiguration const& channelConfig
      : boardConfig.channels
    ) {
      
      ++nConfigured;
      
      raw::Channel_t const channelID = channelConfig.channelID;
      if (channelID == NoChannelID) continue;
      
      auto const channelSlot = static_cast<std::size_t>(channelID);
      if (baselines.size() <= channelSlot)
        baselines.resize(channelSlot + 1U, NoBaseline);
      
      Baseline_t const baseline = channelConfig.baseline;
      if (baseline == NoBaseline) {
        // if this exception is thrown, it means we need a different way than
        // a special baseline value (`NoBaseline`) to track the missing channels
        // (e.g. a map or a std::vector<std::optional<...>>)
        throw art::Exception(art::errors::LogicError)
          << "A configured baseline actually matches the special value "
          << "`NoBaseline` (" << NoBaseline << ")\n";
      } // if bad luck
      
      baselines[channelSlot] = baseline;
      
    } // for all channels on a board
  } // for all boards
  
  return { nConfigured, std::move(baselines) };
} // icarus::PMTWaveformBaselinesFromReadoutConfiguration::extractBaselinesFromConfiguration()



//------------------------------------------------------------------------------
unsigned int
icarus::PMTWaveformBaselinesFromReadoutConfiguration::nChannelsWithBaseline
  () const
{
  return fBaselines.size()
    - std::count(fBaselines.begin(), fBaselines.end(), NoBaseline);
}

//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromReadoutConfiguration::printBaselines()
  const
{
  mf::LogInfo log { fLogCategory };
  
  unsigned int const nChannels = nChannelsWithBaseline();
  assert(nChannels <= fConfigured);
  unsigned int const nMissing = fConfigured - nChannels;
  
  log << "Configuration has " << nChannels << " channels";
  if (nMissing > 0U)
    log << " (plus " << nMissing << " not associated to any channel ID)";
  log << ":";
  for (auto [ channelID, baseline ]: util::enumerate(fBaselines)) {
    if (baseline == NoBaseline) continue;
    log << "\n  channel " << channelID << ": " << baseline;
  } // for all channels
  
} // icarus::PMTWaveformBaselinesFromReadoutConfiguration::printBaselines()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::PMTWaveformBaselinesFromReadoutConfiguration)


//------------------------------------------------------------------------------
