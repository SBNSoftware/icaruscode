/**
 * @file   PMTWaveformBaselines_module.cc
 * @brief  Extracts and writes PMT waveform baselines.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 14, 2020
 */

// ICARUS libraries
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"
// #include "icaruscode/Utilities/DataProductPointerMap.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/Utilities/StatCollector.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
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
#include "TGraph.h"
#include "TProfile.h"

// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::partial_sort_copy()
#include <iterator> // std::distance()
#include <memory> // std::make_unique(), std::allocator
#include <string>
#include <cmath> // std::ceil()
#include <cassert>

//------------------------------------------------------------------------------
namespace icarus { class PMTWaveformBaselines; }

/**
 * @brief Extracts the baseline of PMT waveforms.
 * 
 * This module produces a baseline data product for each optical detector
 * waveform.
 * 
 * The waveforms on the same channels are currently treated as independent
 * (which is less than ideal).
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
 * * `Baselines`: baseline distribution, per channel; average baseline [ADC] per
 *     event per channel; all waveforms on the same channels in a single event
 *     contribute to the average, and channels with no waveforms in an event do
 *     not contribute an entry for that event.
 * 
 * 
 * 
 * Input data products
 * ====================
 * 
 * * `std::vector<raw::OpDetWaveform>`: a single waveform for each recorded
 *      optical detector activity; the activity belongs to a single channel, but
 *      there may be multiple waveforms on the same channel.
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
 * `lar --print-description PMTWaveformBaselines`.
 * 
 * * `OpticalWaveforms` (input tag, mandatory): the data product containing all
 *   optical detector waveforms;
 * * `OutputCategory` (string, default: `"PMTWaveformBaselines"`): label
 *   for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration;
 * * `PlotBaselines` (flag, default: `true`): whether to produce distributions
 *   of the extracted baselines.
 * * `BaselineTimeAverage` (real number, default: `600.0`): binning of the
 *   baseline profile vs. time, in seconds. Requires `PlotBaselines` to be set.
 * 
 */
class icarus::PMTWaveformBaselines: public art::EDProducer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> OpticalWaveforms{
      Name("OpticalWaveforms"),
      Comment("label of input digitized optical waveform data product")
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name("OutputCategory"),
      Comment("tag of the module output to console via message facility"),
      "PMTWaveformBaselines"
      };
    
    fhicl::Atom<bool> PlotBaselines {
      Name("PlotBaselines"),
      Comment("produce plots on the extracted baseline"),
      true
      };
    
    fhicl::Atom<double> BaselineTimeAverage {
      Name("BaselineTimeAverage"),
      Comment("binning of the baseline profile vs. time [s]"),
      [this](){ return PlotBaselines(); }, // enabled if `PlotBaselines` is set
      600.0
      };
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit PMTWaveformBaselines(Parameters const& config);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  /// Prepares the plots to be filled.
  virtual void beginJob() override;
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  /// Remove empty plots.
  virtual void endJob() override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fOpDetWaveformTag; ///< Input optical waveform tag.
  
  bool fPlotBaselines; ///< Whether to produce plots.
  
  /// Width of baseline time profile binning [s]
  double const fBaselineTimeAverage { 0.0 };
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  // --- END Algorithms --------------------------------------------------------
  
  std::size_t fNPlotChannels = 0U; ///< Number of plotted channels
  TH2* fHBaselines = nullptr; ///< All baselines, per channel.
  
  /// For each channel, all event times and their baselines.
  std::vector<std::vector<std::pair<double, double>>> fBaselinesVsTime;
  
  
  /// Creates all the plots to be filled by the module.
  void setupPlots();
  
  /// Removes the empty plots.
  void buildBaselineGraphs();
  
  /// Extracts a baseline as median from a single waveform.
  icarus::WaveformBaseline baselineFromMedian
    (raw::OpDetWaveform const& waveform) const;
  
}; // icarus::PMTWaveformBaselines



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
//--- icarus::PMTWaveformBaselines
//------------------------------------------------------------------------------
icarus::PMTWaveformBaselines::PMTWaveformBaselines
  (Parameters const& config)
  : art::EDProducer(config)
  // configuration
  , fOpDetWaveformTag(config().OpticalWaveforms())
  , fPlotBaselines(config().PlotBaselines())
  , fBaselineTimeAverage(config().BaselineTimeAverage())
  , fLogCategory(config().OutputCategory())
{
  //
  // optional configuration parameters
  //
  
  //
  // configuration report (currently, more like a placeholder)
  //
  mf::LogInfo(fLogCategory)
    << "Using the standard (median) algorithm, waveform by waveform, on '"
    << fOpDetWaveformTag.encode() << "'.";
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  
  //
  // declaration of output
  //
  produces<std::vector<icarus::WaveformBaseline>>();
  produces<art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform>>();
  
} // icarus::PMTWaveformBaselines::PMTWaveformBaselines()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselines::beginJob() {
  
  //
  // set up the plots, if needed
  //
  if (fPlotBaselines) setupPlots();
  
} // icarus::PMTWaveformBaselines::beginJob()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselines::produce(art::Event& event) {
  
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
  
  std::vector<lar::util::StatCollector<double>> averages;
  if (fHBaselines || !fBaselinesVsTime.empty())
    averages.resize(fNPlotChannels);
  
  art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform> baselineToWaveforms;
  
  art::PtrMaker<icarus::WaveformBaseline> const makeBaselinePtr(event);
  
  for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
    assert(iWaveform == baselines.size());
    
    icarus::WaveformBaseline const baseline { baselineFromMedian(waveform) };
    
    if (!averages.empty())
      averages[waveform.ChannelNumber()].add(baseline.baseline());
    baselines.push_back(baseline);
    
    baselineToWaveforms.addSingle(
      makeBaselinePtr(iWaveform),
      art::Ptr<raw::OpDetWaveform>(waveformHandle, iWaveform)
      );
    
  } // for waveforms
  
  //
  // plot filling
  //
  if (!averages.empty()) {
    
    double const eventTime = static_cast<double>(event.time().timeHigh())
      + static_cast<double>(event.time().timeHigh()) * 1e-9;
    
    for (auto const& [ channel, stat ]: util::enumerate(averages)) {
      if (stat.N() == 0) continue;
      
      double const aveBaseline = stat.Average();
      
      fHBaselines->Fill(double(channel), aveBaseline);
      
      if (channel < fBaselinesVsTime.size())
        fBaselinesVsTime[channel].emplace_back(eventTime, aveBaseline);
      
    } // for baselines
    
  } // if plots
  
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
  
  
} // icarus::PMTWaveformBaselines::produce()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselines::endJob() {
  
  if (fPlotBaselines) buildBaselineGraphs();
  
} // icarus::PMTWaveformBaselines::endJob()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselines::setupPlots() {
  
  auto const& tfs = *(art::ServiceHandle<art::TFileService>());
  auto const& geom = *(lar::providerFrom<geo::Geometry const>());
  
  fNPlotChannels = geom.NOpChannels();
  
  fHBaselines = tfs.make<TH2F>(
    "Baselines",
    "PMT baseline;channel;baseline per channel [ / 8 ADC ]",
    fNPlotChannels, 0.0, double(fNPlotChannels),
    256, 13312.0, 15360.0
    );
  
  // these are graphs, and it is more convenient to carry around their data
  // in a vector than carrying around the graphs themselves;
  // `buildBaselineGraphs()` will turn that data into graph at end of the job;
  // here we just declare which channels we are going to plot.
  fBaselinesVsTime.resize(fNPlotChannels);
  
} // icarus::PMTWaveformBaselines::setupPlots()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselines::buildBaselineGraphs() {
  
  auto& tfs = *(art::ServiceHandle<art::TFileService>());
  
  art::TFileDirectory graphDir = tfs.mkdir("graphs", "Baseline vs. time");
  art::TFileDirectory profileDir
    = tfs.mkdir("profiles", "Baseline profiles vs. time");
  
  for (auto const channel: util::counter(fBaselinesVsTime.size())) {
    
    auto& timeAndBaselines = fBaselinesVsTime[channel];
    if (timeAndBaselines.empty()) continue;
    
    // sort by time (entries with the same time would be sorted by baseline,
    // but that does not really happen nor it mattered if it happened)
    std::sort(timeAndBaselines.begin(), timeAndBaselines.end());
    
    // graph, one point per event
    auto* const graph = graphDir.makeAndRegister<TGraph>(
      "BaselineCh" + std::to_string(channel),
      "PMT channel #" + std::to_string(channel) + ": baseline vs. run time",
      timeAndBaselines.size()
      );
    assert(graph->GetXaxis());
    assert(graph->GetYaxis());
    graph->GetXaxis()->SetTitle("event time");
    graph->GetYaxis()->SetTitle("baseline  [ ADC ]");
    
    
    // profile, one point every 10 minutes (or best offer)
    double const startTime = timeAndBaselines.front().first;
    double const totalTime = timeAndBaselines.back().first - startTime;
    auto const nBins = std::max(1U,
      static_cast<unsigned int>(std::ceil(totalTime / fBaselineTimeAverage))
      );
    double const endTime = startTime + nBins * fBaselineTimeAverage;
    
    auto* const profile = profileDir.make<TProfile>(
      ("BaselineCh" + std::to_string(channel) + "profile").c_str(),
      ("PMT channel #" + std::to_string(channel) + ": baseline vs. run time")
        .c_str(),
      nBins, startTime, endTime
      );
    assert(profile->GetXaxis());
    assert(profile->GetYaxis());
    profile->GetXaxis()->SetTitle(
      ("event time  [ / " + std::to_string(fBaselineTimeAverage) + " s ]")
      .c_str()
      );
    profile->GetYaxis()->SetTitle("average baseline  [ ADC ]");
    
    for (auto const& [ i, data ]: util::enumerate(timeAndBaselines)) {
      auto const [ time, baseline ] = data;
      graph->SetPoint(i, time, baseline);
      profile->Fill(time, baseline);
    } // for
    
  } // for
  
} // icarus::PMTWaveformBaselines::buildBaselineGraphs()


//------------------------------------------------------------------------------
icarus::WaveformBaseline icarus::PMTWaveformBaselines::baselineFromMedian
  (raw::OpDetWaveform const& waveform) const
{
  
  return icarus::WaveformBaseline
    { waveform.empty()? 0.0f: median(waveform.begin(), waveform.end()) };
  
} // icarus::PMTWaveformBaselines::baselineFromMedian()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::PMTWaveformBaselines)


//------------------------------------------------------------------------------
