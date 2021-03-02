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
#include "lardataalg/Utilities/StatCollector.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "art_root_io/TFileService.h"
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

// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::partial_sort_copy()
#include <iterator> // std::distance()
#include <memory> // std::make_unique(), std::allocator
#include <string>
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
    
  }; // struct Config
  
  using Parameters = art::EDProducer::Table<Config>;
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  explicit PMTWaveformBaselines(Parameters const& config);
  
  // Plugins should not be copied or assigned.
  PMTWaveformBaselines(PMTWaveformBaselines const&) = delete;
  PMTWaveformBaselines(PMTWaveformBaselines&&) = delete;
  PMTWaveformBaselines& operator=(PMTWaveformBaselines const&) = delete;
  PMTWaveformBaselines& operator=(PMTWaveformBaselines&&) = delete;
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  /// Prepares the plots to be filled.
  virtual void beginJob() override;
  
  /// Creates the data products.
  virtual void produce(art::Event& event) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fOpDetWaveformTag; ///< Input optical waveform tag.
  bool fPlotBaselines; ///< Whether to produce plots.
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  // --- END Algorithms --------------------------------------------------------
  
  TH2* fHBaselines = nullptr; ///< All baselines, per channel.
  
  
  /// Creates all the plots to be filled by the module.
  void setupPlots();
  
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
  if (fHBaselines) averages.resize(fHBaselines->GetNbinsX());
  
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
    
    for (auto const& [ channel, stat ]: util::enumerate(averages))
      if (stat.N() > 0) fHBaselines->Fill(double(channel), stat.Average());
    
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
void icarus::PMTWaveformBaselines::setupPlots() {
  
  auto const& tfs = *(art::ServiceHandle<art::TFileService>());
  auto const& geom = *(lar::providerFrom<geo::Geometry const>());
  
  unsigned int const nChannels = geom.NOpChannels();
  
  fHBaselines = tfs.make<TH2F>(
    "Baselines",
    "PMT baseline;channel;baseline per channel [ / 8 ADC ]",
    nChannels, 0.0, double(nChannels),
    256, 13312.0, 15360.0
    );
  
} // icarus::PMTWaveformBaselines::setupPlots()


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
