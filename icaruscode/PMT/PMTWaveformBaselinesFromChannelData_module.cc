/**
 * @file   PMTWaveformBaselinesFromChannelData_module.cc
 * @brief  Extracts and writes PMT waveform baselines.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 14, 2020
 */

// ICARUS libraries
#include "icaruscode/PMT/Data/WaveformRMS.h"
#include "icarusalg/PMT/Algorithms/SharedWaveformBaseline.h"
#include "icarusalg/Utilities/GroupByIndex.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"
#include "sbnobj/Common/PMT/Data/PMTconfiguration.h"
#include "sbnobj/Common/PMT/Data/V1730Configuration.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // electronics_time
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataalg/Utilities/StatCollector.h" // lar::util::MinMaxCollector
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ProcessingFrame.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"

// ROOT libraries
#include "TH2F.h"
#include "TGraph.h"
#include "TProfile.h"

// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::remove_if(), std::sort()
#include <iterator> // std::distance()
#include <memory> // std::make_unique()
#include <utility> // std::pair, std::move()
#include <string>
#include <optional>
#include <limits>
#include <cmath> // std::round(), std::max(), std::ceil()
#include <cstdint> // std::size_t
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus { class PMTWaveformBaselinesFromChannelData; }

/**
 * @class icarus::PMTWaveformBaselinesFromChannelData
 * @brief Extracts a baseline from PMT waveforms.
 * 
 * This module produces a baseline data product for each optical detector
 * waveform.
 * The algorithm extracts a baseline per channel per event, considering all the
 * waveforms from one channel together, under the assumptions that:
 * 
 * 1. the baseline is not going to change during the few milliseconds of readout
 * 2. the beginning of the waveform is on its baseline level
 * 
 * 
 * The algorithm
 * ==============
 * 
 * The core of the algorithm is described in its own class,
 * `opdet::SharedWaveformBaseline`.
 * 
 * On each event independently, the waveforms are grouped by channel.
 * For each waveform, `opdet::SharedWaveformBaseline` considers their first part
 * for baseline calculation. The number of samples of that part is specified as
 * a fraction of the pre-trigger buffer. The pre-trigger buffer includes the
 * samples that were collected before the physics activity that causes the data
 * acquisition happened, and as such is expected to be almost completely free
 * of physics signal and to be made of just electronics noise. The size of this
 * prebuffer is determined in data by the readout configuration, which can be
 * accessed by specifying the parameter `PMTconfigurationTag`, and in simulation
 * by a digitization module parameter, which can be replicated here by the
 * parameter `PretriggerBufferSize`. The fraction of such buffer used for the
 * baseline estimation is also specified by a configuration parameter
 * (`PretriggerBufferFractionForBaseline`). Depending on the value of the
 * configuration parameter `ExcludeSpillTimeIfMoreThan`, the one waveform that
 * covers the trigger time may be excluded and not used; the rationale is that
 * there is a conspicuous number of events triggered by the late light of a
 * cosmic ray happening before the beam spill, in which case the activity may
 * contaminate the pre-trigger buffer and could bias the estimation of the
 * baseline.
 * 
 * The result of the `opdet::SharedWaveformBaseline` is currently used directly
 * as the baseline for all the waveforms on the channel on that event.
 * 
 * @todo Add run-level information and checks.
 * 
 * 
 * Output data products
 * =====================
 * 
 * * a data product of type `std::vector<icarus::WaveformBaseline>`, with one
 *   baseline per input waveform; the baselines are guaranteed to be in the same
 *   order as the waveforms in the input collection;
 * * a data product of type `std::vector<icarus::WaveformRMS>`, with one
 *   baseline RMS per input waveform; the baselines are guaranteed to be in the
 *   same order as the waveforms in the input collection;
 * * two _art_ associations, with one baseline and one RMS per input optical
 *   waveform; normally this is not needed as long as the original PMT waveform
 *   data product is available; the type of the associations are
 *   `art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform>` and
 *   `art::Assns<icarus::WaveformRMS, raw::OpDetWaveform>`.
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
 * * `DetectorClocksService` for the determination of the optical tick duration
 *     and the trigger time (if the exclusion of waveform with trigger is
 *     enabled).
 * * `TFileService` and `Geometry` if `PlotBaselines` is enabled.
 * 
 * 
 * Configuration parameters
 * =========================
 * 
 * A terse description of the parameters is printed by running
 * `lar --print-description PMTWaveformBaselinesFromChannelData`.
 * 
 * * `OpticalWaveforms` (input tag, mandatory): the data product containing all
 *   optical detector waveforms;
 * * `PretriggerBufferSize` (positive integer, optional): the number of samples
 *     collected by the PMT readout before the PMT signal. This number is used
 *     as base for the size of the initial part of the waveform to be used in
 *     baseline calculation. This and `PMTconfigurationTag` parameters are
 *     exclusive. This one is expected to be preferred for simulated samples.
 * * `PMTconfigurationTag` (input tag): the (run) data product containing the
 *     PMT readout configuration. This is extracted from the run configuration
 *     for data events. The pre-trigger buffer size is read from this data
 *     product, and used as base for the size of the initial part of the
 *     waveform to be used in baseline calculation. This and
 *     `PretriggerBufferSize` parameters are exclusive. This one is expected
 *     to be preferred for data samples.
 * * `PretriggerBufferFractionForBaseline` (real, default: `0.5`): the fraction
 *     of the pre-trigger buffer (see `PretriggerBufferSize` and
 *     `PMTconfigurationTag` parameters) to be used to extract the baseline.
 * * `ExcludeSpillTimeIfMoreThan` (integer, default: disable): the minimum
 *     number of PMT waveforms on a single channel of a single event, in order
 *     for the exclusion of the waveform containing the global trigger to be
 *     enabled. By default (huge number) the feature is disabled.
 * * `AlgoParams` (table): configuration of the core algorithm extracting the
 *     baseline from a set of waveforms on the same channel. The configuration
 *     is as follows:
 *     * `AcceptedSampleRangeRMS` (real): for a waveform to be considered for
 *         the baseline, the values of the samples must stay within a certain
 *         range, which is defined as this number of RMS away in either
 *         directions from a central value that is a rougher estimation of the
 *         baseline.
 *     * `ExcessSampleLimit` (integer): for a waveform to be considered for
 *         the baseline, there must be less than this number of samples in a
 *         row that are outside of the range defined by
 *         `AcceptedSampleRangeRMS` parameter.
 * * `PlotBaselines` (flag, default: `true`): whether to produce distributions
 *   of the extracted baselines.
 * * `BaselineTimeAverage` (real number, default: `600.0`): binning of the
 *   baseline profile vs. time, in seconds. Requires `PlotBaselines` to be set.
 * * `OutputCategory` (string, default: `"PMTWaveformBaselinesFromChannelData"`):
 *   label for the category of messages in the console output; this is the label
 *   that can be used for filtering messages via MessageFacility service
 *   configuration.
 *
 */
class icarus::PMTWaveformBaselinesFromChannelData: public art::SharedProducer {
  
    public:
  
  // --- BEGIN Configuration ---------------------------------------------------
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    /// Configuration of the algorithm parametes.
    struct AlgoConfig {
      
      fhicl::Atom<double> AcceptedSampleRangeRMS {
        Name{ "AcceptedSampleRangeRMS" },
        Comment{
          "The range of samples accepted for baseline, expressed in number of"
          " RMS units from the first baseline estimation"
          }
        };
      
      fhicl::Atom<unsigned int> ExcessSampleLimit {
        Name{ "ExcessSampleLimit" },
        Comment{
          "If there is this number of samples outside the RMS limit,"
          " the waveform is not used"
          }
        };
      
    }; // AlgoConfig
    
    
    fhicl::Atom<art::InputTag> OpticalWaveforms {
      Name{ "OpticalWaveforms" },
      Comment{ "label of input digitized optical waveform data product" }
      };
    
    fhicl::OptionalAtom<unsigned int> PretriggerBufferSize {
      Name{ "PretriggerBufferSize" },
      Comment{ "the number of samples in the pre-trigger readout buffer" }
      };
    
    fhicl::OptionalAtom<art::InputTag> PMTconfigurationTag {
      Name{ "PMTconfigurationTag" },
      Comment{ "PMT readout configuration object" }
      };
    
    fhicl::Atom<float> PretriggerBufferFractionForBaseline {
      Name{ "PretriggerBufferFractionForBaseline" },
      Comment
        { "fraction of the pretrigger buffer used for baseline determination" },
      0.5f
      };
    
    fhicl::Atom<unsigned int> ExcludeSpillTimeIfMoreThan {
      Name{ "ExcludeSpillTimeIfMoreThan" },
      Comment{
        "do not include the waveform at spill time"
        " if there are at least this number of off-spill PMT waveforms"
      },
      std::numeric_limits<unsigned int>::max()
      };
    
    fhicl::TableAs<opdet::SharedWaveformBaseline::Params_t, AlgoConfig>
    AlgoParams
      {
        Name{ "AlgoParams" },
        Comment{ "baseline algorithm parameters" }
      };
    
    fhicl::Atom<bool> PlotBaselines {
      Name{ "PlotBaselines" },
      Comment{ "produce plots on the extracted baseline" },
      true
      };
    
    fhicl::Atom<double> BaselineTimeAverage {
      Name{ "BaselineTimeAverage" },
      Comment{ "binning of the baseline profile vs. time [s]" },
      [this](){ return PlotBaselines(); }, // enabled if `PlotBaselines` is set
      600.0
      };
    
    fhicl::Atom<std::string> OutputCategory {
      Name{ "OutputCategory" },
      Comment{ "tag of the module output to console via message facility" },
      "PMTWaveformBaselinesFromChannelData"
      };
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  // --- END Configuration -----------------------------------------------------
  
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  explicit PMTWaveformBaselinesFromChannelData
    (Parameters const& config, art::ProcessingFrame const&);
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Framework hooks -------------------------------------------------
  /// Prepares the plots to be filled.
  virtual void beginJob(art::ProcessingFrame const& frame) override;
  
  /// Reads needed PMT configuration.
  virtual void beginRun(art::Run& run, art::ProcessingFrame const&) override;
  
  /// Creates the data products.
  virtual void produce
    (art::Event& event, art::ProcessingFrame const&) override;
  
  /// Remove empty plots.
  virtual void endJob(art::ProcessingFrame const& frame) override;
  
  // --- END Framework hooks ---------------------------------------------------
  
  
    private:
  
  // --- BEGIN Configuration variables -----------------------------------------
  
  art::InputTag const fOpDetWaveformTag; ///< Input optical waveform tag.
  
  /// Minimum number of waveforms when excluding trigger one is allowed.
  unsigned int fExcludeSpillTimeIfMoreThan;
  
  // this may be updated during the job:
  unsigned int fPretriggerSamples; ///< Samples in the pre-trigger buffer.
  
  /// Tag for PMT configuration data product.
  std::optional<art::InputTag> fPMTconfigTag;
  
  float const fSampleFraction; ///< Fraction of pretrigger buffer to use.
  
  // this may be updated during the job:
  /// Parameters for the baseline algorithm.
  opdet::SharedWaveformBaseline::Params_t fAlgoParams;
  
  bool const fPlotBaselines; ///< Whether to produce plots.
  
  /// Width of baseline time profile binning [s]
  double const fBaselineTimeAverage { 0.0 };
  
  std::string const fLogCategory; ///< Category name for the console output stream.
  
  // --- END Configuration variables -------------------------------------------
  
  
  // --- BEGIN Service variables -----------------------------------------------
  
  /// Duration of a PMT digitizer tick.
  util::quantities::intervals::nanoseconds fOpticalTick;
  
  // --- END Service variables -------------------------------------------------
  
  
  // --- BEGIN Algorithms ------------------------------------------------------
  
  // --- END Algorithms --------------------------------------------------------
  
  /// Record of baseline information to be written.
  struct BaselineInfo_t {
    icarus::WaveformBaseline baseline;
    icarus::WaveformRMS RMS;
  }; // BaselineInfo_t
  
  std::size_t fNPlotChannels = 0U; ///< Number of plotted channels
  TH2* fHBaselines = nullptr; ///< All baselines, per channel.
  
  /// For each channel, all event times and their baselines.
  std::vector<std::vector<std::pair<double, double>>> fBaselinesVsTime;
  
  /// Removes `waveforms` containing `time`, retuning how many were removed.
  unsigned int removeWaveformsAround
    (std::vector<raw::OpDetWaveform const*>& waveforms, double time) const;
  
  /// Returns the smallest pre-trigger buffer size available among all boards.
  unsigned int getPretriggerBuffer
    (sbn::PMTconfiguration const& PMTconfig) const;
  
  /// Creates all the plots to be filled by the module.
  void setupPlots(art::ProcessingFrame const& frame);
  
  /// Removes the empty plots.
  void buildBaselineGraphs(art::ProcessingFrame const& frame);
  
}; // icarus::PMTWaveformBaselinesFromChannelData



//------------------------------------------------------------------------------
//--- Implementation
//------------------------------------------------------------------------------
namespace icarus {
  
  opdet::SharedWaveformBaseline::Params_t convert
    (icarus::PMTWaveformBaselinesFromChannelData::Config::AlgoConfig const& config)
  {
    return {
        0U                               // nSample (to be changed run by run)
      , config.AcceptedSampleRangeRMS()  // nRMS
      , config.ExcessSampleLimit()       // nExcessSamples
      };
  } // convert(icarus::PMTWaveformBaselinesFromChannelData::Config::AlgoConfig)
  
} // namespace icarus


//------------------------------------------------------------------------------
namespace {
  
  template <typename Vect, typename Index = typename Vect::size_type>
  class VectorExpander {
      public:
    using Vector_t = Vect; ///< Type of bound vector.
    using Index_t = Index; ///< Type used for indexing the vector.
    using value_type = typename Vector_t::value_type;
    using reference = typename Vector_t::reference;
    
    /// Constructor: binds to the specified vector (forever).
    VectorExpander(Vector_t& v, value_type defaultValue = value_type{})
      : fVector{ &v }, fDefValue{ std::move(defaultValue) } {}
    
    /// Returns a reference to the element `index` (created as needed).
    reference operator() (Index_t index) const
      {
        auto const ix = static_cast<typename Vector_t::size_type>(index);
        if (fVector->size() <= ix) fVector->resize(ix + 1, fDefValue);
        return (*fVector)[ix];
      }
    
    
      private:
    Vector_t* fVector; ///< Bound vector.
    value_type fDefValue; ///< Value used when resizing.
    
  }; // VectorExpander
  
} // local namespace


//------------------------------------------------------------------------------
//--- icarus::PMTWaveformBaselinesFromChannelData
//------------------------------------------------------------------------------

icarus::PMTWaveformBaselinesFromChannelData::PMTWaveformBaselinesFromChannelData
  (Parameters const& config, art::ProcessingFrame const& frame)
  : art::SharedProducer(config)
  // configuration
  , fOpDetWaveformTag(config().OpticalWaveforms())
  , fExcludeSpillTimeIfMoreThan(config().ExcludeSpillTimeIfMoreThan())
  , fPretriggerSamples(config().PretriggerBufferSize().value_or(0U))
  , fPMTconfigTag(config().PMTconfigurationTag())
  , fSampleFraction(config().PretriggerBufferFractionForBaseline())
  , fAlgoParams(config().AlgoParams())
  , fPlotBaselines(config().PlotBaselines())
  , fBaselineTimeAverage(config().BaselineTimeAverage())
  , fLogCategory(config().OutputCategory())
  // service caching
  , fOpticalTick(
      detinfo::DetectorTimings
        { frame.serviceHandle<detinfo::DetectorClocksService>()->DataForJob() }
        .OpticalClockPeriod()
    )
  // algorithms
{
  
  if (fPlotBaselines)
//     serialize<art::InEvent>(art::TFileService::resource_name()); // TODO isn't art supposed to provide this method?
      serializeExternal<art::InEvent>(std::string{ "TFileService" });
  else
    async<art::InEvent>();
  
  
  //
  // optional configuration parameters
  //
  if (fPMTconfigTag.has_value() == (fPretriggerSamples != 0)) {
    throw art::Exception(art::errors::Configuration)
      << "Exactly one between parameters '"
      << config().PMTconfigurationTag.name() << "' and '"
      << config().PretriggerBufferFractionForBaseline.name()
      << "' (with a positive value) must be specified!\n";
  }
  
  //
  // configuration report
  //
  {
    mf::LogInfo log{ fLogCategory };
    log << "Using the standard (median) algorithm, waveform by waveform, on '"
      << fOpDetWaveformTag.encode() << "'";
  }
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  if (fPMTconfigTag) consumes<sbn::PMTconfiguration>(*fPMTconfigTag);
  
  //
  // declaration of output
  //
  produces<std::vector<icarus::WaveformBaseline>>();
  produces<std::vector<icarus::WaveformRMS>>();
  produces<art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform>>();
  produces<art::Assns<icarus::WaveformRMS, raw::OpDetWaveform>>();
  
} // icarus::PMTWaveformBaselinesFromChannelData::PMTWaveformBaselinesFromChannelData()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromChannelData::beginJob
  (art::ProcessingFrame const& frame)
{
  
  //
  // set up the plots, if needed
  //
  if (fPlotBaselines) setupPlots(frame);
  
} // icarus::PMTWaveformBaselinesFromChannelData::beginJob()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromChannelData::beginRun
  (art::Run& run, art::ProcessingFrame const&)
{
  
  if (fPMTconfigTag) { // update the pretrigger samples number
    fPretriggerSamples = getPretriggerBuffer
      (run.getProduct<sbn::PMTconfiguration>(*fPMTconfigTag));
  }
  assert(fPretriggerSamples > 0U);
  
  fAlgoParams.nSample = std::max(
    1U,
    static_cast<unsigned int>(std::round(fSampleFraction * fPretriggerSamples))
    );
  
  {
    // not clear whether this is debug or info level
    mf::LogInfo log{ fLogCategory };
    log << run.id() << ": baseline algorithm configuration:\n";
    fAlgoParams.dump(log, " - ");
  }
  
} // icarus::PMTWaveformBaselinesFromChannelData::beginRun()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromChannelData::produce
  (art::Event& event, art::ProcessingFrame const& frame)
{
  
  detinfo::DetectorTimings const detTimings{
    frame.serviceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    };
  
  detinfo::timescales::electronics_time const triggerTime
     = detTimings.TriggerTime();
  
  mf::LogDebug{ fLogCategory }
    << event.id() << " trigger time: " << triggerTime;
  
  //
  // fetch input
  //
  auto const& waveformHandle
    = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
  auto const& waveforms = *waveformHandle;
  
  //
  // compute all the baselines (and plot them)
  //
  opdet::SharedWaveformBaseline const sharedWaveformBaselineAlgo
    { fAlgoParams, fLogCategory };
  
  double const eventTime = static_cast<double>(event.time().timeHigh())
    + static_cast<double>(event.time().timeHigh()) * 1e-9;
  
  std::vector<BaselineInfo_t> channelBaselines;
  VectorExpander baselineForChannel { channelBaselines };
  
  icarus::ns::util::GroupByIndex<raw::OpDetWaveform> waveformsByChannel
    { waveforms, std::mem_fn(&raw::OpDetWaveform::ChannelNumber) };
  
  for
    (auto const& [ channel, allWaveforms ]: util::enumerate(waveformsByChannel))
  {
    if (allWaveforms.empty()) continue;
    
    mf::LogTrace{ fLogCategory }
      << "Processing " << allWaveforms.size() << " waveforms for channel "
      << channel;
    
    //
    // remove global trigger waveform
    //
    std::vector<raw::OpDetWaveform const*> waveforms{ allWaveforms };
    if (waveforms.size() >= fExcludeSpillTimeIfMoreThan) {
      
      unsigned int const nExcluded
        = removeWaveformsAround(waveforms, triggerTime.value());
      if (nExcluded > 0U) {
        mf::LogTrace{ fLogCategory }
          << "Removed " << nExcluded << "/" << (waveforms.size() + nExcluded)
          << " waveforms at trigger time " << triggerTime;
      }
      
    } // if many waveforms
    
    //
    // extract baseline
    //
    opdet::SharedWaveformBaseline::BaselineInfo_t const baseline
      = sharedWaveformBaselineAlgo(waveforms);

    mf::LogTrace{ fLogCategory }
      << "Channel " << channel << ": baseline " << baseline.baseline
      << " ADC# from " << baseline.nSamples << " samples in "
      << baseline.nWaveforms << "/" << waveforms.size()
      << " waveforms; found RMS=" << baseline.RMS << " ADC#";
    
    auto const chIndex = static_cast<std::size_t>(channel);
    baselineForChannel(chIndex) = { baseline.baseline, baseline.RMS };
    
    if (fHBaselines) fHBaselines->Fill(double(channel), baseline.baseline);
    if (chIndex < fBaselinesVsTime.size())
      fBaselinesVsTime[chIndex].emplace_back(eventTime, baseline.baseline);
    
  } // for grouped waveforms
  
  //
  // assign baselines to waveforms
  //
  
  std::vector<icarus::WaveformBaseline> baselines;
  baselines.reserve(waveforms.size());
  art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform> baselineToWaveforms;
  art::PtrMaker<icarus::WaveformBaseline> const makeBaselinePtr(event);
  
  std::vector<icarus::WaveformRMS> RMSs;
  RMSs.reserve(waveforms.size());
  art::Assns<icarus::WaveformRMS, raw::OpDetWaveform> RMStoWaveforms;
  art::PtrMaker<icarus::WaveformRMS> const makeRMSPtr(event);
  
  for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {
    
    BaselineInfo_t baselineInfo
      = channelBaselines[static_cast<std::size_t>(waveform.ChannelNumber())];
    
    art::Ptr<raw::OpDetWaveform> const waveformPtr{ waveformHandle, iWaveform };
    
    baselines.push_back(baselineInfo.baseline);
    RMSs.push_back(baselineInfo.RMS);
    baselineToWaveforms.addSingle(makeBaselinePtr(iWaveform), waveformPtr);
    RMStoWaveforms.addSingle(makeRMSPtr(iWaveform), waveformPtr);
    
  } // for
  
  //
  // output
  //
  event.put(
    std::make_unique<std::vector<icarus::WaveformBaseline>>
      (std::move(baselines))
    );
  event.put
    (std::make_unique<std::vector<icarus::WaveformRMS>>(std::move(RMSs)));
  event.put(
    std::make_unique<art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform>>
      (std::move(baselineToWaveforms))
    );
  event.put(
    std::make_unique<art::Assns<icarus::WaveformRMS, raw::OpDetWaveform>>
      (std::move(RMStoWaveforms))
    );
  
  
} // icarus::PMTWaveformBaselinesFromChannelData::produce()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromChannelData::endJob
  (art::ProcessingFrame const& frame)
{
  
  if (fPlotBaselines) buildBaselineGraphs(frame);
  
} // icarus::PMTWaveformBaselinesFromChannelData::endJob()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromChannelData::setupPlots
  (art::ProcessingFrame const& frame)
{
  
  auto const& tfs = *(frame.serviceHandle<art::TFileService>());
  auto const& geom = *(frame.serviceHandle<geo::Geometry>()->provider());
  
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
  
} // icarus::PMTWaveformBaselinesFromChannelData::setupPlots()


//------------------------------------------------------------------------------
unsigned int icarus::PMTWaveformBaselinesFromChannelData::getPretriggerBuffer
  (sbn::PMTconfiguration const& PMTconfig) const
{
  lar::util::MinMaxCollector<unsigned int> prebuffer;
  for (sbn::V1730Configuration const& board: PMTconfig.boards)
    prebuffer.add(board.preTriggerTicks());
  if (!prebuffer.has_data()) {
    throw art::Exception{ art::errors::Unknown }
      << "No boards found in the PMT configuration!\n";
  }
  else if (prebuffer.max() > prebuffer.min()) {
    throw art::Exception{ art::errors::Unknown }
      << "Found different pre-trigger readout buffer sizes between "
      << prebuffer.min() << " and " << prebuffer.max() << " ticks.\n";
  }
  return prebuffer.min();
} // icarus::PMTWaveformBaselinesFromChannelData::getPretriggerBuffer()


//------------------------------------------------------------------------------
unsigned int icarus::PMTWaveformBaselinesFromChannelData::removeWaveformsAround
  (std::vector<raw::OpDetWaveform const*>& waveforms, double time) const
{
  
  // tick duration in the same unit as the electronics time scale times
  // (that is microseconds, actually)
  double const tickDuration
    = detinfo::timescales::electronics_time::interval_t{ fOpticalTick }.value();
  
  auto const hasTime = [tickDuration,time](raw::OpDetWaveform const* waveform)
    {
      double const relTime = time - waveform->TimeStamp();
      return (relTime >= 0.0) && (relTime < (waveform->size() * tickDuration));
    };
  
  auto const iLast
    = std::remove_if(waveforms.begin(), waveforms.end(), hasTime);
  
  unsigned int const nRemoved = std::distance(iLast, waveforms.end());
  waveforms.erase(iLast, waveforms.end());
  
  return nRemoved;
  
} // icarus::PMTWaveformBaselinesFromChannelData::removeWaveformsAround()


//------------------------------------------------------------------------------
void icarus::PMTWaveformBaselinesFromChannelData::buildBaselineGraphs
  (art::ProcessingFrame const& frame)
{
  
  auto& tfs = *(frame.serviceHandle<art::TFileService>());
  
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
  
} // icarus::PMTWaveformBaselinesFromChannelData::buildBaselineGraphs()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::PMTWaveformBaselinesFromChannelData)


//------------------------------------------------------------------------------
