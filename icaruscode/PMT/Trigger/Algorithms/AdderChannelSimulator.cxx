/**
* @file   icaruscode/PMT/Trigger/Algorithms/AdderChannelSimulator.cxx
* @brief  Simulation of adder board output.
* @author Gianluca Petrillo (petrillo@slac.stanford.edu)
* @date   August 25, 2025
* @see    icaruscode/PMT/Trigger/Algorithms/AdderChannelSimulator.h
* 
*/

#undef NDEBUG // FIXME

// library header
#include "icaruscode/PMT/Trigger/Algorithms/AdderChannelSimulator.h"

// ICARUS and framework libraries
#include "icaruscode/PMT/Trigger/Algorithms/FilteringCircuit.h"
#include "icaruscode/PMT/Trigger/Algorithms/WaveformFilterAlg.h"
#include "sbnalg/Utilities/ToStr.h" // util::toStr
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"
#include "cetlib_except/exception.h"
#include "TDirectory.h"
#include "TGraph.h"

// C/C++ standard libraries
#include <cstddef> // std::ptrdiff_t
#include <cmath> // std::isnormal()
#include <iomanip>
#include <limits>
#include <string>
#include <utility> // std::move()
#include <vector>


// -----------------------------------------------------------------------------
using detinfo::timescales::electronics_time;
using util::quantities::intervals::microseconds;


// -----------------------------------------------------------------------------
namespace {
  
  /// Use `out << dumpCollection(v)` to print all elements in one line...
  template <typename Coll>
  struct dumpCollection {
    static constexpr auto OhSoMany = std::numeric_limits<std::size_t>::max();
    Coll const* coll = nullptr;
    std::size_t perLine = OhSoMany;
    std::size_t maxElement = OhSoMany;
    
    dumpCollection(
      Coll const& coll, std::size_t perLine = 0, std::size_t maxElement = OhSoMany
      )
      : coll{ &coll }, perLine{ perLine }, maxElement{ maxElement } {}
    
    template <typename Stream>
    void operator() (Stream& out) const {
      using std::size, std::begin, std::end;
      out << "(" << size(*coll) << " elements) [";
      std::size_t leftInLine = 0;
      for (auto const& [ i, e ]: util::enumerate(*coll)) {
        if ((perLine > 0) && (leftInLine-- == 0)) {
          leftInLine = perLine;
          out << "\n  [" << i << "] ";
        }
        if (i > maxElement) {
          out << " ... and " << (size(*coll) - i) << " more";
          break;
        }
        out << " " << e;
      }
      out << " ]";
    }
  }; // struct dumpCollection
  
  template <typename Stream, typename Coll>
  Stream& operator<< (Stream&& out, dumpCollection<Coll> const& collDumper) {
    collDumper(out);
    return out;
  }
  
} // local namespace


// -----------------------------------------------------------------------------
// ---  icarus::trigger::AdderChannelSimulator
// -----------------------------------------------------------------------------
raw::OpDetWaveform
icarus::trigger::AdderChannelSimulator::AdderSignal_t::makeOpDetWaveform
  (ADCsettings_t const& ADCsettings, Voltage_t baseline /* = 0 */) const
{
  using std::begin, std::end;
  raw::OpDetWaveform waveform{ startTime, raw::Channel_t{ channel } };
  WaveformSamples_t const ADCs = samples * ADCsettings.mV2ADC() + baseline;
  waveform.Waveform()
    = std::vector<raw::ADC_Count_t>{ begin(ADCs), end(ADCs) };
  return waveform;
} // icarus::trigger::AdderChannelSimulator::AdderSignal_t::makeOpDetWaveform()


// -----------------------------------------------------------------------------
icarus::trigger::AdderChannelSimulator::AdderChannelSimulator(
  AdderChannelID channel,
  ADCsettings_t ADCsettings,
  double splitFraction,
  FilteringCircuit const* reshapingFilter,
  AdderCalibrationDatabase::RunCalibration const* calibration,
  TDirectory* plotDir /* = nullptr */,
  std::string logCategory /* = "AdderChannelSimulator" */
)
  : icarus::ns::util::mfLoggingClass{ std::move(logCategory) }
  , fChannel         { channel }
  , fADCsettings     { std::move(ADCsettings) }
  , fPMTsplitterScale{ splitFraction / (1.0 - splitFraction) }
  , fPlotDir         { plotDir }
  , fReshapingFilter { reshapingFilter }
  , fCalibration     { calibration }
  , fFilterAlg       {
      fReshapingFilter
        ? std::make_unique<WaveformFilterAlg>(fADCsettings.samplingRate())
        : nullptr
      }
{
  //
  // parameter check
  //
  if (splitFraction <= 0 || splitFraction >= 1) {
    throw cet::exception{ "AdderChannelSimulator" }
      << "Splitter fraction must be in ] 0, 1 [ (specified: " << splitFraction
      << ").\n";
  }
  
} // icarus::trigger::AdderChannelSimulator::AdderChannelSimulator()


// -----------------------------------------------------------------------------
// need to postpone destructor definition to a place where `WaveformFilterAlg`
// is fully defined and `std::unique_ptr<WaveformFilterAlg> fFilterAlg`
// can be destroyed.
icarus::trigger::AdderChannelSimulator::~AdderChannelSimulator() = default;


// -----------------------------------------------------------------------------
auto icarus::trigger::AdderChannelSimulator::simulate(
  TimeInterval_t const& timeInterval,
  std::vector<WaveformWithBaseline const*> const& waveforms
) const -> std::pair<AdderSignal_t, std::vector<WaveformWithBaseline const*>>
{
  mfLogTrace() << "Simulating adder channel " << fChannel
    << " with " << waveforms.size() << " waveforms in " << timeInterval;
  
  unsigned int nClippedInput = 0; // will receive the number of clipped channels
  
  // addWaveformsInInterval() converts early to floating point before summing,
  //   preventing ADC count data type (16-bit) overflow
  
  // std::valarray<float> PMTsum
  // std::vector<WaveformWithBaseline const*> contributingWaveforms
  auto [ PMTsum, contributingWaveforms ]
    = addWaveformsInInterval(waveforms, timeInterval, &nClippedInput);
  
  PMTsum *= fPMTsplitterScale;
  
  if (fPlotDir) {
    std::unique_ptr<TGraph> graph
      = makeWaveformPlot(timeInterval.start, PMTsum, "PMTsum", "#sum PMT", false);
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);
    writePlot(graph.get());
  }
  
  mfLogTrace()
    << "∑ PMT waveform with " << PMTsum.size() << " samples created from "
    << contributingWaveforms.size() << " waveforms; " << nClippedInput << "/"
    << waveforms.size() << " were found to be clipped.";
    
  WaveformSamples_t reshapedPMTsum = reshapeWaveform(PMTsum);
  assert(reshapedPMTsum.size() == PMTsum.size());
  if (fPlotDir) {
    std::unique_ptr<TGraph> graph = makeWaveformPlot
      (timeInterval.start, reshapedPMTsum, "ReshapedPMTsum", "reshaped #sum PMT", false);
    graph->SetLineColor(kOrange);
    graph->SetLineWidth(2);
    writePlot(graph.get());
  }
  
  // we extract the time offset first, while `reshapedPMTsum` is still valid
  microseconds const timeOffset = fCalibration
    ? fCalibration->timeOffset(fChannel, reshapedPMTsum): microseconds{ 0 };
  
  WaveformSamples_t adderWaveformSamples
    = calibrateReshapedPMTsum(std::move(reshapedPMTsum), fChannel);
  
  electronics_time const startTime = timeInterval.start + timeOffset;
  
  assert(adderWaveformSamples.size() == PMTsum.size());
  if (fPlotDir) {
    std::unique_ptr<TGraph> graph = makeWaveformPlot
      (startTime, adderWaveformSamples, "SimulatedAdder", "simulated adder", false);
    graph->SetLineColor(kRed);
    graph->SetLineWidth(2);
    writePlot(graph.get());
  }
  
  // we need to convert the result into a standard waveform sample type;
  // this initialization pattern should use move-assignment
  AdderSignal_t adderWaveform{
      fChannel                        // channel
    , std::move(adderWaveformSamples) // samples
    , startTime.value()               // startTime
    , nClippedInput                   // nClipped
    };
    
  return { std::move(adderWaveform), std::move(contributingWaveforms) };
} // icarus::trigger::AdderChannelSimulator::simulate()


//------------------------------------------------------------------------------
auto icarus::trigger::AdderChannelSimulator::addWaveformsInInterval(
  std::vector<WaveformWithBaseline const*> const& waveforms,
  TimeInterval_t const& timeInterval,
  unsigned int* pNClipped /* = nullptr */
) const
  -> std::pair<WaveformSamples_t, std::vector<WaveformWithBaseline const*>>
{
  using elapsed_time_t = electronics_time::interval_t;
  
  class TimeToWaveformTick {
    
    electronics_time fStartTime; ///< Waveform start time.
    elapsed_time_t fTickDuration; ///< Duration of a tick.
    
    // truncates, but if value is closer than `tol` to an integer, rounds to it
    static constexpr std::ptrdiff_t truncate(double v, double tol = 0.999)
      {
        std::ptrdiff_t tv = static_cast<std::ptrdiff_t>(std::trunc(v));
        if (v - tv >= tol) ++tv;
        return tv;
      }
    
      public:
    
    TimeToWaveformTick
      (raw::OpDetWaveform const& waveform, elapsed_time_t tickDuration)
      : TimeToWaveformTick{ waveform.TimeStamp(), tickDuration }
      {}
    TimeToWaveformTick
      (raw::TimeStamp_t startTime, elapsed_time_t tickDuration)
      : fStartTime{ startTime }, fTickDuration{ tickDuration }
      {}
    
    std::ptrdiff_t operator() (electronics_time time) const
      { return (*this)(time - fStartTime); }
    std::ptrdiff_t operator() (elapsed_time_t delta) const
      { return toTicks(delta, fTickDuration); }
    
    electronics_time startTime() const { return fStartTime; }
    
    // using `round()` to work around rounding errors
    static std::ptrdiff_t toTicks
      (elapsed_time_t delta, elapsed_time_t tickDuration)
      { return truncate(delta / tickDuration); };
    
  }; // TimeToWaveformTick
  
  
  if (pNClipped) *pNClipped = 0;
  
  std::ptrdiff_t const nSamples = TimeToWaveformTick::toTicks
    (timeInterval.duration(), fADCsettings.samplingTime);
  assert(nSamples >= 0);
  
  if (waveforms.empty() || nSamples == 0) return {};
  
  Voltage_t baselineSum = 0.0;
  for (WaveformWithBaseline const* wb: waveforms)
    baselineSum += wb->hasBaseline()? wb->baseline().baseline(): Voltage_t{ 0 };
  
  // start with the total baseline, we then subtract all signals.
  WaveformSamples_t samples(baselineSum, nSamples);
  
  auto log = mfLogTrace();
  log << "Adding " << waveforms.size() << " waveforms in " << timeInterval
    << " (" << nSamples << " samples)";
  for (WaveformWithBaseline const* wb: waveforms) {
    
    raw::OpDetWaveform const& waveform = wb->waveform();
    
    TimeToWaveformTick const toTick{ waveform, fADCsettings.samplingTime };
    
    std::ptrdiff_t const beginTick = toTick(timeInterval.start);
    assert((beginTick >= 0) && (beginTick < std::ptrdiff_t(waveform.size())));
    std::ptrdiff_t const endTick = toTick(timeInterval.stop);
    assert((endTick > 0) && (endTick <= std::ptrdiff_t(waveform.size())));
    assert(endTick - beginTick == nSamples);
    
    bool clipped = false;
    for (std::ptrdiff_t i = beginTick; i < endTick; ++i) {
      samples[i] -= waveform[i];
      if (waveform[i] == 0) clipped = true;
    }
    log << "\n  [CH=" << waveform.ChannelNumber() << " at "
      << toTick.startTime() << "] samples [#" << beginTick
      << "; #" << endTick << "[ (" << (endTick - beginTick) << " ticks)";
    
    if (pNClipped && clipped) {
      ++*pNClipped;
      log << " (clipped)";
    }
  } // for waveforms
  
  // convert to mV
  samples *= fADCsettings.ADC2mV();
  
  return { std::move(samples), waveforms };
} // icarus::trigger::AdderChannelSimulator::::addWaveformsInInterval()


//------------------------------------------------------------------------------
auto icarus::trigger::AdderChannelSimulator::reshapeWaveform
  (WaveformSamples_t PMTsum) const -> WaveformSamples_t
{
  
  if (!fReshapingFilter) return PMTsum;
  
  assert(fFilterAlg);
  
  auto spectrum = fFilterAlg->waveformToSpectrum(PMTsum);
  
  double const samplingTime
    = fFilterAlg->samplingPeriod().convertInto<util::quantities::second>().value();
  auto reshapedSpectrum
    = fReshapingFilter->apply(samplingTime, std::move(spectrum));
  
  auto reshapedWaveform
    = fFilterAlg->spectrumToWaveform(std::move(reshapedSpectrum));
  
  return reshapedWaveform;
  
} // icarus::trigger::AdderChannelSimulator::::reshapeWaveform()


//------------------------------------------------------------------------------
auto icarus::trigger::AdderChannelSimulator::calibrateReshapedPMTsum
  (WaveformSamples_t PMTsum, AdderChannelID channel) const -> WaveformSamples_t
{
  
  if (!fCalibration) return PMTsum;
  
  double const c = fCalibration->calibrationFactor(channel, PMTsum);
  if (!std::isnormal(c)) {
    mfLogError() << "Adder calibration factor is invalid (" << c << ")!";
    return PMTsum;
  }
  mfLogTrace() << "  applying calibration factor " << c << " to adder channel "
    << channel;
  return PMTsum * c;
  
} // icarus::trigger::AdderChannelSimulator::calibrateReshapedPMTsum()


// -----------------------------------------------------------------------------
template <typename Samples>
std::unique_ptr<TGraph> icarus::trigger::AdderChannelSimulator::makeWaveformPlot(
  electronics_time startTime, Samples const& samples,
  std::string const& typeName /* = "Waveform" */,
  std::string const& descr /* = "" */,
  bool write /* = true  */
) const {
  
  if (!fPlotDir) return nullptr;
  
  using std::size, std::cbegin;
  using util::toStr;
  
  using value_type = typename Samples::value_type;
  
  
  std::size_t const nSamples = size(samples);
  value_type const* sampleValues = &*cbegin(samples);
  assert(sampleValues != nullptr);
  
  std::vector<value_type> times(nSamples);
  for (auto const i: util::counter(nSamples))
    times[i] = (startTime + i * fADCsettings.samplingTime).value();
  
  // if this is changed, update the class documentation
  std::string const name = toStr << typeName << "_" << fPlotDir->GetName()
    << "_T" << std::fixed << std::setprecision(3) << std::setfill('0') << std::setw(8)
    << startTime.value();
  
  std::string title = toStr << fPlotDir->GetTitle() << ": " << descr;
  if (title.find(';') == std::string::npos) {
    title += toStr << ";electronics time  [ #mus ];voltage  [ mV / "
      << fADCsettings.samplingTime << " ]";
  }
  
  // in principle we don't need to go into the destination directory
  // when creating a graph (we should when creating a histogram or tree)
  TDirectory::TContext ctx{ fPlotDir };
  auto graph = std::make_unique<TGraph>((int) nSamples, times.data(), sampleValues);
  graph->SetNameTitle(name.c_str(), title.c_str());
  
  if (write) writePlot(graph.get());
  
  return graph;
} // icarus::trigger::AdderChannelSimulator::makeWaveformPlot()


//------------------------------------------------------------------------------
int icarus::trigger::AdderChannelSimulator::writePlot(TObject* obj) const {
  if (!obj) return 0;
  TDirectory::TContext ctx{ fPlotDir };
  return obj->Write();
}


//------------------------------------------------------------------------------

