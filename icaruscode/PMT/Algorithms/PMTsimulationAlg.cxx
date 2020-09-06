/**
 * @file    icaruscode/PMT/Algorithms/PMTsimulationAlg.cxx
 * @brief   Algorithms for the simulation of ICARUS PMT channels.
 * @date    October 16, 2018
 * @see     `icaruscode/PMT/Algorithms/PMTsimulationAlg.h`
 *
 * Implementation file.
 */

// this library header
#include "icaruscode/PMT/Algorithms/PMTsimulationAlg.h"

// ICARUS libraries
#include "icaruscode/PMT/Algorithms/AsymGaussPulseFunction.h"
#include "icaruscode/Utilities/WaveformOperations.h"


// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::begin(), util::end()

// framework libraries
#include "cetlib_except/exception.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandExponential.h"

// C++ standard libaries
#include <chrono> // std::chrono::high_resolution_clock
#include <unordered_map>
#include <algorithm> // std::accumulate()
#include <utility> // std::move(), std::cref(), ...
#include <limits> // std::numeric_limits
#include <cmath> // std::signbit(), std::pow()


// -----------------------------------------------------------------------------
// ---  icarus::opdet::PMTsimulationAlg
// -----------------------------------------------------------------------------

// The reason for this being static is feeble... in short: because it can be;
// the memory it takes is already "static" anyway, so no gain in there.
// At least, this clarifies the absolute nature of this object.
util::FastAndPoorGauss<32768U, float> const
  icarus::opdet::PMTsimulationAlg::fFastGauss;


// -----------------------------------------------------------------------------
double
icarus::opdet::PMTsimulationAlg::ConfigurationParameters_t::PMTspecs_t::
  multiplicationStageGain(unsigned int i /* = 1 */) const
{
  // if all stages were born equal:
  // return std::pow(gain, 1.0 / nDynodes());
  double const k = dynodeK;
  double const mu = gain;
  unsigned const N = nDynodes();
  double prodRho = 1.0;
  for (double rho: voltageDistribution) prodRho *= rho;
  double const aVk
    = std::pow(mu / std::pow(prodRho, k), 1.0/static_cast<double>(N));
  return aVk * std::pow(voltageDistribution.at(i - 1), k);
} // icarus::...::PMTsimulationAlg::...::PMTspecs_t::multiplicationStageGain()


// -----------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::ConfigurationParameters_t::PMTspecs_t::
  setVoltageDistribution
  (std::vector<double>&& Rs)
{
  voltageDistribution = std::move(Rs);
  double total = std::accumulate
    (voltageDistribution.begin(), voltageDistribution.end(), 0.0);
  for (auto& R: voltageDistribution) R /= total;
} // icarus::...::PMTsimulationAlg::...::PMTspecs_t::setPMTvoltageDistribution()


// -----------------------------------------------------------------------------
icarus::opdet::PMTsimulationAlg::PMTsimulationAlg
  (ConfigurationParameters_t const& config)
  : fParams(config)
  , fQE(fParams.QEbase / fParams.larProp->ScintPreScale())
  , fSampling(fParams.clockData->OpticalClock().Frequency())
  , fNsamples(fParams.readoutEnablePeriod * fSampling) // us * MHz cancels out
  , wsp(
    *(fParams.pulseFunction),
    fSampling,
    fParams.pulseSubsamples, // tick subsampling
    1.0e-4_ADCf // stop sampling when ADC counts are below this value
    )
  , fNoiseAdder(fParams.useFastElectronicsNoise
      ? &icarus::opdet::PMTsimulationAlg::AddNoise_faster
      : &icarus::opdet::PMTsimulationAlg::AddNoise
    )
{
  using namespace util::quantities::electronics_literals;

  //  mf::LogDebug("PMTsimulationAlg") << "Sampling = " << fSampling << std::endl;

  // shape of single pulse

  // Correction due to scalling factor applied during simulation if any
  mf::LogDebug("PMTsimulationAlg") << "PMT corrected efficiency = " << fQE;

  if (fQE >= 1.0001) {
    mf::LogWarning("PMTsimulationAlg")
      << "WARNING: Quantum efficiency set in the configuration (QE="
        << fParams.QEbase << ") seems to be too large!"
      "\nThe photon visibility library is assumed to already include a"
        " quantum efficiency of " << fParams.larProp->ScintPreScale() <<
        " (`ScintPreScale` setting of `LArProperties` service)"
        " and here we are requesting a higher one."
      "\nThis configuration is not supported by `PMTsimulationAlg`,"
        " and this simulation will effectively apply a quantum efficiency of "
        << fParams.larProp->ScintPreScale();
      ;
  }

  // check that the sampled waveform has a sufficiently large range, so that
  // tails are below 10^-3 ADC counts (absolute value);
  // if this test fails, it's better to reduce the threshold in wsp constructor
  wsp.checkRange(1.0e-3_ADCf, "PMTsimulationAlg");

} // icarus::opdet::PMTsimulationAlg::PMTsimulationAlg()


// -----------------------------------------------------------------------------
std::tuple<std::vector<raw::OpDetWaveform>, std::optional<sim::SimPhotons>>
  icarus::opdet::PMTsimulationAlg::simulate(sim::SimPhotons const& photons)
{
  std::optional<sim::SimPhotons> photons_used;

  Waveform_t const waveform = CreateFullWaveform(photons, photons_used);

  return {
    CreateFixedSizeOpDetWaveforms(photons.OpChannel(), waveform),
    std::move(photons_used)
    };
  
} // icarus::opdet::PMTsimulationAlg::simulate()


//------------------------------------------------------------------------------
auto icarus::opdet::PMTsimulationAlg::makeGainFluctuator() const {

  using Fluctuator_t = GainFluctuator<CLHEP::RandPoisson>;

  if (fParams.doGainFluctuations) {
    double const refGain = fParams.PMTspecs.firstStageGain();
    return Fluctuator_t
      { refGain, CLHEP::RandPoisson{ *fParams.gainRandomEngine, refGain } };
  }
  else return Fluctuator_t{}; // default-constructed does not fluctuate anything

} // icarus::opdet::PMTsimulationAlg::makeGainFluctuator()


//------------------------------------------------------------------------------
auto icarus::opdet::PMTsimulationAlg::CreateFullWaveform
  (sim::SimPhotons const& photons, std::optional<sim::SimPhotons>& photons_used)
  const -> Waveform_t
{

    using namespace util::quantities::time_literals;
    using namespace util::quantities::frequency_literals;
    using namespace util::quantities::electronics_literals;
    using namespace detinfo::timescales;

    detinfo::DetectorTimings const& timings
      = detinfo::makeDetectorTimings(fParams.clockData);

    tick const endSample = tick::castFrom(fNsamples);

    //
    // collect the amount of photoelectrons arriving at each subtick;
    // the waveform is split in groups of photons at the same relative subtick
    // (i.e. first all the photons on the first subtick of a tick, then
    // all the photons on the second subtick of a tick, and so on);
    // storage is by subtick group (vector index is the subtick), then by
    // tick (unordered map index is the tick number).
    //
    std::vector<std::unordered_map<tick, unsigned int>> peMaps
      (wsp.nSubsamples());

    // returns tick and relative subtick number
    TimeToTickAndSubtickConverter const toTickAndSubtick(peMaps.size());

//     auto start = std::chrono::high_resolution_clock::now();
    
    if (photons_used) {
      photons_used->clear();
      photons_used->SetChannel(photons.OpChannel());
    }
    for(auto const& ph : photons) {
      if (!KicksPhotoelectron()) continue;

      if (photons_used) photons_used->push_back(ph); // copy

      simulation_time const photonTime { ph.Time };

      trigger_time const mytime
        = timings.toTriggerTime(photonTime)
        - fParams.triggerOffsetPMT
        ;
      if ((mytime < 0.0_us) || (mytime >= fParams.readoutEnablePeriod)) continue;

      auto const [ tick, subtick ]
        = toTickAndSubtick(mytime.quantity() * fSampling);
      /*
      mf::LogTrace("PMTsimulationAlg")
        << "Photon at " << photonTime << ", optical time " << mytime
        << " => tick " << tick_d
        << " => sample " << tick << " subsample " << subtick
        ;
      */
      if (tick >= endSample) continue;
      ++peMaps[subtick][tick];
    } // for photons

//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> diff = end-start;
//     std::cout << "\tcollected pes... " << photons.OpChannel() << " " << diff.count() << std::endl;
//     start=std::chrono::high_resolution_clock::now();

    //
    // add the collected photoelectrons to the waveform
    //
    Waveform_t waveform(fNsamples, fParams.baseline);
    
    unsigned int nTotalPE [[gnu::unused]] = 0U; // unused if not in `debug` mode
    double nTotalEffectivePE [[gnu::unused]] = 0U; // unused if not in `debug` mode

    auto gainFluctuation = makeGainFluctuator();

    // go though all subsamples (starting each at a fraction of a tick)
    for (auto const& [ iSubsample, peMap ]: util::enumerate(peMaps)) {

      // this is the waveform sampling for the selected subsample:
      auto const& subsample = wsp.subsample(iSubsample);

      for (auto const& [ startTick, nPE ]: peMap) {
        nTotalPE += nPE;

        double const nEffectivePE = gainFluctuation(nPE);
        nTotalEffectivePE += nEffectivePE;

        AddPhotoelectrons(
          subsample, waveform, startTick,
          static_cast<WaveformValue_t>(nEffectivePE)
          );

      } // for sample
    } // for subsamples
    MF_LOG_TRACE("PMTsimulationAlg")
      << nTotalPE << " photoelectrons at "
      << std::accumulate(
           peMaps.begin(), peMaps.end(), 0U,
           [](unsigned int n, auto const& map){ return n + map.size(); }
         )
      << " times in channel " << photons.OpChannel()
      ;

//       end=std::chrono::high_resolution_clock::now(); diff = end-start;
//       std::cout << "\tadded pes... " << photons.OpChannel() << " " << diff.count() << std::endl;
//       start=std::chrono::high_resolution_clock::now();

      if(fParams.ampNoise > 0.0_ADCf) (this->*fNoiseAdder)(waveform);
      if(fParams.darkNoiseRate > 0.0_Hz) AddDarkNoise(waveform);

//       end=std::chrono::high_resolution_clock::now(); diff = end-start;
//       std::cout << "\tadded noise... " << photons.OpChannel() << " " << diff.count() << std::endl;
//       start=std::chrono::high_resolution_clock::now();

      // Implementing saturation effects;
      // waveform is negative, and saturation is a minimum ADC count
      // TODO use waveform_operations (what is polarity is different?)
      auto const saturationLevel
        = fParams.baseline + fParams.saturation*wsp.peakAmplitude();
      std::replace_if(waveform.begin(),waveform.end(),
		      [saturationLevel](auto s) -> bool{return s < saturationLevel;},
		      saturationLevel);

//       end=std::chrono::high_resolution_clock::now(); diff = end-start;
//       std::cout << "\tadded saturation... " << photons.OpChannel() << " " << diff.count() << std::endl;
    
    return waveform;
  } // CreateFullWaveform()

  auto icarus::opdet::PMTsimulationAlg::CreateBeamGateTriggers() const
    -> std::vector<optical_tick>
  {
    using namespace util::quantities::time_literals;
    using detinfo::timescales::trigger_time;

    detinfo::DetectorTimings const& timings
      = detinfo::makeDetectorTimings(fParams.clockData);

    std::vector<optical_tick> trigger_locations;
    trigger_locations.reserve(fParams.beamGateTriggerNReps);
    trigger_time const beamTime = timings.toTriggerTime(timings.BeamGateTime());

    for(unsigned int i_trig=0; i_trig<fParams.beamGateTriggerNReps; ++i_trig) {
      trigger_time const trig_time
        = beamTime
        + i_trig * fParams.beamGateTriggerRepPeriod
        - fParams.triggerOffsetPMT
        ;
      if (trig_time < 0_us) continue;
      if (trig_time > fParams.readoutEnablePeriod) break;
      trigger_locations.push_back
        (optical_tick::castFrom(trig_time.quantity()*fSampling));
    }
    return trigger_locations;
  }

  auto icarus::opdet::PMTsimulationAlg::FindTriggers(Waveform_t const& wvfm) const
    -> std::vector<optical_tick>
  {
    std::vector<optical_tick> trigger_locations;

    // find all ticks at which we would trigger readout
    bool above_thresh=false;
    for(size_t i_t=0; i_t<wvfm.size(); ++i_t){

      auto const val { fParams.pulsePolarity* (wvfm[i_t]-fParams.baseline) };

      if(!above_thresh && val>=fParams.thresholdADC){
	above_thresh=true;
	trigger_locations.push_back(optical_tick::castFrom(i_t));
      }
      else if(above_thresh && val<fParams.thresholdADC){
	above_thresh=false;
      }

    }//end loop over waveform

    // next, add the triggers injected at beam gate time
    if (fParams.createBeamGateTriggers) {
      auto beamGateTriggers = CreateBeamGateTriggers();

      // insert the new triggers and sort them
      trigger_locations.insert(trigger_locations.end(),
        beamGateTriggers.begin(), beamGateTriggers.end());
      std::inplace_merge(
        trigger_locations.begin(),
        trigger_locations.end() - beamGateTriggers.size(),
        trigger_locations.end()
        );
    }

    return trigger_locations;
  }

//------------------------------------------------------------------------------
std::vector<raw::OpDetWaveform>
icarus::opdet::PMTsimulationAlg::CreateFixedSizeOpDetWaveforms
  (raw::Channel_t opChannel, Waveform_t const& waveform) const
{
  /*
   * Plan:
   * 
   * 1. set up
   * 2. get the trigger points of the waveform
   * 3. define the size of data around each trigger to commit to waveforms
   *    (also merge contiguous and overlapping intervals)
   * 4. create the actual `raw::OpDetWaveform` objects
   * 
   */
  
  //
  // parameters check and setup
  //
  
  // not a big deal if this assertion fails, but a bit more care needs to be
  // taken in comparisons and subtractions
  static_assert(
    std::is_signed_v<optical_tick::value_t>,
    "This algorithm requires tick type to be signed."
    );
  
  using namespace detinfo::timescales; // electronics_time, time_interval, ...

  auto const pretrigSize = optical_time_ticks::castFrom(fParams.pretrigSize());
  auto const posttrigSize = optical_time_ticks::castFrom(fParams.posttrigSize());
  
  detinfo::DetectorTimings const& timings
    = detinfo::makeDetectorTimings(fParams.clockData);

  // first viable tick number: since this is the item index in `wvfm`, it's 0
  optical_tick const firstTick { 0 };

  // use hardware trigger time plus the configured offset as waveform start time
  OpDetWaveformMaker_t createOpDetWaveform {
    waveform,
    timings.TriggerTime() + time_interval{ fParams.triggerOffsetPMT },
    1.0 / fSampling
    };
  
  //
  // get the PMT channel triggers to consider
  //
  
  // prepare the set of triggers
  std::vector<optical_tick> const trigger_locations = FindTriggers(waveform);
  auto const tend = trigger_locations.end();
  
  // find the first viable trigger
  auto tooEarlyTrigger = [earliest = firstTick + pretrigSize](optical_tick t)
    { return t < earliest; };
  auto iNextTrigger
    = std::find_if_not(trigger_locations.begin(), tend, tooEarlyTrigger);
  
  //
  // collect all buffer ranges and merge them
  //
  using BufferRange_t = OpDetWaveformMaker_t::BufferRange_t;
  auto makeBuffer
    = [pretrigSize, posttrigSize](optical_tick triggerTime) -> BufferRange_t
    { return { triggerTime - pretrigSize, triggerTime + posttrigSize }; }
    ;
  
  std::vector<BufferRange_t> buffers;
  buffers.reserve(std::distance(iNextTrigger, tend)); // worst case
  
  auto earliestBufferStart { firstTick };
  while (iNextTrigger != tend) {
    
    BufferRange_t const buffer = makeBuffer(*iNextTrigger);
    
    if (buffer.first < earliestBufferStart) { // extend the previous buffer
      assert(!buffers.empty()); // guaranteed because we skipped early triggers
      buffers.back().second = buffer.second;
    }
    else buffers.emplace_back(buffer);
    
    earliestBufferStart = buffer.second;
    
    ++iNextTrigger;
  } // while
  
  //
  // turn each buffer into a waveform
  //
  MF_LOG_TRACE("PMTsimulationAlg")
    << "Channel #" << opChannel << ": " << buffers.size() << " waveforms for "
    << trigger_locations.size() << " triggers"
    ;
  std::vector<raw::OpDetWaveform> output_opdets;
  for (BufferRange_t const& buffer: buffers) {
    
    output_opdets.push_back(createOpDetWaveform(opChannel, buffer));
    
  } // for buffers
  
  return output_opdets;
} // icarus::opdet::PMTsimulationAlg::CreateFixedSizeOpDetWaveforms()


// -----------------------------------------------------------------------------
bool icarus::opdet::PMTsimulationAlg::KicksPhotoelectron() const
  { return CLHEP::RandFlat::shoot(fParams.randomEngine) < fQE; }


// -----------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::AddPhotoelectrons(
  PulseSampling_t const& pulse, Waveform_t& wave, tick const time_bin,
  WaveformValue_t const n
) const {

  if (n == 0.0) return;

  if (n == 1.0) {
    // simple addition
    AddPulseShape(pulse, wave, time_bin, std::plus<ADCcount>());
  }
  else {
    // multiply each `pulse` sample by `n`:
    AddPulseShape(pulse, wave, time_bin, [n](auto a, auto b) { return a+n*b; });
  }

} // icarus::opdet::PMTsimulationAlg::AddPhotoelectrons()


// -----------------------------------------------------------------------------
template <typename Combine>
void icarus::opdet::PMTsimulationAlg::AddPulseShape(
  PulseSampling_t const& pulse, Waveform_t& wave, tick const time_bin,
  Combine combination
  ) const
{
  std::size_t const min = time_bin.value();
  std::size_t const max = std::min(min + pulse.size(), fNsamples);
  if (min >= max) return;

  std::transform(
    wave.begin() + min, wave.begin() + max,
    pulse.begin(),
    wave.begin() + min,
    combination
    );

} // icarus::opdet::PMTsimulationAlg::AddPulseShape()


// -----------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::AddNoise(Waveform_t& wave) const {

  CLHEP::RandGaussQ random
    (*fParams.elecNoiseRandomEngine, 0.0, fParams.ampNoise.value());
  for(auto& sample: wave) {
    ADCcount const noise { static_cast<float>(random.fire()) }; // Gaussian noise
    sample += noise;
  } // for sample

} // PMTsimulationAlg::AddNoise()


// -----------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::AddNoise_faster(Waveform_t& wave) const {

  /*
    * Compared to AddNoise(), we use a somehow faster random generator;
    * to squeeze the CPU cycles, we avoid the CLHEP interface as much as
    * possible; the random number from the engine is immediately converted
    * to single precision, and the rest of the math happens in there as well.
    * No virtual interfaces nor indirection is involved within this function
    * (except for CLHEP random engine). We generate a normal variable _z_
    * (standard deviation 1, mean 0) and we just scale it to the desired
    * standard deviation, not bothering to add the mean offset of 0.
    * Note that unless the random engine is multi-thread safe, this function
    * won't gain anything from multi-threading.
    */
  auto& engine = *fParams.elecNoiseRandomEngine;

  for(auto& sample: wave) {
    sample += fParams.ampNoise * fFastGauss(engine.flat()); // Gaussian noise
  } // for sample

} // PMTsimulationAlg::AddNoise_faster()


// -----------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::AddDarkNoise(Waveform_t& wave) const {
  /*
   * We assume leakage current ("dark noise") is completely stochastic and
   * distributed uniformly in time with a fixed and known rate.
   *
   * In these condition, the time between two consecutive events follows an
   * exponential distribution with as decay constant the same rate.
   * We extact all the leakage events (including the first one) according to
   * that distribution.
   *
   * We follow the "standard" approach of subsampling of tick as for the
   * photoelectron.
   *
   */
  using namespace util::quantities::frequency_literals;

  if (fParams.darkNoiseRate <= 0.0_Hz) return; // no dark noise

  // CLHEP random objects do not understand quantities, so we use scalars;
  // we choose to work with nanosecond
  CLHEP::RandExponential random(*(fParams.darkNoiseRandomEngine),
    (1.0 / fParams.darkNoiseRate).convertInto<nanoseconds>().value());

  // time to stop at: full duration of the waveform
  nanoseconds const maxTime = static_cast<double>(wave.size()) / fSampling;

  // the time of first leakage event:
  nanoseconds darkNoiseTime { random.fire() };

  TimeToTickAndSubtickConverter const toTickAndSubtick(wsp.nSubsamples());

  auto gainFluctuation = makeGainFluctuator();

  MF_LOG_TRACE("PMTsimulationAlg")
    << "Adding dark noise (" << fParams.darkNoiseRate << ") up to " << maxTime;

  while (darkNoiseTime < maxTime) {

    auto const [ tick, subtick ] = toTickAndSubtick(darkNoiseTime * fSampling);

    double const n = gainFluctuation(1.0); // leakage is one photoelectron
    MF_LOG_TRACE("PMTsimulationAlg")
      << " * at " << darkNoiseTime << " (" << tick << ", subsample " << subtick
      << ") x" << n;

    AddPhotoelectrons
      (wsp.subsample(subtick), wave, tick, static_cast<WaveformValue_t>(n));

    // time of the next leakage event:
    darkNoiseTime += nanoseconds{ random.fire() };

  } // while

} // icarus::opdet::PMTsimulationAlg::AddDarkNoise()



// -----------------------------------------------------------------------------
auto icarus::opdet::PMTsimulationAlg::TimeToTickAndSubtickConverter::operator()
  (double const tick_d) const -> std::tuple<tick, SubsampleIndex_t>
{
  double const tickNumber_d = std::floor(tick_d);
  double const subtick = std::floor((tick_d - tickNumber_d) * fNSubsamples);
  return {
    tick::castFrom(tickNumber_d),
    static_cast<SubsampleIndex_t>(subtick)
    };
} // icarus::opdet::PMTsimulationAlg::TimeToTickAndSubtickConverter::operator()


// -----------------------------------------------------------------------------
template <typename Rand>
double icarus::opdet::PMTsimulationAlg::GainFluctuator<Rand>::operator()
  (double const n)
  { return fRandomGain? (n * fRandomGain->fire() / fReferenceGain): n; }


// -----------------------------------------------------------------------------
// ---  icarus::opdet::PMTsimulationAlgMaker
// -----------------------------------------------------------------------------
icarus::opdet::PMTsimulationAlgMaker::PMTsimulationAlgMaker
  (Config const& config)
{
  using util::quantities::microsecond;
  using util::quantities::nanosecond;
  using util::quantities::hertz;
  using util::quantities::megahertz;
  using util::quantities::picocoulomb;
  using ADCcount = icarus::opdet::PMTsimulationAlg::ADCcount;

  //
  // readout settings
  //
  fBaseConfig.readoutEnablePeriod      = config.ReadoutEnablePeriod();
  fBaseConfig.readoutWindowSize        = config.ReadoutWindowSize();
  fBaseConfig.baseline                 = ADCcount(config.Baseline());
  fBaseConfig.pulsePolarity            = config.PulsePolarity();
  fBaseConfig.pretrigFraction          = config.PreTrigFraction();
  fBaseConfig.triggerOffsetPMT         = config.TriggerOffsetPMT();

  //
  // PMT settings
  //
  auto const& PMTspecs = config.PMTspecs();
  fBaseConfig.saturation               = config.Saturation();
  fBaseConfig.QEbase                   = config.QE();
  fBaseConfig.PMTspecs.dynodeK         = PMTspecs.DynodeK();
  fBaseConfig.PMTspecs.setVoltageDistribution
                                        (PMTspecs.VoltageDistribution());
  fBaseConfig.PMTspecs.gain            = PMTspecs.Gain();
  fBaseConfig.doGainFluctuations       = config.FluctuateGain();

  //
  // single photoelectron response
  //
  fBaseConfig.pulseSubsamples          = config.PulseSubsamples();

  //
  // dark noise
  //
  fBaseConfig.darkNoiseRate            = hertz(config.DarkNoiseRate());

  //
  // electronics noise
  //
  fBaseConfig.ampNoise                 = ADCcount(config.AmpNoise());
  fBaseConfig.useFastElectronicsNoise  = config.FastElectronicsNoise();

  //
  // trigger
  //
  fBaseConfig.thresholdADC             = ADCcount(config.ThresholdADC());
  fBaseConfig.createBeamGateTriggers   = config.CreateBeamGateTriggers();
  fBaseConfig.beamGateTriggerRepPeriod = microsecond(config.BeamGateTriggerRepPeriod());
  fBaseConfig.beamGateTriggerNReps     = config.BeamGateTriggerNReps();

  //
  // parameter checks
  //
  if (std::abs(fBaseConfig.pulsePolarity) != 1.0) {
    throw cet::exception("PMTsimulationAlg")
      << "Pulse polarity settings can be only +1.0 or -1.0 (got: "
      << fBaseConfig.pulsePolarity << ")\n"
      ;
  } // check pulse polarity

  if (fBaseConfig.doGainFluctuations) {
    double const mu0 = fBaseConfig.PMTspecs.firstStageGain();
    if (!std::isnormal(mu0) || (mu0 < 0.0)) {
      cet::exception e("PMTsimulationAlg");
      e << "PMT gain " << fBaseConfig.PMTspecs.gain
        << ", dynode resistance values {";
      for (double rho: fBaseConfig.PMTspecs.voltageDistribution)
        e << " " << rho;
      e << " } and dynode k constant " << fBaseConfig.PMTspecs.dynodeK
        << " resulted in an invalid gain " << mu0
        << " for the first multiplication stage!\n";
      throw e;
    }
  } // if invalid gain fluctuation

} // icarus::opdet::PMTsimulationAlgMaker::PMTsimulationAlgMaker()


//-----------------------------------------------------------------------------
std::unique_ptr<icarus::opdet::PMTsimulationAlg>
icarus::opdet::PMTsimulationAlgMaker::operator()(
  detinfo::LArProperties const& larProp,
  detinfo::DetectorClocksData const& clockData,
  SinglePhotonResponseFunc_t const& SPRfunction,
  CLHEP::HepRandomEngine& mainRandomEngine,
  CLHEP::HepRandomEngine& darkNoiseRandomEngine,
  CLHEP::HepRandomEngine& elecNoiseRandomEngine,
  bool trackSelectedPhotons /* = false */
  ) const
{
  return std::make_unique<PMTsimulationAlg>(makeParams(
    larProp, clockData,
    SPRfunction,
    mainRandomEngine, darkNoiseRandomEngine, elecNoiseRandomEngine,
    trackSelectedPhotons
    ));

} // icarus::opdet::PMTsimulationAlgMaker::operator()


//-----------------------------------------------------------------------------
auto icarus::opdet::PMTsimulationAlgMaker::makeParams(
  detinfo::LArProperties const& larProp,
  detinfo::DetectorClocksData const& clockData,
  SinglePhotonResponseFunc_t const& SPRfunction,
  CLHEP::HepRandomEngine& mainRandomEngine,
  CLHEP::HepRandomEngine& darkNoiseRandomEngine,
  CLHEP::HepRandomEngine& elecNoiseRandomEngine,
  bool trackSelectedPhotons /* = false */
  ) const -> PMTsimulationAlg::ConfigurationParameters_t
{
  using namespace util::quantities::electronics_literals;
  
  //
  // set the configuration
  //
  auto params = fBaseConfig;

  //
  // set up parameters
  //
  params.larProp = &larProp;
  params.clockData = &clockData;

  params.pulseFunction = &SPRfunction;

  params.randomEngine = &mainRandomEngine;
  params.gainRandomEngine = params.randomEngine;
  params.darkNoiseRandomEngine = &darkNoiseRandomEngine;
  params.elecNoiseRandomEngine = &elecNoiseRandomEngine;
  
  params.trackSelectedPhotons = trackSelectedPhotons;
  
  //
  // setup checks
  //
  bool const expectedNegativePolarity
    = (SPRfunction.peakAmplitude() < 0.0_ADCf);
  
  if (std::signbit(params.pulsePolarity) != expectedNegativePolarity) {
    throw cet::exception("PMTsimulationAlg")
      << "Inconsistent settings: pulse polarity declared "
      << params.pulsePolarity << ", but photoelectron waveform amplitude is "
      << SPRfunction.peakAmplitude()
      << "\n"
      ;
  } // check polarity consistency

  return params;
  
} // icarus::opdet::PMTsimulationAlgMaker::create()


//-----------------------------------------------------------------------------
