/**
 * @file    icaruscode/Light/Algorithms/PMTsimulationAlg.cxx
 * @brief   Algorithms for the simulation of ICARUS PMT channels.
 * @date    October 16, 2018
 * @see     `icaruscode/Light/Algorithms/PMTsimulationAlg.h`
 *
 * Implementation file.
 */

// this library header
#include "icaruscode/PMT/Algorithms/PMTsimulationAlg.h"

// LArSoft libraries
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "larcorealg/CoreUtils/counter.h"

// framework libraries
#include "cetlib_except/exception.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandExponential.h"

// C++ standard libaries
#include <utility> // std::move()


// -----------------------------------------------------------------------------
// ---  icarus::opdet::DiscretePhotoelectronPulse
// -----------------------------------------------------------------------------
auto icarus::opdet::DiscretePhotoelectronPulse::sampleShape(
  PulseFunction_t const& pulseShape, gigahertz samplingFreq, double rightSigmas,
  unsigned int nSubsamples
) -> SampledFunction_t
{
  using namespace util::quantities::time_literals;

  nanoseconds const pulseDuration
    = pulseShape.peakTime() + rightSigmas * pulseShape.rightSigma();
  std::size_t const pulseSize = samplingFreq * pulseDuration;
  nanoseconds const roundPulseDuration = pulseSize / samplingFreq;

  return SampledFunction_t{
    pulseShape,
    0.0_ns, roundPulseDuration, // time range of the pulse shape
    static_cast<gsl::index>(pulseSize), // number of samples
    static_cast<gsl::index>(nSubsamples) // how many subsamples per tick
    };

} // icarus::opdet::DiscretePhotoelectronPulse::sampleShape()


// -----------------------------------------------------------------------------
bool icarus::opdet::DiscretePhotoelectronPulse::checkRange
  (ADCcount limit, std::string const& outputCat /* = "" */) const
{
  assert(pulseLength() > 0);
  auto const low = *(fSampledShape.subsample(0).begin());
  auto const high
    = *(fSampledShape.subsample(fSampledShape.nSubsamples() - 1).rbegin());

  bool const bLowOk = (low.abs() < limit);
  bool const bHighOk = (high.abs() < limit);
  if (bLowOk && bHighOk) return true;
  if (!outputCat.empty()) {
    mf::LogWarning log(outputCat);
    log << "Check on sampled photoelectron waveform template failed!";
    if (!bLowOk) {
      log
        << "\n => low tail after " << (shape().peakTime() / shape().leftSigma())
          << " standard deviations is still at " << low
        ;
    }
    if (!bHighOk) {
      log
        << "\n => high tail after "
          << ((duration() - shape().peakTime()) / shape().rightSigma())
          << " standard deviations is still at " << high
        ;
    }
    log << "\nShape parameters:";
    shape().dump(log, "  ", "");
  } // if writing a message on failure
  return false;
} // icarus::opdet::DiscretePhotoelectronPulse::checkRange()



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
  , fSampling(fParams.timeService->OpticalClock().Frequency())
  , fNsamples(fParams.readoutEnablePeriod * fSampling) // us * MHz cancels out
  , wsp(
    { // PhotoelectronPulseWaveform
      // amplitude is a charge, so we have to wrap the arm of the constructor to
      // accept it as ADC count (`value()` makes `meanAmplitude` lose its unit)
      // `ADCcount` conversion is redundant but left for clarity
      ADCcount(fParams.ADC * fParams.meanAmplitude.value()), // amplitude
      fParams.transitTime,                 // peak time
      riseTimeToRMS(fParams.riseTime),     // sigma left
      riseTimeToRMS(fParams.fallTime)      // sigma right
    },
    fSampling,
    6.0,                    // 6 std. dev. of tail should suffice
    fParams.pulseSubsamples // tick subsampling
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

  printConfiguration
    (mf::LogDebug("PMTsimulationAlg") << "PMT simulation configuration:\n");

  // check that the sampled waveform has a sufficiently large range, so that
  // tails are below 10^-3 ADC counts (absolute value)
  wsp.checkRange(1e-3_ADCf, "PMTsimulationAlg");

} // icarus::opdet::PMTsimulationAlg::PMTsimulationAlg()


// -----------------------------------------------------------------------------
std::vector<raw::OpDetWaveform> icarus::opdet::PMTsimulationAlg::simulate
(sim::SimPhotons const& photons, sim::SimPhotons &photons_used)
{
  std::vector<raw::OpDetWaveform> waveforms; // storage of the results

  Waveform_t waveform;
  CreateFullWaveform(waveform, photons, photons_used);
  CreateOpDetWaveforms(photons.OpChannel(), waveform, waveforms);
  return waveforms;

} // icarus::opdet::PMTsimulationAlg::simulate()


//------------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::CreateFullWaveform(Waveform_t & waveform,
							 sim::SimPhotons const& photons,
							 sim::SimPhotons& photons_used)
{

    using namespace util::quantities::time_literals;
    using namespace util::quantities::frequency_literals;
    using namespace util::quantities::electronics_literals;
    using namespace detinfo::timescales;

    detinfo::DetectorTimings const& timings
      = detinfo::makeDetectorTimings(fParams.timeService);

    waveform.resize(fNsamples,fParams.baseline);
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

    auto start = std::chrono::high_resolution_clock::now();

    photons_used.clear();
    photons_used.SetChannel(photons.OpChannel());
    for(auto const& ph : photons) {
      if (!KicksPhotoelectron()) continue;

      photons_used.push_back(ph);

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

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    //std::cout << "\tcollected pes... " << photons.OpChannel() << " " << diff.count() << std::endl;
    start=std::chrono::high_resolution_clock::now();

    //
    // add the collected photoelectrons to the waveform
    //
    unsigned int nTotalPE [[gnu::unused]] = 0U; // unused if not in `debug` mode
    double nTotalEffectivePE [[gnu::unused]] = 0U; // unused if not in `debug` mode

    double const refGain = fParams.PMTspecs.firstStageGain();
    CLHEP::RandPoisson randomGainFluctuation
      (*fParams.gainRandomEngine, refGain);

    // go though all subsamples (starting each at a fraction of a tick)
    for (auto const& [ iSubsample, peMap ]: util::enumerate(peMaps)) {
      
      // this is the waveform sampling for the selected subsample:
      auto const& subsample = wsp.subsample(iSubsample);
      
      for (auto const& [ startTick, nPE ]: peMap) {
        nTotalPE += nPE;

        // add gain fluctuations in the conversion
        double nEffectivePE = nPE;
        if (fParams.doGainFluctuations) {
          nEffectivePE *= randomGainFluctuation.fire() / refGain;
        }
        nTotalEffectivePE += nEffectivePE;

        if (nEffectivePE == 0.0) continue;
        if (nEffectivePE == 1.0) AddSPE(subsample, waveform, startTick);
        else {
          AddPhotoelectrons(
            subsample, waveform, startTick,
            static_cast<WaveformValue_t>(nEffectivePE)
            );
        }
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

      end=std::chrono::high_resolution_clock::now(); diff = end-start;
      //std::cout << "\tadded pes... " << photons.OpChannel() << " " << diff.count() << std::endl;
      start=std::chrono::high_resolution_clock::now();

      if(fParams.ampNoise > 0.0_ADCf) (this->*fNoiseAdder)(waveform);
      if(fParams.darkNoiseRate > 0.0_Hz) AddDarkNoise(waveform);

      end=std::chrono::high_resolution_clock::now(); diff = end-start;
      //std::cout << "\tadded noise... " << photons.OpChannel() << " " << diff.count() << std::endl;
      start=std::chrono::high_resolution_clock::now();

      // Implementing saturation effects;
      // waveform is negative, and saturation is a minimum ADC count
      auto const saturationLevel
        = fParams.baseline + fParams.saturation*wsp.shape().amplitude();
      std::replace_if(waveform.begin(),waveform.end(),
		      [saturationLevel](auto s) -> bool{return s < saturationLevel;},
		      saturationLevel);

      end=std::chrono::high_resolution_clock::now(); diff = end-start;
      //std::cout << "\tadded saturation... " << photons.OpChannel() << " " << diff.count() << std::endl;

  } // CreateFullWaveform()

  auto icarus::opdet::PMTsimulationAlg::CreateBeamGateTriggers() const
    -> std::vector<optical_tick>
  {
    using namespace util::quantities::time_literals;
    using detinfo::timescales::trigger_time;

    detinfo::DetectorTimings const& timings
      = detinfo::makeDetectorTimings(fParams.timeService);

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
void icarus::opdet::PMTsimulationAlg::CreateOpDetWaveforms(raw::Channel_t const& opch,
					  Waveform_t const& wvfm,
					  std::vector<raw::OpDetWaveform>& output_opdets)
  {
    //std::cout << "Finding triggers in " << opch << std::endl;

    using namespace detinfo::timescales; // electronics_time, time_interval

    detinfo::DetectorTimings const& timings
      = detinfo::makeDetectorTimings(fParams.timeService);

    electronics_time const PMTstartTime
      = timings.TriggerTime() + time_interval{ fParams.triggerOffsetPMT };


    std::vector<optical_tick> trigger_locations = FindTriggers(wvfm);
    auto iNextTrigger = trigger_locations.begin();
    optical_tick nextTrigger
      = trigger_locations.empty()
      ? std::numeric_limits<optical_tick>::max()
      : *iNextTrigger
      ;

    //std::cout << "Creating opdet waveforms in " << opch << std::endl;

    bool in_pulse=false;
    size_t trig_start=0,trig_stop=wvfm.size();

    auto const pretrigSize = fParams.pretrigSize();
    auto const posttrigSize = fParams.posttrigSize();
    MF_LOG_TRACE("PMTsimulationAlg")
      << "Channel #" << opch << ": " << trigger_locations.size() << " triggers";
    for (std::size_t const i_t: util::counter(wvfm.size())) {

      auto const thisTick = optical_tick::castFrom(i_t);

      //if we are at a trigger point, open the window
      if (thisTick == nextTrigger) {

        // update the next trigger
        nextTrigger
          = (++iNextTrigger == trigger_locations.end())
          ? std::numeric_limits<optical_tick>::max()
          : *iNextTrigger
          ;


	//if not already in a pulse
	if(!in_pulse){
	  in_pulse=true;
	  trig_start = i_t>pretrigSize ? i_t-pretrigSize : 0;
	  trig_stop  = (wvfm.size()-1-i_t)>posttrigSize ? i_t+posttrigSize : wvfm.size();

	}
	//else, if we are already in a pulse, extend it
	else if(in_pulse){
	  trig_stop  = (wvfm.size()-1-i_t)>posttrigSize ? i_t+posttrigSize : wvfm.size();
	}
      }

      //ok, now, if we are in a pulse but have reached its end, store the waveform
      if(in_pulse && i_t==trig_stop-1){
	// time should be absolute (on electronics time scale)
	output_opdets.emplace_back(
	     // start of the waveform (tick #0) in the full optical reading
	  raw::TimeStamp_t{ PMTstartTime + trig_start/fSampling },
				    opch,
				    trig_stop-trig_start );
	auto& outputWaveform = output_opdets.back().Waveform();
	outputWaveform.reserve(trig_stop - trig_start);
	std::transform(wvfm.begin()+trig_start,wvfm.begin()+trig_stop,
	  std::back_inserter(outputWaveform),
	  [](auto ADC){ return ADC.value(); }
	  );
	in_pulse=false;
      }

    } // for i_t (loop over waveform)
  } // icarus::opdet::PMTsimulationAlg::CreateOpDetWaveforms()

  bool icarus::opdet::PMTsimulationAlg::KicksPhotoelectron() const
    { return CLHEP::RandFlat::shoot(fParams.randomEngine) < fQE; }


  void icarus::opdet::PMTsimulationAlg::AddPhotoelectrons(
    PulseSampling_t const& pulse, Waveform_t& wave, tick const time_bin,
    WaveformValue_t const n
  ) const {


    std::size_t const min = time_bin.value();
    std::size_t const max = std::min(min + pulse.size(), fNsamples);
    if (min >= max) return;

    std::transform(
      wave.begin() + min, wave.begin() + max,
      pulse.begin(),
      wave.begin() + min,
      [n](auto a, auto b) { return a+n*b; }
      );

  } // PMTsimulationAlg::AddPhotoelectrons()

  void icarus::opdet::PMTsimulationAlg::AddSPE
    (PulseSampling_t const& pulse, Waveform_t& wave, tick const time_bin) const
  {
    std::size_t const min = time_bin.value();
    std::size_t const max = std::min(min + pulse.size(), fNsamples);
    if (min >= max) return;

    std::transform(
      wave.begin() + min, wave.begin() + max,
      pulse.begin(),
      wave.begin() + min,
      std::plus<ADCcount>()
      );

  }


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

  MF_LOG_TRACE("PMTsimulationAlg")
    << "Adding dark noise (" << fParams.darkNoiseRate << ") up to " << maxTime;
  
  while (darkNoiseTime < maxTime) {
    
    auto const [ tick, subtick ] = toTickAndSubtick(darkNoiseTime * fSampling);
    
    double const n = 1.0;
    MF_LOG_TRACE("PMTsimulationAlg")
      << " * at " << darkNoiseTime << " (" << tick << ", subsample " << subtick
      << ") x" << n;
    
    // TODO: add gain fluctuation
    AddSPE(wsp.subsample(subtick), wave, tick);
    
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
  fBaseConfig.ADC                      = config.ADC();
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
  fBaseConfig.transitTime              = config.TransitTime();
  fBaseConfig.riseTime                 = config.RiseTime();
  fBaseConfig.fallTime                 = config.FallTime();
  fBaseConfig.meanAmplitude            = config.MeanAmplitude();
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

  if (fBaseConfig.pulsePolarity != fBaseConfig.expectedPulsePolarity())
  {
    throw cet::exception("PMTsimulationAlg")
      << "Inconsistent settings: pulse polarity (" << fBaseConfig.pulsePolarity
      << "), photoelectron waveform amplitude (" << fBaseConfig.meanAmplitude
      << ") and ADC-per-charge calibration factor (" << fBaseConfig.ADC
      << ")\n"
      ;
  } // check polarity consistency

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
  detinfo::DetectorClocks const& detClocks,
  CLHEP::HepRandomEngine& mainRandomEngine,
  CLHEP::HepRandomEngine& darkNoiseRandomEngine,
  CLHEP::HepRandomEngine& elecNoiseRandomEngine
  ) const
{
  // set the configuration
  auto params = fBaseConfig;

  // set up parameters
  params.larProp = &larProp;
  params.timeService = &detClocks;

  params.randomEngine = &mainRandomEngine;
  params.gainRandomEngine = params.randomEngine;
  params.darkNoiseRandomEngine = &darkNoiseRandomEngine;
  params.elecNoiseRandomEngine = &elecNoiseRandomEngine;

  return std::make_unique<PMTsimulationAlg>(params);

} // icarus::opdet::PMTsimulationAlgMaker::create()

