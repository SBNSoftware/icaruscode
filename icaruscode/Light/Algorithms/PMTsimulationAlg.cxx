/**
 * @file    icaruscode/Light/Algorithms/PMTsimulationAlg.cxx
 * @brief   Algorithms for the simulation of ICARUS PMT channels.
 * @date    October 16, 2018
 * @see     `icaruscode/Light/Algorithms/PMTsimulationAlg.h`
 *
 * Implementation file.
 */

// this library header
#include "icaruscode/Light/Algorithms/PMTsimulationAlg.h"

// framework libraries
#include "cetlib_except/exception.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"

// C++ standard libaries
#include <utility> // std::move()


// -----------------------------------------------------------------------------
// ---  icarus::opdet::DiscretePhotoelectronPulse
// -----------------------------------------------------------------------------
auto icarus::opdet::DiscretePhotoelectronPulse::sampleShape(
  PulseFunction_t const& pulseShape, gigahertz samplingFreq, double rightSigmas
) -> std::vector<ADCcount>
{
  std::size_t const pulseSize = samplingFreq
    * (pulseShape.peakTime() + rightSigmas * pulseShape.rightSigma());
  std::vector<ADCcount> samples(pulseSize);
  for (std::size_t i = 0; i < pulseSize; ++i)
    samples[i] = pulseShape(static_cast<double>(i)/samplingFreq);
  return samples;
} // icarus::opdet::DiscretePhotoelectronPulse::sampleShape()


// -----------------------------------------------------------------------------
bool icarus::opdet::DiscretePhotoelectronPulse::checkRange
  (ADCcount limit, std::string const& outputCat /* = "" */) const
{
  assert(pulseLength() > 0);
  bool const bLowOk = (fSampledShape.front().abs() < limit);
  bool const bHighOk = (fSampledShape.back().abs() < limit);
  if (bLowOk && bHighOk) return true;
  if (!outputCat.empty()) {
    mf::LogWarning log(outputCat);
    log << "Check on sampled photoelectron waveform template failed!";
    if (!bLowOk) {
      log
        << "\n => low tail after " << (shape().peakTime() / shape().leftSigma())
          << " standard deviations is still at " << fSampledShape.front()
        ;
    }
    if (!bHighOk) {
      log
        << "\n => high tail after "
          << ((duration() - shape().peakTime()) / shape().rightSigma())
          << " standard deviations is still at " << fSampledShape.back()
        ;
    }
    log << "\nShape parameters:";
    shape().dump(log, "  ", "");
  } // if writing a message on failure
  return false;
} // icarus::opdet::DiscretePhotoelectronPulse::checkRange()


// -----------------------------------------------------------------------------   

// -----------------------------------------------------------------------------   
// ---  icarus::opdet::PMTsimulationAlg
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
    6.0        // 6 std. dev. of tail should suffice
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
  
  printConfiguration(mf::LogDebug("PMTsimulationAlg") << "PMT simulation configuration:\n");
  // check that the sampled waveform has a sufficiently large range, so that
  // tails are below 10^-3 ADC counts (absolute value)
  wsp.checkRange(1e-3_ADCf, "PMTsimulationAlg");
   
} // icarus::opdet::PMTsimulationAlg::setup()


// -----------------------------------------------------------------------------   
std::vector<raw::OpDetWaveform> icarus::opdet::PMTsimulationAlg::simulate
  (sim::SimPhotons const& photons)
{
  std::vector<raw::OpDetWaveform> waveforms; // storage of the results
  
  Waveform_t waveform;
  CreateFullWaveform(waveform, photons);
  CreateOpDetWaveforms(photons.OpChannel(), waveform, waveforms);
  return waveforms;
  
} // icarus::opdet::PMTsimulationAlg::simulate()


//------------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::CreateFullWaveform(Waveform_t & waveform,
					sim::SimPhotons const& photons){
    using util::quantities::tick;
    using namespace util::quantities::time_literals;
    
    waveform.resize(fNsamples,fParams.baseline);
    
    // collect the amount of photoelectrons arriving at each tick
    std::unordered_map<tick,unsigned int> peMap;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for(auto const& ph : photons) {
      if (!KicksPhotoelectron()) continue;
      
      auto const mytime = microsecond(
        fParams.timeService->G4ToElecTime(nanosecond(ph.Time)+fParams.transitTime)
        - fParams.timeService->TriggerTime()
        )
        - fParams.triggerOffsetPMT
        ;
      if ((mytime < 0.0_us) || (mytime >= fParams.readoutEnablePeriod)) continue;
      
      tick const iSample { static_cast<std::size_t>(mytime * fSampling) };
      if (iSample >= fNsamples) continue;
      ++peMap[iSample];
    } // for photons
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    //std::cout << "\tcollected pes... " << photons.OpChannel() << " " << diff.count() << std::endl;
    start=std::chrono::high_resolution_clock::now();

      // add the collected photoelectrons to the waveform
    unsigned int nTotalPE [[gnu::unused]] = 0U; // unused if not in `debug` mode
    for(auto const& pe : peMap){
      auto const nPE = pe.second;
      nTotalPE += nPE;
      if (nPE == 0) continue;
      if (nPE == 1) AddSPE(pe.first,waveform); // faster if n = 1
      else AddPhotoelectrons(pe.first, nPE, waveform);
    }
    LOG_TRACE("PMTsimulationAlg") 
      << nTotalPE << " photoelectrons at " << peMap.size()
      << " times in channel " << photons.OpChannel();

      end=std::chrono::high_resolution_clock::now(); diff = end-start;
      //std::cout << "\tadded pes... " << photons.OpChannel() << " " << diff.count() << std::endl;
      start=std::chrono::high_resolution_clock::now();

      if(fParams.ampNoise>0.) AddNoise(waveform);
      if(fParams.darkNoiseRate>0.) AddDarkNoise(waveform);

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

  std::set<size_t> icarus::opdet::PMTsimulationAlg::CreateBeamGateTriggers() const
  {
    using namespace util::quantities::time_literals;
    
    std::set<size_t> trigger_locations;
    microsecond const deltaBeamTime
      { fParams.timeService->BeamGateTime() - fParams.timeService->TriggerTime() };

    for(unsigned int i_trig=0; i_trig<fParams.beamGateTriggerNReps; ++i_trig) {
      microsecond const trig_time
        = deltaBeamTime
        + i_trig * fParams.beamGateTriggerRepPeriod
        - fParams.triggerOffsetPMT
        ;
      if (trig_time < 0_us || trig_time > fParams.readoutEnablePeriod) continue;
      trigger_locations.insert(size_t(trig_time*fSampling));
    }
    return trigger_locations;
  }

  std::set<size_t> icarus::opdet::PMTsimulationAlg::FindTriggers(Waveform_t const& wvfm) const
  {
    std::set<size_t> trigger_locations;
    if (fParams.createBeamGateTriggers) trigger_locations = CreateBeamGateTriggers();
    
    bool above_thresh=false;

    //next, find all ticks at which we would trigger readout
    for(size_t i_t=0; i_t<wvfm.size(); ++i_t){
      
      auto const val { fParams.pulsePolarity* (wvfm[i_t]-fParams.baseline) };

      if(!above_thresh && val>=fParams.thresholdADC){
	above_thresh=true;
	trigger_locations.insert(i_t);
      }
      else if(above_thresh && val<fParams.thresholdADC){
	above_thresh=false;
      } 
      
    }//end loop over waveform   
    return trigger_locations;
  }

//------------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::CreateOpDetWaveforms(raw::Channel_t const& opch,
					  Waveform_t const& wvfm,
					  std::vector<raw::OpDetWaveform>& output_opdets)
  {
    //std::cout << "Finding triggers in " << opch << std::endl;

    std::set<size_t> trigger_locations = FindTriggers(wvfm);

    //std::cout << "Creating opdet waveforms in " << opch << std::endl;

    bool in_pulse=false;
    size_t trig_start=0,trig_stop=wvfm.size();

    auto const pretrigSize = fParams.pretrigSize();
    auto const posttrigSize = fParams.posttrigSize();
    LOG_TRACE("PMTsimulationAlg")
      << "Channel #" << opch << ": " << trigger_locations.size() << " triggers";
    for(size_t i_t=0; i_t<wvfm.size(); ++i_t){

      //if we are at a trigger point, open the window
      if(trigger_locations.count(i_t)==1){

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
	output_opdets.emplace_back( raw::TimeStamp_t((trig_start/fSampling + fParams.triggerOffsetPMT)),
				    opch,
				    trig_stop-trig_start );
	output_opdets.back().Waveform().assign(wvfm.begin()+trig_start,wvfm.begin()+trig_stop);
	in_pulse=false;
      }
      
    }//end loop over waveform
  }

  bool icarus::opdet::PMTsimulationAlg::KicksPhotoelectron() const
    { return CLHEP::RandFlat::shoot(fParams.randomEngine) < fQE; } 

  
  void icarus::opdet::PMTsimulationAlg::AddPhotoelectrons(size_t time_bin, unsigned int n, Waveform_t& wave) const {
    
    if (time_bin >= fNsamples) return;
    
    std::size_t const min = time_bin;
    std::size_t const max = std::min(time_bin + wsp.pulseLength(), fNsamples);

    std::transform(wave.begin()+min,wave.begin()+max,wsp.begin(),wave.begin()+min,
		   [n](auto a, auto b) { return a+n*b; });
		   //addmultiple<n>());

  } // PMTsimulationAlg::AddPhotoelectrons()
  
  void icarus::opdet::PMTsimulationAlg::AddSPE(size_t time_bin, Waveform_t& wave){
     
    if (time_bin >= fNsamples) return;
    
    std::size_t const min = time_bin;
    std::size_t const max = std::min(time_bin + wsp.pulseLength(), fNsamples);

    std::transform(wave.begin()+min,wave.begin()+max,wsp.begin(),wave.begin()+min,std::plus<ADCcount>());
    
  }
  
  void icarus::opdet::PMTsimulationAlg::AddNoise(Waveform_t& wave){
    
    CLHEP::RandGauss random(*fParams.elecNoiseRandomEngine, 0.0, fParams.ampNoise);
    for(auto& sample: wave) {
      ADCcount const noise { static_cast<float>(random.fire()) }; //gaussian noise
      sample += noise;
    } // for sample
    
  } // PMTsimulationAlg::AddNoise()
  
  void icarus::opdet::PMTsimulationAlg::AddDarkNoise(Waveform_t& wave)
  {
    using util::quantities::nanosecond;
    if (fParams.darkNoiseRate <= 0.0) return; // no dark noise
    size_t timeBin=0;
    CLHEP::RandExponential random(*(fParams.darkNoiseRandomEngine),
      (1.0/fParams.darkNoiseRate).convertInto<nanosecond>());
    double darkNoiseTime = random.fire();
    while (darkNoiseTime < wave.size()){
      timeBin = (darkNoiseTime);
      AddSPE(timeBin,wave);
      // Find next time to add dark noise
      darkNoiseTime += random.fire();
    }
  }
  


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
  fBaseConfig.readoutEnablePeriod      = microsecond(config.ReadoutEnablePeriod());
  fBaseConfig.readoutWindowSize        = config.ReadoutWindowSize();
  fBaseConfig.ADC                      = config.ADC();
  fBaseConfig.baseline                 = ADCcount(config.Baseline());
  fBaseConfig.pulsePolarity            = config.PulsePolarity();
  fBaseConfig.pretrigFraction          = config.PreTrigFraction();
  fBaseConfig.triggerOffsetPMT         = microsecond(config.TriggerOffsetPMT());
  
  //
  // PMT settings
  //
  fBaseConfig.saturation               = config.Saturation();
  fBaseConfig.QEbase                   = config.QE();
  
  //
  // single photoelectron response
  //
  fBaseConfig.transitTime              = nanosecond(config.TransitTime());
  fBaseConfig.riseTime                 = nanosecond(config.RiseTime());
  fBaseConfig.fallTime                 = nanosecond(config.FallTime());
  fBaseConfig.meanAmplitude            = picocoulomb(config.MeanAmplitude());
  
  //
  // dark noise
  //
  fBaseConfig.darkNoiseRate            = hertz(config.DarkNoiseRate());
  
  //
  // electronics noise
  //
  fBaseConfig.ampNoise                 = ADCcount(config.AmpNoise());
  
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
  params.darkNoiseRandomEngine = &darkNoiseRandomEngine;
  params.elecNoiseRandomEngine = &elecNoiseRandomEngine;

  return std::make_unique<PMTsimulationAlg>(params);
   
} // icarus::opdet::PMTsimulationAlgMaker::create()

