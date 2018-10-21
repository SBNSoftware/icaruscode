/**
 * @file    icaruscode/Light/Algorithms/PMTsimulationAlg.cxx
 * @brief   Algorithms for the simulation of ICARUS PMT channels.
 * @date    October 16, 2018
 * @see     icaruscode/Light/Algorithms/PMTsimulationAlg.h
 *
 * Implementation file.
 */

// this library header
#include "icaruscode/Light/Algorithms/PMTsimulationAlg.h"


// -----------------------------------------------------------------------------   
// icarus::opdet::PMTsimulationAlg::PMTsimulationAlg(Config const& config)
//   : fReadoutEnablePeriod(config.ReadoutEnablePeriod())
//   {}


// -----------------------------------------------------------------------------   
icarus::opdet::PMTsimulationAlg::PMTsimulationAlg
  (fhicl::ParameterSet const& p)
{
//     fInputModuleName = p.get< art::InputTag >("InputModule" );
    fTransitTime     = p.get< double >("TransitTime"  ); //ns
    fADC             = p.get< double >("ADC"          ); //voltage to ADC factor
    fBaseline        = p.get< double >("Baseline"     ); //in ADC counts (may be fractional)
    fFallTime        = p.get< double >("FallTime"     ); //in ns
    fRiseTime        = p.get< double >("RiseTime"     ); //in ns
    fMeanAmplitude   = p.get< double >("MeanAmplitude"); //in pC
    fAmpNoise        = p.get< double >("AmpNoise"     ); //in ADC
    fDarkNoiseRate   = p.get< double >("DarkNoiseRate"); //in Hz
    
    fReadoutWindowSize   = p.get<size_t>("ReadoutWindowSize"); ///ReadoutWindowSize
    fPretrigFraction     = p.get<float>("PreTrigFraction");    ///Fraction of window size to be before "trigger"
    fThresholdADC        = p.get<float>("ThresholdADC");       ///ADC Threshold for self-triggered readout
    fPulsePolarity       = p.get<int>("PulsePolarity");        ///Pulse polarity (=1 for positive, =-1 for negative)
//     TODO is this given by DetectorClocks? should it?
    fTriggerOffsetPMT    = p.get<double>("TriggerOffsetPMT");   ///Time (us) relative to trigger that readout begins
    fReadoutEnablePeriod = p.get<double>("ReadoutEnablePeriod"); ///Time (us) for which pmt readout is enabled

    fCreateBeamGateTriggers = p.get<bool>("CreateBeamGateTriggers"); ///Option to create unbiased readout around beam spill
    fBeamGateTriggerRepPeriod = p.get<double>("BeamGateTriggerRepPeriod"); ///Repetition Period (us) for BeamGateTriggers
    fBeamGateTriggerNReps = p.get<size_t>("BeamGateTriggerNReps"); ///Number of beamgate trigger reps to produce
    
    fSaturation      = p.get< double >("Saturation"   ); //in number of p.e.
    fQEbase         = p.get< double >("QE"           ); //PMT quantum efficiency
// 
    fPretrigSize = fPretrigFraction*fReadoutWindowSize;
    fPosttrigSize = fReadoutWindowSize-fPretrigSize;
   
} // PMTsimulationAlg::PMTsimulationAlg()


//-----------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::setup(
  detinfo::LArProperties const& larProp,
  detinfo::DetectorClocks const& detClocks,
  CLHEP::HepRandomEngine& mainRandomEngine, 
  CLHEP::HepRandomEngine& darkNoiseRandomEngine, 
  CLHEP::HepRandomEngine& elecNoiseRandomEngine 
  )
{
  fTimeService = &detClocks;
  
  fRandomEngine = &mainRandomEngine;
  fDarkNoiseRandomEngine = &darkNoiseRandomEngine;
  fElecNoiseRandomEngine = &elecNoiseRandomEngine;

  fSampling = fTimeService->OpticalClock().Frequency();
  fNsamples = fReadoutEnablePeriod * fSampling ; // us * MHz cancels out
  //  mf::LogDebug("PMTsimulationAlg") << "Sampling = " << fSampling << " MHz." << std::endl;
  
  // shape of single pulse
  sigma1 = fRiseTime / 1.687;
  sigma2 = fFallTime / 1.687;
    
  pulsesize = (int)((6*sigma2+fTransitTime)*1.e3*fSampling);
  wsp.resize(pulsesize);
  for (int i = 0; i<pulsesize; i++)
    wsp[i] = Pulse1PE(static_cast<double>(i)*1.e3/fSampling);
    
  
  // Correction due to scalling factor applied during simulation if any
  fQE = fQEbase / larProp.ScintPreScale();
  mf::LogDebug("SimPMTICARUS") << "PMT corrected efficiency = " << fQE;
    
  if (fQE >= 1.0001) {
    mf::LogWarning("SimPMTICARUS") << "WARNING: Quantum efficiency set in fhicl file "
				   << fQEbase
				   << " seems to be too large! Final QE must be equal"
				   << " or smaller than the scintillation pre scale applied"
				   << " at simulation time. Please check this number (ScintPreScale): "
				   << larProp.ScintPreScale();
  }
  
  printConfiguration(mf::LogDebug("PMTsimulationAlg") << "PMT simulation configuration:\n");
   
} // icarus::opdet::PMTsimulationAlg::setup()


//-----------------------------------------------------------------------------
std::vector<raw::OpDetWaveform> icarus::opdet::PMTsimulationAlg::simulate
  (sim::SimPhotons const& photons)
{
  std::vector<raw::OpDetWaveform> waveforms; // storage of the results
  
  Waveform_t waveform;
  std::vector<unsigned int> PhotoelectronsPerSample;
  CreateFullWaveform(waveform, PhotoelectronsPerSample, photons);
  CreateOpDetWaveforms(photons.OpChannel(), waveform, waveforms);
  return waveforms;
  
} // icarus::opdet::PMTsimulationAlg::simulate()


//------------------------------------------------------------------------------
void icarus::opdet::PMTsimulationAlg::CreateFullWaveform(Waveform_t & waveform,
					std::vector<unsigned int> & PhotoelectronsPerSample,
					sim::SimPhotons const& photons){

    //auto& waveform = fFullWaveforms[opch];
    waveform.resize(fNsamples,fBaseline);
    PhotoelectronsPerSample.resize(fNsamples,0U);
    std::unordered_map<unsigned int,unsigned int> peMap;
    // collect the amount of photoelectrons arriving at each tick
    //std::vector<unsigned int> PhotoelectronsPerSample(fNsamples, 0U);
    
    //for(auto const& photons : pmtVector){
    //if(raw::Channel_t(photons.OpChannel())!=opch) continue;
    
    //std::cout << "Creating waveform " << photons.OpChannel() << " " << photons.size() << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    //auto& waveform = fFullWaveforms[photons.OpChannel()];
    //waveform.resize(fNsamples,fBaseline);
    
    for(auto const& ph : photons) {
      if (!KicksPhotoelectron()) continue;
      
      double const mytime = fTimeService->G4ToElecTime(ph.Time+fTransitTime)-fTimeService->TriggerTime()-fTriggerOffsetPMT;
      if ((mytime < 0.0) || (mytime >= fReadoutEnablePeriod)) continue;
      
      std::size_t const iSample = static_cast<std::size_t>(mytime * fSampling);
      if (iSample >= fNsamples) continue;
      //++PhotoelectronsPerSample[iSample];
      ++peMap[iSample];
    } // for photons
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    //std::cout << "\tcollected pes... " << photons.OpChannel() << " " << diff.count() << std::endl;
    start=std::chrono::high_resolution_clock::now();

      // add the collected photoelectrons to the waveform
      //for (std::size_t iSample = 0; iSample < fNsamples; ++iSample) {
    unsigned int nTotalPE = 0U;
    for(auto const& pe : peMap){
      auto const nPE = pe.second;//PhotoelectronsPerSample[iSample];
      nTotalPE += nPE;
      if (nPE == 0) continue;
      if (nPE == 1) AddSPE(pe.first,waveform);//AddSPE(iSample, waveform); // faster if n = 1
      else AddPhotoelectrons(pe.first, nPE, waveform);//AddPhotoelectrons(iSample, nPE, waveform);
    }
    std::cout << nTotalPE << " photoelectrons at " << peMap.size() << " times in channel " << photons.OpChannel() << std::endl;

      end=std::chrono::high_resolution_clock::now(); diff = end-start;
      //std::cout << "\tadded pes... " << photons.OpChannel() << " " << diff.count() << std::endl;
      start=std::chrono::high_resolution_clock::now();

      if(fAmpNoise>0.) AddNoise(waveform);
      if(fDarkNoiseRate>0.) AddDarkNoise(waveform);

      end=std::chrono::high_resolution_clock::now(); diff = end-start;
      //std::cout << "\tadded noise... " << photons.OpChannel() << " " << diff.count() << std::endl;
      start=std::chrono::high_resolution_clock::now();

      // Implementing saturation effects;
      // waveform is negative, and saturation is a minimum ADC count
      auto const saturationLevel = fBaseline + fSaturation*fADC*fMeanAmplitude;
      std::replace_if(waveform.begin(),waveform.end(),
		      [saturationLevel](float s) -> bool{return s < saturationLevel;},
		      saturationLevel);
      //for (auto& sample: waveform) {
      //if (sample < saturationLevel) sample = saturationLevel;
      //}
      //}
      
      end=std::chrono::high_resolution_clock::now(); diff = end-start;
      //std::cout << "\tadded saturation... " << photons.OpChannel() << " " << diff.count() << std::endl;

  } // CreateFullWaveform()

  std::set<size_t> icarus::opdet::PMTsimulationAlg::CreateBeamGateTriggers() const
  {
    double trig_time;
    std::set<size_t> trigger_locations;
    for(size_t i_trig=0; i_trig<fBeamGateTriggerNReps; ++i_trig){
      trig_time = (fTimeService->BeamGateTime()-fTimeService->TriggerTime())+fBeamGateTriggerRepPeriod*i_trig-fTriggerOffsetPMT;
      if(trig_time<0 || trig_time>fReadoutEnablePeriod) continue;
      trigger_locations.insert(size_t(trig_time*fSampling));
    }
    return trigger_locations;
  }

  std::set<size_t> icarus::opdet::PMTsimulationAlg::FindTriggers(Waveform_t const& wvfm) const
  {
    std::set<size_t> trigger_locations;
    if (fCreateBeamGateTriggers) trigger_locations = CreateBeamGateTriggers();
    
    short val;
    bool above_thresh=false;

    //next, find all ticks at which we would trigger readout
    for(size_t i_t=0; i_t<wvfm.size(); ++i_t){
      
      val = fPulsePolarity*(short)(wvfm[i_t]-fBaseline);

      if(!above_thresh && val>=fThresholdADC){
	above_thresh=true;
	trigger_locations.insert(i_t);
      }
      else if(above_thresh && val<fThresholdADC){
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

    std::cout << "Channel #" << opch << ": " << trigger_locations.size() << " triggers" << std::endl;
    for(size_t i_t=0; i_t<wvfm.size(); ++i_t){

      //if we are at a trigger point, open the window
      if(trigger_locations.count(i_t)==1){

	//if not already in a pulse
	if(!in_pulse){
	  in_pulse=true;
	  trig_start = i_t>fPretrigSize ? i_t-fPretrigSize : 0;
	  trig_stop  = (wvfm.size()-1-i_t)>fPosttrigSize ? i_t+fPosttrigSize : wvfm.size();
	  
	}
	//else, if we are already in a pulse, extend it
	else if(in_pulse){
	  trig_stop  = (wvfm.size()-1-i_t)>fPosttrigSize ? i_t+fPosttrigSize : wvfm.size();
	}
      }

      //ok, now, if we are in a pulse but have reached its end, store the waveform
      if(in_pulse && i_t==trig_stop-1){
	output_opdets.emplace_back( raw::TimeStamp_t((trig_start/fSampling + fTriggerOffsetPMT)),
				    opch,
				    trig_stop-trig_start );
	output_opdets.back().Waveform().assign(wvfm.begin()+trig_start,wvfm.begin()+trig_stop);
	in_pulse=false;
      }
      
    }//end loop over waveform
  }

  bool icarus::opdet::PMTsimulationAlg::KicksPhotoelectron() const
    { return CLHEP::RandFlat::shoot(fRandomEngine) < fQE; } 

  
  double icarus::opdet::PMTsimulationAlg::Pulse1PE(double time) const//single pulse waveform
  {
    double const sigma = (time < fTransitTime)? sigma1: sigma2;
    return (fADC*fMeanAmplitude*std::exp(-sqr(time - fTransitTime)/(2.0*sqr(sigma))));
  }

  void icarus::opdet::PMTsimulationAlg::AddPhotoelectrons(size_t time_bin, unsigned int n, Waveform_t& wave) const {
    
    if (time_bin >= fNsamples) return;
    
    std::size_t const min = time_bin;
    std::size_t const max = std::min(time_bin + pulsesize, fNsamples);

    std::transform(wave.begin()+min,wave.begin()+max,wsp.begin(),wave.begin()+min,
		   [n](float a, float b) -> float{return a+n*b;});
		   //addmultiple<n>());

    //for (std::size_t i = min; i < max; ++i) {
    //wave[i] += n * wsp[i-min];
    //}
      
  } // PMTsimulationAlg::AddPhotoelectrons()
  
  void icarus::opdet::PMTsimulationAlg::AddSPE(size_t time_bin, Waveform_t& wave){
     
    if (time_bin >= fNsamples) return;
    
    std::size_t const min = time_bin;
    std::size_t const max = std::min(time_bin + pulsesize, fNsamples);

    std::transform(wave.begin()+min,wave.begin()+max,wsp.begin(),wave.begin()+min,std::plus<float>());
    //for (std::size_t i = min; i < max; ++i) {
    //wave[i] += wsp[i-min];
    //}
    
  }
  
  void icarus::opdet::PMTsimulationAlg::AddNoise(Waveform_t& wave){
    
    CLHEP::RandGauss random(*fElecNoiseRandomEngine, 0.0, fAmpNoise);
    for(auto& sample: wave) {
      double const noise = random.fire(); //gaussian noise
      sample += noise;
    } // for sample
    
  } // PMTsimulationAlg::AddNoise()
  
  void icarus::opdet::PMTsimulationAlg::AddDarkNoise(Waveform_t& wave)
  {
    if (fDarkNoiseRate <= 0.0) return; // no dark noise
    size_t timeBin=0;
    CLHEP::RandExponential random(*fDarkNoiseRandomEngine, (1.0/fDarkNoiseRate)*1e9);
    // Multiply by 10^9 since fDarkNoiseRate is in Hz (conversion from s to ns)
    double darkNoiseTime = random.fire();
    while (darkNoiseTime < wave.size()){
      timeBin = (darkNoiseTime);
      AddSPE(timeBin,wave);
      // Find next time to add dark noise
      darkNoiseTime += random.fire();
    }
  }
  

