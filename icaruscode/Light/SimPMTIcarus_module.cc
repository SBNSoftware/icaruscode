////////////////////////////////////////////////////////////////////////
// Class:       SimPMTIcarus
// Plugin Type: producer (art v2_09_06)
// File:        SimPMTIcarus_module.cc
//
// Generated at Wed Feb  7 15:06:56 2018 by Andrea Falcone using cetskelgen
// from cetlib version v3_01_03.

//Based on SimPMTSBND_module.cc by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <chrono>
#include <algorithm>
#include <functional>

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "nutools/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"


namespace {
  template <typename T>
  T sqr(T v) { return v*v; }
} // local namespace

namespace opdet{
  
  class SimPMTIcarus : public art::EDProducer {
  public:
    explicit SimPMTIcarus(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    SimPMTIcarus(SimPMTIcarus const &) = delete;
    SimPMTIcarus(SimPMTIcarus &&) = delete;
    SimPMTIcarus & operator = (SimPMTIcarus const &) = delete;
    SimPMTIcarus & operator = (SimPMTIcarus &&) = delete;
    
    // Required functions.
    void produce(art::Event & e) override;
    
  private:
    /// Type internally used for storing waveforms.
    using Waveform_t = std::vector<float>;
    
    // Declare member data here.
    art::InputTag fInputModuleName;
    
    double fSampling;       //wave sampling frequency (GHz)
    std::size_t fNsamples; //Samples per waveform
    double fQE;             //PMT quantum efficiency
    
    size_t fReadoutWindowSize;     ///ReadoutWindowSize in samples
    float  fPretrigFraction;       ///Fraction of window size to be before "trigger"
    float  fThresholdADC;          ///ADC Threshold for self-triggered readout
    int    fPulsePolarity;         ///Pulse polarity (=1 for positive, =-1 for negative)
    double  fTriggerOffsetPMT;      ///Time (us) relative to trigger when pmt readout starts
    double  fReadoutEnablePeriod;  ///Time (us) for which pmt readout is enabled

    size_t fPretrigSize;
    size_t fPosttrigSize;

    bool fCreateBeamGateTriggers; ///Option to create unbiased readout around beam spill
    double fBeamGateTriggerRepPeriod; ///Repetition Period (us) for BeamGateTriggers
    size_t fBeamGateTriggerNReps; ///Number of beamgate trigger reps to produce

    //Single PE parameters
    double fFallTime;       //fall time of 1PE in ns
    double fRiseTime;      //rise time in ns
    double fTransitTime;   //to be added to pulse minimum time
    double sigma1;
    double sigma2;
    double fMeanAmplitude;  //in pC
    
    void AddSPE(size_t time_bin, Waveform_t& wave); // add single pulse to auxiliary waveform
    /// Add `n` standard pulses starting at the specified `time_bin` of `wave`.
    void AddPhotoelectrons(size_t time_bin, unsigned int n, Waveform_t& wave) const;
    double Pulse1PE(double time) const;
    
    Waveform_t wsp; //single photon pulse vector
    
    int pulsesize; //size of 1PE waveform
    
    double fADC;      //charge to ADC convertion scale
    double fBaseline; //waveform baseline
    double fAmpNoise; //amplitude of gaussian noise
    double fDarkNoiseRate; //in Hz
    double fSaturation; //equivalent to the number of p.e. that saturates the electronic signal	
    
    
    void AddNoise(Waveform_t& wave); //add noise to baseline
    void AddDarkNoise(Waveform_t& wave); //add noise to baseline
    void AddPhoton(sim::OnePhoton const& ph, Waveform_t& wvfm);
    void CreateFullWaveforms(std::vector<sim::SimPhotons> const& pmtVector);
    void CreateFullWaveform(Waveform_t&,
			    std::vector<unsigned int>&,
			    sim::SimPhotons const&);
    std::set<size_t> CreateBeamGateTriggers() const;
    std::set<size_t> FindTriggers(Waveform_t const& wvfm) const;
    void CreateOpDetWaveforms(raw::Channel_t const& opch,
			      Waveform_t const& wvfm,
			      std::vector<raw::OpDetWaveform>& output_opdets);

    std::unordered_map< raw::Channel_t, Waveform_t > fFullWaveforms;
    detinfo::DetectorClocks const* fTimeService = nullptr; ///< DetectorClocks service provider.
    CLHEP::HepRandomEngine* fRandomEngine = nullptr; ///< Main random stream engine.
    
    /// Returns a random response whether a photon generates a photoelectron.
    bool KicksPhotoelectron() const;
    
  };
  
  
  SimPMTIcarus::SimPMTIcarus(fhicl::ParameterSet const & p)
  // :
  // Initialize member data here.
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    
    fInputModuleName = p.get< art::InputTag >("InputModule" );
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
    fTriggerOffsetPMT    = p.get<double>("TriggerOffsetPMT");   ///Time (us) relative to trigger that readout begins
    fReadoutEnablePeriod = p.get<double>("ReadoutEnablePeriod"); ///Time (us) for which pmt readout is enabled

    fCreateBeamGateTriggers = p.get<bool>("CreateBeamGateTriggers"); ///Option to create unbiased readout around beam spill
    fBeamGateTriggerRepPeriod = p.get<double>("BeamGateTriggerRepPeriod"); ///Repetition Period (us) for BeamGateTriggers
    fBeamGateTriggerNReps = p.get<size_t>("BeamGateTriggerNReps"); ///Number of beamgate trigger reps to produce
    
    fSaturation      = p.get< double >("Saturation"   ); //in number of p.e.
    double temp_fQE  = p.get< double >("QE"           ); //PMT quantum efficiency

    fPretrigSize = fPretrigFraction*fReadoutWindowSize;
    fPosttrigSize = fReadoutWindowSize-fPretrigSize;
    
    //Correction due to scalling factor applied during simulation if any
    auto const* LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
    fQE = temp_fQE/(LarProp->ScintPreScale());
    
    fTimeService = lar::providerFrom< detinfo::DetectorClocksService >();
    fSampling = fTimeService->OpticalClock().Frequency();
    fNsamples = fReadoutEnablePeriod * fSampling ; //us * MHz cancels out
    
    mf::LogDebug("SimPMTICARUS") << "PMT corrected efficiency = " << fQE << std::endl;
    
    if(fQE>1.0001)
      mf::LogDebug("SimPMTICARUS") << "WARNING: Quantum efficiency set in fhicl file "
				   << temp_fQE
				   << " seems to be too large! Final QE must be equal"
				   << " or smaller than the scintillation pre scale applied"
				   << " at simulation time. Please check this number (ScintPreScale): "
				   << LarProp->ScintPreScale() << std::endl;
    
    mf::LogDebug("SimPMTICARUS") << "Sampling = " << fSampling << " MHz." << std::endl;


    //shape of single pulse
    sigma1 = fRiseTime/1.687;
    sigma2 = fFallTime/1.687;
    
    pulsesize=(int)((6*sigma2+fTransitTime)*1.e3*fSampling);
    wsp.resize(pulsesize);
    
    for(int i=0; i<pulsesize; i++){
      wsp[i]=(Pulse1PE(static_cast< double >(i)*1.e3/fSampling));
    }
    
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed";
    // currently using only one stream, but it is recommended that one stream is used
    // per task, so that a change in one task does not cascade on all others.
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, p, "Seed");
    
  }

  void SimPMTIcarus::CreateFullWaveforms(std::vector<sim::SimPhotons> const& pmtVector){
    std::vector<unsigned int> PhotoelectronsPerSample(fNsamples, 0U);
    for(auto const& photons : pmtVector)
      CreateFullWaveform(fFullWaveforms[photons.OpChannel()],PhotoelectronsPerSample,photons);
  }

  void SimPMTIcarus::CreateFullWaveform(Waveform_t & waveform,
					std::vector<unsigned int> & PhotoelectronsPerSample,
					sim::SimPhotons const& photons){

    //auto& waveform = fFullWaveforms[opch];
    waveform.clear();
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
    for(auto const& pe : peMap){
      auto const nPE = pe.second;//PhotoelectronsPerSample[iSample];
      if (nPE == 0) continue;
      if (nPE == 1) AddSPE(pe.first,waveform);//AddSPE(iSample, waveform); // faster if n = 1
      else AddPhotoelectrons(pe.first, nPE, waveform);//AddPhotoelectrons(iSample, nPE, waveform);
    }

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

  }//end CreateFullWaveform

  std::set<size_t> SimPMTIcarus::CreateBeamGateTriggers() const
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
  
  std::set<size_t> SimPMTIcarus::FindTriggers(Waveform_t const& wvfm) const
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
  
  void SimPMTIcarus::CreateOpDetWaveforms(raw::Channel_t const& opch,
					  Waveform_t const& wvfm,
					  std::vector<raw::OpDetWaveform>& output_opdets)
  {
    //std::cout << "Finding triggers in " << opch << std::endl;

    std::set<size_t> trigger_locations = FindTriggers(wvfm);

    //std::cout << "Creating opdet waveforms in " << opch << std::endl;

    bool in_pulse=false;
    size_t trig_start=0,trig_stop=wvfm.size();

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
  
  bool SimPMTIcarus::KicksPhotoelectron() const
    { return CLHEP::RandFlat::shoot(fRandomEngine) < fQE; } 

  void SimPMTIcarus::AddPhoton(sim::OnePhoton const& ph, Waveform_t& wvfm){
    if(CLHEP::RandFlat::shoot(fRandomEngine) >= fQE) return; // hit&miss random selection 
    
    double const mytime = fTimeService->G4ToElecTime(ph.Time+fTransitTime)-fTimeService->TriggerTime()-fTriggerOffsetPMT;

    if(mytime<0 || mytime>fReadoutEnablePeriod)
      return;
    AddSPE(size_t(mytime*fSampling),wvfm);
  }
  
  void SimPMTIcarus::produce(art::Event & e)
  {
    mf::LogDebug("SimPMTICARUS") << e.id() << std::endl;

    // update the pointer to the random engine (probably not necessary...)
    fRandomEngine = &(art::ServiceHandle<art::RandomNumberGenerator>()->getEngine());
    fTimeService  = lar::providerFrom< detinfo::DetectorClocksService >();

    auto pulseVecPtr = std::make_unique< std::vector< raw::OpDetWaveform > > ();
    
    auto const& pmtVector = *(e.getValidHandle< std::vector<sim::SimPhotons> >(fInputModuleName));
    
    //CreateFullWaveforms(pmtVector);
    Waveform_t waveform;
    std::vector<unsigned int> PhotoelectronsPerSample;
    for(auto const& photons : pmtVector){
      CreateFullWaveform(waveform,PhotoelectronsPerSample,photons);
      CreateOpDetWaveforms(photons.OpChannel(),waveform,*pulseVecPtr);
    }
    e.put(std::move(pulseVecPtr));
    
  }
  
  DEFINE_ART_MODULE(SimPMTIcarus)

  double SimPMTIcarus::Pulse1PE(double time) const//single pulse waveform
  {
    double const sigma = (time < fTransitTime)? sigma1: sigma2;
    return (fADC*fMeanAmplitude*std::exp(-sqr(time - fTransitTime)/(2.0*sqr(sigma))));
  }

  //template <unsigned int N> float addmultiple (const float& x, const float& y) {return x+N*y;}
  
  void SimPMTIcarus::AddPhotoelectrons(size_t time_bin, unsigned int n, Waveform_t& wave) const {
    
    if (time_bin >= fNsamples) return;
    
    std::size_t const min = time_bin;
    std::size_t const max = std::min(time_bin + pulsesize, fNsamples);

    std::transform(wave.begin()+min,wave.begin()+max,wsp.begin(),wave.begin()+min,
		   [n](float a, float b) -> float{return a+n*b;});
		   //addmultiple<n>());

    //for (std::size_t i = min; i < max; ++i) {
    //wave[i] += n * wsp[i-min];
    //}
      
  } // SimPMTIcarus::AddPhotoelectrons()
  
  void SimPMTIcarus::AddSPE(size_t time_bin, Waveform_t& wave){
     
    if (time_bin >= fNsamples) return;
    
    std::size_t const min = time_bin;
    std::size_t const max = std::min(time_bin + pulsesize, fNsamples);

    std::transform(wave.begin()+min,wave.begin()+max,wsp.begin(),wave.begin()+min,std::plus<float>());
    //for (std::size_t i = min; i < max; ++i) {
    //wave[i] += wsp[i-min];
    //}
    
  }
  
  void SimPMTIcarus::AddNoise(Waveform_t& wave){
    
    CLHEP::RandGauss random(*fRandomEngine, 0.0, fAmpNoise);
    for(auto& sample: wave) {
      double const noise = random.fire(); //gaussian noise
      sample += noise;
    } // for sample
    
  } // SimPMTIcarus::AddNoise()
  
  void SimPMTIcarus::AddDarkNoise(Waveform_t& wave)
  {
    if (fDarkNoiseRate <= 0.0) return; // no dark noise
    size_t timeBin=0;
    CLHEP::RandExponential random(*fRandomEngine, (1.0/fDarkNoiseRate)*1e9);
    // Multiply by 10^9 since fDarkNoiseRate is in Hz (conversion from s to ns)
    double darkNoiseTime = random.fire();
    while (darkNoiseTime < wave.size()){
      timeBin = (darkNoiseTime);
      AddSPE(timeBin,wave);
      // Find next time to add dark noise
      darkNoiseTime += random.fire();
    }
  }
  
}

