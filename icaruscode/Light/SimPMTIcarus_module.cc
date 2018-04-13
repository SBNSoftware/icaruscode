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
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <sstream>
#include <fstream>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TF1.h"

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
    
    // Declare member data here.
    art::InputTag fInputModuleName;
    
    double fSampling;       //wave sampling frequency (GHz)
    double fReadoutWindow;  //waveform time interval (ns)
    unsigned int fNsamples; //Samples per waveform
    double fPreTrigger;     //(ns)
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
    double fTransitSpread; //transit time spread
    double sigma1;
    double sigma2;
    double fMeanAmplitude;  //in pC
    
    void AddSPE(size_t time_bin, std::vector<double>& wave); // add single pulse to auxiliary waveform
    double Pulse1PE(double time) const;
    
    std::vector<double> wsp; //single photon pulse vector
    
    int pulsesize; //size of 1PE waveform
    
    double fADC;      //charge to ADC convertion scale
    double fBaseline; //waveform baseline
    double fAmpNoise; //amplitude of gaussian noise
    double fDarkNoiseRate; //in Hz
    double fSaturation; //equivalent to the number of p.e. that saturates the electronic signal	
    
    
    void AddNoise(std::vector<double>& wave); //add noise to baseline
    void AddDarkNoise(std::vector<double>& wave); //add noise to baseline
    void AddPhoton(sim::OnePhoton const& ph, std::vector<double>& wvfm);
    void CreateFullWaveforms(std::vector<sim::SimPhotons> const& pmtVector);
    void CreateBeamGateTriggers(std::set<size_t>& trigger_locations);
    void FindTriggers(std::vector<double> const& wvfm,
		      std::set<size_t>& trigger_locations);
    void CreateOpDetWaveforms(raw::Channel_t const& opch,
			      std::vector<double> const& wvfm,
			      std::vector<raw::OpDetWaveform>& output_opdets);

    std::unordered_map< raw::Channel_t,std::vector<double> > fFullWaveforms;
    
  };
  
  
  SimPMTIcarus::SimPMTIcarus(fhicl::ParameterSet const & p)
  // :
  // Initialize member data here.
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    
    fInputModuleName = p.get< std::string >("InputModule" );
    fTransitTime     = p.get< double >("TransitTime"  ); //ns
    fADC             = p.get< double >("ADC"          ); //voltage to ADC factor
    fBaseline        = p.get<uint16_t>("Baseline"     ); //in ADC
    fFallTime        = p.get< double >("FallTime"     ); //in ns
    fRiseTime        = p.get< double >("RiseTime"     ); //in ns
    fMeanAmplitude   = p.get< double >("MeanAmplitude"); //in pC
    fAmpNoise        = p.get< double >("AmpNoise"     ); //in ADC
    fDarkNoiseRate   = p.get< double >("DarkNoiseRate"); //in Hz
    
    fReadoutWindowSize   = p.get<size_t>("ReadoutWindowSize"); ///ReadoutWindowSize
    fPretrigFraction     = p.get<float>("PretrigFraction");    ///Fraction of window size to be before "trigger"
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
    
    auto const* TimeService  = lar::providerFrom< detinfo::DetectorClocksService >();
    fSampling = TimeService->OpticalClock().Frequency();
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


  //Random number engine initialization
    int seed = time(NULL);
    gRandom = new TRandom3(seed);
    
    //shape of single pulse
    sigma1 = fRiseTime/1.687;
    sigma2 = fFallTime/1.687;
    
    pulsesize=(int)((6*sigma2+fTransitTime)*1.e3*fSampling);
    wsp.resize(pulsesize);
    
    for(int i=0; i<pulsesize; i++){
      wsp[i]=(Pulse1PE(static_cast< double >(i)*1.e3/fSampling));
    }
  }

  void SimPMTIcarus::CreateFullWaveforms(std::vector<sim::SimPhotons> const& pmtVector){

    fFullWaveforms.clear();

    for(auto const& photons : pmtVector){

      fFullWaveforms[photons.OpChannel()].resize(fNsamples,fBaseline);
      for(auto const& ph : photons)
	AddPhoton(ph,fFullWaveforms[photons.OpChannel()]);
      if(fAmpNoise>0.) AddNoise(fFullWaveforms[photons.OpChannel()]);
      AddDarkNoise(fFullWaveforms[photons.OpChannel()]);

      for(size_t k=0; k<fNsamples; k++){ //Implementing saturation effects
	if(fFullWaveforms[photons.OpChannel()][k]<(fBaseline+fSaturation*fADC*fMeanAmplitude))	
	  fFullWaveforms[photons.OpChannel()][k]=fBaseline+fSaturation*fADC*fMeanAmplitude;
      }
    }

  }//end CreateFullWaveforms

  void SimPMTIcarus::CreateBeamGateTriggers(std::set<size_t>& trigger_locations)
  {
    double trig_time;
    auto const* TimeService  = lar::providerFrom< detinfo::DetectorClocksService >();
    
    for(size_t i_trig=0; i_trig<fBeamGateTriggerNReps; ++i_trig){
      trig_time = (TimeService->BeamGateTime()-TimeService->TriggerTime())+fBeamGateTriggerRepPeriod*i_trig-fTriggerOffsetPMT;
      if(trig_time<0 || trig_time>fReadoutEnablePeriod) continue;
      trigger_locations.insert(size_t(trig_time*fSampling));
    }
    
  }
  
  void SimPMTIcarus::FindTriggers(std::vector<double> const& wvfm,
				  std::set<size_t>& trigger_locations)
  {
    trigger_locations.clear();
    if(fCreateBeamGateTriggers) CreateBeamGateTriggers(trigger_locations);
    
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
    
  }
  
  void SimPMTIcarus::CreateOpDetWaveforms(raw::Channel_t const& opch,
					  std::vector<double> const& wvfm,
					  std::vector<raw::OpDetWaveform>& output_opdets)
  {
    std::set<size_t> trigger_locations;
    FindTriggers(wvfm,trigger_locations);

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
  
  void SimPMTIcarus::AddPhoton(sim::OnePhoton const& ph, std::vector<double>& wvfm){
    auto const* TimeService  = lar::providerFrom< detinfo::DetectorClocksService >();
    double mytime = TimeService->G4ToElecTime(ph.Time+fTransitTime)-TimeService->TriggerTime()-fTriggerOffsetPMT;

    if(mytime<0 || mytime>fReadoutEnablePeriod)
      return;
    if((gRandom->Uniform(1.0))<(fQE)) AddSPE(size_t(mytime*fSampling),wvfm);
  }
  
  void SimPMTIcarus::produce(art::Event & e)
  {
    mf::LogDebug("SimPMTICARUS") <<"Event: " << e.id().event() << std::endl;

    std::unique_ptr< std::vector< raw::OpDetWaveform > > pulseVecPtr(std::make_unique< std::vector< raw::OpDetWaveform > > ());
    
    art::Handle< std::vector<sim::SimPhotons> > pmtHandle;
    e.getByLabel(fInputModuleName, pmtHandle);
    
    if(!pmtHandle.isValid()){
      mf::LogError("SimPMTICARUS") << "Did not find any G4 photons from producer." << std::endl;
    }
    
    auto const& pmtVector(*pmtHandle);
    
    CreateFullWaveforms(pmtVector);
    for(auto const& full_wvfm : fFullWaveforms)
      CreateOpDetWaveforms(full_wvfm.first,full_wvfm.second,*pulseVecPtr);

    e.put(std::move(pulseVecPtr));
    
  }
  
  DEFINE_ART_MODULE(SimPMTIcarus)

  double SimPMTIcarus::Pulse1PE(double time) const//single pulse waveform
  {
    if (time < fTransitTime) return (fADC*fMeanAmplitude*std::exp(-1.0*pow(time - fTransitTime,2.0)/(2.0*pow(sigma1,2.0))));
    else return (fADC*fMeanAmplitude*std::exp(-1.0*pow(time - fTransitTime,2.0)/(2.0*pow(sigma2,2.0))));
    
  }
  
  void SimPMTIcarus::AddSPE(size_t time_bin, std::vector<double>& wave){
    
    size_t min=0;
    size_t max=0;
    
    if(time_bin<fNsamples){
      min=time_bin;
      max=time_bin+pulsesize < fNsamples ? time_bin+pulsesize : fNsamples;

      for(size_t i = min; i<= max; i++){
	wave[i]+= wsp[i-min];
      }
      
    }
  }
  
  void SimPMTIcarus::AddNoise(std::vector<double>& wave){
    
    double noise = 0.0;
    
    for(size_t i = 0; i<wave.size(); i++){
      noise = gRandom->Gaus(0,fAmpNoise); //gaussian noise
      wave[i] += noise;
    }
    
  }
  
  void SimPMTIcarus::AddDarkNoise(std::vector< double >& wave)
  {
    size_t timeBin=0;

    // Multiply by 10^9 since fDarkNoiseRate is in Hz (conversion from s to ns)
    double darkNoiseTime = static_cast< double >(gRandom->Exp((1.0/fDarkNoiseRate)*1000000000.0));
    while (darkNoiseTime < wave.size()){
      timeBin = (darkNoiseTime);
      AddSPE(timeBin,wave);
      // Find next time to add dark noise
      darkNoiseTime += static_cast< double >(gRandom->Exp((1.0/fDarkNoiseRate)*1000000000.0));
    }
  }
  
}

