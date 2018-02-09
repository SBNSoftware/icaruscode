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

class SimPMTIcarus;


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
  std::string fInputModuleName;
  double fSampling;       //wave sampling frequency (GHz)
  double fReadoutWindow;  //waveform time interval (ns)
  unsigned int fNsamples; //Samples per waveform
  double fPreTrigger;     //(ns)
  double fQE;             //PMT quantum efficiency

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
  
  void AddNoise(std::vector<double>& wave); //add noise to baseline
  void AddDarkNoise(std::vector<double>& wave); //add noise to baseline


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
  fBaseline        = p.get< double >("Baseline"     ); //in ADC
  fFallTime        = p.get< double >("FallTime"     ); //in ns
  fRiseTime        = p.get< double >("RiseTime"     ); //in ns
  fMeanAmplitude   = p.get< double >("MeanAmplitude"); //in pC
  fAmpNoise        = p.get< double >("AmpNoise"     ); //in ADC
  fDarkNoiseRate   = p.get< double >("DarkNoiseRate"); //in Hz
  fReadoutWindow   = p.get< double >("ReadoutWindow"); //in ns
  fPreTrigger      = p.get< double >("PreTrigger"   ); //in ns
  double temp_fQE  = p.get< double >("QE"           ); //PMT quantum efficiency

//Correction due to scalling factor applied during simulation if any
  auto const *LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
  fQE = temp_fQE/(LarProp->ScintPreScale());
  
  std::cout << "PMT corrected efficiency = " << fQE << std::endl;

  if(fQE>1.0001)
	std::cout << "WARNING: Quantum efficiency set in fhicl file " << temp_fQE << " seems to be too large! Final QE must be equal or smaller than the scintillation pre scale applied at simulation time. Please check this number (ScintPreScale): " << LarProp->ScintPreScale() << std::endl;

  if(p.get <double>("Sampling")==-1){
    auto const *timeService = lar::providerFrom< detinfo::DetectorClocksService >();
    fSampling = timeService->OpticalClock().Frequency();
  }else	
    fSampling = p.get <double>("Sampling");
  
  std::cout << "Sampling = " << fSampling << " GHz." << std::endl;
  
  fNsamples = (int)((fPreTrigger+fReadoutWindow)*fSampling);
  
//Random number engine initialization
  int seed = time(NULL);
  gRandom = new TRandom3(seed);

  //shape of single pulse
  sigma1 = fRiseTime/1.687;
  sigma2 = fFallTime/1.687;

  pulsesize=(int)((6*sigma2+fTransitTime)*fSampling);
  wsp.resize(pulsesize);

  for(int i=0; i<pulsesize; i++){
	wsp[i]=(Pulse1PE(static_cast< double >(i)/fSampling));
	//	std::cout << wsp[i] << std::endl;
  }


}

void SimPMTIcarus::produce(art::Event & e)
{
  std::unique_ptr< std::vector< raw::OpDetWaveform > > pulseVecPtr(std::make_unique< std::vector< raw::OpDetWaveform > > ());

  // Implementation of required member function here.
  std::cout <<"Event: " << e.id().event() << std::endl;

  art::Handle< std::vector<sim::SimPhotons> > pmtHandle;
  e.getByLabel(fInputModuleName, pmtHandle);

  if(!pmtHandle.isValid()){
    std::cout <<Form("Did not find any G4 photons from a producer: %s", "largeant") << std::endl;
  }
  
  std::cout << "Number of photon channels: " << pmtHandle->size() << std::endl;
  unsigned int nChannels = pmtHandle->size();

  std::vector<std::vector<short unsigned int>> waveforms(nChannels,std::vector<short unsigned int> (fNsamples,0));

  std::vector<std::vector<double>> waves(nChannels,std::vector<double>(fNsamples,fBaseline));

  int ch;
  double t_min = 1e15;

  for (auto const& simphotons : (*pmtHandle)){
	ch = simphotons.OpChannel();
  	std::cout <<"Channel: " << ch << std::endl;
        t_min = 1e15;

	for(size_t i=0; i<simphotons.size(); i++){
		if(simphotons[i].Time<t_min) t_min = simphotons[i].Time;
	}

	for(size_t i=0; i<simphotons.size(); i++){
	  if((gRandom->Uniform(1.0))<(fQE)) AddSPE((fPreTrigger+simphotons[i].Time-t_min)*fSampling,waves[ch]);
	}

	if(fAmpNoise>0.0) AddNoise(waves[ch]);
	if (fDarkNoiseRate > 0.0) AddDarkNoise(waves[ch]);

	//  for(size_t k=0; k<fNsamples; k++)std::cout << waves[ch][k] << std::endl;
	waveforms[ch] = std::vector<short unsigned int> (waves[ch].begin(), waves[ch].end());

	raw::OpDetWaveform adcVec(t_min, (unsigned int)ch, waveforms[ch]);//including pre trigger window and transit time
	pulseVecPtr->emplace_back(std::move(adcVec));
	
  }

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

