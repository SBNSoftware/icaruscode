////////////////////////////////////////////////////////////////////////
// Class:       SimPMTIcarus
// Plugin Type: producer (art v2_09_06)
// File:        SimPMTIcarus_module.cc
//
// Generated at Wed Feb  7 15:06:56 2018 by Andrea Falcone using cetskelgen
// from cetlib version v3_01_03.

//Based on SimPMTSBND_module.cc by L. Paulucci and F. Marinho
////////////////////////////////////////////////////////////////////////

// ICARUS libraries
#include "icaruscode/Light/Algorithms/PMTsimulationAlg.h"

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
//    using Waveform_t = std::vector<float>;
    
    // Declare member data here.
    art::InputTag fInputModuleName;
    
//    double fSampling;       //wave sampling frequency (GHz)
//    std::size_t fNsamples; //Samples per waveform
//    double fQE;             //PMT quantum efficiency
//    
//    size_t fReadoutWindowSize;     ///ReadoutWindowSize in samples
//    float  fPretrigFraction;       ///Fraction of window size to be before "trigger"
//    float  fThresholdADC;          ///ADC Threshold for self-triggered readout
//    int    fPulsePolarity;         ///Pulse polarity (=1 for positive, =-1 for negative)
//    double  fTriggerOffsetPMT;      ///Time (us) relative to trigger when pmt readout starts
//    double  fReadoutEnablePeriod;  ///Time (us) for which pmt readout is enabled
//
//    size_t fPretrigSize;
//    size_t fPosttrigSize;
//
//    bool fCreateBeamGateTriggers; ///Option to create unbiased readout around beam spill
//    double fBeamGateTriggerRepPeriod; ///Repetition Period (us) for BeamGateTriggers
//    size_t fBeamGateTriggerNReps; ///Number of beamgate trigger reps to produce
//
//    //Single PE parameters
//    double fFallTime;       //fall time of 1PE in ns
//    double fRiseTime;      //rise time in ns
//    double fTransitTime;   //to be added to pulse minimum time
//    double sigma1;
//    double sigma2;
//    double fMeanAmplitude;  //in pC
    
//     Waveform_t wsp; //single photon pulse vector
//     
//     int pulsesize; //size of 1PE waveform
//     
//     double fADC;      //charge to ADC convertion scale
//     double fBaseline; //waveform baseline
//     double fAmpNoise; //amplitude of gaussian noise
//     double fDarkNoiseRate; //in Hz
//     double fSaturation; //equivalent to the number of p.e. that saturates the electronic signal	
    
    
//    std::unordered_map< raw::Channel_t, Waveform_t > fFullWaveforms;
    
    /// The actual simulation algorithm.
    icarus::opdet::PMTsimulationAlgMaker makePMTsimulator;
    
  };
  
  
  SimPMTIcarus::SimPMTIcarus(fhicl::ParameterSet const & p)
    : makePMTsimulator(p)
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    
//    icarus::opdet::PMTsimulationAlg::Config config;
//    PMTsimulator = std::make_unique<icarus::opdet::PMTsimulationAlg>(config);
    
    // create three random engines for three independent tasks;
    // obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed";
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "Efficiencies", p, "Seed");
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "DarkNoise", p, "DarkNoiseSeed");
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "ElectronicsNoise", p, "ElectronicsNoiseSeed");
    
  }

  void SimPMTIcarus::produce(art::Event & e)
  {
    mf::LogDebug("SimPMTICARUS") << e.id() << std::endl;
    
    //
    // fetch the input
    //
    auto pulseVecPtr = std::make_unique< std::vector< raw::OpDetWaveform > > ();
    
    //
    // prepare the output
    //
    auto const& pmtVector = *(e.getValidHandle< std::vector<sim::SimPhotons> >(fInputModuleName));
    
    //
    // prepare the algorithm
    //
    auto PMTsimulator = makePMTsimulator(
      *(lar::providerFrom<detinfo::LArPropertiesService>()),
      *(lar::providerFrom<detinfo::DetectorClocksService>()),
      art::ServiceHandle<art::RandomNumberGenerator>()->getEngine("Efficiencies"),
      art::ServiceHandle<art::RandomNumberGenerator>()->getEngine("DarkNoise"),
      art::ServiceHandle<art::RandomNumberGenerator>()->getEngine("ElectronicsNoise")
      );
    
    //
    // run the algorithm
    //
    for(auto const& photons : pmtVector) {
      auto const& channelWaveforms = PMTsimulator->simulate(photons);
      std::move(channelWaveforms.cbegin(), channelWaveforms.cend(), std::back_inserter(*pulseVecPtr));
    } // for

    mf::LogInfo("SimPMTIcarus") << "Generated " << pulseVecPtr->size()
      << " waveforms out of " << pmtVector.size() << " optical channels.";
    
    //
    // save the result
    //
    e.put(std::move(pulseVecPtr));
    
  } // SimPMTIcarus::produce()
  
  DEFINE_ART_MODULE(SimPMTIcarus)

}

