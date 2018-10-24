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

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "nutools/RandomUtils/NuRandomService.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/TableFragment.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C/C++ standard library
#include <vector>
#include <iterator> // std::back_inserter()
#include <memory> // std::make_unique()
#include <utility> // std::move()


namespace opdet{
  
  class SimPMTIcarus : public art::EDProducer {
  public:
    
    struct Config {
      
      using Comment = fhicl::Comment;
      using Name = fhicl::Name;
      
      fhicl::Atom<art::InputTag> inputModule {
        Name("InputModule"),
        Comment("simulated photons to be digitised (sim::SimPhotons)")
        };
      
      fhicl::TableFragment<icarus::opdet::PMTsimulationAlgMaker::Config>
        algoConfig;
      
    }; // struct Config
    
    using Parameters = art::EDProducer::Table<Config>;
    
    explicit SimPMTIcarus(fhicl::ParameterSet const & p);
  //  explicit SimPMTIcarus(Parameters const& config);
    
    
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
    
    /// The actual simulation algorithm.
    icarus::opdet::PMTsimulationAlgMaker makePMTsimulator;
    
  }; // class SimPMTIcarus
  
  
  // ---------------------------------------------------------------------------
  // --- SimPMTIcarus implementation
  // ---------------------------------------------------------------------------
#if 0
  SimPMTIcarus::SimPMTIcarus(Parameters const& config)
    : fInputModuleName(config().inputModule())
    , makePMTsimulator(config().algoConfig())
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    
    // create three random engines for three independent tasks;
    // obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed";
    art::ServiceHandle<rndm::NuRandomService>()->createEngine
      (*this, "HepJamesRandom", "Efficiencies", p, "Seed");
    art::ServiceHandle<rndm::NuRandomService>()->createEngine
      (*this, "HepJamesRandom", "DarkNoise", p, "DarkNoiseSeed");
    art::ServiceHandle<rndm::NuRandomService>()->createEngine
      (*this, "HepJamesRandom", "ElectronicsNoise", p, "ElectronicsNoiseSeed");
    
  } // SimPMTIcarus::SimPMTIcarus()
  
#endif // 0
  // ---------------------------------------------------------------------------
  SimPMTIcarus::SimPMTIcarus(fhicl::ParameterSet const & p)
    : fInputModuleName(p.get< art::InputTag >("InputModule"))
    , makePMTsimulator(p)
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    
    // create three random engines for three independent tasks;
    // obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed";
    art::ServiceHandle<rndm::NuRandomService>()->createEngine
      (*this, "HepJamesRandom", "Efficiencies", p, "Seed");
    art::ServiceHandle<rndm::NuRandomService>()->createEngine
      (*this, "HepJamesRandom", "DarkNoise", p, "DarkNoiseSeed");
    art::ServiceHandle<rndm::NuRandomService>()->createEngine
      (*this, "HepJamesRandom", "ElectronicsNoise", p, "ElectronicsNoiseSeed");
    
  } // SimPMTIcarus::SimPMTIcarus()
  
  
  // ---------------------------------------------------------------------------
  void SimPMTIcarus::produce(art::Event & e)
  {
    mf::LogDebug("SimPMTIcarus") << e.id();
    
    //
    // fetch the input
    //
    auto pulseVecPtr = std::make_unique< std::vector< raw::OpDetWaveform > > ();
    
    //
    // prepare the output
    //
    auto const& pmtVector
      = *(e.getValidHandle< std::vector<sim::SimPhotons> >(fInputModuleName));
    
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
      std::move(
        channelWaveforms.cbegin(), channelWaveforms.cend(),
        std::back_inserter(*pulseVecPtr)
        );
    } // for

    mf::LogInfo("SimPMTIcarus") << "Generated " << pulseVecPtr->size()
      << " waveforms out of " << pmtVector.size() << " optical channels.";
    
    //
    // save the result
    //
    e.put(std::move(pulseVecPtr));
    
  } // SimPMTIcarus::produce()
  
  DEFINE_ART_MODULE(SimPMTIcarus)

} // namespace icarus

