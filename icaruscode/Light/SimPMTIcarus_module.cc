/**
 * @file   icaruscode/Light/SimPMTIcarus_module.cc
 * @see    `icarus::opdet::PMTsimulationAlg`
 * 
 * Based on `SimPMTSBND_module.cc` by L. Paulucci and F. Marinho.
 */


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
  
  /**
   * @brief Simulates the digitization of ICARUS PMT response and trigger.
   * 
   * The module is a simple interface to the simulation algorithm,
   * `icarus::opdet::PMTsimulationAlg`.
   * 
   * 
   * Configuration
   * ==============
   * 
   * Apart from the input collection of propagated photons, all the
   * configuration parameters are passed directly to the
   * `icarus::opdet::PMTsimulationAlg` algorithm.
   * 
   * The module also utilizes three random number engines.
   * Currently, no configuration interface is provided to directly control their
   * seeds, which is delegated to `rndm::NuRandomService` service.
   * 
   * 
   * Input
   * ======
   * 
   * The module utilizes as input a collection of `sim::SimPhotons`, each
   * containing the photons propagated to a single optical detector channel.
   * 
   * 
   * Output
   * =======
   * 
   * A collection of optical detector waveforms
   * (`std::vector<raw::OpDetWaveform>`) is produced.
   * See `icarus::opdet::PMTsimulationAlg` algorithm documentation for details.
   * 
   * 
   * Requirements
   * =============
   * 
   * This module currently requires LArSoft services:
   * * `DetectorClocksService` for timing conversions and settings
   * * `LArPropertiesService` for the scintillation yield(s)
   * 
   * Three random streams are also used.
   * 
   */
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
    
    explicit SimPMTIcarus(Parameters const& config);
    
    
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

    CLHEP::HepRandomEngine&  fEfficiencyEngine;
    CLHEP::HepRandomEngine&  fDarkNoiseEngine;
    CLHEP::HepRandomEngine&  fElectronicsNoiseEngine;
    
  }; // class SimPMTIcarus
  
  
  // ---------------------------------------------------------------------------
  // --- SimPMTIcarus implementation
  // ---------------------------------------------------------------------------
  SimPMTIcarus::SimPMTIcarus(Parameters const& config)
    : EDProducer{config}
    , fInputModuleName(config().inputModule())
    , makePMTsimulator(config().algoConfig())
    , fEfficiencyEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "Efficiencies", config, "SeedEfficinecy"))
    , fDarkNoiseEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "DarkNoise", config, "SeedDarkNoise"))
    , fElectronicsNoiseEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "ElectronicsNoise", config, "SeedElectronicsNoise"))
 {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    
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
      fEfficiencyEngine,
      fDarkNoiseEngine,
      fElectronicsNoiseEngine
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
  
  
  // ---------------------------------------------------------------------------
  DEFINE_ART_MODULE(SimPMTIcarus)
  
  
  // ---------------------------------------------------------------------------

} // namespace icarus
