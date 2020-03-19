/**
 * @file   icaruscode/Light/SimPMTIcarus_module.cc
 * @see    `icarus::opdet::PMTsimulationAlg`
 * 
 * Based on `SimPMTSBND_module.cc` by L. Paulucci and F. Marinho.
 */


// ICARUS libraries
#include "icaruscode/PMT/SinglePhotonPulseFunctionTool.h"
#include "icaruscode/PMT/Algorithms/PMTsimulationAlg.h"
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanosecond
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/TableFragment.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C/C++ standard library
#include <vector>
#include <atomic> // std::atomic_flag
#include <iterator> // std::back_inserter()
#include <memory> // std::make_unique()
#include <utility> // std::move()
#include <optional>


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
   * All the configuration parameters are passed directly to the
   * `icarus::opdet::PMTsimulationAlg` algorithm, with the following exceptions:
   * * **InputModule** (input tag): tag of the data product containing the
   *   `sim::SimPhotons` collection to be digitized;
   * * **EfficiencySeed**, **DarkNoiseSeed** and **ElectronicsNoiseSeed**
   *   (integers, optional): if specified, each number is used to seed the
   *   pertaining random engine; otherwise, the seed is assigned by
   *   `NuRandomService`;
   * * **ElectronicsNoiseRandomEngine** and **DarkNoiseRandomEngine** (strings,
   *   default: `HepJamesRandom`): name of the random number generation
   *   algorithm to use; valid values are all the ones supported by
   *   `art::RandomNumberGenerator` (which match the random engine classes
   *   derived from `clhep::HepRandomEngine` in CLHEP 2.3 except
   *   `NonRandomEngine` and `RandEngine`);
   * * **WritePhotons** (boolean, default: `false`): saves in an additional
   *   `sim::SimPhotons` collection the photons effectively contributing to
   *   the waveforms; currently, no selection ever happens and all photons are
   *   contributing, making this collection the same as the input one.
   * 
   * See the @ref ICARUS_PMTSimulationAlg_RandomEngines "documentation" of
   * `icarus::PMTsimulationAlg` for the purpose of the three random number
   * engines.
   * 
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
   * If `WritePhotons` configuration parameter is set `true`, a collection of
   * the scintillation photons (`std::vector<sim::SimPhotons>`) which
   * contributed to the waveforms is also produced.
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
      
    struct Config
    {
      using Comment = fhicl::Comment;
      using Name = fhicl::Name;
      
      fhicl::Atom<art::InputTag> inputModuleLabel {
          Name("InputModule"),
          Comment("simulated photons to be digitised (sim::SimPhotons)")
      };
      
      fhicl::DelegatedParameter SinglePhotonResponse {
        fhicl::Name("SinglePhotonResponse"),
        fhicl::Comment(
          "parameters describing the single photon response"
          " (SinglePhotonPulseFunctionTool tool)"
          )
        };

      fhicl::TableFragment<icarus::opdet::PMTsimulationAlgMaker::Config> algoConfig;
      
      fhicl::Atom<bool> writePhotons {
          Name("WritePhotons"),
          Comment
            ("writes the scintillation photon contributing to the waveforms"),
          false
      };
      
      rndm::SeedAtom EfficiencySeed {
        Name("EfficiencySeed"),
        Comment("fix the seed for stocastic photon detection efficiency")
        };
      
      rndm::SeedAtom DarkNoiseSeed {
        Name("DarkNoiseSeed"),
        Comment("fix the seed for stocastic dark noise generation")
        };
      
      rndm::SeedAtom ElectronicsNoiseSeed {
        Name("ElectronicsNoiseSeed"),
        Comment("fix the seed for stocastic electronics noise generation")
        };
      
      fhicl::Atom<std::string> electronicsNoiseRandomEngine {
          Name("ElectronicsNoiseRandomEngine"),
          Comment("type of random engine to use for electronics noise"),
          "HepJamesRandom"
      };

      fhicl::Atom<std::string> darkNoiseRandomEngine {
          Name("DarkNoiseRandomEngine"),
          Comment("type of random engine to use for dark noise"),
          "HepJamesRandom"
      };

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
    
    /// Type of single photoelectron response function.
    using SinglePhotonResponseFunc_t
      = icarus::opdet::SinglePhotonResponseFunc_t const;
    
    /// Input tag for simulated scintillation photons (or photoelectrons).
    art::InputTag fInputModuleName;
    
    bool fWritePhotons { false }; ///< Whether to save contributing photons.
    
    /// Single photoelectron response function.
    std::unique_ptr<SinglePhotonResponseFunc_t> const fSinglePhotonResponseFunc;
    
    /// The actual simulation algorithm.
    icarus::opdet::PMTsimulationAlgMaker makePMTsimulator;

    CLHEP::HepRandomEngine&  fEfficiencyEngine;
    CLHEP::HepRandomEngine&  fDarkNoiseEngine;
    CLHEP::HepRandomEngine&  fElectronicsNoiseEngine;
    
    
    /// True if `firstTime()` has already been called.
    std::atomic_flag fNotFirstTime;
    
    /// Returns whether no other event has been processed yet.
    bool firstTime() { return !fNotFirstTime.test_and_set(); }
    
  }; // class SimPMTIcarus
  
  
  // ---------------------------------------------------------------------------
  // --- SimPMTIcarus implementation
  // ---------------------------------------------------------------------------
SimPMTIcarus::SimPMTIcarus(Parameters const& config)
    : EDProducer{config}
    , fInputModuleName(config().inputModuleLabel())
    , fWritePhotons(config().writePhotons())
    , fSinglePhotonResponseFunc{
        art::make_tool<icarus::opdet::SinglePhotonPulseFunctionTool>
          (config().SinglePhotonResponse.get<fhicl::ParameterSet>())
          ->getPulseFunction()
      }
    , makePMTsimulator(config().algoConfig())
    , fEfficiencyEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine
        (*this, "HepJamesRandom", "Efficiencies", config().EfficiencySeed)
      )
    , fDarkNoiseEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(
        *this,
        config().darkNoiseRandomEngine(),
        "DarkNoise",
        config().DarkNoiseSeed
      ))
    , fElectronicsNoiseEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(
        *this,
        config().electronicsNoiseRandomEngine(),
        "ElectronicsNoise",
        config().ElectronicsNoiseSeed
      ))
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    if (fWritePhotons) produces<std::vector<sim::SimPhotons> >();
  } // SimPMTIcarus::SimPMTIcarus()
  
  
  // ---------------------------------------------------------------------------
  void SimPMTIcarus::produce(art::Event & e)
  {
    mf::LogDebug("SimPMTIcarus") << e.id();
    
    //
    // fetch the input
    //
    auto pulseVecPtr = std::make_unique< std::vector< raw::OpDetWaveform > > ();
    
    std::unique_ptr<std::vector<sim::SimPhotons>> simphVecPtr;
    if (fWritePhotons)
      simphVecPtr = std::make_unique< std::vector< sim::SimPhotons > > ();
    
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
      *fSinglePhotonResponseFunc,
      fEfficiencyEngine,
      fDarkNoiseEngine,
      fElectronicsNoiseEngine,
      fWritePhotons
      );
    
    if (firstTime()) {
      mf::LogDebug log { "SimPMTIcarus" };
      log << "PMT simulation configuration (first event):\n";
      PMTsimulator->printConfiguration(log);
    } // if first time
    
    //
    // run the algorithm
    //
    for(auto const& photons : pmtVector) {
      
      auto const& [ channelWaveforms, photons_used ]
        = PMTsimulator->simulate(photons);
      std::move(
        channelWaveforms.cbegin(), channelWaveforms.cend(),
        std::back_inserter(*pulseVecPtr)
        );
      if (simphVecPtr && photons_used)
        simphVecPtr->emplace_back(std::move(photons_used.value()));
      
    } // for

    mf::LogInfo("SimPMTIcarus") << "Generated " << pulseVecPtr->size()
      << " waveforms out of " << pmtVector.size() << " optical channels.";
    
    //
    // save the result
    //
    e.put(std::move(pulseVecPtr));
    if (simphVecPtr) e.put(std::move(simphVecPtr));
  } // SimPMTIcarus::produce()
  
  
// ---------------------------------------------------------------------------
  DEFINE_ART_MODULE(SimPMTIcarus)
  
  
  // ---------------------------------------------------------------------------

} // namespace icarus
