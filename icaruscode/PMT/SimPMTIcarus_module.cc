/**
 * @file   icaruscode/PMT/SimPMTIcarus_module.cc
 * @see    `icarus::opdet::PMTsimulationAlg`
 * 
 * Based on `SimPMTSBND_module.cc` by L. Paulucci and F. Marinho.
 */


// ICARUS libraries
#include "icaruscode/PMT/PMTpedestalGeneratorTool.h"
#include "icaruscode/PMT/PMTnoiseGeneratorTool.h"
#include "icaruscode/PMT/SinglePhotonPulseFunctionTool.h"
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker
#include "icaruscode/PMT/Algorithms/PMTsimulationAlg.h"
#include "icaruscode/PMT/Algorithms/PedestalGeneratorAlg.h"
#include "icaruscode/PMT/Algorithms/NoiseGeneratorAlg.h"
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanosecond
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "nurandom/RandomUtils/NuRandomService.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
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
#include <tuple> // std::tie()
#include <atomic> // std::atomic_flag
#include <iterator> // std::back_inserter()
#include <memory> // std::make_unique()
#include <utility> // std::move()
#include <optional>


namespace icarus::opdet {
  
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
   * * **SinglePhotonResponse** (configuration table): configuration of the
   *   _art_ tool delivering the single photon response function (see below);
   * * **Pedestal** (configuration table): configuration of the _art_ tool
   *   generating the waveform pedestal and electronics noise;
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
   * * **MakeMetadata** (boolean, default: `true`): produces waveform metadata
   *   objects.
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
   * If `MakeMetadata` configuration parameter is set `true`, a collection of
   * `std::vector<sbn::OpDetWaveformMeta>` is produced, one per waveform, in the
   * same order as the waveform data product (satisfying the
   * @ref LArSoftProxyDefinitionParallelData "parallel data product"
   * requirement). A regular
   * `art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>` association is
   * also produced. The on-trigger/on-beam flags are assigned according to the
   * times reported by `DetectorClocks` service provider; the trigger in
   * particular is unlikely to be meaningful, given that to emulate it the
   * waveforms are typically needed.
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
   * 
   * Single photon response function tool
   * -------------------------------------
   * 
   * The single photon response function passed to the underlying algorithm,
   * `icarus::opdet::PMTsimulationAlg`, is extracted from an _art_ tool.
   * The specifications of the tool are:
   * 
   * * the _art_ tool class must implement the
   *   `icarus::opdet::SinglePhotonPulseFunctionTool` interface;
   * * the resulting function describes the signal as digitized by the PMT
   *   readout;
   * * the function must therefore express the level of the response as
   *   digitized ADC count above the baseline; this level is negative if the
   *   signal polarity is negative;
   * * the time scale of this function has `0` as the photoelectron production
   *   time.
   * 
   * These requirements imply that the tool (or the underlying function) must
   * include electronics effects like the the time it takes to form the signal,
   * the signal shape at PMT output, the delay and distortion of the cables and
   * the ones of the digitization electronics, but the baseline is assumed flat
   * at 0.
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
          Comment("simulated photons to be digitized (sim::SimPhotons)")
      };
      
      fhicl::DelegatedParameter SinglePhotonResponse {
        fhicl::Name("SinglePhotonResponse"),
        fhicl::Comment(
          "parameters describing the single photon response"
          " (SinglePhotonPulseFunctionTool tool)"
          )
        };

      fhicl::DelegatedParameter Pedestal {
        fhicl::Name("Pedestal"),
        fhicl::Comment(
          "parameters describing the pedestal generation algorithm"
          " (PMTnoiseGeneratorTool tool), including the electronics noise"
          )
        };

      fhicl::TableFragment<icarus::opdet::PMTsimulationAlgMaker::Config> algoConfig;
      
      fhicl::Atom<bool> MakeMetadata {
        Name("MakeMetadata"),
        Comment("writes a metadata object for each waveform"),
        true
      };
      
      fhicl::Atom<bool> writePhotons {
          Name("WritePhotons"),
          Comment
            ("writes the scintillation photon contributing to the waveforms"),
          false
      };
      
      rndm::SeedAtom EfficiencySeed {
        Name("EfficiencySeed"),
        Comment("fix the seed for stochastic photon detection efficiency")
        };
      
      rndm::SeedAtom DarkNoiseSeed {
        Name("DarkNoiseSeed"),
        Comment("fix the seed for stochastic dark noise generation")
        };
      
      rndm::SeedAtom ElectronicsNoiseSeed {
        Name("ElectronicsNoiseSeed"),
        Comment("fix the seed for stochastic electronics noise generation")
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
    
    /// Type of the pedestal generator algorithm.
    using PedestalGenerator_t
      = icarus::opdet::PMTpedestalGeneratorTool::Generator_t;
    
    /// Input tag for simulated scintillation photons (or photoelectrons).
    art::InputTag fInputModuleName;
    
    bool fMakeMetadata; ///< Whether to produce waveform metadata.
    bool fWritePhotons { false }; ///< Whether to save contributing photons.
    
    CLHEP::HepRandomEngine&  fEfficiencyEngine;
    CLHEP::HepRandomEngine&  fDarkNoiseEngine;
    CLHEP::HepRandomEngine&  fElectronicsNoiseEngine;
    
    /// Single photoelectron response function.
    std::unique_ptr<SinglePhotonResponseFunc_t> const fSinglePhotonResponseFunc;
    
    /// Pedestal generation algorithm (including electronics noise).
    std::unique_ptr<PedestalGenerator_t> const fPedestalGen;
    
    /// The actual simulation algorithm.
    icarus::opdet::PMTsimulationAlgMaker makePMTsimulator;

    
    /// True if `firstTime()` has already been called.
    std::atomic_flag fNotFirstTime;
    
    /// Returns the metadata of `waveforms` and their associations.
    std::pair<
      std::vector<sbn::OpDetWaveformMeta>,
      art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>
      >
      makeMetadata(
        art::Event const& event,
        std::vector<raw::OpDetWaveform> const& waveforms,
        detinfo::DetectorTimings const& detTimings
      ) const;
    
    /// Returns whether no other event has been processed yet.
    bool firstTime() { return !fNotFirstTime.test_and_set(); }
    
  }; // class SimPMTIcarus
  
  
  // ---------------------------------------------------------------------------
  // --- SimPMTIcarus implementation
  // ---------------------------------------------------------------------------
SimPMTIcarus::SimPMTIcarus(Parameters const& config)
    : EDProducer{config}
    // configuration
    , fInputModuleName(config().inputModuleLabel())
    , fMakeMetadata(config().MakeMetadata())
    , fWritePhotons(config().writePhotons())
    // random engines
    , fEfficiencyEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
          createEngine(0, "HepJamesRandom", "Efficiencies"),
          "HepJamesRandom",
          "Efficiencies",
          config().EfficiencySeed
      ))
    , fDarkNoiseEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
        createEngine(0, config().darkNoiseRandomEngine(), "DarkNoise"),
        config().darkNoiseRandomEngine(),
        "DarkNoise",
        config().DarkNoiseSeed
      ))
    , fElectronicsNoiseEngine(art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
        createEngine(0, config().electronicsNoiseRandomEngine(), "ElectronicsNoise"),
        config().electronicsNoiseRandomEngine(),
        "ElectronicsNoise",
        config().ElectronicsNoiseSeed
      ))
    // algorithms
    , fSinglePhotonResponseFunc{
        art::make_tool<icarus::opdet::SinglePhotonPulseFunctionTool>
          (config().SinglePhotonResponse.get<fhicl::ParameterSet>())
          ->getPulseFunction()
      }
    , fPedestalGen{
        art::make_tool<icarus::opdet::PMTpedestalGeneratorTool>
          (config().Pedestal.get<fhicl::ParameterSet>())
          ->makeGenerator(fElectronicsNoiseEngine)
      }
    , makePMTsimulator(config().algoConfig())
  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    if (fMakeMetadata) {
      produces<std::vector<sbn::OpDetWaveformMeta>>();
      produces<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>();
    }
    if (fWritePhotons) produces<std::vector<sim::SimPhotons> >();
    
    fNotFirstTime.clear(); // superfluous in C++20
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
    art::Handle<std::vector<sim::SimPhotons> > pmtVector;
    art::Handle<std::vector<sim::SimPhotonsLite> > pmtLiteVector;
    pmtVector = e.getHandle< std::vector<sim::SimPhotons> >(fInputModuleName);
    if(!pmtVector.isValid())
      pmtLiteVector = e.getHandle< std::vector<sim::SimPhotonsLite> >(fInputModuleName);
    
    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    //
    // prepare the algorithm
    //
    auto PMTsimulator = makePMTsimulator(
      e.time().value(), // using the event generation time as beam time stamp
      *(lar::providerFrom<detinfo::LArPropertiesService>()),
      clockData,
      *fSinglePhotonResponseFunc,
      *fPedestalGen,
      fEfficiencyEngine,
      fDarkNoiseEngine,
      fElectronicsNoiseEngine,
      fWritePhotons
      );
    
    if (firstTime()) {
      mf::LogInfo log { "SimPMTIcarus" };
      log << "PMT simulation configuration (first event):\n";
      PMTsimulator->printConfiguration(log);
    } // if first time
    
    //
    // run the algorithm
    //
    unsigned int nopch = 0;
    if(pmtVector.isValid()) {
      nopch = pmtVector->size();
      for(auto const& photons : *pmtVector) {
      
        // Make an empty SimPhotonsLite with the same channel number.

        sim::SimPhotonsLite lite_photons(photons.OpChannel());

        auto const& [ channelWaveforms, photons_used ]
          = PMTsimulator->simulate(photons, lite_photons);
        std::move(
          channelWaveforms.cbegin(), channelWaveforms.cend(),
          std::back_inserter(*pulseVecPtr)
          );
        if (simphVecPtr && photons_used)
          simphVecPtr->emplace_back(std::move(photons_used.value()));

      } // for
    }
    else if(pmtLiteVector.isValid()) {
      nopch = pmtLiteVector->size();
      for(auto const& lite_photons : *pmtLiteVector) {

        // Make an empty SimPhotons with the same channel number.

        sim::SimPhotons photons(lite_photons.OpChannel);
      
        auto const& [ channelWaveforms, photons_used ]
          = PMTsimulator->simulate(photons, lite_photons);
        std::move(
          channelWaveforms.cbegin(), channelWaveforms.cend(),
          std::back_inserter(*pulseVecPtr)
          );

      } // for
    }

    mf::LogInfo("SimPMTIcarus") << "Generated " << pulseVecPtr->size()
      << " waveforms out of " << nopch << " optical channels.";
    
    // waveform metadata
    std::unique_ptr<std::vector<sbn::OpDetWaveformMeta>> metadataVec;
    std::unique_ptr<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
      metadataAssns;
    if (fMakeMetadata) {
      metadataVec = std::make_unique<std::vector<sbn::OpDetWaveformMeta>>();
      metadataAssns =
        std::make_unique<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>
        ();
      
      auto const detTimings = detinfo::makeDetectorTimings
        (art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e))
        ;
      std::tie(*metadataVec, *metadataAssns)
        = makeMetadata(e, *pulseVecPtr, detTimings);
    }
    
    //
    // save the result
    //
    e.put(std::move(pulseVecPtr));
    if (fMakeMetadata) {
      e.put(std::move(metadataVec));
      e.put(std::move(metadataAssns));
    }
    if (simphVecPtr) e.put(std::move(simphVecPtr));
  } // SimPMTIcarus::produce()
  
  
  // ---------------------------------------------------------------------------
  std::pair<
    std::vector<sbn::OpDetWaveformMeta>,
    art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>
    >
  SimPMTIcarus::makeMetadata(
    art::Event const& event,
    std::vector<raw::OpDetWaveform> const& waveforms,
    detinfo::DetectorTimings const& detTimings
  ) const {
    
    std::vector<sbn::OpDetWaveformMeta> meta;
    art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta> assns;
    
    art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr{ event };
    art::PtrMaker<sbn::OpDetWaveformMeta> const makeMetaPtr{ event };
    
    sbn::OpDetWaveformMetaMaker const makeOpDetWaveformMeta{ detTimings };
    
    for (auto const& [iWaveform, waveform ]: util::enumerate(waveforms)) {
      meta.push_back(makeOpDetWaveformMeta(waveform));
      assns.addSingle(makeWaveformPtr(iWaveform), makeMetaPtr(iWaveform));
    } // for waveforms
    
    return { std::move(meta), std::move(assns) };
    
  } // SimPMTIcarus::makeMetadata()
  
  
  // ---------------------------------------------------------------------------
  DEFINE_ART_MODULE(SimPMTIcarus)
  
  
  // ---------------------------------------------------------------------------

} // namespace icarus::opdet
