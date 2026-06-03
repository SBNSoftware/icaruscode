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
#include "icaruscode/PMT/Calibration/ICARUSPhotonCalibratorServiceFromDB.h"
#include "icaruscode/PMT/Status/IPMTChannelStatusService.h"
#include "icaruscode/PMT/Algorithms/OpDetWaveformMetaUtils.h" // OpDetWaveformMetaMaker
#include "icaruscode/PMT/Algorithms/PMTsimulationAlg.h"
#include "icaruscode/PMT/Algorithms/PedestalGeneratorAlg.h"
#include "icaruscode/PMT/Algorithms/NoiseGeneratorAlg.h"
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"
#include "icaruscode/IcarusObj/OpDetWaveformMeta.h"
#include "icaruscode/Timing/IPMTTimingCorrectionService.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

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
#include "art_root_io/TFileService.h"
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

// ROOT libraries
#include "TTree.h"

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
   * If the `sim::SimPhotons` collection is not available, a collection of
   * `sim::SimPhotonsLite` is used instead if found.
   * 
   * 
   * Output
   * =======
   * 
   * A collection of optical detector waveforms
   * (`std::vector<raw::OpDetWaveform>`) is produced.
   * See `icarus::opdet::PMTsimulationAlg` algorithm documentation for details.
   * 
   * A collection of baselines (`std::vector<icarus::WaveformBaseline>`) and
   * associations to their waveform is also produced, using the configured
   * baseline as value for all the waveforms.
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
      
      fhicl::Atom<bool> UseChannelStatusDatabase {
        Name("UseChannelStatusDatabase"),
        Comment("skip simulation of channels flagged in the PMT channel status DB"),
        false
      };

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

      fhicl::Atom<bool> DebugTree {
        Name("DebugTree"),
        Comment("enable debug tree output"),
        false
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
    virtual void beginJob() override;
    void produce(art::Event & e) override;
    
  private:
    
    /// Type of single photoelectron response function.
    using SinglePhotonResponseFunc_t
      = icarus::opdet::SinglePhotonResponseFunc_t const;
    
    /// Type of the pedestal generator algorithm.
    using PedestalGenerator_t
      = icarus::opdet::PMTpedestalGeneratorTool::Generator_t;
    
    /// Debug tree data (one entry per optical channel per event)
    struct DebugInfo_t {
      // debug tree branch variables (one entry per optical channel per event)
      // event identification
      int run = -1;    ///< Run number.
      int event = -1;  ///< Event number.
      int channel = -1; ///< Optical channel number.
      // per-channel timing parameters
      float timeDelay_us = 0;          ///< Per-channel timing delay applied [us]: laser+cosmics DB corrections.
      float triggerOffsetPMT_us = 0;   ///< Global PMT readout start offset relative to trigger [us]
      // waveform geometry
      int nSamples = 0;    ///< number of ticks.
      int nSubsamples = 0; ///< number of subsamples (sub-tick interpolation steps)

      // full-window waveform (before zero-suppression splitting)
      std::vector<float> waveform_adc; ///< full waveform ADC values [ADC counts].

      // per-photon info
      int nDetectedPhotons = 0; ///< Number of photons accepted after QE cut.
      int nInputPhotons = 0;    ///< Total number of photons in sim::SimPhotons.
      std::vector<float> photon_simTime_ns;    ///< G4 simulation time of each detected photon [ns].
      std::vector<float> photon_waveformTime_us; ///< Photon time in waveform coordinates [us]: toTriggerTime - triggerOffsetPMT + timeDelay.
      std::vector<float> photon_start_x;      ///< X position of the photon emission point [cm].
      std::vector<float> photon_start_y;      ///< Y position of the photon emission point [cm].
      std::vector<float> photon_start_z;      ///< Z position of the photon emission point [cm].
      std::vector<int>   photon_tick;         ///< Waveform tick index where the photon was placed.
      std::vector<uint16_t> photon_subtick;   ///< Sub-tick index (subsample bin) where the photon was placed.

      // per-PE-deposit info
      int nPEDeposits = 0;                      ///< Number of distinct PE deposit entries.
      std::vector<int>      pedeposit_tick;     ///< Waveform tick of the PE deposit.
      std::vector<uint16_t> pedeposit_subtick;  ///< Sub-tick index of the PE deposit.
      std::vector<uint16_t> pedeposit_nPE;      ///< Integer number of photoelectrons at this deposit (after QE, before gain fluctuation).
      std::vector<float>    pedeposit_nEffectivePE; ///< Effective PE count after gain fluctuation (what was actually added to the waveform).
      std::vector<float>    pedeposit_gainFactor;   ///< Ratio nEffectivePE/nPE: encodes the per-deposit gain fluctuation (DB-calibrated or Poisson).
    }; // DebugInfo_t


    /// Input tag for simulated scintillation photons (or photoelectrons).
    art::InputTag fInputModuleName;
    
    bool fMakeMetadata; ///< Whether to produce waveform metadata.
    bool fWritePhotons { false }; ///< Whether to save contributing photons.
    bool fDoTimingDelays; ///< Whether timing delay corrections are applied.
    bool fUseGainCalibDB; ///< Whether per-channel gain calibration from DB is applied.
    bool fUseChannelStatusDB; ///< Whether to skip non-ON channels using the status DB.

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

    bool fDoDebugTree { false }; ///< Whether to enable debug tree output.
    TTree* fDebugTree { nullptr }; ///< Debug tree pointer (owned by TFileService).
    DebugInfo_t fDebugInfo; ///< Collected debug information.

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
    , fDoTimingDelays(config().algoConfig().ApplyTimingDelays())
    , fUseGainCalibDB(config().algoConfig().UseGainDatabase())
    , fUseChannelStatusDB(config().UseChannelStatusDatabase())
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
    , fDoDebugTree(config().DebugTree())

  {
    // Call appropriate produces<>() functions here.
    produces<std::vector<raw::OpDetWaveform>>();
    produces<std::vector<icarus::WaveformBaseline>>();
    produces<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>();
    if (fMakeMetadata) {
      produces<std::vector<sbn::OpDetWaveformMeta>>();
      produces<art::Assns<raw::OpDetWaveform, sbn::OpDetWaveformMeta>>();
    }
    if (fWritePhotons) produces<std::vector<sim::SimPhotons> >();
    
    fNotFirstTime.clear(); // superfluous in C++20
  } // SimPMTIcarus::SimPMTIcarus()
  
  // ---------------------------------------------------------------------------
  void SimPMTIcarus::beginJob()
  {
    if (fDoDebugTree) {
      art::ServiceHandle<art::TFileService> tfs;
      fDebugTree = tfs->make<TTree>("SimPMTDebug","Debug tree for SimPMTIcarus");
      fDebugTree->Branch("run", &fDebugInfo.run, "run/I");
      fDebugTree->Branch("event", &fDebugInfo.event, "event/I");
      fDebugTree->Branch("channel", &fDebugInfo.channel, "channel/I");
      fDebugTree->Branch("nSamples", &fDebugInfo.nSamples, "nSamples/I");
      fDebugTree->Branch("nSubsamples", &fDebugInfo.nSubsamples, "nSubsamples/I");
      fDebugTree->Branch("timeDelay_us", &fDebugInfo.timeDelay_us, "timeDelay_us/F");
      fDebugTree->Branch("triggerOffsetPMT_us", &fDebugInfo.triggerOffsetPMT_us, "triggerOffsetPMT_us/F");
      fDebugTree->Branch("waveform_adc", &fDebugInfo.waveform_adc);
      fDebugTree->Branch("nInputPhotons", &fDebugInfo.nInputPhotons, "nInputPhotons/I");
      fDebugTree->Branch("nDetectedPhotons", &fDebugInfo.nDetectedPhotons, "nDetectedPhotons/I");
      fDebugTree->Branch("photon_start_x", &fDebugInfo.photon_start_x);
      fDebugTree->Branch("photon_start_y", &fDebugInfo.photon_start_y);
      fDebugTree->Branch("photon_start_z", &fDebugInfo.photon_start_z);
      fDebugTree->Branch("photon_simTime_ns", &fDebugInfo.photon_simTime_ns);
      fDebugTree->Branch("photon_waveformTime_us", &fDebugInfo.photon_waveformTime_us);
      fDebugTree->Branch("photon_tick", &fDebugInfo.photon_tick);
      fDebugTree->Branch("photon_subtick", &fDebugInfo.photon_subtick);
      fDebugTree->Branch("nPEDeposits", &fDebugInfo.nPEDeposits, "nPEDeposits/I");
      fDebugTree->Branch("pedeposit_tick", &fDebugInfo.pedeposit_tick);
      fDebugTree->Branch("pedeposit_subtick", &fDebugInfo.pedeposit_subtick);
      fDebugTree->Branch("pedeposit_nPE", &fDebugInfo.pedeposit_nPE);
      fDebugTree->Branch("pedeposit_nEffectivePE", &fDebugInfo.pedeposit_nEffectivePE);
      fDebugTree->Branch("pedeposit_gainFactor", &fDebugInfo.pedeposit_gainFactor);
    }
  } // SimPMTIcarus::beginJob()

  
  // ---------------------------------------------------------------------------
  void SimPMTIcarus::produce(art::Event & e)
  {
    mf::LogDebug("SimPMTIcarus") << e.id();
    
    //
    // fetch the input
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
      fDoTimingDelays ? lar::providerFrom<icarusDB::IPMTTimingCorrectionService>() : nullptr,
      fUseGainCalibDB ? art::ServiceHandle<calib::ICARUSPhotonCalibratorServiceFromDB>()->provider() : nullptr,
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
      log << "useChannelStatusDB:  " << std::boolalpha << fUseChannelStatusDB << "\n";
      PMTsimulator->printConfiguration(log);
    } // if first time
    
    //
    // run the algorithm
    //

    // Prefer SimPhotons if available.
    // Make sure that there are parallel inputs for both formats;
    bool const useLitePhotons = !pmtVector.isValid();
    
    // storage for the photons that are not used (but still required)
    std::vector<sim::SimPhotons> fakePhotons;
    std::vector<sim::SimPhotonsLite> fakeLitePhotons;
    
    // these are the data collections passed to the algorithm
    std::vector<sim::SimPhotons> const* pInputPhotons = &fakePhotons;
    std::vector<sim::SimPhotonsLite> const* pInputLitePhotons = &fakeLitePhotons;
    
    // fill the relevant fake input vectors, and prepare the pointers
    if (useLitePhotons) {
      for (sim::SimPhotonsLite const& ph: *pmtLiteVector)
        fakePhotons.emplace_back(ph.OpChannel);
      pInputLitePhotons = pmtLiteVector.product();
    }
    else {
      for (sim::SimPhotons const& ph: *pmtVector)
        fakeLitePhotons.emplace_back(ph.OpChannel());
      pInputPhotons = pmtVector.product();
    }
    
    assert(pInputLitePhotons->size() == pInputPhotons->size());
    
    auto simphVecPtr = fWritePhotons? std::make_unique<std::vector<sim::SimPhotons>>(): nullptr;
    auto pulseVecPtr = std::make_unique<std::vector<raw::OpDetWaveform>>();
    auto baselineVecPtr = std::make_unique<std::vector<icarus::WaveformBaseline>>();
    
    for(auto const& [ photons, litePhotons ]
      : util::zip(*pInputPhotons, *pInputLitePhotons)
    ) {
      
      int const channel = photons.OpChannel();
      assert(channel == litePhotons.OpChannel);
      
      if (fUseChannelStatusDB
        && !lar::providerFrom<icarusDB::IPMTChannelStatusService>()->isGood(channel))
        continue;

      auto const& [ channelWaveforms, channelPedestals, photons_used, debug ]
        = PMTsimulator->simulate(photons, litePhotons, fDoDebugTree);
      assert(channelWaveforms.size() == channelPedestals.size());
      std::move(
        channelWaveforms.cbegin(), channelWaveforms.cend(),
        std::back_inserter(*pulseVecPtr)
        );
      std::move(
        channelPedestals.cbegin(), channelPedestals.cend(),
        std::back_inserter(*baselineVecPtr)
        );
      assert(pulseVecPtr->size() == baselineVecPtr->size());
      
      if (!useLitePhotons && simphVecPtr && photons_used)
        simphVecPtr->emplace_back(std::move(photons_used.value()));
      
      if(fDoDebugTree && debug && !debug->photons.empty()) {
        
        fDebugInfo = DebugInfo_t{}; // reset all information

        // fill debug tree variables
        fDebugInfo.run = e.id().run();
        fDebugInfo.event = e.id().event();
        fDebugInfo.channel = debug->opChannel;
        fDebugInfo.timeDelay_us = debug->timeDelay_us;
        fDebugInfo.triggerOffsetPMT_us = debug->triggerOffsetPMT_us;
        fDebugInfo.nSamples = debug->nSamples;
        fDebugInfo.nSubsamples = debug->nSubsamples;
        // fill full-window waveform
        fDebugInfo.waveform_adc = debug->waveform;
        // fill per photon info
        fDebugInfo.nInputPhotons = photons.size();
        fDebugInfo.nDetectedPhotons = debug->photons.size();
        for(auto const& phot : debug->photons) {
          fDebugInfo.photon_start_x.push_back(phot.startX);
          fDebugInfo.photon_start_y.push_back(phot.startY);
          fDebugInfo.photon_start_z.push_back(phot.startZ);
          fDebugInfo.photon_simTime_ns.push_back(phot.simTime_ns);
          fDebugInfo.photon_waveformTime_us.push_back(phot.trigTime_us);
          fDebugInfo.photon_tick.push_back(phot.tick);
          fDebugInfo.photon_subtick.push_back(phot.subtick);
        }
        // fill per PE deposit info
        fDebugInfo.nPEDeposits = debug->peDeposits.size();
        for(auto const& pedep : debug->peDeposits) {
          fDebugInfo.pedeposit_tick.push_back(pedep.tick);
          fDebugInfo.pedeposit_subtick.push_back(pedep.subtick);
          fDebugInfo.pedeposit_nPE.push_back(pedep.nPE);
          fDebugInfo.pedeposit_nEffectivePE.push_back(pedep.nEffectivePE);
          fDebugInfo.pedeposit_gainFactor.push_back(pedep.gainFactor());
        }

        fDebugTree->Fill();
      } // if debug tree
      
    } // for scintillation photons
    
    mf::LogInfo("SimPMTIcarus") << "Generated " << pulseVecPtr->size()
      << " waveforms out of " << pInputPhotons->size() << " optical channels.";
    
    // waveform baselines
    auto baselineAssns
      = std::make_unique<art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>>();
    art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr{ e };
    art::PtrMaker<icarus::WaveformBaseline> const makeBaselinePtr{ e };
    for (std::size_t const iWaveform: util::counter(pulseVecPtr->size())) {
      baselineAssns->addSingle
        (makeWaveformPtr(iWaveform), makeBaselinePtr(iWaveform));
    } // for
    
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
    e.put(std::move(baselineVecPtr));
    e.put(std::move(baselineAssns));
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
