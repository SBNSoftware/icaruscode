/**
 * @file   icaruscode/PMT/ICARUSOpHitFinder_module.cc
 * @brief  LArSoft module for PMT hit finding.
 * @date   May 6, 2022
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * 
 * Offers `opdet::ICARUSOpHitFinder` (see in its documentation why).
 */

// ICARUS libraries
#include "icaruscode/PMT/OpReco/Algorithms/PedAlgoFixed.h"
#include "icaruscode/PMT/OpReco/Algorithms/OpRecoFactoryStuff.h"
#include "icaruscode/PMT/Data/WaveformRMS.h"
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"
// LArSoft libraries
#include "larana/OpticalDetector/OpHitFinder/AlgoCFD.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoFixedWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSlidingWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoRollingMean.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h" // pmtana::Waveform_t
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larreco/Calibrator/IPhotonCalibratorService.h"
#include "larreco/Calibrator/PhotonCalibratorStandard.h"

// framework libraries
#include "art/Framework/Core/ReplicatedProducer.h"
#include "art/Framework/Core/ProcessingFrame.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"

// C++ standard libraries
#include <algorithm> // std::copy_if(), std::binary_search()
#include <vector>
#include <variant>
#include <memory> // std::unique_ptr
#include <string>
#include <functional> // std::mem_fn()
#include <utility> // std::move()
#include <cassert>


// -----------------------------------------------------------------------------
namespace {
  
  using opdet::factory::Decl;
  
  /// Types of _art_ interface.
  using ArtTraits
    = opdet::factory::FWInterfaceTraits<art::Event, art::ConsumesCollector>;
  
  // ===========================================================================
  // BEGIN ===  DECLARE NEW ALGORITHMS HERE
  // ===========================================================================
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // optical hit algorithm factory
  opdet::factory::AlgorithmFactory<pmtana::PMTPulseRecoBase> const
  HitAlgoFactory {
      "Name" // find the name of the algorithm under "Name" in its configuration
    , opdet::factory::Decl<pmtana::AlgoThreshold    >{"Threshold"    }
    , opdet::factory::Decl<pmtana::AlgoSiPM         >{"SiPM"         }
    , opdet::factory::Decl<pmtana::AlgoSlidingWindow>{"SlidingWindow"}
    , opdet::factory::Decl<pmtana::AlgoFixedWindow  >{"FixedWindow"  }
    , opdet::factory::Decl<pmtana::AlgoCFD          >{"CFD"          }
    };

  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // pedestal algorithm factory
  template <typename Algo>
  using FWIFPedAlgo
    = opdet::factory::FWInterfaced<Algo, pmtana::PMTPedestalBase, ArtTraits>;
  
  opdet::factory::AlgorithmFactory
    <opdet::factory::FWInterfacedIF<pmtana::PMTPedestalBase, ArtTraits>> const
  PedAlgoFactory {
      "Name" // find the name of the algorithm under "Name" in its configuration
    , Decl<FWIFPedAlgo<pmtana::PedAlgoEdges      >>{"Edges"      }
    , Decl<FWIFPedAlgo<pmtana::PedAlgoRollingMean>>{"RollingMean"}
    , Decl<FWIFPedAlgo<pmtana::PedAlgoUB         >>{"UB"         }
    , Decl<FWIFPedAlgo<pmtana::PedAlgoFixed      >>{"Fixed"      }
    };
  
  // ===========================================================================
  // END =====  DECLARE NEW ALGORITHMS HERE
  // ===========================================================================
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
} // local namespace


// -----------------------------------------------------------------------------
namespace opdet { class ICARUSOpHitFinder; }
/**
 * @class opdet::ICARUSOpHitFinder
 * @brief Extracts PMT activity as optical hits (`recob::OpHit`).
 * 
 * This module runs hit-finding algorithms, and is tailored to allow for an
 * update of the hit finding parameters event by event.
 * 
 * 
 * Motivation
 * -----------
 * 
 * This module is a branch of the standard LArSoft `opdet::OpHitFinder` by
 * Gleb Sinev, Ben Jones and others.
 * That module itself is an interface to the hit finding mini-framework also
 * in `larana`, authored by Kazuhiro Terao.
 * 
 * That framework independently analyzes each waveform and performs hit finding
 * using general timing and geometry information (from `detinfo::DetectorClocks`
 * and `geo::GeometryCore` service providers), a threshold parameter and
 * configurable algorithms for baseline evaluation, hit finding and calibration.
 * 
 * It is believed that the baseline of ICARUS PMT is remarkably stable in the
 * short time, and it is conceivable to use a single steady value for many
 * waveforms on each channel. Unfortunately the existing module is not suitable
 * for a baseline that may change event by event or in time: the mini-framework
 * does not explicitly allow this, and the module is not equipped with the
 * necessary workarounds.
 * 
 * This module is a replica of LArSoft's, with the additional workarounds needed
 * to allow for a baseline learnt from an external source and changing event by
 * event.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * The configuration parameters are intentionally left the same as in
 * `opdet::OpHitFinder`, with a few exceptions and additions where needed.
 * 
 * * `InputModule` (input tag): data product with the optical waveforms to
 *   analyze.
 * * `GenModule` (input tag, mandatory): data product with the beam gate (either
 *   generated or from trigger information).
 * * `UseStartTime` (flag, default: `false`): store the start time in the hit
 *   instead of the peak time.
 * * `ChannelMasks` (list of channel numbers, default: empty): skip waveforms
 *   on the channels specified in this list.
 * * `HitAlgoPset` (table, mandatory): configuration of the hit finding
 *   algorithm; its content depends on the algorithm itself, but the following
 *   elements are nonetheless mandatory:
 *     * `Name` (text): algorithm name; the supported algorithms are hard-coded
 *       by name.
 * * `PedAlgoPset` (table, mandatory): configuration of the pedestal
 *   algorithm; its content depends on the algorithm itself, but the following
 *   elements are nonetheless mandatory:
 *     * `Name` (text): algorithm name; the supported algorithms are hard-coded
 *       by name.
 *     
 * * `HitThreshold` (real): hit threshold [ADC counts]
 * * `UseCalibrator` (flag, default: `false`): if set, use the calibration
 *   service configured in the job; otherwise, use the simpler calibration
 *   configured in this module (see the following parameters).
 * * `AreaToPE` (flag): whether the conversion factor goes from hit area (total
 *   ADC) or from amplitude (ADC), to photoelectrons.
 * * `SPEArea` (real): area or amplitude (depending on `AreaToPE`) of the
 *   response signal to a single photoelectron.
 * * `SPEShift` (real, default: `0`)
 * 
 * Noticeable differences with LArSoft's `opdet::OpHitFinder` configuration:
 * 
 * * multiple input collections are not supported, so the `InputModule` is a
 *   complete input tag and no list of instance names can be specified.
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * In order to support the algorithms that require an event-specific
 * configuration, the module is not shared. While the algorithms could be
 * written to keep track of different event configurations, at the time the
 * algorithms are executed on a waveform they have no mean to know which event
 * the waveforms belong to. Rather than writing complex workarounds on this
 * problem, the module is replicated so that there are multiple algorithms
 * (and multiple managers) but each of them sees an event at a time, in a way
 * that the module (replica) can predict.
 * 
 */
class opdet::ICARUSOpHitFinder: public art::ReplicatedProducer {
    public:
  
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> InputModule {
      Name{ "InputModule" },
      Comment{ "Data product with the optical waveforms to process" }
      };
    
    fhicl::OptionalAtom<art::InputTag> GenModule {
      Name{ "GenModule" },
      Comment{ "Data product with the beam gates" }
      };
    
    fhicl::Sequence<raw::Channel_t> ChannelMasks {
      Name{ "ChannelMasks" },
      Comment{ "List of channels to skip [none by default]" },
      std::vector<raw::Channel_t>{}
      };
    
    fhicl::Atom<float> HitThreshold {
      Name{ "HitThreshold" },
      Comment{ "Hit reconstruction threshold [ADC#]" }
      };
    
    fhicl::Atom<bool> UseStartTime {
      Name{ "UseStartTime" },
      Comment{ "Store the start time instead of the peak time" },
      false
      };
    
    fhicl::DelegatedParameter HitAlgoPset {
      Name{ "HitAlgoPset" },
      Comment{
        "parameters of the hit finding algorithm."
        " The parameters must include the algorithm \"Name\", one of: "
        + HitAlgoFactory.names(", ") + "."
        }
      };
    
    fhicl::DelegatedParameter PedAlgoPset {
      Name{ "PedAlgoPset" },
      Comment{ 
        "parameters of the pedestal extraction algorithm."
        " The parameters must include the algorithm \"Name\", one of: "
        + PedAlgoFactory.names(", ") + "."
        }
      };
    
    fhicl::Atom<bool> UseCalibrator {
      Name{ "UseCalibrator" },
      Comment{ "Use the photon calibration service configured in the job" },
      false
      };
    
    fhicl::Atom<bool> AreaToPE {
      Name{ "AreaToPE" },
      Comment{ "Whether the `SPEArea` parameter refers to area or amplitude" },
      [this](){ return !UseCalibrator(); }
      };
    
    fhicl::Atom<float> SPEArea {
      Name{ "SPEArea" },
      Comment{ "area or amplitude of PMT response to single photoelectron" },
      [this](){ return !UseCalibrator(); }
      };
    
    fhicl::Atom<float> SPEShift {
      Name{ "SPEShift" },
      Comment{ "shift on the single photoelectron response" },
      [this](){ return !UseCalibrator(); }
      };
    
  }; // Config
  
  using Parameters = art::ReplicatedProducer::Table<Config>;
  
  explicit ICARUSOpHitFinder
    (Parameters const& config, art::ProcessingFrame const& frame);

  virtual void produce(art::Event&, art::ProcessingFrame const&) override;

    private:
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  art::InputTag const fWaveformTags; ///< Input PMT waveform data product tags.
  art::InputTag const fBeamGateTag; ///< Data product with beam gates.
  /// Sorted list of channels to skip.
  std::vector<unsigned int> const fChannelMasks;
  float const fHitThreshold; ///< Hit finding threshold.
  bool const fUseStartTime; ///< Whether to store start instead of peak time.
  // --- END ---- Configuration parameters -------------------------------------
  
  // --- BEGIN -- Cached service values ----------------------------------------
  /// Storage for our own calibration algorithm, if any.
  std::unique_ptr<calib::IPhotonCalibrator> fMyCalib;
  
  unsigned int const fMaxOpChannel; ///< Number of channels in the detector.
  
  /// The calibration algorithm to be used (not owned).
  calib::IPhotonCalibrator const* fCalib = nullptr;
  
  // --- END ---- Cached service values ----------------------------------------
  
  
  // --- BEGIN -- Algorithms ---------------------------------------------------
  using FWInterfacedPedAlgo
    = opdet::factory::FWInterfacedIF<pmtana::PMTPedestalBase, ArtTraits>;
  
  pmtana::PulseRecoManager fPulseRecoMgr;
  std::unique_ptr<pmtana::PMTPulseRecoBase> const fThreshAlg;
  std::unique_ptr<FWInterfacedPedAlgo> const fPedAlg;
  
  // --- END ---- Algorithms ---------------------------------------------------
  
  /// Optionally reads the beam gates from `fBeamGateTag`, empty if none.
  std::vector<sim::BeamGateInfo const*> fetchBeamGates
    (art::Event const& event) const;

  /// Returns a vector with copies of only the waveforms not in masked channels.
  std::vector<raw::OpDetWaveform> selectWaveforms
    (std::vector<raw::OpDetWaveform> const& waveforms) const;


}; // opdet::ICARUSOpHitFinder


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
// ---  specializations
// -----------------------------------------------------------------------------
/// Framework interface to `pmtana::PedAlgoFixed`.
template <>
class opdet::factory::FWInterfaced
  <pmtana::PedAlgoFixed, pmtana::PMTPedestalBase, ArtTraits>
  : public opdet::factory::FWInterfacedBase
      <pmtana::PedAlgoFixed, pmtana::PMTPedestalBase, ArtTraits>
{
  using IFBase_t = opdet::factory::FWInterfacedBase
    <pmtana::PedAlgoFixed, pmtana::PMTPedestalBase, ArtTraits>;
  
    public:
  FWInterfaced(fhicl::ParameterSet const& pset): IFBase_t{ pset } {}
  
  
    private:
  
  /// Declares the consumables.
  virtual void doInitialize(typename IFBase_t::Module_t& module)
    {
      pmtana::PedAlgoFixed const& algo = *(IFBase_t::getAlgo());
      module.template consumes<std::vector<raw::OpDetWaveform>>
        (art::InputTag{ algo.waveformSourceName() });
      module.template consumes<std::vector<icarus::WaveformBaseline>>
        (art::InputTag{ algo.pedestalSourceName() });
      module.template consumes<std::vector<icarus::WaveformRMS>>
        (art::InputTag{ algo.pedestalRMSName() });
    }
  
  /// Reads and reports to the algorithm the pedestals for this event.
  virtual void doBeginEvent(typename IFBase_t::Event_t const& event)
    {
      pmtana::PedAlgoFixed& algo = *(IFBase_t::getAlgo());
      
      //
      // collect the data
      //
      auto const& waveforms
         = event.template getProduct<std::vector<raw::OpDetWaveform>>
         (art::InputTag{ algo.waveformSourceName() });
      auto const& pedSrc
         = event.template getProduct<std::vector<icarus::WaveformBaseline>>
         (art::InputTag{ algo.pedestalSourceName() });
      auto const& rmsSrc
         = event.template getProduct<std::vector<icarus::WaveformRMS>>
         (art::InputTag{ algo.pedestalRMSName() });
      assert(waveforms.size() == pedSrc.size());
      assert(waveforms.size() == rmsSrc.size());
      
      //
      // translate and store
      //
      
      // event number as ID:
      pmtana::PedAlgoFixed::InputSet_t inputSet { event.event() };
      
      inputSet.waveforms.reserve(waveforms.size());
      std::transform(
        waveforms.begin(), waveforms.end(), back_inserter(inputSet.waveforms),
        [](auto& obj){ return &obj; }
        );
      
      inputSet.pedestals.reserve(pedSrc.size());
      std::transform(
        pedSrc.begin(), pedSrc.end(), back_inserter(inputSet.pedestals),
        std::mem_fn(&icarus::WaveformBaseline::baseline)
        );
      
      inputSet.RMSs.reserve(rmsSrc.size());
      std::transform(
        rmsSrc.begin(), rmsSrc.end(), back_inserter(inputSet.RMSs),
        std::mem_fn(&icarus::WaveformRMS::RMS)
        );

      algo.setParameters(std::move(inputSet));
      
    } // doBeginEvent()
  
  /// Clear the event information (just "safety").
  virtual void doEndEvent(typename IFBase_t::Event_t const&)
    { IFBase_t::getAlgo()->clearParameters(); }
  
  
}; // FWInterfaced<pmtana::PedAlgoFixed, pmtana::PMTPedestalBase>


// -----------------------------------------------------------------------------
namespace {

  /// Sorts in place and returns a vector.
  template <typename T>
  std::vector<T>& sortVector(std::vector<T>& v)
    { std::sort(v.begin(), v.end()); return v; }
  
  /// Sorts and returns a copy of the vector.
  template <typename T>
  std::vector<T> sortedVector(std::vector<T> v)
    { sortVector(v); return v; }
  
  
  /**
   * @brief A "smart" pointer that may or may not own its data.
   * @tparam T type of the data to point to
   * 
   * This pointer can be used when depending on a run-time condition a data set
   * may come from an immutable, owning source (like `art::Event`) or from a
   * local processing (e.g. a selection of the above).
   * 
   * If initialized with a pointer to the data, this pointer will refer to it
   * without owning; if initialized with a r-value (reference) to data, the data
   * is stolen and owned from now on.
   */
  template <typename T>
  class MaybeOwnerPtr {
      public:
    using element_type = T;
    using pointer = T*;
    
    /// Points to the given pointer without owning it.
    MaybeOwnerPtr(pointer ptr) noexcept: fStorage{ ptr } {}
    
    /// Steals the data and owns it.
    MaybeOwnerPtr(element_type&& data): fStorage{ std::move(data) } {}
    
    pointer get() const { return std::visit(PointerGetter{}, fStorage); }
    
    pointer operator->() const { return get(); }
    auto& operator*() const { return *(get()); }
    
      private:
    
    struct PointerGetter {
      pointer operator() (pointer ptr) const noexcept { return ptr; }
      pointer operator() (element_type& item) const noexcept { return &item; }
    }; // PointerGetter
    
    using store_t = std::variant<element_type, pointer>;
    
    store_t fStorage; ///< Contains the data or pointer to it.
    
  }; // MaybeOwnerPtr
  
  
} // local namespace


// -----------------------------------------------------------------------------
opdet::ICARUSOpHitFinder::ICARUSOpHitFinder
  (Parameters const& params, art::ProcessingFrame const& frame)
  : art::ReplicatedProducer{ params, frame }
  // configuration
  , fWaveformTags{ params().InputModule() }
  , fBeamGateTag{ params().GenModule().value_or(art::InputTag{}) }
  , fChannelMasks{ sortedVector(params().ChannelMasks()) }
  , fHitThreshold{ params().HitThreshold() }
  , fUseStartTime{ params().UseStartTime() }
  // caches
  , fMaxOpChannel{ frame.serviceHandle<geo::Geometry>()->MaxOpChannel() }
  // algorithms
  , fPulseRecoMgr{}
  , fThreshAlg
    { HitAlgoFactory.create(params().HitAlgoPset.get<fhicl::ParameterSet>()) }
  , fPedAlg
    { PedAlgoFactory.create(params().PedAlgoPset.get<fhicl::ParameterSet>()) }
{
  
  //
  // declaration of input
  //
  consumes<std::vector<raw::OpDetWaveform>>(fWaveformTags);
  if (!fBeamGateTag.empty())
    mayConsume<std::vector<sim::BeamGateInfo>>(fBeamGateTag);
  
  //
  // calibration initialization
  //
  if (params().UseCalibrator()) {
    fCalib = frame.serviceHandle<calib::IPhotonCalibratorService>()->provider();
  }
  else {
    fMyCalib = std::make_unique<calib::PhotonCalibratorStandard>(
        params().AreaToPE()
      , params().SPEArea()
      , params().SPEShift()
      );
    fCalib = fMyCalib.get();
  } // if ... else
  
  //
  // register the algorithms in the manager
  //
  fPulseRecoMgr.AddRecoAlgo(fThreshAlg.get());
  fPulseRecoMgr.SetDefaultPedAlgo(&(fPedAlg->algo()));
  
  // framework hooks to the algorithms
  if (!fChannelMasks.empty() && (fPedAlg->algo().Name() == "Fixed")) {
    /* TODO
     * The reason why this is not working is complicate.
     * The interface of ophit mini-framework is quite minimal and tight,
     * it accepts only `std::vector<short>` input and nothing more.
     * The whole reason of this module is to overcome that limitation;
     * since there is no channel information, we have to guess which channel
     * this is the data of, and the way we do it is by address of the data;
     * `raw::OpDetWaveform` derives from `std::vector<short>`, so when we use
     * the original data product, the waveform data `std::vector<short>` is
     * actually the same object and we can compare with the address of the data
     * product, which we have because of the framework hooks here introduced;
     * but if the `ChannelMasks` is not empty, we are forced to copy the
     * original data in a new location (excluding the masked channels) and the
     * new data has a different address. In principle we could still try to
     * look for which `raw::OpDetWaveform` data product has the same content as
     * the `std::vector<short>` we were given, which may be a bit tedious
     * (thousand of vector comparisons to do, although all except the target one
     * will be short because the noise pattern are different).
     */
    throw art::Exception{ art::errors::Configuration }
      << "Pedestal algorithm \"Fixed\" will not work"
      " since a channel mask list is specified.\n";
  }
  fPedAlg->initialize(consumesCollector());
  
  //
  // declare output products
  //
  produces<std::vector<recob::OpHit>>();

} // opdet::ICARUSOpHitFinder::ICARUSOpHitFinder()


//----------------------------------------------------------------------------
void opdet::ICARUSOpHitFinder::produce
  (art::Event& event, art::ProcessingFrame const& frame)
{

  //
  // read and select the waveforms
  //
  auto const& allWaveforms
    = event.getProduct<std::vector<raw::OpDetWaveform>>(fWaveformTags);
  
  // this pointer owns the waveforms only if they come from `selectWaveforms()`
//   MaybeOwnerPtr<std::vector<raw::OpDetWaveform> const> waveforms
//     = fChannelMasks.empty()? &allWaveforms: selectWaveforms(allWaveforms);
  using WaveformsPtr_t = MaybeOwnerPtr<std::vector<raw::OpDetWaveform> const>;
  WaveformsPtr_t waveforms = fChannelMasks.empty()
    ? WaveformsPtr_t{ &allWaveforms }
    : WaveformsPtr_t{ selectWaveforms(allWaveforms) }
    ;
  
  if (fChannelMasks.empty() && (&allWaveforms != &*waveforms)) {
    // if this happens, contact the author
    throw art::Exception{ art::errors::LogicError }
      << "Bug in MaybeOwnerPtr!\n";
  }
  assert(fChannelMasks.empty() || (&allWaveforms == &*waveforms));
  
  //
  // run the algorithm
  //
  
  // framework hooks to the algorithms
  fPedAlg->beginEvent(event);
  
  std::vector<sim::BeamGateInfo const*> const beamGateArray
    = fetchBeamGates(event);

  auto const& geom = *(frame.serviceHandle<geo::Geometry>()->provider());
  auto const clockData
    = frame.serviceHandle<detinfo::DetectorClocksService const>()
      ->DataFor(event)
    ;
  
  std::vector<recob::OpHit> opHits;

  RunHitFinder(
    *waveforms, opHits,
    fPulseRecoMgr, *fThreshAlg,
    geom,
    fHitThreshold,
    clockData,
    *fCalib,
    fUseStartTime
    );
  
  mf::LogInfo{ "ICARUSOpHitFinder" }
    << "Found " << opHits.size() << " hits from " << waveforms->size()
    << " waveforms.";
  
  // framework hooks to the algorithms
  fPedAlg->endEvent(event);
  
  //
  // store results into the event
  //
  event.put(std::make_unique<std::vector<recob::OpHit>>(std::move(opHits)));
  
} // opdet::ICARUSOpHitFinder::produce()


//------------------------------------------------------------------------------
std::vector<sim::BeamGateInfo const*> opdet::ICARUSOpHitFinder::fetchBeamGates
  (art::Event const& event) const
{
  std::vector<const sim::BeamGateInfo*> beamGateArray;
  if (fBeamGateTag.empty()) return beamGateArray;
  try {
    event.getView(fBeamGateTag, beamGateArray);
  }
  catch (art::Exception const& err) {
    if (err.categoryCode() != art::errors::ProductNotFound) throw;
  }
  return beamGateArray;
} // opdet::ICARUSOpHitFinder::fetchBeamGates()


//----------------------------------------------------------------------------
std::vector<raw::OpDetWaveform> opdet::ICARUSOpHitFinder::selectWaveforms
  (std::vector<raw::OpDetWaveform> const& waveforms) const
{
  std::vector<raw::OpDetWaveform> selected;
  
  auto const isNotMasked = [this](raw::OpDetWaveform const& waveform)
    {
      return !std::binary_search
        (fChannelMasks.begin(), fChannelMasks.end(), waveform.ChannelNumber());
    };
  
  std::copy_if
    (waveforms.begin(), waveforms.end(), back_inserter(selected), isNotMasked);
  
  return selected;
} // selectWaveforms()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(opdet::ICARUSOpHitFinder)

// -----------------------------------------------------------------------------
