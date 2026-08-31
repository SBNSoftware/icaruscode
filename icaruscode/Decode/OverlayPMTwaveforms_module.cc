/**
 * @file   OverlayPMTwaveforms_module.cc
 * @brief  Overlays simulated PMT waveforms on top of data PMT waveforms.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 23, 2026
 */

// ICARUS libraries
#include "sbncode/Utilities/AssnsUtils.h" // sbn::RebindAssociatedProducts()
#include "icarusalg/PMT/Algorithms/OverlayPMTwaveformAlg.h"
#include "icarusalg/Utilities/TimeInterval.h"
#include "icarusalg/Utilities/TimeIntervalConfig.h" // TimeIntervalOptionalTable
#include "sbnobj/ICARUS/PMT/Data/WaveformBaseline.h"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // detinfo::timescales::electronics_time
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds, ...
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"

// framework libraries
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"

// C/C++ standard libraries
#include <cassert>
#include <memory> // std::make_unique()
#include <optional>
#include <string>
#include <utility> // std::move()
#include <vector>


//------------------------------------------------------------------------------
using namespace util::quantities::time_literals;


//------------------------------------------------------------------------------
namespace sbn { class OverlayPMTwaveforms; }

//------------------------------------------------------------------------------
/**
 * @brief Overlays simulated PMT waveforms on top of data PMT waveforms.
 *
 * This module reads two collections of `raw::OpDetWaveform`:
 * - a data collection (the output segmentation matches these waveforms)
 * - a simulated collection (added on top of data where time overlaps)
 *
 * The overlay is performed by `sbn::opdet::OverlayPMTwaveformAlg`.
 * 
 * The overlay is based on the data waveforms: if simulation waveforms overlap
 * them, their overlapping samples are added to the data waveforms; the
 * simulation samples not in the data waveform coverage intervals are ignored:
 * data waveforms are never extended nor shrunk.
 * See `sbn::opdet::OverlayPMTwaveformAlg` algorithm documentation for more
 * details.
 * 
 * It is possible to limit the overlay to a (single) time interval
 * (`LimitToTimeInterval` configuration parameter), in which case the data
 * waveforms which do not overlap that interval will be copied in output
 * unchanged. The ones overlapping that interval undergo full overlay in their
 * whole range, even if wider than the limit interval.
 *
 *
 * Input
 * -----
 * 
 * * `std::vector<raw::OpDetWaveform>` (`DataWaveformTag`): input data
 *   waveforms.
 * * `std::vector<raw::OpDetWaveform>` (`SimWaveformTag`): simulation waveforms
 *   to be overlaid to the data ones. They are required not to overlap within
 *   one channel.
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>` (`DataBaselineAssns`)
 *   associations between data waveforms and their baseline value.
 *   These associations are rebased to the new waveforms. Because the actual
 *   baseline values are not used, the associated baseline data product does
 *   not need to actually be available.
 * * `art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>` (`SimBaselineAssns`)
 *   associations between simulated waveforms and their baseline value.
 *   The associated baseline data product must also be available.
 *   The recommended baseline is the one actually used for the generation of
 *   the simulated waveforms.
 * 
 * 
 * Output
 * ------
 * 
 * * `std::vector<raw::OpDetWaveform>`: overlaid data + simulation waveforms.
 *   There is always one overlay waveform for each data waveform, in the same
 *   order as the input. Channel, time interval coverage and baseline are the
 *   same for the overlay waveform as in its original data waveform.
 * * `std::vector<icarus::WaveformBaseline>` (if `DuplicateBaselines` is set):
 *   a new baseline collection, one-to-one
 *   @ref LArSoftProxyDefinitionParallelData "parallel data product", with the
 *   baseline values identical to the ones from the input waveforms.
 * * `art::Assns<icarus::WaveformBaseline, raw::OpDetWaveform>>`: association
 *   between the new overlaid waveforms and their baseline value. Because the
 *   overlaying does not change the waveform baseline, the baseline values
 *   are the same as the input one. If `DuplicateBaselines` is not set, the
 *   baseline objects associated to the input data waveforms are directly
 *   reassociated to the overlaid ones, otherwise the ones from the new data
 *   product are used.
 * 
 *
 * Requirements
 * ------------
 *
 * ### Noise disabled and check
 * 
 * Simulation of electronics noise must be disabled. Noise is already present
 * in the data waveforms and must not be double-counted.
 * 
 * The module checks that the simulated waveforms begin with a constant baseline
 * region with no noise (`BaselineCheckLength` configuration parameter):
 * if noise is disabled, the first samples of the waveform are expected to be
 * exactly equal to the baseline.
 * However, depending on the simulation settings, some waveforms may fail the
 * check (false positive). For this reason, a fraction of failure can be set
 * (`BaselineCheckForgivenessFraction`) that will be ignored.
 * The check is considered failed only if more than this fraction of waveforms
 * fail. For example, `BaselineCheckForgivenessFraction` set to `0.5` means that
 * up to half of the waveform can fail and still pass the check; a value of
 * `0.0` means that a single failing waveform will cause the check to fail.
 * 
 * This check can be completely disabled by setting the check length to zero.
 *
 *
 * Configuration
 * -------------
 *
 * * `DataWaveformTag` (input tag, mandatory): the data waveforms to be used.
 *   This data product drives the overlay.
 * * `DataBaselineAssns` (input tag, default: as `DataWaveformTag`): the
 *   associations between data waveforms and their baselines.
 * * `SimWaveformTag` (input tag, mandatory): the simulated waveforms to be
 *   overlaid on top of the data.
 * * `SimBaselineAssns` (input tag, default: as `SimWaveformTag`): the
 *   associations between data waveforms and their baselines.
 * * `LimitToTimeInterval` (time interval table, optional): if specified, only
 *   the data waveforms overlapping this interval are overlaid, while the others
 *   are copied verbatim in the output. Times are expressed on the electronics
 *   time scale (the same as the waveform timestamps), and the table format is
 *   documented in `icarus::ns::fhicl::TimeIntervalConfig`.
 * * `BaselineCheckLength` (microseconds, default: `"0.5 us"`): interval for
 *   the validation of the requirement of no noise in simulation waveforms.
 *   A zero value disables the check.
 * * `BaselineCheckForgivenessFraction` (real number, default: `0.5`): 
 * * `DuplicateBaselines` (flag, default: `false`): if set to true, a new
 *   baseline collection is produced instead of reusing the existing one.
 *   This is useful for modules that use it directly rather than via
 *   associations.
 * * `LogCategory` (string, default: `"OverlayPMTwaveforms"`): name of the
 *   message facility stream for module console output.
 *
 */
class sbn::OverlayPMTwaveforms: public art::SharedProducer {

  public:

  // --- BEGIN Configuration -----------------------------------------------------
  struct Config {

    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    using nanoseconds = util::quantities::intervals::nanoseconds;
    using microseconds = util::quantities::intervals::microseconds;
    using electronics_time = detinfo::timescales::electronics_time;
    
    struct BaselineCheckConfig {
      
      fhicl::Atom<microseconds> Length {
        Name{ "Length" },
        Comment{
          "expected duration of noiseless initial baseline in simulated waveforms (0 disables)"
          },
        0.5_us
        };

      fhicl::Atom<float> ForgivenessFraction {
        Name{ "ForgivenessFraction" },
        Comment{ "tolerated fraction of waveform failing the baseline check" },
        0.0
        };

    }; // BaselineCheckConfig
    

    fhicl::Atom<art::InputTag> DataWaveformTag {
      Name{ "DataWaveformTag" },
      Comment{ "tag of input data PMT waveforms" }
      };

    fhicl::OptionalAtom<art::InputTag> DataBaselineAssns {
      Name{ "DataBaselineAssns" },
      Comment{ "tag of baseline association for data waveforms (as DataWaveformTag if omitted)" }
      };

    fhicl::Atom<art::InputTag> SimWaveformTag {
      Name{ "SimWaveformTag" },
      Comment{ "tag of input simulated PMT waveforms" }
      };

    fhicl::OptionalAtom<art::InputTag> SimBaselineAssns {
      Name{ "SimBaselineAssns" },
      Comment{ "tag of baseline association for simulated waveforms (as SimWaveformTag if omitted)" }
      };

    icarus::ns::fhicl::TimeIntervalOptionalTable<electronics_time> LimitToTimeInterval {
      Name{ "LimitToTimeInterval" },
      Comment{
        "optional time interval (electronics time scale) where overlay is applied"
        }
      };

    fhicl::Atom<bool> DuplicateBaselines {
      Name{ "DuplicateBaselines" },
      Comment
        { "create a new baseline collection instead of using the input one" },
      false
      };

    fhicl::Table<BaselineCheckConfig> BaselineCheck {
      Name{ "BaselineCheck" },
      Comment{ "baseline check configuration" }
      };


    fhicl::Atom<std::string> LogCategory {
      Name{ "LogCategory" },
      Comment{ "MessageFacility category used for output messages" },
      "OverlayPMTwaveforms"
      };

  }; // struct Config

  using Parameters = art::SharedProducer::Table<Config>;
  // --- END Configuration -------------------------------------------------------

  explicit OverlayPMTwaveforms(Parameters const& config, art::ProcessingFrame const&);

  void produce(art::Event& event, art::ProcessingFrame const&) override;

  private:

  // aliases
  using nanoseconds = util::quantities::intervals::nanoseconds;
  using electronics_time = detinfo::timescales::electronics_time;

  using Interval_t = icarus::ns::util::TimeInterval<electronics_time>;
  
  using BaselineAssns_t = art::Assns<raw::OpDetWaveform, icarus::WaveformBaseline>;

  /// Exception thrown when a pointer was not valid.
  struct PtrProductNotFoundError;

  // --- BEGIN ---  Configuration  ---------------------------------------------
  
  /// Input data waveform tag.
  art::InputTag const fDataWaveformsTag;
  
  /// Input data baseline-to-waveform association tag.
  art::InputTag const fDataBaselineAssnsTag;
  
  /// Input simulated waveform tag.
  art::InputTag const fSimWaveformsTag;
  
  /// Input simulated baseline-to-waveform association tag.
  art::InputTag const fSimBaselineAssnsTag;
  
  /// Limit to data waveforms overlapping this interval.
  std::optional<Interval_t> const fLimit;
  
  /// Whether to produce a baseline collection data product.
  bool const fDuplicateBaselines;
  
  std::string const fLogCategory; ///< Name of console message stream.
  
  // ---  END  ---  Configuration  ---------------------------------------------

  nanoseconds const fOpticalTick; ///< Cached size of sampling tick.

  sbn::opdet::OverlayPMTwaveformAlg fAlgo; ///< Overlay algorithm.


  /// Prepares a waveform input collection for the algorithm.
  std::vector<sbn::opdet::OverlayPMTwaveformAlg::InputWaveform_t>
  prepareWaveformInput(
    std::vector<raw::OpDetWaveform> const& waveforms,
    art::FindOneP<icarus::WaveformBaseline> const& toBaseline
    ) const;

  /// Makes the message of `notFound` more meaningful with info from `event`.
  static art::Exception dressPtrNotFoundException(
    std::string const& message,
    PtrProductNotFoundError const& notFound,
    art::Event const& event
    );

}; // class sbn::OverlayPMTwaveforms


//------------------------------------------------------------------------------
namespace {
  
  /// Moves `data` to a unique pointer, ready for `art::Event::put()`.
  template <typename T>
  [[nodiscard]] std::unique_ptr<T> moveProduct(T& data)
    { return std::make_unique<T>(std::move(data)); }
  
} // local namespace


//------------------------------------------------------------------------------
//---  sbn::OverlayPMTwaveforms
//------------------------------------------------------------------------------
struct sbn::OverlayPMTwaveforms::PtrProductNotFoundError: public art::Exception
{
  
  art::ProductID id;
  std::size_t key = std::numeric_limits<std::size_t>::max();
  
  PtrProductNotFoundError(
    art::ProductID const& id = art::ProductID{},
    std::size_t key = std::numeric_limits<std::size_t>::max()
    )
    : art::Exception{ art::errors::ProductNotFound }
    , id{ id }, key{ key }
    {}
  
  template <typename T>
  explicit PtrProductNotFoundError(art::Ptr<T> const& ptr)
    : PtrProductNotFoundError{ ptr.id(), ptr.key() }
    {}
  
}; // sbn::OverlayPMTwaveforms::PtrProductNotFoundError


//------------------------------------------------------------------------------
sbn::OverlayPMTwaveforms::OverlayPMTwaveforms
  (Parameters const& config, art::ProcessingFrame const&)
  : art::SharedProducer{ config }
  // configuration
  , fDataWaveformsTag{ config().DataWaveformTag() }
  , fDataBaselineAssnsTag{ config().DataBaselineAssns().value_or(fDataWaveformsTag) }
  , fSimWaveformsTag{ config().SimWaveformTag() }
  , fSimBaselineAssnsTag{ config().SimBaselineAssns().value_or(fSimWaveformsTag) }
  , fLimit{ config().LimitToTimeInterval() }
  , fDuplicateBaselines{ config().DuplicateBaselines() }
  , fLogCategory{ config().LogCategory() }
  // caches
  , fOpticalTick{
      detinfo::makeDetectorTimings(
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob()
        ).OpticalClockPeriod()
    }
  , fAlgo{
      sbn::opdet::OverlayPMTwaveformAlg::Config{
        fOpticalTick,
        {
          config().BaselineCheck().Length(),
          config().BaselineCheck().ForgivenessFraction()
        },
        fLogCategory
      }
    }
{
  async<art::InEvent>();

  consumes<std::vector<raw::OpDetWaveform>>(fDataWaveformsTag);
  consumes<BaselineAssns_t>(fDataBaselineAssnsTag);

  consumes<std::vector<raw::OpDetWaveform>>(fSimWaveformsTag);
  consumes<BaselineAssns_t>(fSimBaselineAssnsTag);

  produces<std::vector<raw::OpDetWaveform>>();
  produces<BaselineAssns_t>();
  if (fDuplicateBaselines)
    produces<std::vector<icarus::WaveformBaseline>>();

  {
    mf::LogInfo log{ fLogCategory };

    log << "OverlayPMTwaveforms configuration:"
      << "\n - data waveforms: '" << fDataWaveformsTag.encode() << "'"
      << "\n - data baselines association: '" << fDataBaselineAssnsTag.encode() << "'"
      << "\n - simulated waveforms: '" << fSimWaveformsTag.encode() << "'"
      << "\n - simulated baselines association: '" << fSimBaselineAssnsTag.encode() << "'"
      ;
    if (fDuplicateBaselines)
      log << "\n - create a new baseline collection data product";
    else
      log << "\n - associate waveforms to the existing baseline data product";

    fAlgo.dumpConfiguration(log, "    ", "\n - ");

  } // configuration report block

} // sbn::OverlayPMTwaveforms::OverlayPMTwaveforms()


//------------------------------------------------------------------------------
void sbn::OverlayPMTwaveforms::produce(art::Event& event, art::ProcessingFrame const&) {

  detinfo::DetectorTimings const detTimings = detinfo::makeDetectorTimings(
    art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)
    );

  assert(detTimings.OpticalClockPeriod() == fOpticalTick);

  //
  // input fetching and preparation
  //
  auto const& dataHandle = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fDataWaveformsTag);
  auto const& simHandle = event.getValidHandle<std::vector<raw::OpDetWaveform>>(fSimWaveformsTag);

  std::vector<raw::OpDetWaveform> const& dataWaveforms = *dataHandle;
  std::vector<raw::OpDetWaveform> const& simWaveforms = *simHandle;

  art::FindOneP<icarus::WaveformBaseline> const findDataBaseline{
    dataHandle, event, fDataBaselineAssnsTag
    };
  if (!findDataBaseline.isValid()) {
    throw art::Exception{ art::errors::ProductNotFound }
      << "OverlayPMTwaveforms: baseline association product for data waveforms '"
      << fDataBaselineAssnsTag.encode() << "' is not available or not valid.\n";
  } // if invalid data baseline finder

  art::FindOneP<icarus::WaveformBaseline> const findSimBaseline{
    simHandle, event, fSimBaselineAssnsTag
    };
  if (!findSimBaseline.isValid()) {
    throw art::Exception{ art::errors::ProductNotFound }
      << "OverlayPMTwaveforms: baseline association product for simulated waveforms '"
      << fSimBaselineAssnsTag.encode() << "' is not available or not valid.\n";
  } // if invalid simulated baseline finder

  std::vector<sbn::opdet::OverlayPMTwaveformAlg::InputWaveform_t> dataIn;
  try {
    dataIn = prepareWaveformInput(dataWaveforms, findDataBaseline);
  }
  catch (PtrProductNotFoundError const& notFoundError) {
    throw dressPtrNotFoundException(
      "Error processing data input from '" + fDataWaveformsTag.encode() + "'.",
      notFoundError, event
      );
  }
  catch (art::Exception const& e) {
    throw art::Exception{ e.categoryCode(), "", e }
      << "Error while processing data input from '"
      << fDataWaveformsTag.encode() << "'.\n";
  }

  std::vector<sbn::opdet::OverlayPMTwaveformAlg::InputWaveform_t> simIn;
  try {
    simIn = prepareWaveformInput(simWaveforms, findSimBaseline);
  }
  catch (PtrProductNotFoundError const& notFoundError) {
    throw dressPtrNotFoundException(
      "Error processing simulation input from '"
        + fSimWaveformsTag.encode() + "'.",
      notFoundError, event
      );
  }
  catch (art::Exception const& e) {
    throw art::Exception{ e.categoryCode(), "", e }
      << "Error while processing simulation input from '"
      << fSimWaveformsTag.encode() << "'.\n";
  }
  
  //
  // overlay
  //
  sbn::opdet::OverlayPMTwaveformAlg::OverlaidWaveforms const algoRes
    = fAlgo.overlay(dataIn, simIn, fLimit);
  assert(algoRes.waveforms.size() == dataWaveforms.size());

  //
  // data products
  //
  std::vector<raw::OpDetWaveform> outWaveforms = std::move(algoRes.waveforms);

  // check that input and output waveforms share everything but sample values
  for ([[maybe_unused]] auto const& [ dataWaveform, overlayWaveform  ]
    : util::zip(dataWaveforms, outWaveforms)
  ) {
    assert(overlayWaveform.ChannelNumber() == dataWaveform.ChannelNumber());
    assert(overlayWaveform.TimeStamp() == dataWaveform.TimeStamp());
    assert(overlayWaveform.size() == dataWaveform.size());
  } // for output waveforms
  
  mf::LogDebug{ fLogCategory }
    << "OverlayPMTwaveforms wrote " << outWaveforms.size() << " waveforms.";

  std::optional<std::vector<icarus::WaveformBaseline>> outBaselines;
  BaselineAssns_t outBaselineAssns;
  
  art::PtrMaker<raw::OpDetWaveform> const makeWaveformPtr{ event };
  
  if (fDuplicateBaselines) {
    outBaselines.emplace();
    
    art::PtrMaker<icarus::WaveformBaseline> const makeBaselinePtr{ event };
    
    for (auto const& [ iWaveform, waveform ]: util::enumerate(outWaveforms)) {
      
      outBaselines->push_back(*findDataBaseline.at(iWaveform)); // copy
      
      outBaselineAssns.addSingle
        (makeWaveformPtr(iWaveform), makeBaselinePtr(iWaveform));
      
    } // for output waveforms
    
  }
  else {
    // associate the new waveforms to the old baselines
    outBaselineAssns = sbn::RebindAssociatedProducts(
      event.getProduct<BaselineAssns_t>(fDataBaselineAssnsTag),
      makeWaveformPtr
      );
  }
  assert(outBaselineAssns.size() == outWaveforms.size());
  
  event.put(moveProduct(outWaveforms));
  if (outBaselines)
    event.put(moveProduct(*outBaselines));
  event.put(moveProduct(outBaselineAssns));

} // sbn::OverlayPMTwaveforms::produce()


//------------------------------------------------------------------------------
auto sbn::OverlayPMTwaveforms::prepareWaveformInput(
  std::vector<raw::OpDetWaveform> const& waveforms,
  art::FindOneP<icarus::WaveformBaseline> const& toBaseline
) const
  -> std::vector<sbn::opdet::OverlayPMTwaveformAlg::InputWaveform_t>
{
  
  std::vector<sbn::opdet::OverlayPMTwaveformAlg::InputWaveform_t> input;
  input.reserve(waveforms.size());

  for (auto const& [ iWaveform, waveform ]: util::enumerate(waveforms)) {

    art::Ptr<icarus::WaveformBaseline> const& baseline = toBaseline.at(iWaveform);
    if (!baseline) {
      throw PtrProductNotFoundError{ baseline }
        << "OverlayPMTwaveforms: missing baseline associated to waveform #"
        << iWaveform << " (CH=" << waveform.ChannelNumber()
        << " TS=" << waveform.TimeStamp() << ").\n";
    } // if no baseline

    input.push_back({ &waveform, baseline.get() });

  } // for waveforms
  return input;
} // sbn::OverlayPMTwaveforms::prepareWaveformInput()


//------------------------------------------------------------------------------
art::Exception sbn::OverlayPMTwaveforms::dressPtrNotFoundException(
  std::string const& message,
  PtrProductNotFoundError const& notFoundError,
  art::Event const& event
) {
  art::Exception e{ notFoundError.categoryCode(), "", notFoundError };
  art::BranchDescription const* descr = (notFoundError.id != art::ProductID{})
    ? event.getProductDescription(notFoundError.id).get(): nullptr;
  e << "Could not retrieve an art::Ptr";
  if (descr) {
    e << "<" << descr->friendlyClassName()
      << "> pointing to entry #" << notFoundError.key << " of '"
      << descr->inputTag().encode() << "'";
  }
  e << ".\n" << message << "\n";
  return e;
} // sbn::OverlayPMTwaveforms::dressPtrNotFoundException()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::OverlayPMTwaveforms)
