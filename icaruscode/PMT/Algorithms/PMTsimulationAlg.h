/**
 * @file    icaruscode/PMT/Algorithms/PMTsimulationAlg.h
 * @brief   Algorithms for the simulation of ICARUS PMT channels.
 * @date    October 16, 2018
 * @see     `icaruscode/PMT/Algorithms/PMTsimulationAlg.cxx`
 *
 * These algoritms were originally extracted from the module
 * `SimPMTIcarus_module.cc`, which was in turn based on
 * `SimPMTSBND_module.cc` by L. Paulucci and F. Marinho.
 * Heavy hands of Wesley Ketchum (ketchum@fnal.gov) and Gianluca Petrillo
 * (petrillo@slac.standord.edu) for the ICARUS customization.
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_PMTSIMULATIONALG_H
#define ICARUSCODE_PMT_ALGORITHMS_PMTSIMULATIONALG_H


// ICARUS libraries
#include "icaruscode/PMT/Algorithms/DiscretePhotoelectronPulse.h"
#include "icaruscode/PMT/Algorithms/PhotoelectronPulseFunction.h"
#include "icarusalg/Utilities/SampledFunction.h"
#include "icarusalg/Utilities/FastAndPoorGauss.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h"
#include "lardataalg/Utilities/quantities_fhicl.h" // microsecond from FHiCL
#include "lardataalg/Utilities/quantities/spacetime.h" // microsecond, ...
#include "lardataalg/Utilities/quantities/frequency.h" // hertz, gigahertz
#include "lardataalg/Utilities/quantities/electronics.h" // tick, counts_f
#include "lardataalg/Utilities/quantities/electromagnetism.h" // picocoulomb

// framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C++ standard library
#include <vector>
#include <string>
#include <tuple>
#include <optional>
#include <ios> // std::boolalpha
#include <utility> // std::forward()
#include <memory> // std::unique_ptr()
#include <functional> // std::plus
#include <cmath> // std::abs(), std::exp()
#include <cstdlib> // std::size_t


// -----------------------------------------------------------------------------
namespace icarus::opdet {
  // in this Issue:

  using namespace util::quantities::electromagnetism_literals;
  using namespace util::quantities::electronics_literals;
  
  /// Type for single photon response shape function: nanosecond -> ADC counts.
  using SinglePhotonResponseFunc_t
    = DiscretePhotoelectronPulse::PulseFunction_t;

  template <typename SampleType> class OpDetWaveformMakerClass;
  
  class PMTsimulationAlg;
  
  class PMTsimulationAlgMaker;
  
} // namespace icarus::opdet


// -----------------------------------------------------------------------------
/// Helper class to cut a `raw::OpDetWaveform` from a longer waveform data.
template <typename SampleType>
class icarus::opdet::OpDetWaveformMakerClass {
  
    public:
  using Sample_t = SampleType;
  
  /// Type of waveform data.
  using WaveformData_t = std::vector<Sample_t>;
  
  using BufferRange_t = std::pair
    <detinfo::timescales::optical_tick, detinfo::timescales::optical_tick>;
  
  WaveformData_t const& fWaveform; ///< Full data from the PMT channel.
  
  /// Time of the first sample in waveform.
  detinfo::timescales::electronics_time const fPMTstartTime;
  
  util::quantities::nanosecond fSamplingPeriod; /// Sampling period.
  
  /// Constructor: waveform data, start time and sampling period
  OpDetWaveformMakerClass(
    WaveformData_t const& waveform,
    detinfo::timescales::electronics_time PMTstartTime,
    util::quantities::nanosecond samplingPeriod
    );
  
  // @{
  /// Returns an `raw::OpDetWaveform` with data at the `bufferRange`.
  raw::OpDetWaveform create
    (raw::Channel_t opChannel, BufferRange_t const& bufferRange) const;
  raw::OpDetWaveform operator()
    (raw::Channel_t opChannel, BufferRange_t const& bufferRange) const
    { return create(opChannel, bufferRange); }
  // @}

}; // class icarus::opdet::OpDetWaveformMakerClass<>


// -----------------------------------------------------------------------------
/**
 * @brief Algorithm class for the full simulation of PMT channels.
 *
 * The algorithm creates simulated PMT waveforms as read out by ICARUS,
 * including the generation of trigger primitives.
 * Contributions to the waveforms include:
 *  * physical photons
 *  * dark noise
 *  * electronics noise
 *
 * The algorithm processes an optical channel at a time, independently
 * and uncorrelated from the other channels.
 * For each channel, multiple waveforms may be generated according to the
 * readout parameters.
 *
 *
 * Activity sources
 * =================
 *
 * Physical photons
 * -----------------
 *
 * Photons are read from `sim::SimPhotons` data objects, each one pertaining
 * a single optical detector channel.
 * Each photon on the channel is assumed to have successfully reached the
 * external surface of the photocathode, with the wavelength shifter.
 * Depending on the upstream simulation, and in particular on the photon
 * visibility library settings, the photon might have also already passed
 * the wavelength shifting and even triggered the conversion to a detectable
 * photoelectron.
 *
 * Quantum efficiency is simulated to determine if each photon converts into
 * a photoelectron on the internal side of the photocathode. The target
 * quantum efficiency is specified via the `QE` configuration parameter.
 * It is assumed that some level of quantum efficiency has already been
 * simulated upstream: more precisely, that the quantum efficiency already
 * applied is in the amount returned by
 * `detinfo::LArProperties::ScintPreScale()`. Therefore:
 *
 * 1. the quantum efficiency applied here is only the residual one to go
 *    from `detinfo::LArProperties::ScintPreScale()` to the value in `QE`
 * 2. there is no implement here to _increase_ quantum efficiency, i.e.
 *    `QE` must not exceed `detinfo::LArProperties::ScintPreScale()`
 * 3. if the configuration specifies a target quantum efficiency `QE` larger
 *    than the one applied upstream
 *    (`detinfo::LArProperties::ScintPreScale()`), a warning message is
 *    printed, and no change to quantum efficiency is performed
 *
 * Note that if the upstream code has not applied any quantum efficiency,
 * the configuration should give a `detinfo::LArProperties::ScintPreScale()`
 * of 1.0.
 *
 * @note If the photon visibility library already includes the probability
 *       of the photon converting to a photoelectron, the quantum efficiency
 *       check here should be skipped by setting the efficiency to 1.
 *
 * For each converting photon, a photoelectron is added to the channel by
 * placing a template waveform shape into the channel waveform.
 *
 * The timestamp of each waveform is based on the same scale as the trigger
 * time, as defined by `detinfo::DetectorClocks::TriggerTime()`.
 * On that scale, the timestamp pins down the time of the first sample of
 * the waveform. Note that this is typically earlier than when the actual
 * signal starts. More precisely, the signal is defined to start at an
 * interest point (see `FindTriggers()` for their definition), and the
 * waveform starts (at tick #0) earlier than that by a fraction
 * `PreTrigFraction` of the readout window size `ReadoutWindowSize`
 * (both are configuration parameters of the algorithm), allowing for that
 * amount of pre-trigger data.
 *
 * The configuration parameter `TriggerOffsetPMT` describes how much earlier
 * than the trigger time the optical readout has started. Note that if
 * an interest point (see above) happens early after optical readout has
 * started, there might be not enough data to fill the pre-trigger data.
 * In such case, the interest point will just be located earlier than usual
 * within the final waveform. This situation may be caused for example by
 * asynchronous physics events like scintillation light from cosmic rays
 * or radioactive decay of the detector materials, or from a fluctuation
 * of the noise.
 *
 *
 * Photoelectrons
 * ---------------
 *
 * The response of the PMT to a single photoelectron is passed to the
 * algorithm as a full blown function of type `SinglePhotonResponseFunc_t`.
 * The function needs to be valid for the lifetime of the algorithm, since
 * the algorithm refers to without owning it, and it is expected not to
 * change during that time. See `icarus::opdet::SimPMTIcarus`
 *
 * To account for gain fluctuations, that shape is considered to correspond
 * to a nominal gain (`PMTspecs.gain` configuration parameter), which is
 * then fluctuated to obtain the effective gain. This feature can be
 * disabled by setting configuration parameter `FluctuateGain` to `false`.
 * The approximation used here is that the fluctuation is entirely due to
 * the first stage of multiplication. The gain on the first stage is
 * described as a random variable with Poisson distribution around the mean
 * gain. The gain on a single photoelectron at the first stage is, in fact,
 * an integral number in each case.
 * The time spread of the signal may be increased by the difference in time
 * of the different branches of the multiplication avalanche. Therefore,
 * increasing or decreasing the number of branches, as it is done by
 * changing the gain of the first stage, the time evolution of the signal
 * will also be likewise affected.
 * At this time we do not take this aspect into account in the simulation.
 * For the nominal settings of a Hamamatsu 5912 photomultiplier
 * (gain 10 ^7^, high multiplication on the first stage) the gain at the
 * first stage is around 20, causing a fluctuation of about 20%.
 * If the multiplication were equally distributed across the stages, that
 * fluctuation would be almost 45%.
 *
 * The first stage gain is computed by
 * `icarus::opdet::PMTsimulationAlg::ConfigurationParameters_t::PMTspecs_t::multiplicationStageGain()`.
 *
 *
 * Dark noise
 * -----------
 *
 * Dark noise, i.e. the noise originating by "spontaneous" emission of a
 * photoelectron in the photocathode without any external stimulation, is
 * simulated by randomly extracting the time such emission happens.
 * Each emission causes a photoelectron template waveform to be added at the
 * extracted time.
 * The rate of dark noise emission is set by configuration with
 * `DarkNoiseRate` parameter.
 *
 *
 * Electronics noise
 * ------------------
 *
 * Electronics noise is described by Gaussian fluctuations of a given
 * standard deviation, controlled by the configuration parameter `AmpNoise`.
 * No noise correlation is simulated neither in time nor in space.
 *
 *
 * Configuration
 * ==============
 *
 * PMT specifications
 * -------------------
 *
 * PMT specifications are used to evaluate the variance of the gain.
 * The details of the calculation are documented in
 * `icarus::opdet::PMTsimulationAlg::ConfigurationParameters_t::PMTspecs_t::multiplicationStageGain()`.
 *
 * The available parameters include:
 *
 * * `gain` (default: `1e7`): the nominal gain of the PMT; this is just a
 *     reference value.
 * * `voltageDistribution` is a sequence of values, one for each stage of
 *     the multiplication chain of the PMT. Each number represents the
 *     relative size of the resistance that determines the fall of the
 *     potential on that stage. Only the stages that contribute to the gain
 *     need to be included. The absolute value of each element is
 *     inconsequential. For example, a 10-stage PMT with the first stage
 *     having twice the resistance of all the other would be represented by
 *     the setting `[ 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]`.
 * * `dynodeK` (default: `0.75`) represents the dependence of the gain of
 *     a stage on the potential drop:
 *     @f$ \mu_{i} \propto (\Delta V_{i})^{k} @f$ (with @f$ \mu_{i} @f$ the
 *     gain for stage @f$ i @f$, @f$ \Delta V_{i} @f$ the drop of potential
 *     of that stage and @f$ k @f$ the parameter set by `dynodeK`.
 *
 *
 * Random number generators
 * -------------------------
 *
 * @anchor ICARUS_PMTSimulationAlg_RandomEngines
 *
 * Three independent random engines are currently used in the simulation:
 *
 * * "main" random engine:
 *     * residual quantum efficiency;
 *     * gain fluctuations;
 * * "dark noise" engine: dark current noise only;
 * * "electronics noise" engine: electronics noise only.
 *
 *
 * Structure of the algorithm
 * ===========================
 *
 * _This section needs completion._
 * 
 * The algorithm is serviceable immediately after construction.
 * Construction relies on a custom parameter data structure.
 * 
 * An utility, `PMTsimulationAlgMaker`, splits the set up in two parts:
 * 
 * 1. configuration, where the full set of parameters is learned;
 * 2. set up, where service providers, random number engines and the external
 *    single photon response function are acquired.
 * 
 * This is supposed to make the creation of the filling of parameter structure
 * and the creation of the algorithm easier.
 * At that point, a single `sim::SimPhotons` can be processed (`simulate()`)
 * at a time, or multiple at the same time (see the multithreading notes
 * below).
 * 
 * The function used to describe the single particle response is customizable
 * and it must in fact be specified by the caller, since there is no default
 * form. The function must implement the `PulseFunction_t` interface.
 * 
 * 
 * Multithreading notes
 * ---------------------
 * 
 * The algorithm processes one channel at a time, and it does not depend on
 * event-level information; therefore, the same algorithm object can be used
 * to process many events in sequence.
 * On the other hand, multithreading is impaired by the random number
 * generation, in the sense that multithreading will break reproducibility
 * if the random engine is not magically thread-resistant.
 * 
 * If the set up is event-dependent, then this object can't be used for
 * multiple events at the same time. There is no global state, so at least
 * different instances of the algorithm can be run at the same time.
 * In fact, the creation of an algorithm is expected to take negligible time
 * compared to its run time on a single event, and it is conceivable to create
 * a new algorithm instance for each event.
 * 
 */
class icarus::opdet::PMTsimulationAlg {

    public:
  using microseconds = util::quantities::microsecond;
  using nanoseconds = util::quantities::nanosecond;
  using hertz = util::quantities::hertz;
  using megahertz = util::quantities::megahertz;
  using picocoulomb = util::quantities::picocoulomb;
  using tick = util::quantities::tick;
  using ADCcount = DiscretePhotoelectronPulse::ADCcount;

  using time_interval = detinfo::timescales::time_interval;
  using optical_tick = detinfo::timescales::optical_tick;

  /// Type holding all configuration parameters for this algorithm.
  struct ConfigurationParameters_t {

    struct PMTspecs_t {

      /// Voltage distribution of the PMT. Each number represents the
      /// relative weight of the resistor between the two arms of a
      /// multiplication stage.
      std::vector<double> voltageDistribution;

      /// Gain from stage with voltage dV is proportional to dV^K.
      double dynodeK;

      double gain; ///< Total typical gain of a PMT.

      /**
       * @brief Returns the gain of the specified multiplication stage.
       * @param i index of multiplication stage (default: first, `1`)
       *
       * The gain is assumed to be the product of gains from each
       * multiplication stage. The stages are supposed to be connected
       * by @f$ N @f$ resistors of known value, whose weight relative to the
       * total (series) resistance is in `PMTvoltageDistribution`.
       * The total gain (known from `gain`) is:
       * @f[ \mu = \prod_{i} \mu_{i} @f]
       * and the gain of each stage @f$ i @f$ is
       * @f[ \mu_{i} = a (\Delta V_{i})^{k} @f]
       * with @f$ k @f$ a known constant (`dynodeK`) and @f$ a @f$ an
       * unknown one. Considered the total applied voltage (cathode to last
       * dynode) to be @f$ \Delta V @f$, the total resistance
       * @f$ R = \sum_{i} R_{i} @f$ and the weight of each stage
       * @f$ \rho_{i} = R_{i} / R @f$ (stored in `PMTvoltageDistribution`),
       * the potential on stage @f$ i @f$ is
       * @f[ \Delta V_{i} = \Delta V \rho_{i} @f]
       * (supporting a circuit current of @f$ \Delta V / R @f$) and
       * therefore
       * @f[ \mu = \prod_{i} a (\Delta V \rho_{i})^{k} @f]
       * that allows to find
       * @f[ a \Delta V^{k} = \sqrt[N]{\frac{\mu}{(\prod_{i} \rho_{i})^{k}}} @f]
       * With this constant known, the gain of each stage is also known:
       * @f[ \mu_{i} = a (\Delta V)^{k} (\rho_{i})^{k} @f]
       *
       * This function returns @f$ \mu_{i} @f$, with `i` starting from `1`
       * to `nDynodes()` included.
       */
      double multiplicationStageGain(unsigned int i = 1) const;

      /// Returns the gain from the first stage of PMT multiplication.
      double firstStageGain() const { return multiplicationStageGain(1U); }

      /// Number of dynodes in the PMTs.
      unsigned int nDynodes() const { return voltageDistribution.size(); }

      /// @}

      /// Sets `voltageDistribution` by stealing and normalizing `Rs`.
      void setVoltageDistribution(std::vector<double>&& Rs);

    }; // struct PMTspecs_t


    /// @{
    /// @name High level configuration parameters.

    double QEbase;         ///< Uncorrected PMT quantum efficiency.

    size_t readoutWindowSize;     ///< ReadoutWindowSize in samples
    float  pretrigFraction;       ///< Fraction of window size to be before "trigger"
    ADCcount thresholdADC; ///< ADC Threshold for self-triggered readout
    int    pulsePolarity;         ///< Pulse polarity (=1 for positive, =-1 for negative)
    microseconds triggerOffsetPMT; ///< Time relative to trigger when PMT readout starts TODO make it a `trigger_time` point

    microseconds readoutEnablePeriod;  ///< Time (us) for which pmt readout is enabled

    bool createBeamGateTriggers; ///< Option to create unbiased readout around beam spill
    microseconds beamGateTriggerRepPeriod; ///< Repetition Period (us) for BeamGateTriggers TODO make this a time_interval
    size_t beamGateTriggerNReps; ///< Number of beamgate trigger reps to produce

    unsigned int pulseSubsamples = 1U; ///< Number of tick subsamples.

    ADCcount baseline; //waveform baseline
    ADCcount ampNoise; //amplitude of gaussian noise
    bool useFastElectronicsNoise; ///< Whether to use fast generator for electronics noise.
    hertz darkNoiseRate;
    float saturation; //equivalent to the number of p.e. that saturates the electronic signal
    PMTspecs_t PMTspecs; ///< PMT specifications.
    bool doGainFluctuations; ///< Whether to simulate fain fluctuations.
    /// @}

    /// @{
    /// @name Setup parameters


    detinfo::LArProperties const* larProp = nullptr; ///< LarProperties service provider.

    detinfo::DetectorClocksData const* clockData = nullptr;

    /// Single photon response function.
    SinglePhotonResponseFunc_t const* pulseFunction;

    /// Main random stream engine.
    CLHEP::HepRandomEngine* randomEngine = nullptr;

    /// Random stream engine for gain fluctuations.
    CLHEP::HepRandomEngine* gainRandomEngine = nullptr;

    /// Dark noise random stream engine.
    CLHEP::HepRandomEngine* darkNoiseRandomEngine = nullptr;

    /// Electronics noise random stream engine.
    CLHEP::HepRandomEngine* elecNoiseRandomEngine = nullptr;

    /// Whether to track the scintillation photons used.
    bool trackSelectedPhotons = false;
    
    /// @}

    /// @{
    /// @name Derivative configuration parameters.

    std::size_t pretrigSize() const { return pretrigFraction * readoutWindowSize; }
    std::size_t posttrigSize() const { return readoutWindowSize - pretrigSize(); }

    /// @}

  }; // ConfigurationParameters_t



  /// Constructor.
  PMTsimulationAlg(ConfigurationParameters_t const& config);


  /**
   * @brief Returns the waveforms originating from simulated photons.
   * @param photons all the photons simulated to land on the channel
   * @return a list of optical waveforms, response to those photons,
   *         and which photons were used (if requested)
   *
   * Due to threshold readout, a single channel may result in multiple
   * waveforms, which are all on the same channel but disjunct in time.
   * 
   * The second element of the return value is optional and filled only
   * if the `trackSelectedPhotons` configuration parameter is set to `true`.
   * In that case, the returned `sim::SimPhotons` contains a copy of each of
   * the `photons` contributing to any of the waveforms.
   */
  std::tuple<std::vector<raw::OpDetWaveform>, std::optional<sim::SimPhotons>>
    simulate(sim::SimPhotons const& photons);

  /// Prints the configuration into the specified output stream.
  template <typename Stream>
  void printConfiguration(Stream&& out, std::string indent = "") const;



    private:
  
  using OpDetWaveformMaker_t
    = icarus::opdet::OpDetWaveformMakerClass<ADCcount>;
  
  /// Type internally used for storing waveforms.
  using Waveform_t = OpDetWaveformMaker_t::WaveformData_t;
  using WaveformValue_t = ADCcount::value_t; ///< Numeric type in waveforms.

  /// Type of sampled pulse shape: sequence of samples, one per tick.
  using PulseSampling_t = DiscretePhotoelectronPulse::Subsample_t;

  /// Type of member function to add electronics noise.
  using NoiseAdderFunc_t = void (PMTsimulationAlg::*)(Waveform_t&) const;


  // --- BEGIN -- Helper functors ----------------------------------------------
  /// Functor to convert tick point into a tick number and a subsample index.
  class TimeToTickAndSubtickConverter {

    double const fNSubsamples; ///< Number of subsamples.

      public:
    using SubsampleIndex_t = DiscretePhotoelectronPulse::SubsampleIndex_t;

    TimeToTickAndSubtickConverter(unsigned int nSubsamples)
      : fNSubsamples(static_cast<double>(nSubsamples)) {}

    /// Converts the `tick_d` in a subsample number and tick number.
    std::tuple<tick, SubsampleIndex_t> operator() (double const tick_d) const;

  }; // TimeToTickAndSubtickConverter


  /// Applies a random gain fluctuation to the specified number of
  /// photoelectrons.
  template <typename Rand>
  class GainFluctuator {

    std::optional<Rand> fRandomGain; ///< Random gain extractor (optional).
    double const fReferenceGain = 0.0; ///< Reference (average) gain.

      public:
    GainFluctuator() = default;
    GainFluctuator(double const refGain, Rand&& randomGain)
      : fRandomGain(std::move(randomGain))
      , fReferenceGain(refGain)
      {}

    /// Returns the new number of photoelectrons after fluctuation from `n`.
    double operator() (double const n);

  }; // GainFluctuator

  /// Returns a configured gain fluctuator object.
  auto makeGainFluctuator() const;

  // --- END -- Helper functors ------------------------------------------------


  ConfigurationParameters_t fParams; ///< Complete algorithm configuration.

  double fQE;            ///< PMT quantum efficiency.
  megahertz fSampling;   ///< Wave sampling frequency [MHz].
  std::size_t fNsamples; ///< Samples per waveform.

  DiscretePhotoelectronPulse wsp; /// Single photon pulse (sampled).

  NoiseAdderFunc_t const fNoiseAdder; ///< Selected electronics noise method.

  ///< Transformation uniform to Gaussian for electronics noise.
  static util::FastAndPoorGauss<32768U, float> const fFastGauss;

  /**
   * @brief Creates `raw::OpDetWaveform` objects from simulated photoelectrons.
   * @param photons the simulated list of photoelectrons
   * @param photons_used (_output_) list of used photoelectrons
   * @return a collection of digitised `raw::OpDetWaveform` objects
   * 
   * This function performs the digitization of a optical detector channel which
   * is collecting the photoelectrons in the `photons` list.
   * 
   * The photoelectrons are already screened for quantum efficiency: all of them
   * are considered for use. It is still possible, if a reduction of the quantum
   * efficiency is requested, that some of them are discarded.
   * 
   * The `photons_used` output argument is constructed and filled only if the
   * configuration of the algorithm requires the creation of a list of used
   * photons.
   * 
   */
  Waveform_t CreateFullWaveform(
    sim::SimPhotons const& photons,
    std::optional<sim::SimPhotons>& photons_used
    ) const;
  
  /**
   * @brief Creates `raw::OpDetWaveform` objects from a waveform data.
   * @param opChannel number of optical detector channel the data belongs to
   * @param waveform the waveform data
   * @return a collection of `raw::OpDetWaveform`
   * 
   * The waveform data is a sequence of samples on the optical detector channel,
   * starting at the beginning of the optical time clock, that is set by the
   * algorithm configuration as the global hardware trigger time (configured
   * in `detinfo::DetectorClocks`) and an offset.
   * It may be obtained from `CreateFullWaveform()`.
   * 
   * Multiple waveforms may be returned for a single channel, since a form of
   * zero suppression is applied by the simulated hardware.
   * The waveforms are guaranteed to be non-overlapping, non-contiguous and
   * sorted with increasing timestamp.
   * 
   * The procedure attempts to mimic the working mode of CAEN V1730B boards
   * as described on page 32 of the manual "UM2792_V1730_V1725_rev2.pdf"
   * in SBN DocDB 15024.
   * 
   * 
   * 
   * Algorithm details
   * ------------------
   * 
   * In the board readout language, the "trigger" is the per-channel
   * information that the signal on the channel has passed the threshold.
   * 
   * It is critical to understand the details of the generation of the triggers.
   * At the time of writing, the algorithm in `FindTriggers()` emits a trigger
   * every time the signal passes from under threshold to beyond threshold.
   * Assuming that the threshold is at less than one photoelectron,
   * this kind of behavior implies that every scintillation photons arriving on
   * the tail of the first (large?) signal will add a new waveform in tail.
   * It is not clear what happens if _two_ more triggers happen while the window
   * from the first trigger is still open.
   * 
   * The assumptions in the code include:
   * 
   * 1. overlapped triggers are not discarded (p. 32);
   * 2. board is set for self-trigger (p. 38);
   * 3. each waveform has a fixed duration (except the ones overlapping the
   *    previous one);
   * 4. when a trigger starts the recording window on a channel, only the last
   *    trigger happening during that window ("overlapping"), if any, is kept.
   * 5. at decoding time, contiguous buffers are merged in a single waveform
   * 
   * The fixed duration of the waveform if the sum of pre- and post-trigger
   * window length.
   * If during this window another trigger is emitted, a new window will follow
   * the current one and it will be long enough to contain all the post-trigger
   * data.
   */
  std::vector<raw::OpDetWaveform> CreateFixedSizeOpDetWaveforms
    (raw::Channel_t opChannel, Waveform_t const& waveform) const;
  
  
  /**
   * @brief Adds a pulse to a waveform, starting at a given tick.
   * @tparam Combine binary operation combining two ADC counts into one
   * @param pulse the sampling to add to the waveform
   * @param wave the waveform the pulse will be added to
   * @param time_bin the tick of the waveform where the pulse starts
   * @param combination how to combine the pulse and the waveform
   *
   * This is the internal implementation of `AddPhotoelectrons()`.
   *
   * The `combination` functor behaves as a binary function taking the
   * existing `wave` sample and the sample from the `pulse` at the same time
   * and returning their combination as a new sample value.
   */
  template <typename Combine>
  void AddPulseShape(
    PulseSampling_t const& pulse, Waveform_t& wave, tick const time_bin,
    Combine combination
    ) const;
  
  /**
   * @brief Adds a number of pulses to a waveform, starting at a given tick.
   * @param pulse the sampling to add, scaled, to the waveform
   * @param wave the waveform the pulses will be added to
   * @param time_bin the tick of the waveform where the pulses start being added
   * @param n the number of pulses added (it may be fractional)
   *
   * All the samples of `pulse` are scaled by the factor `n` and then _added_
   * to the sampling waveform `wave`, starting from the `time_bin` sample of
   * this waveform.
   *
   * The `pulse` samples are a sequence of ADC counts describing the single
   * photoelectron pulse shape. The waveform is also a sequence of samples
   * representing a optical detector channel digitized waveform, starting at
   * tick #0.
   */
  void AddPhotoelectrons(
    PulseSampling_t const& pulse, Waveform_t& wave, tick const time_bin,
    WaveformValue_t const n
    ) const;
  
  
  void AddNoise(Waveform_t& wave) const; //add noise to baseline
  /// Same as `AddNoise()` but using an alternative generator.
  void AddNoise_faster(Waveform_t& wave) const;
  // Add "dark" noise to baseline.
  void AddDarkNoise(Waveform_t& wave) const;
  
  /**
   * @brief Ticks in the specified waveform where some signal activity starts.
   * @return a collection of ticks with interesting activity, sorted
   * @see `CreateBeamGateTriggers()`
   *
   * We define an "interest point" a time when some activity in the
   * waveform is considered interesting enough to be recorded.
   * This method returns a list of interest points, in the form of the
   * index they are located at in the waveform `wvfm`.
   *
   * In general (but with the exception noted below), a time becomes an
   * interest point if the sample recorded at that time is above the threshold
   * set by the configuration parameter `ThresholdADC`.
   *
   * These interest points are local readout triggers that drive zero
   * suppression on the optical readout channel and that are not necessarily
   * causing any level of event trigger.
   *
   * This method also adds the mandatory beam gate interest points as
   * explained in `CreateBeamGateTriggers()` documentation. These are
   * additional interest points that are added independently of whether there
   * is actual interesting activity in there.
   */
  std::vector<optical_tick> FindTriggers(Waveform_t const& wvfm) const;
  
  /**
   * @brief Generate periodic interest points regardless the actual activity.
   * @return a collection of ticks where we pretend interesting activity to be
   * @see `FindTriggers()`
   *
   * This methods emits a list of interest points according to the algorithm
   * configuration. More precisely, if `CreateBeamGateTriggers` is configured
   * `true`, `BeamGateTriggerNReps` interest points are generated at
   * `BeamGateTriggerRepPeriod` intervals, starting from the beam gate time
   * as defined by `detinfo::DetectorClocks::BeamGateTime()`.
   *
   * See `FindTriggers()` for the meaning of "interest point".
   *
   * @note It is assumed that tick `0` happens at a time defined by
   *       `triggerOffsetPMT` configuration parameter _after_ the trigger
   *       (but since the value of that parameter is expected to be negative,
   *       tick `0` effectively happens _before_ the trigger).
   */
  std::vector<optical_tick> CreateBeamGateTriggers() const;

  /// Returns a random response whether a photon generates a photoelectron.
  bool KicksPhotoelectron() const;
}; // class PMTsimulationAlg


// -----------------------------------------------------------------------------

/// Returns a new `PMTsimulationAlg` with an updated configuration.
class icarus::opdet::PMTsimulationAlgMaker {

     public:
  using microseconds = util::quantities::microsecond;
  using nanoseconds = util::quantities::nanosecond;
  using hertz = util::quantities::hertz;
  using picocoulomb = util::quantities::picocoulomb;

  struct PMTspecConfig {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<double> DynodeK {
      Name("DynodeK"),
      Comment("exponent to the voltage in multiplication gain expression"),
      0.75 // middle of Hamamatsu 5912 range [ 0.7 -- 0.8 ]
      };
    fhicl::Sequence<double> VoltageDistribution {
      Name("VoltageDistribution"),
      Comment("voltage distribution (relative resistor value)"),
      { 17.4, 3.4, 5.0, 3.33, 1.67, 1.0, 1.2, 1.5, 2.2, 3.0, 2.4 }
      // Hamamatsu 5912
      };
    fhicl::Atom<double> Gain {
      Name("Gain"),
      Comment("average total gain (from one photoelectron to full signal)"),
      1.0e7
      };

  }; // struct PMTspecConfig

  
  /// Main algorithm FHiCL configuration.
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    //
    // readout settings
    //
    fhicl::Atom<microseconds> ReadoutEnablePeriod {
      Name("ReadoutEnablePeriod"),
      Comment("Time for which PMT readout is enabled [us]")
      // mandatory
      };
    fhicl::Atom<double> ReadoutWindowSize {
      Name("ReadoutWindowSize"),
      Comment
        ("Duration of a single PMT readout acquisition window [samples]")
      // mandatory
      };
    fhicl::Atom<float> Baseline {
      Name("Baseline"),
      Comment("Waveform baseline (may be fractional) [ADC]")
      // mandatory
      };
    fhicl::Atom<int> PulsePolarity {
      Name("PulsePolarity"),
      Comment("Pulse polarity: 1 for positive, -1 for negative")
      // mandatory
      };
    fhicl::Atom<double> PreTrigFraction {
      Name("PreTrigFraction"),
      Comment("fraction of the readout window located earlier than the readout trigger")
      // mandatory
      };

    //
    // PMT settings
    //
    fhicl::Atom<float> Saturation {
      Name("Saturation"),
      Comment("photomultiplier saturation (as number of photoelectrons)")
      // mandatory
      };
    fhicl::Atom<double> QE {
      Name("QE"),
      Comment("total photoelectron quantum efficiency")
      // mandatory
      };
    fhicl::Table<PMTspecConfig> PMTspecs {
      Name("PMTspecs"),
      Comment("collection of PMT characteristics"),
      };
    fhicl::Atom<bool> FluctuateGain {
      Name("FluctuateGain"),
      Comment("include gain fluctuation in the photoelectron response"),
      true
      };

    //
    // single photoelectron response
    //
    fhicl::Atom<unsigned int> PulseSubsamples {
      Name("PulseSubsamples"),
      Comment
        ("split each tick by this many subsamples to increase PMT timing simulation"),
      1U
      };

    //
    // dark noise
    //
    fhicl::Atom<hertz> DarkNoiseRate {
      Name("DarkNoiseRate"),
      Comment("Frequency of \"spontaneous\" emission of a dark noise photoelectron [Hz]")
      // mandatory
      };

    //
    // electronics noise
    //
    fhicl::Atom<double> AmpNoise {
      Name("AmpNoise"),
      Comment("RMS of the electronics noise fluctuations [ADC counts]")
      // mandatory
      };
    fhicl::Atom<bool> FastElectronicsNoise {
      Name("FastElectronicsNoise"),
      Comment
        ("use an approximate and faster random generator for electronics noise"),
      true
      };

    //
    // trigger
    //
    fhicl::Atom<float> ThresholdADC {
      Name("ThresholdADC"),
      Comment("Threshold for self-triggered readout [ADC counts]")
      // mandatory
      };
    fhicl::Atom<bool> CreateBeamGateTriggers {
      Name("CreateBeamGateTriggers"),
      Comment("Whether to create unbiased readout trigger at beam spill")
      // mandatory
      };
    fhicl::Atom<microseconds> BeamGateTriggerRepPeriod {
      Name("BeamGateTriggerRepPeriod"),
      Comment("Repetition period for beam gate generated readout triggers [us]")
      // mandatory
      };
    fhicl::Atom<std::size_t> BeamGateTriggerNReps {
      Name("BeamGateTriggerNReps"),
      Comment("Number of beam gate readout triggers to generate")
      // mandatory
      };
    fhicl::Atom<microseconds> TriggerOffsetPMT {
      Name("TriggerOffsetPMT"),
      Comment("Time  when readout begins, relative to readout trigger [us]")
      // mandatory
      };


  }; // struct Config


  /// Constructor.
  PMTsimulationAlgMaker(Config const& config);

  /**
   * @brief Creates and returns a new algorithm instance.
   * @param larProp instance of `detinfo::LArProperties` to be used
   * @param detClocks instance of `detinfo::DetectorClocks` to be used
   * @param SPRfunction function to use for the single photon response
   * @param mainRandomEngine main random engine (quantum efficiency, etc.)
   * @param darkNoiseRandomEngine random engine for dark noise simulation
   * @param elecNoiseRandomEngine random engine for electronics noise simulation
   * @param trackSelectedPhotons (default: `false`) keep track and return
   *                             a copy of the scintillation photons used
   *
   * All random engines are required in this interface, even if the
   * configuration disabled noise simulation.
   */
  std::unique_ptr<PMTsimulationAlg> operator()(
    detinfo::LArProperties const& larProp,
    detinfo::DetectorClocksData const& detClocks,
    SinglePhotonResponseFunc_t const& SPRfunction,
    CLHEP::HepRandomEngine& mainRandomEngine,
    CLHEP::HepRandomEngine& darkNoiseRandomEngine,
    CLHEP::HepRandomEngine& elecNoiseRandomEngine,
    bool trackSelectedPhotons = false
    ) const;

  /**
   * @brief Returns a data structure to construct the algorithm.
   * @param larProp instance of `detinfo::LArProperties` to be used
   * @param detClocks instance of `detinfo::DetectorClocks` to be used
   * @param SPRfunction function to use for the single photon response
   * @param mainRandomEngine main random engine (quantum efficiency, etc.)
   * @param darkNoiseRandomEngine random engine for dark noise simulation
   * @param elecNoiseRandomEngine random engine for electronics noise simulation
   *
   * Returns a data structure ready to be used to construct a
   * `PMTsimulationAlg` algorithm object, based on the configuration passed
   * at construction time and the arguments specified in this function call.
   * 
   * All random engines are required in this interface, even if the
   * configuration disabled noise simulation.
   */
  PMTsimulationAlg::ConfigurationParameters_t makeParams(
    detinfo::LArProperties const& larProp,
    detinfo::DetectorClocksData const& clockData,
    SinglePhotonResponseFunc_t const& SPRfunction,
    CLHEP::HepRandomEngine& mainRandomEngine,
    CLHEP::HepRandomEngine& darkNoiseRandomEngine,
    CLHEP::HepRandomEngine& elecNoiseRandomEngine,
    bool trackSelectedPhotons = false
    ) const;

    private:
  /// Part of the configuration learned from configuration files.
  PMTsimulationAlg::ConfigurationParameters_t fBaseConfig;

}; // class PMTsimulationAlgMaker



//-----------------------------------------------------------------------------
//--- template implementation
//-----------------------------------------------------------------------------
//--- icarus::opdet::OpDetWaveformMakerClass
// -----------------------------------------------------------------------------
template <typename SampleType>
icarus::opdet::OpDetWaveformMakerClass<SampleType>::OpDetWaveformMakerClass(
  WaveformData_t const& waveform,
  detinfo::timescales::electronics_time PMTstartTime,
  util::quantities::nanosecond samplingPeriod
  )
  : fWaveform(waveform)
  , fPMTstartTime(PMTstartTime)
  , fSamplingPeriod(samplingPeriod)
{}


// -----------------------------------------------------------------------------
template <typename SampleType>
raw::OpDetWaveform icarus::opdet::OpDetWaveformMakerClass<SampleType>::create
  (raw::Channel_t opChannel, BufferRange_t const& bufferRange) const
{
  std::size_t const start
    = std::min(std::size_t(bufferRange.first.value()), fWaveform.size());
  std::size_t const end
    = std::min(std::size_t(bufferRange.second.value()), fWaveform.size());
  assert(start <= end);
  
  // start of the waveform (tick #0) in the full optical reading
  raw::TimeStamp_t const timeStamp { fPMTstartTime + start * fSamplingPeriod };
  
  // create a new waveform preallocating enough room for the full buffer
  raw::OpDetWaveform outputWaveform(timeStamp, opChannel, end - start);
  
  // copy the buffer (need to unwrap the ADCcount value)
  std::transform(
    fWaveform.begin() + start,
    fWaveform.begin() + end,
    std::back_inserter(outputWaveform),
    [](auto sample){ return sample.value(); }
    );
  
  return outputWaveform;
} // icarus::opdet::OpDetWaveformMakerClass<>::create()


//-----------------------------------------------------------------------------
//--- icarus::opdet::PMTsimulationAlg
//-----------------------------------------------------------------------------
template <typename Stream>
void icarus::opdet::PMTsimulationAlg::printConfiguration
  (Stream&& out, std::string indent /* = "" */) const
{
  using namespace util::quantities::electronics_literals;

  out
            << indent << "Baseline:            " << fParams.baseline
    << '\n' << indent << "ReadoutWindowSize:   " << fParams.readoutWindowSize << " ticks"
    << '\n' << indent << "PreTrigFraction:     " << fParams.pretrigFraction
    << '\n' << indent << "ThresholdADC:        " << fParams.thresholdADC
    << '\n' << indent << "Saturation:          " << fParams.saturation << " p.e."
    << '\n' << indent << "doGainFluctuations:  "
      << std::boolalpha << fParams.doGainFluctuations
    << '\n' << indent << "PulsePolarity:       " << ((fParams.pulsePolarity == 1)? "positive": "negative") << " (=" << fParams.pulsePolarity << ")"
    << '\n' << indent << "Sampling:            " << fSampling;
  if (fParams.pulseSubsamples > 1U)
    out << " (subsampling: x" << fParams.pulseSubsamples << ")";
  out
    << '\n' << indent << "Samples/waveform:    " << fNsamples << " ticks"
    << '\n' << indent << "Gain at first stage: " << fParams.PMTspecs.firstStageGain()
    ;

  out << '\n' << indent << "Electronics noise:   ";
  if (fParams.ampNoise > 0_ADCf) {
    out << fParams.ampNoise << " RMS ("
      << (fParams.useFastElectronicsNoise? "faster": "slower") << " algorithm)";
  }
  else out << "none";

  if (fParams.createBeamGateTriggers) {
    out << '\n' << indent << "Create " << fParams.beamGateTriggerNReps
      << " beam gate triggers, one every " << fParams.beamGateTriggerRepPeriod << ".";
  }
  else out << '\n' << indent << "Do not create beam gate triggers.";

  out << '\n' << indent << "... and more.";

  out << '\n' << indent << "Template photoelectron waveform settings:"
    << '\n';
  wsp.dump(std::forward<Stream>(out), indent + "  ");
  
  out << '\n' << indent << "Track used photons:  "
    << std::boolalpha << fParams.trackSelectedPhotons
    << '\n';
} // icarus::opdet::PMTsimulationAlg::printConfiguration()


//-----------------------------------------------------------------------------


#endif // ICARUSCODE_PMT_ALGORITHMS_PMTSIMULATIONALG_H
