/**
 * @file    icaruscode/PMT/Algorithms/PMTsimulationAlg.h
 * @brief   Algorithms for the simulation of ICARUS PMT channels.
 * @date    October 16, 2018
 * @see     `icaruscode/PMT/Algorithms/PMTsimulationAlg.cxx`
 *
 * These algoritms were originally extracted from the module
 * `SimPMTICARUS_module.cc`, which was in turnb based on 
 * `SimPMTSBND_module.cc` by L. Paulucci and F. Marinho
 */

#ifndef ICARUSCODE_PMT_ALGORITHMS_PMTSIMULATIONALG_H
#define ICARUSCODE_PMT_ALGORITHMS_PMTSIMULATIONALG_H


// ICARUS libraries

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorTimingTypes.h" // picocoulomb
#include "lardataalg/Utilities/quantities_fhicl.h" // microsecond from FHiCL
#include "lardataalg/Utilities/quantities/spacetime.h" // microsecond, ...
#include "lardataalg/Utilities/quantities/frequency.h" // hertz, gigahertz
#include "lardataalg/Utilities/quantities/electronics.h" // tick, counts_f
#include "lardataalg/Utilities/quantities/electromagnetism.h" // picocoulomb

// framework libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// CLHEP libraries
#include "CLHEP/Random/RandEngine.h" // CLHEP::HepRandomEngine

// C++ standard library
#include <chrono> // std::chrono::high_resolution_clock
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <algorithm> // std::transform()
#include <utility> // std::forward()
#include <memory> // std::unique_ptr()
#include <functional> // std::plus
#include <cmath> // std::abs(), std::exp()
#include <cstdlib> // std::size_t


namespace icarus::opdet {


  using namespace util::quantities::electromagnetism_literals;


  // -------------------------------------------------------------------------
  /**
   * @brief Describes the waveform from a single photoelectron.
   * @tparam Time type of time unit to be used
   *
   * This functor (class behaving like a function) describes the shape of the
   * response to a single photoelectron.
   *
   * It is currently implemented as an asymmetric Gaussian shape.
   */
  template <typename T>
  class PhotoelectronPulseWaveform {

      public:
   /// Type for ADC counts (floating point).
   using ADCcount = util::quantities::counts_f;

   using Time = T; ///< Type of time being used.

    /**
     * @brief Constructor: assigns the parameters of the shape.
     * @brief amplitude the maximum amplitudes of the shape (at transition)
     * @brief peakTime the time of the maximum amplitude of the shape
     * @brief sigmaLeft the standard deviation of the shape before transition
     * @brief sigmaRight the standard deviation of the shape after transition
     *
     * The time parameters (`peakTime`, `sigmaLeft` and `sigmaRight`) must be
     * measured in same unit. The `peakTime` defined the position of the shape
     * with respect to time 0.
     *
     */
    PhotoelectronPulseWaveform(
      ADCcount amplitude,
      Time peakTime,
      Time sigmaLeft,
      Time sigmaRight
      )
      : fAmplitude(amplitude)
      , fTransitTime(peakTime)
      , fSigmaL(sigmaLeft)
      , fSigmaR(sigmaRight)
      {}

    // @{
    // @name Parameter accessors.

    Time peakTime() const { return fTransitTime; }
    Time leftSigma() const { return fSigmaL; }
    Time rightSigma() const { return fSigmaR; }
    ADCcount amplitude() const { return fAmplitude; }

    /// @}

    /**
     * @brief Evaluates the pulse at the given time.
     * @param time time to evaluate the shape at
     * @see `PhotoelectronPulseWaveform()`
     *
     * The unit and scale of the time are defined by the transition time passed
     * at construction.
     */
    ADCcount evaluateAt(Time time) const
      {
        return Gaussian(time,
          peakTime(), ((time < peakTime())? leftSigma(): rightSigma()),
          amplitude()
          );
      } // evaluateAt()

    /// Alias of `evaluateAt()`.
    ADCcount operator() (Time time) const { return evaluateAt(time); }

    // @{
    /**
     * @brief Prints on stream the parameters of this shape.
     * @tparam Stream type of stream to write into
     * @param out the stream to write into
     * @param indent indentation string, prepended to all lines except first
     * @param indentFirst indentation string prepended to the first line
     */
    template <typename Stream>
    void dump(Stream&& out,
      std::string const& indent, std::string const& firstIndent
      ) const;
    template <typename Stream>
    void dump(Stream&& out, std::string const& indent = "") const
      { dump(std::forward<Stream>(out), indent, indent); }
    // @}

    /// Returns the value of normal distribution at specified point.
    static ADCcount Gaussian(Time x, Time mean, Time sigma, ADCcount amplitude)
      { return amplitude * std::exp(-sqr((x - mean)/sigma)/2.0); }

      private:
    ADCcount fAmplitude; ///< Amplitude of the Gaussian shapes at peak (transition).
    Time fTransitTime; ///< Time of transition between the two forms of shape.
    Time fSigmaL; ///< RMS parameter of the shape before transition.
    Time fSigmaR; ///< RMS parameter of the shape after transition.

    template <typename V>
    static V sqr(V value) { return value * value; }

  }; // class PhotoelectronPulseWaveform<>


  // -------------------------------------------------------------------------
  /**
   * @brief Precomputed digitised shape of a given function.
   *
   * The sampling happens internally in double precision.
   */
  class DiscretePhotoelectronPulse {
      public:
    using gigahertz = util::quantities::gigahertz;
    using nanoseconds = util::quantities::nanosecond;

    /// Type of shape (times are in nanoseconds).
    using PulseFunction_t = PhotoelectronPulseWaveform<nanoseconds>;
    using ADCcount = PulseFunction_t::ADCcount;

    using Time_t = nanoseconds;
    using Tick_t = util::quantities::tick;

    static_assert(!std::is_same<Time_t, Tick_t>(),
      "Time and tick must be different types!");

    /**
     * @brief Constructor: samples the pulse.
     * @param pulseShape the shape to be pulsed; times in nanoseconds
     * @param samplingFreq frequency of samples [gigahertz]
     * @param rightSigmas sample until this standard deviations after peak
     *
     * Samples start from time 0 (as defined by the pulse shape).
     */
    DiscretePhotoelectronPulse
      (PulseFunction_t&& pulseShape, gigahertz samplingFreq, double rightSigmas);

    /// Returns the length of the sampled pulse in ticks.
    std::size_t pulseLength() const { return fSampledShape.size(); }

    /// Returns the value sampled for the specified `tick` (the first is `0`).
    ADCcount operator[] (std::size_t tick) const
      { return fSampledShape[tick]; }

    /// Returns the value sampled for the specified `tick` (the first is `0`).
    ADCcount operator() (Tick_t tick) const
      { return this->operator[](tick.value()); }

    /// Evaluates the shape at the specified time.
    ADCcount operator() (Time_t time) const { return fShape(time); }

    /// @{
    /// @name Iterator interface.

    auto cbegin() const { return fSampledShape.cbegin(); }
    auto cend() const { return fSampledShape.cend(); }
    auto begin() const { return fSampledShape.begin(); }
    auto end() const { return fSampledShape.end(); }

    /// @}

    /// Returns the function which was sampled.
    PulseFunction_t const& shape() const { return fShape; }

    /// Returns the sampling frequency (same units as entered).
    gigahertz samplingFrequency() const { return fSamplingFreq; }

    /// Returns the sampling period (inverse of frequency).
    nanoseconds samplingPeriod() const { return 1.0 / samplingFrequency(); }

    /// Returns the duration of the waveform in time units.
    /// @see `pulseLength()`
    nanoseconds duration() const { return pulseLength() * samplingPeriod(); }

    // @{
    /**
     * @brief Prints on stream the parameters of this shape.
     * @tparam Stream type of stream to write into
     * @param out the stream to write into
     * @param indent indentation string, prepended to all lines except first
     * @param indentFirst indentation string prepended to the first line
     */
    template <typename Stream>
    void dump(Stream&& out,
      std::string const& indent, std::string const& firstIndent
      ) const;
    template <typename Stream>
    void dump(Stream&& out, std::string const& indent = "") const
      { dump(std::forward<Stream>(out), indent, indent); }
    // @}

    /**
     * @brief Checks that the waveform tails not sampled are negligible.
     * @param limit threshold below which the waveform is considered negligible
     * @param outputCat _(default: empty)_ message facility output category
     *        to use for messages
     * @return whether the two tails are negligible
     *
     * If `outputCat` is empty (default) no message is printed.
     * Otherwise, in case of failure a message is sent to message facility
     * (under category `outputCat`) describing the failure(s).
     */
    bool checkRange(ADCcount limit, std::string const& outputCat) const;

      private:
    PulseFunction_t fShape;  ///< Analytical shape of the pules.
    gigahertz fSamplingFreq; ///< Sampling frequency.

    std::vector<ADCcount> const fSampledShape; ///< Pulse shape, discretized.

    /// Builds the sampling cache.
    static std::vector<ADCcount> sampleShape(
      PulseFunction_t const& pulseShape,
      gigahertz samplingFreq, double rightSigmas
      );

  }; // class DiscretePhotoelectronPulse<>


  // -------------------------------------------------------------------------

  /** ************************************************************************
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
   * The response of the PMT to a single photoelectron is a fixed shape
   * described in `icarus:opdet::PhotoelectronPulseWaveform`.
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
   *
   *
   */
  class PMTsimulationAlg {

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

      float ADC;      ///< charge to ADC conversion scale
      nanoseconds transitTime; ///< to be added to pulse minimum time
      picocoulomb meanAmplitude;
      nanoseconds fallTime;
      nanoseconds riseTime;

      ADCcount baseline; //waveform baseline
      ADCcount ampNoise; //amplitude of gaussian noise
      hertz darkNoiseRate;
      float saturation; //equivalent to the number of p.e. that saturates the electronic signal
      PMTspecs_t PMTspecs; ///< PMT specifications.
      bool doGainFluctuations; ///< Whether to simulate fain fluctuations.
      /// @}

      /// @{
      /// @name Setup parameters


      detinfo::LArProperties const* larProp = nullptr; ///< LarProperties service provider.
      detinfo::DetectorClocks const* timeService = nullptr; ///< DetectorClocks service provider.

      ///< Main random stream engine.
      CLHEP::HepRandomEngine* randomEngine = nullptr;

      ///< Random stream engine for gain fluctuations.
      CLHEP::HepRandomEngine* gainRandomEngine = nullptr;

      ///< Dark noise random stream engine.
      CLHEP::HepRandomEngine* darkNoiseRandomEngine = nullptr;

      ///< Electronics noise random stream engine.
      CLHEP::HepRandomEngine* elecNoiseRandomEngine = nullptr;

      /// @}

      /// @{
      /// @name Derivative configuration parameters.
      std::size_t pretrigSize() const { return pretrigFraction * readoutWindowSize; }
      std::size_t posttrigSize() const { return readoutWindowSize - pretrigSize(); }
      int expectedPulsePolarity() const
        { return ((ADC * meanAmplitude) < 0.0_pC)? -1: +1; }

      /// @}

    }; // ConfigurationParameters_t



    /// Constructor.
    PMTsimulationAlg(ConfigurationParameters_t const& config);


    /**
     * @brief Returns the waveforms originating from simulated photons.
     * @param photons all the photons simulated to land on the channel
     * @return a list of optical waveforms, response to those photons
     *
     * Due to threshold readout, a single channel may result in multiple
     * waveforms, which are all on the same channel but disjunct in time.
     */
    std::vector<raw::OpDetWaveform> simulate
    (sim::SimPhotons const& photons, sim::SimPhotons &photons_used);

    /// Prints the configuration into the specified output stream.
    template <typename Stream>
    void printConfiguration(Stream&& out, std::string indent = "") const;


    /// Converts rise time (10% to 90%) into a RMS under Gaussian hypothesis.
    template <typename T>
    static T riseTimeToRMS(T riseTime)
      {
        return riseTime / (
          std::sqrt(2.0)
          * (std::sqrt(-std::log(0.1)) - std::sqrt(-std::log(0.9)))
          );
      }


      private:
    /// Type internally used for storing waveforms.
    using Waveform_t = std::vector<ADCcount>;
    using WaveformValue_t = ADCcount::value_t; ///< Numeric type in waveforms.

    ConfigurationParameters_t fParams; ///< Complete algorithm configuration.

    double fQE;            ///< PMT quantum efficiency.
    megahertz fSampling;   ///< Wave sampling frequency [MHz].
    std::size_t fNsamples; ///< Samples per waveform.

    DiscretePhotoelectronPulse wsp; /// Single photon pulse (sampled).

  void CreateFullWaveform
  (Waveform_t&, sim::SimPhotons const&, sim::SimPhotons &);

  void CreateOpDetWaveforms(raw::Channel_t const& opch,
          Waveform_t const& wvfm,
          std::vector<raw::OpDetWaveform>& output_opdets);


  void AddSPE(tick time_bin, Waveform_t& wave); // add single pulse to auxiliary waveform
  /// Add `n` standard pulses starting at the specified `time_bin` of `wave`.
  void AddPhotoelectrons
    (tick time_bin, WaveformValue_t n, Waveform_t& wave) const;


  void AddNoise(Waveform_t& wave); //add noise to baseline
  void AddDarkNoise(Waveform_t& wave); //add noise to baseline

  /**
   * @brief Ticks in the specified waveform where some signal activity starts.
   * @return a collection of ticks with interesting activity
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
  std::set<optical_tick> FindTriggers(Waveform_t const& wvfm) const;

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
  std::set<optical_tick> CreateBeamGateTriggers() const;

    /// Returns a random response whether a photon generates a photoelectron.
    bool KicksPhotoelectron() const;
  }; // class PMTsimulationAlg


  // -------------------------------------------------------------------------

  /// Returns a new `PMTsimulationAlg` with an updated configuration.
  class PMTsimulationAlgMaker {

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
      fhicl::Atom<float> ADC {
        Name("ADC"),
        Comment("Voltage to ADC conversion factor")
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
      fhicl::Atom<nanoseconds> TransitTime {
        Name("TransitTime"),
        Comment("Single photoelectron: peak time from the beginning of the waveform [ns]")
        // mandatory
        };
      fhicl::Atom<picocoulomb> MeanAmplitude {
        Name("MeanAmplitude"),
        Comment("Single photoelectron: signal amplitude at peak [pC]")
        // mandatory
        };
      fhicl::Atom<nanoseconds> RiseTime {
        Name("RiseTime"),
        Comment("Single photoelectron: rise time (10% to 90%, sigma * ~1.687) [ns]")
        // mandatory
        };
      fhicl::Atom<nanoseconds> FallTime {
        Name("FallTime"),
        Comment("Single photoelectron: fall time (90% to 10%, sigma * ~1.687) [ns]")
        // mandatory
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
     * @param mainRandomEngine main random engine (quantum efficiency, etc.)
     * @param darkNoiseRandomEngine random engine for dark noise simulation
     * @param elecNoiseRandomEngine random engine for electronics noise simulation
     *
     * All random engines are required in this interface, even if the
     * configuration disabled noise simulation.
     */
    std::unique_ptr<PMTsimulationAlg> operator()(
      detinfo::LArProperties const& larProp,
      detinfo::DetectorClocks const& detClocks,
      CLHEP::HepRandomEngine& mainRandomEngine,
      CLHEP::HepRandomEngine& darkNoiseRandomEngine,
      CLHEP::HepRandomEngine& elecNoiseRandomEngine
      ) const;

      private:
    /// Part of the configuration learned from configuration files.
    PMTsimulationAlg::ConfigurationParameters_t fBaseConfig;

  }; // class PMTsimulationAlgMaker


  // -------------------------------------------------------------------------

} // namespace icarus::opdet


//-----------------------------------------------------------------------------
//--- template implementation
//-----------------------------------------------------------------------------
// --- icarus::opdet::PhotoelectronPulseWaveform
// -----------------------------------------------------------------------------
template <typename T>
template <typename Stream>
void icarus::opdet::PhotoelectronPulseWaveform<T>::dump(Stream&& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  out
       << firstIndent << "Pulse shape: asymmetric Gaussian with peak at "
          << peakTime() << " and amplitude " << amplitude() << ":"
    << '\n' << indent << "  (t <  " << peakTime() << "): sigma " << leftSigma()
    << '\n' << indent << "  (t >= " << peakTime() << "): sigma " << rightSigma()
    << '\n';
} // icarus::opdet::PhotoelectronPulseWaveform::dump()


// -----------------------------------------------------------------------------
// --- icarus::opdet::DiscretePhotoelectronPulse
// -----------------------------------------------------------------------------
inline icarus::opdet::DiscretePhotoelectronPulse::DiscretePhotoelectronPulse
  (PulseFunction_t&& pulseShape, gigahertz samplingFreq, double rightSigmas)
  : fShape(std::move(pulseShape))
  , fSamplingFreq(samplingFreq)
  , fSampledShape(sampleShape(fShape, fSamplingFreq, rightSigmas))
  {}


//-----------------------------------------------------------------------------
template <typename Stream>
void icarus::opdet::DiscretePhotoelectronPulse::dump(Stream&& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  out << firstIndent << "Sampled pulse waveform " << pulseLength()
    << " samples long (" << duration()
    << " long, sampled at " << samplingFrequency()
    << "); pulse shape:"
    << "\n" << indent;
  shape().dump(std::forward<Stream>(out), indent + "  ", "");
} // icarus::opdet::DiscretePhotoelectronPulse::dump()


//-----------------------------------------------------------------------------
//--- icarus::opdet::PMTsimulationAlg
//-----------------------------------------------------------------------------
template <typename Stream>
void icarus::opdet::PMTsimulationAlg::printConfiguration
  (Stream&& out, std::string indent /* = "" */) const
{
  out
            << indent << "Baseline:            " << fParams.baseline
    << '\n' << indent << "ReadoutWindowSize:   " << fParams.readoutWindowSize << " ticks"
    << '\n' << indent << "PreTrigFraction:     " << fParams.pretrigFraction
    << '\n' << indent << "ThresholdADC:        " << fParams.thresholdADC
    << '\n' << indent << "Saturation:          " << fParams.saturation << " p.e."
    << '\n' << indent << "doGainFluctuations:  "
      << std::boolalpha << fParams.doGainFluctuations
    << '\n' << indent << "PulsePolarity:       " << ((fParams.pulsePolarity == 1)? "positive": "negative") << " (=" << fParams.pulsePolarity << ")"
    << '\n' << indent << "Sampling:            " << fSampling
    << '\n' << indent << "Samples/waveform:    " << fNsamples << " ticks"
    << '\n' << indent << "Gain at first stage: " << fParams.PMTspecs.firstStageGain()
    ;
  if (fParams.createBeamGateTriggers) {
    out << '\n' << indent << "Create " << fParams.beamGateTriggerNReps
      << " beam gate triggers, one every " << fParams.beamGateTriggerRepPeriod << ".";
  }
  else out << '\n' << indent << "Do not create beam gate triggers.";
  
  out << '\n' << indent << "... and more.";
  
  out << '\n' << indent << "Template photoelectron waveform settings:"
    << '\n';
  wsp.dump(std::forward<Stream>(out), indent + "  ");
  out << '\n';
} // icarus::opdet::PMTsimulationAlg::printConfiguration()


//-----------------------------------------------------------------------------

 
#endif // ICARUSCODE_PMT_ALGORITHMS_PMTSIMULATIONALG_H

