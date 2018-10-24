/**
 * @file    icaruscode/Light/Algorithms/PMTsimulationAlg.h
 * @brief   Algorithms for the simulation of ICARUS PMT channels.
 * @date    October 16, 2018
 * @see     icaruscode/Light/Algorithms/PMTsimulationAlg.cxx
 *
 * These algoritms were originally extracted from the module
 * `SimPMTICARUS_module.cc`, which was in turnb based on 
 * `SimPMTSBND_module.cc` by L. Paulucci and F. Marinho
 */

#ifndef ICARUSCODE_LIGHT_ALGORITHMS_PMTSIMULATIONALG_H
#define ICARUSCODE_LIGHT_ALGORITHMS_PMTSIMULATIONALG_H


#include <vector> // FIXME remove with larsoft v07_07_02/v07_08_00 (bug in DetectorClocks.h)
#include <string> // FIXME remove with larsoft v07_07_02/v07_08_00 (bug in DetectorClocks.h)

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataalg/DetectorInfo/LArProperties.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
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


namespace icarus {
  namespace opdet {
    
    // -------------------------------------------------------------------------
    /// Utility class to make a number become a unique type.
    template <typename T>
    struct Value {
      using Datum_t = T;
      
      constexpr Value(Datum_t v): fValue(v) {}
      
      operator Datum_t() const { return fValue; }
      
        private:
      Datum_t fValue;
    }; // class Value<>
    
    
    // -------------------------------------------------------------------------
    /**
     * @brief Describes the waveform from a single photoelectron.
     * 
     * This functor (class behaving like a function) describes the shape of the response
     * to a single photoelectron.
     *
     * It is currently implemented as an asymmetric Gaussian shape.
     */
    class PhotoelectronPulseWaveform {
      
        public:
     
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
        double amplitude,
        double peakTime,
        double sigmaLeft,
        double sigmaRight
        )
        : fAmplitude(amplitude)
        , fTransitTime(peakTime)
        , fSigmaL(sigmaLeft)
        , fSigmaR(sigmaRight)
        {}
      
      // @{
      // @name Parameter accessors.

      double peakTime() const { return fTransitTime; }
      double leftSigma() const { return fSigmaL; }
      double rightSigma() const { return fSigmaR; }
      double amplitude() const { return fAmplitude; }
      
      /// @}
       
      /**
       * @brief Evaluates the pulse at the given time.
       * @param time time to evaluate the shape at
       * @see `PhotoelectronPulseWaveform()`
       *
       * The unit and scale of the time are defined by the transition time passed
       * at construction.
       */
      double evaluateAt(double time) const
        {
          return Gaussian(time,
            peakTime(), ((time < peakTime())? leftSigma(): rightSigma()),
            amplitude()
            );
        } // evaluateAt()
      
      /// Alias of `evaluateAt()`.
      double operator() (double time) const { return evaluateAt(time); }
      
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
      static double Gaussian(double x, double mean, double sigma, double amplitude)
        { return amplitude * std::exp(-sqr((x - mean)/sigma)/2.0); }
      
        private: 
      double fAmplitude; ///< Amplitude of the Gaussian shapes at peak (transition).
      double fTransitTime; ///< Time of transition between the two forms of shape.
      double fSigmaL; ///< RMS parameter of the shape before transition.
      double fSigmaR; ///< RMS parameter of the shape after transition.
      
      template <typename T>
      static constexpr T sqr(T value) { return value * value; }
      
    }; // class PhotoelectronPulseWaveform
 
    
    // -------------------------------------------------------------------------
    /**
     * @brief Precomputed digitised shape of a given function.
     * @tparam T type of the samples being stored
     *
     * The sampling happens internally in double precision.
     */
    template <typename T>
    class DiscretePhotoelectronPulse {
        public:
      using PulseFunction_t = PhotoelectronPulseWaveform; ///< Type of shape.
      
      struct Time_t: public Value<double> {};
      struct Tick_t: public Value<std::size_t> {};
      
      /**
       * @brief Constructor: samples the pulse.
       * @param pulseShape the shape to be pulsed
       * @param samplingFreq frequency of samples
       *                     [inverse of time unit in `pulseSape`]
       * @param rightSigmas sample until this standard deviations after peak
       *
       * Samples start from time 0 (as defined by the pulse shape).
       */
      DiscretePhotoelectronPulse
        (PulseFunction_t&& pulseShape, double samplingFreq, double rightSigmas);
      
      /// Returns the length of the sampled pulse in ticks.
      std::size_t pulseLength() const { return fSampledShape.size(); }
      
      /// Returns the value sampled for the specified `tick` (the first is `0`).
      T operator[] (std::size_t tick) const { return fSampledShape[tick]; }
      
      /// Returns the value sampled for the specified `tick` (the first is `0`).
      T operator() (Tick_t tick) const { return this->operator[](tick); }
      
      /// Evaluates the shape at the specified time.
      T operator() (Time_t time) const { return fShape(time); }
      
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
      double samplingFrequency() const { return fSamplingFreq; }
      
      /// Returns the sampling period (inverse of frequency).
      double samplingPeriod() const { return 1.0 / samplingFrequency(); }
      
      /// Returns the duration of the waveform in time units.
      /// @see `pulseLength()`
      double duration() const { return pulseLength() / samplingFrequency(); }

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
      
      
        private:
      PulseFunction_t fShape; ///< Analytical shape of the pules.
      double fSamplingFreq;   ///< Sampling frequency.
      
      std::vector<T> const fSampledShape; ///< Pulse shape, discretised.
      
      /// Builds the sampling cache.
      static std::vector<T> sampleShape(
        PulseFunction_t const& pulseShape,
        double samplingFreq, double rightSigmas
        );
      
    }; // class DiscretePhotoelectronPulse<>
    

    // -------------------------------------------------------------------------

    /** ************************************************************************
     * @brief Algorithm class for the full simulation of PMT channels.
     *
     *
     */
    class PMTsimulationAlg {
      
        public:
    
      /// Type holding all configuration parameters for this algorithm.
      struct ConfigurationParameters_t {    
        
        /// @{
        /// @name High level configuration parameters.
    
        double QEbase;         ///< Uncorrected PMT quantum efficiency.
        
        size_t readoutWindowSize;     ///ReadoutWindowSize in samples
        float  pretrigFraction;       ///Fraction of window size to be before "trigger"
        float  thresholdADC;          ///ADC Threshold for self-triggered readout
        int    pulsePolarity;         ///Pulse polarity (=1 for positive, =-1 for negative)
        double triggerOffsetPMT;      ///Time (us) relative to trigger when pmt readout starts
            
        double  readoutEnablePeriod;  ///Time (us) for which pmt readout is enabled
    
        bool createBeamGateTriggers; ///Option to create unbiased readout around beam spill
        double beamGateTriggerRepPeriod; ///Repetition Period (us) for BeamGateTriggers
        size_t beamGateTriggerNReps; ///Number of beamgate trigger reps to produce
    
        double ADC;      //charge to ADC convertion scale
        double transitTime;   //to be added to pulse minimum time
        double meanAmplitude;  //in pC
        double fallTime; //in ns
        double riseTime; //in ns
    
        double baseline; //waveform baseline    double fAmpNoise; //amplitude of gaussian noise
        double ampNoise; //amplitude of gaussian noise
        double darkNoiseRate; //in Hz
        double saturation; //equivalent to the number of p.e. that saturates the electronic signal	
        /// @}
        
        /// @{
        /// @name Setup parameters

        
        detinfo::LArProperties const* larProp = nullptr; ///< LarProperties service provider.
        detinfo::DetectorClocks const* timeService = nullptr; ///< DetectorClocks service provider.
        CLHEP::HepRandomEngine* randomEngine = nullptr; ///< Main random stream engine.
        CLHEP::HepRandomEngine* darkNoiseRandomEngine = nullptr; ///< Dark noise random stream engine.
        CLHEP::HepRandomEngine* elecNoiseRandomEngine = nullptr; ///< Electronics noise random stream engine.
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
       * @return a list of optical waveforms, response to those photons
       *
       * Due to threshold readout, a single channel may result in multiple
       * waveforms, which are all on the same channel but disjunct in time.
       */
      std::vector<raw::OpDetWaveform> simulate
        (sim::SimPhotons const& photons);
      
      /// Prints the configuration into the specified output stream.
      template <typename Stream>
      void printConfiguration(Stream&& out, std::string indent = "") const;
      
        private:
      /// Type internally used for storing waveforms.
      using Waveform_t = std::vector<float>;
      
      ConfigurationParameters_t fParams; ///< Complete algorithm configuration.
    
      double fQE;            ///< PMT quantum efficiency.
      double fSampling;      ///< Wave sampling frequency [MHz].
      std::size_t fNsamples; ///< Samples per waveform.

      DiscretePhotoelectronPulse<float> wsp; /// Single photon pulse (sampled).
    
    void CreateFullWaveform
      (Waveform_t&, std::vector<unsigned int>&, sim::SimPhotons const&);
    
    void CreateOpDetWaveforms(raw::Channel_t const& opch,
			      Waveform_t const& wvfm,
			      std::vector<raw::OpDetWaveform>& output_opdets);


    void AddSPE(size_t time_bin, Waveform_t& wave); // add single pulse to auxiliary waveform
    /// Add `n` standard pulses starting at the specified `time_bin` of `wave`.
    void AddPhotoelectrons(size_t time_bin, unsigned int n, Waveform_t& wave) const;

    
    
    void AddNoise(Waveform_t& wave); //add noise to baseline
    void AddDarkNoise(Waveform_t& wave); //add noise to baseline
    
    std::set<size_t> CreateBeamGateTriggers() const;
    std::set<size_t> FindTriggers(Waveform_t const& wvfm) const;

    
      /// Returns a random response whether a photon generates a photoelectron.
      bool KicksPhotoelectron() const;
    }; // class PMTsimulationAlg
    
    
    // -------------------------------------------------------------------------
    
    /// Returns a new `PMTsimulationAlg` with an updated configuration.
    class PMTsimulationAlgMaker {
      
         public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        
        fhicl::Atom<double> ReadoutEnablePeriod {
          Name("ReadoutEnablePeriod"),
          Comment("Time for which PMT readout is enabled [us]")
          };
        
      }; // struct Config
      
      /// Constructor.
      PMTsimulationAlgMaker(Config const& config); 
      PMTsimulationAlgMaker(fhicl::ParameterSet const& pset); 
      
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
    
  
  } // namespace opdet
} // namespace icarus


//-----------------------------------------------------------------------------
//--- template implementation
//-----------------------------------------------------------------------------
// --- icarus::opdet::PhotoelectronPulseWaveform
// -----------------------------------------------------------------------------
template <typename Stream>
void icarus::opdet::PhotoelectronPulseWaveform::dump(Stream&& out,
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
template <typename T>
icarus::opdet::DiscretePhotoelectronPulse<T>::DiscretePhotoelectronPulse
  (PulseFunction_t&& pulseShape, double samplingFreq, double rightSigmas)
  : fShape(std::move(pulseShape))
  , fSamplingFreq(samplingFreq)
  , fSampledShape(sampleShape(fShape, fSamplingFreq, rightSigmas))
  {}


// -----------------------------------------------------------------------------
template <typename T>
std::vector<T> icarus::opdet::DiscretePhotoelectronPulse<T>::sampleShape
  (PulseFunction_t const& pulseShape, double samplingFreq, double rightSigmas)
{
  std::size_t const pulseSize = samplingFreq
    * (pulseShape.peakTime() + rightSigmas * pulseShape.rightSigma());
  std::vector<T> samples(pulseSize);
  for (std::size_t i = 0; i < pulseSize; ++i)
    samples[i] = pulseShape(static_cast<double>(i)/samplingFreq);
  return samples;
} // icarus::opdet::DiscretePhotoelectronPulse<T>::sampleShape()


//-----------------------------------------------------------------------------
template <typename T>
template <typename Stream>
void icarus::opdet::DiscretePhotoelectronPulse<T>::dump(Stream&& out,
  std::string const& indent, std::string const& firstIndent
  ) const
{
  out << firstIndent << "Sampled pulse waveform " << pulseLength()
    << " samples long (" << duration()
    << " time units long, sampled at " << samplingFrequency()
    << "); pulse shape:"
    << "\n" << indent;
  shape().dump(std::forward<Stream>(out), indent + "  ", "");
} // icarus::opdet::DiscretePhotoelectronPulse<T>::dump()


//-----------------------------------------------------------------------------
//--- icarus::opdet::PMTsimulationAlg
//-----------------------------------------------------------------------------
template <typename Stream>
void icarus::opdet::PMTsimulationAlg::printConfiguration
  (Stream&& out, std::string indent /* = "" */) const
{
  out
            << indent << "Baseline:          " << fParams.baseline << " ADC"
    << '\n' << indent << "ReadoutWindowSize: " << fParams.readoutWindowSize << " ticks"
    << '\n' << indent << "PreTrigFraction:   " << fParams.pretrigFraction
    << '\n' << indent << "ThresholdADC:      " << fParams.thresholdADC << " ADC"
    << '\n' << indent << "Saturation:        " << fParams.saturation << " ADC"
    << '\n' << indent << "PulsePolarity:     " << ((fParams.pulsePolarity == 1)? "positive": "negative") << " (=" << fParams.pulsePolarity << ")"
    << '\n' << indent << "Sampling:          " << fSampling << " MHz"
    << '\n' << indent << "Samples/waveform:  " << fNsamples << " ticks"
    ;
  if (fParams.createBeamGateTriggers) {
    out << '\n' << indent << "Create " << fParams.beamGateTriggerNReps
      << " beam gate triggers, one every " << fParams.beamGateTriggerRepPeriod << " us.";
  }
  else out << '\n' << indent << "Do not create beam gate triggers.";
  
  out << '\n' << indent << "... and more.";
  
  out << '\n' << indent << "Template photoelectron waveform settings:"
    << '\n';
  wsp.dump(std::forward<Stream>(out), indent + "  ");
  out << '\n';
} // icarus::opdet::PMTsimulationAlg::printConfiguration()


//-----------------------------------------------------------------------------

 
#endif // ICARUSCODE_LIGHT_ALGORITHMS_PMTSIMULATIONALG_H

