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
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"

// C++ standard library
#include <chrono> // std::chrono::high_resolution_clock
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm> // std::transform()
#include <functional> // std::plus
#include <cstdlib> // std::size_t


namespace {
  template <typename T>
  T sqr(T v) { return v*v; }
} // local namespace

namespace icarus {
  namespace opdet {
    
    /**
     * @brief Algorithm class for the full simulation of PMT channels.
     *
     *
     */
    class PMTsimulationAlg {
      
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
      PMTsimulationAlg(Config const& config); 
      PMTsimulationAlg(fhicl::ParameterSet const& pset); 
      
      /**
       * @brief Set up for the next event.
       * @param larProp instance of `detinfo::LArProperties` to be used
       * @param detClocks instance of `detinfo::DetectorClocks` to be used
       * @param mainRandomEngine main random engine (quantum efficiency, etc.)
       * @param darkNoiseRandomEngine random engine for dark noise simulation
       * @param elecNoiseRandomEngine random engine for electronics noise simulation
       *
       * All random engines are required in this interface, even if the
       * configuration disabled noise simulation.
       */
      void setup(
        detinfo::LArProperties const& larProp,
        detinfo::DetectorClocks const& detClocks,
        CLHEP::HepRandomEngine& mainRandomEngine, 
        CLHEP::HepRandomEngine& darkNoiseRandomEngine, 
        CLHEP::HepRandomEngine& elecNoiseRandomEngine 
        );    
      
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
      
    
    void CreateFullWaveform
      (Waveform_t&, std::vector<unsigned int>&, sim::SimPhotons const&);
    
    void CreateOpDetWaveforms(raw::Channel_t const& opch,
			      Waveform_t const& wvfm,
			      std::vector<raw::OpDetWaveform>& output_opdets);
    
    double fSampling;       //wave sampling frequency (GHz)
    std::size_t fNsamples; //Samples per waveform

    double fQEbase;         ///< Uncorrected PMT quantum efficiency.
    double fQE;             //PMT quantum efficiency
    
    size_t fReadoutWindowSize;     ///ReadoutWindowSize in samples
    float  fPretrigFraction;       ///Fraction of window size to be before "trigger"
    float  fThresholdADC;          ///ADC Threshold for self-triggered readout
    int    fPulsePolarity;         ///Pulse polarity (=1 for positive, =-1 for negative)
    double  fTriggerOffsetPMT;      ///Time (us) relative to trigger when pmt readout starts
        
    double  fReadoutEnablePeriod;  ///Time (us) for which pmt readout is enabled

    size_t fPretrigSize;
    size_t fPosttrigSize;
// 
    bool fCreateBeamGateTriggers; ///Option to create unbiased readout around beam spill
    double fBeamGateTriggerRepPeriod; ///Repetition Period (us) for BeamGateTriggers
    size_t fBeamGateTriggerNReps; ///Number of beamgate trigger reps to produce
// 
//     //Single PE parameters
//     double fFallTime;       //fall time of 1PE in ns
//     double fRiseTime;      //rise time in ns
    double fTransitTime;   //to be added to pulse minimum time
    double sigma1;
    double sigma2;
    double fMeanAmplitude;  //in pC
    double fFallTime; //in ns
    double fRiseTime; //in ns

    void AddSPE(size_t time_bin, Waveform_t& wave); // add single pulse to auxiliary waveform
    /// Add `n` standard pulses starting at the specified `time_bin` of `wave`.
    void AddPhotoelectrons(size_t time_bin, unsigned int n, Waveform_t& wave) const;

    double Pulse1PE(double time) const;
    
    Waveform_t wsp; //single photon pulse vector
    
    int pulsesize; //size of 1PE waveform
    
    double fADC;      //charge to ADC convertion scale
    double fBaseline; //waveform baseline
    double fAmpNoise; //amplitude of gaussian noise
    double fDarkNoiseRate; //in Hz
    double fSaturation; //equivalent to the number of p.e. that saturates the electronic signal	
    
    void AddNoise(Waveform_t& wave); //add noise to baseline
    void AddDarkNoise(Waveform_t& wave); //add noise to baseline
    
    std::set<size_t> CreateBeamGateTriggers() const;
    std::set<size_t> FindTriggers(Waveform_t const& wvfm) const;

    detinfo::DetectorClocks const* fTimeService = nullptr; ///< DetectorClocks service provider.
    CLHEP::HepRandomEngine* fRandomEngine = nullptr; ///< Main random stream engine.
    CLHEP::HepRandomEngine* fDarkNoiseRandomEngine = nullptr; ///< Dark noise random stream engine.
    CLHEP::HepRandomEngine* fElecNoiseRandomEngine = nullptr; ///< Electronics noise random stream engine.
    
      /// Returns a random response whether a photon generates a photoelectron.
      bool KicksPhotoelectron() const;
    }; // class PMTsimulationAlg
  
  } // namespace opdet
} // namespace icarus


//-----------------------------------------------------------------------------
template <typename Stream>
void icarus::opdet::PMTsimulationAlg::printConfiguration
  (Stream&& out, std::string indent /* = "" */) const
{
  out
            << indent << "Baseline:          " << fBaseline << " ADC"
    << '\n' << indent << "ReadoutWindowSize: " << fReadoutWindowSize << " ticks"
    << '\n' << indent << "PreTrigFraction:   " << fPretrigFraction
    << '\n' << indent << "ThresholdADC:      " << fThresholdADC << " ADC"
    << '\n' << indent << "Saturation:        " << fSaturation << " ADC"
    << '\n' << indent << "PulsePolarity:     " << ((fPulsePolarity == 1)? "positive": "negative") << " (=" << fPulsePolarity << ")"
    << '\n' << indent << "Sampling:          " << fSampling << " MHz"
    << '\n' << indent << "Samples/waveform:  " << fNsamples << " ticks"
    ;
  if (fCreateBeamGateTriggers) 
    out << '\n' << indent << "Create " << fBeamGateTriggerNReps << " beam gate triggers, one every " << fBeamGateTriggerRepPeriod << " us.";
  else out << '\n' << indent << "Do not create beam gate triggers.";
  
  
  out << '\n' << indent << "... and more.";
  out << '\n';
} // icarus::opdet::PMTsimulationAlg::printConfiguration()


//-----------------------------------------------------------------------------

 
#endif // ICARUSCODE_LIGHT_ALGORITHMS_PMTSIMULATIONALG_H

