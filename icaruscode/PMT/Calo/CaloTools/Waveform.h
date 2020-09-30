////////////////////////////////////////////////////////////////////////////////
// Waveform class definition
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef  __WAVEFORM_H
#define __WAVEFORM_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <complex>

#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "TH1D.h"
#include "TMath.h"
#include "TF1.h"
#include "TSpectrum.h"

using namespace std;


namespace pmtcalo
{

  // Fit function for the pulse
  class PulseShapeFunction_ExpGaus {
    public:
      // use constructor to customize the function object
      double operator() (double* x, double * par){
        double t0 = par[0];
        //double mu = par[1];
        double w = par[1];
        double c = par[2];
        double a = par[3];

        double t = x[0];

        return a*c/2.0*TMath::Exp(c*c*w*w/2.0)*TMath::Exp(-1.0*c*(t-t0))
                                     * TMath::Erfc( 1.0/1.414* (c*w-(t-t0)/w) );
      }
  };

  class Waveform
  {

      public:

        struct Pulse
        {
          //Base pulse quantities
          //double start_time = 0; // Start time calculated from the amplitude
          double end_time=0;
          double time_peak = 0;
          double width = 0;
          double amplitude = 0;
          double integral = 0;

          // fitted quantities
          double fit_start_time = 0; // Start time calculated with the fit
          double error_start_time = 0;
          double fit_sigma = 0;
          double error_sigma = 0;
          double fit_mu = 0;
          double error_mu = 0;
          double fit_amplitude = 0;
          double error_amplitude = 0;
          double chi2 = -1;
          double ndf = -1;
          double fitstatus = -1; // O:good, >0: bad,  < 0: not working
        };

        typedef raw::OpDetWaveform Rawdigits_t;
        typedef vector<double> Waveform_t;
        typedef vector<complex<double>> Complex_t;

        Waveform();
        Waveform(fhicl::ParameterSet const& pset);
        ~Waveform();

        // Import data
        void loadData( Rawdigits_t raw_waveform );

        // Getters
        Rawdigits_t getRawWaveform(){ return m_raw_waveform; }
        Waveform_t getWaveform(){ return m_waveform; }
        double getBaselineMean(){return m_baseline_mean;};
        double getBaselineWidth(){return m_baseline_width;}
        double findMedian(std::vector<double> a, size_t n);
        void calculateWaveformMeanAndRMS( double &mean, double &width  );


        // Helpers
        void removeBaseline();
        bool isValidWaveform();
        void clean();
        TH1D* makeWaveformHist();

        // Noise removal
        Complex_t doFFT(Waveform_t m_time_domain);
        Waveform_t doIFFT(Waveform::Complex_t m_frequency_domain);
        void filterNoise();

        // Pulse and charge analysis
        Pulse getLaserPulse();
        Pulse getIntegral();
        void resetPulse(Pulse &pulse);
        float getTotalCharge();

      private:

        size_t m_nsamples;
        std::vector<double> m_trigger_time;
        double adc_to_mV = 0.122;
        double adc_to_pC = adc_to_mV*2.0*0.02; 
        double m_sampling_period = 2.0; // in ns

        Rawdigits_t m_raw_waveform;
        Waveform_t m_waveform;

        // Laser pulse
        int m_startbin;
        int m_nbins;

        bool m_dofit;
        double m_pulsethreshold; // in mV
        std::vector<double> m_fitrange;

        // Noise RMS
        size_t n_sample_baseline;
        double m_baseline_mean;
        double m_baseline_width;

        // Noise filter
        size_t window_size;
        bool reverse;
        float threshold;

      };
}

#endif //__WAVEFORM_H
