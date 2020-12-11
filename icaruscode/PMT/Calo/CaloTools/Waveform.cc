////////////////////////////////////////////////////////////////////////////////
// Waveform methods implementation
//
// mailto:ascarpell@bnl.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef WAVEFORM_CXX
#define WAVEFORM_CXX

#include "Waveform.h"

#include "TH1D.h"
#include "TMath.h"
#include "TF1.h"
#include "TSpectrum.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

namespace pmtcalo
{

  Waveform::Waveform(){}


  //----------------------------------------------------------------------------


  Waveform::Waveform(fhicl::ParameterSet const& pset)
  {

    // Baseline removal
    n_sample_baseline = pset.get<unsigned int>("NSamplesBaseline");

    // Noise removal
    window_size = pset.get<int>("NoiseWindowSize");
    reverse = pset.get<bool>("NoiseReverse");
    threshold = pset.get<double>("NoiseThreshold");

    // Charge integration and start time
    m_startbin = pset.get<int>("IntegralStartBin");
    m_nbins = pset.get<int>("IntegrationWindow");

    // Fit
    m_dofit = pset.get<bool>("DoFit");
    m_pulsethreshold = pset.get<double>( "PulseThreshold");
    m_fitrange = pset.get<std::vector<double>>("FitRange");

  }


  //------------------------------------------------------------------------------


  Waveform::~Waveform(){}


  //------------------------------------------------------------------------------


  void Waveform::loadData( Rawdigits_t raw_waveform )
  {

    m_nsamples= raw_waveform.size();
    m_raw_waveform = raw_waveform;

    for( auto w : raw_waveform ) {

      double value = -1*double(w); // Reverse polarity
      m_waveform.push_back( value );

    }

    removeBaseline();

  }


  //----------------------------------------------------------------------------


  double Waveform::findMedian(std::vector<double> a, size_t n) {
      // First we sort the array
      sort(a.begin(), a.end());

      // check for even case
      if (n % 2 != 0)
          return (double)a[n / 2];

      return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0;
  }


  void Waveform::calculateWaveformMeanAndRMS( double &mean, double &width  )
  {

    // Calculate the baseline as the mean values
    double tmpmean=0;
    for(size_t t=0; t<n_sample_baseline; t++){
      tmpmean += m_waveform.at(t);
    }
    tmpmean /= n_sample_baseline;

    // Calculate the stdev of the baseline
    for(size_t t=0; t<n_sample_baseline; t++){
      width += pow(m_waveform.at(t)-tmpmean, 2);
    }
    width = sqrt( width / (n_sample_baseline-1) );


    // we use the median for this
    mean = findMedian( m_waveform, m_waveform.size() );

    return;

  }


  //----------------------------------------------------------------------------


  void Waveform::removeBaseline()
  {

    m_baseline_mean=0; m_baseline_width=0;
    this->calculateWaveformMeanAndRMS(m_baseline_mean, m_baseline_width);

    // Subtract the baseline from the signal waveform
    std::transform(m_waveform.begin(), m_waveform.end(), m_waveform.begin(),
                              [ & ] (double x) { return x - m_baseline_mean; });

    }


    //--------------------------------------------------------------------------


    bool Waveform::isValidWaveform()
    {
      bool isValid;

      if(m_raw_waveform.size()==0){
        isValid=false;
      }
      else {
        isValid=true;
      }

      return isValid;
    }


    //------------------------------------------------------------------------------


    Waveform::Complex_t Waveform::doFFT(Waveform::Waveform_t m_time_domain)
    {
        Eigen::FFT<double> fft;
        Waveform::Complex_t  m_frequency_domain;
        fft.fwd(m_frequency_domain, m_time_domain);
        return m_frequency_domain;
    }


    Waveform::Waveform_t Waveform::doIFFT(Waveform::Complex_t m_frequency_domain)
    {
        Eigen::FFT<double> fft;
        Waveform::Waveform_t  m_time_domain;
        fft.inv(m_time_domain, m_frequency_domain);
        return m_time_domain;
    }

    //------------------------------------------------------------------------------


    void Waveform::filterNoise()
    {

      //Noise filter algorithm. It will produce a new waveform object after noise
      //mitigation. Noise patterns to be mitigated are selected on a given window at
      //the beginning or at end of the waveform.
      //Arguments:
      //  window_size: set the window to produce the nosie model
      //  reverse: if true, the noise window is calculated at the end of the waveform


      Waveform::Waveform_t tmp_noise(window_size);
      copy(m_waveform.begin(), m_waveform.begin()+window_size, tmp_noise.begin());

      if( reverse ){
        copy(m_waveform.end()-window_size, m_waveform.end(), tmp_noise.begin());
      }

      // Get the noise spectra and get rid of some of the most nasty frequencies
      // which survive above a threshold
      vector<complex<double>> spec = doFFT(tmp_noise);
      vector<complex<double>> tmp_spec(spec.size());
      for(size_t i=0; i<spec.size(); i++ )
      {
        double pwr = sqrt( pow(spec[i].real(), 2) + pow(spec[i].imag(), 2) );
        if( pwr  > threshold ){
          tmp_spec[i] = spec[i];
        }
      }

      // Produce the filtered waveform in time domain
      Waveform::Waveform_t tmp_filter = doIFFT(tmp_spec);

      // Now we mirror the tmp_filter waveform to match the same length of the
      // original waveform. The mirroring is reasonable in this case since the
      // base noise we would like to remove is periodic.
      tmp_filter.resize(m_nsamples);

      size_t groups = ceil(float(m_nsamples)/float(window_size));
      for(size_t group=0; group<groups+1; group++ )
      {
        for( size_t i=0; i<window_size; i++ ){
          if( int(group*window_size+i) < int(m_nsamples) ){

            m_waveform[group*window_size+i] -= tmp_filter[i];

          }
          else{ break; }
        }
      }
    }

    //--------------------------------------------------------------------------


      void Waveform::resetPulse(Waveform::Pulse &pulse)
      {
        // Basic pulse quantities
        //pulse.start_time = 0; // Start time calculated from the amplitude
        pulse.end_time=0;
        pulse.time_peak = 0;
        pulse.width = 0;
        pulse.amplitude = 0;
        pulse.integral = 0;

        // Fitted quantities
        pulse.fit_start_time = 0; // Start time calculated with the fit
        pulse.error_start_time = 0;
        pulse.fit_sigma = 0;
        pulse.error_sigma = 0;
        pulse.fit_mu = 0;
        pulse.error_mu = 0;
        pulse.fit_amplitude = 0;
        pulse.error_amplitude = 0;
        pulse.chi2 = -1;
        pulse.ndf = -1;
        pulse.fitstatus = -1; // O:good, >0: bad,  < 0: not working
      }


    //--------------------------------------------------------------------------


    Waveform::Pulse Waveform::getLaserPulse() {

      // Signal integral over the fixed time window given by m_startbin and
      // lasting for m_nbins;
      //
      // Of this pulse it is always caluclated the charge integrated over
      // the given m_nbins interval, the maximum, the position at maximum
      //
      // if m_dofit is set to true:
      //  1) An exp gauss function is fitted to the pulse
      //  2) A limit on the peak amplutude is configurable to skip waveform under
      //     the selected threshold
      //  3) the fit is perfomed in a rage around the maximum, controlled my
      //     the m_fitrange parameter
      //  4) Fit paramters, errors, chi2, ndf, and covergence flag are returned

      // Now we create the pulse and we fill in the basic information
      Pulse temp_pulse = this->getIntegral();

      // Now we perform the fit and add the new information
      if( m_dofit && (m_pulsethreshold < temp_pulse.amplitude) ) {

        double t_min = temp_pulse.time_peak - m_fitrange[0]; // in ns
        double t_max = temp_pulse.time_peak + m_fitrange[1]; // in ns

        if( t_min < m_startbin*m_sampling_period
                            || t_max > (m_startbin+m_nbins)*m_sampling_period ){
            t_min = m_startbin*m_sampling_period;
            t_max = (m_startbin+m_nbins)*m_sampling_period;
        }

        // Crop the area around the pulse of the waveform and return a TH1D
        TH1D* h_pulse = this->makeWaveformHist();

        char funcname[100]; sprintf(funcname, "func_pulse_direct");

        PulseShapeFunction_ExpGaus function_obj;
        TF1* func = new TF1(funcname, function_obj, t_min, t_max, 4);
        func->SetParNames("t0","#mu","#sigma","a");
        func->SetParameters(temp_pulse.time_peak - 5, 2, 0.1, 2.0*h_pulse->Integral() );
        int status = h_pulse->Fit(funcname,"RQN","",t_min, t_max);

        // Refit two more times with the latest version of the parameters
        double par[5];
        for(int k=0; k<2; k++){
          for(int j=0; j<4; j++){
            par[j] = func->GetParameter(j);
            func->SetParameter(j,par[j]);
          }
          status = h_pulse->Fit(funcname,"RQN","",t_min, t_max);
        }

        //If the fit status is ok we calculated the rising time as the time
        // when the fitted function has value 10% of its max
        double first_spe_time = -999;
        double dt=-999;

        // TODO -> this should be transformed into a function
        if( status ==0 ) {
	  
	        double t_start = m_startbin*m_sampling_period;
	        double t_end = (m_startbin+m_nbins)*m_sampling_period;

          int npoints = 5000; // should grant decent resolution
          dt = (t_end-t_start)/npoints;
          double max=0.0;
          for(int i=0; i<npoints; i++){
            double t=t_start + dt*i;
            if(max < func->Eval(t) ){ max = func->Eval(t); };
          }

          double startval = 0.1 * max;
          for(int i=0; i<npoints; i++){
            double t=t_start + dt*i;
            if(startval < func->Eval(t) ){ first_spe_time = t; break; };
          }
        }

        // Save the fit paramteters to the pulse object
        temp_pulse.fit_start_time = first_spe_time;
        temp_pulse.error_start_time = dt; // TODO: need a correct error propagation
        temp_pulse.fit_mu = func->GetParameter(1);
        temp_pulse.error_mu = func->GetParError(1);
        temp_pulse.fit_sigma = func->GetParameter(2);
        temp_pulse.error_sigma = func->GetParError(2);
        temp_pulse.fit_amplitude = func->GetParameter(3);
        temp_pulse.error_amplitude = func->GetParError(3);
        temp_pulse.chi2 = func->GetChisquare();
        temp_pulse.ndf = func->GetNDF();
        temp_pulse.fitstatus = status;

        delete h_pulse;

      }

      return temp_pulse;

    }


      //------------------------------------------------------------------------

      Waveform::Pulse Waveform::getIntegral() {

        /*
        Just look in the integration window and extract some basic info
        */

        double charge=0, t_peak=0, amp=0;
        for( int t=m_startbin; t<m_startbin+m_nbins; t++ ) {

          charge += m_waveform.at(t);

          if( m_waveform.at(t) > amp ) {
            amp = m_waveform.at(t);
            t_peak = t*m_sampling_period;
          }

        }


        Pulse temp_pulse;

        temp_pulse.time_peak = t_peak;
        temp_pulse.amplitude = amp * adc_to_mV ; // in mV
        temp_pulse.integral = charge * adc_to_pC;

        return temp_pulse;

      }



      //------------------------------------------------------------------------

      float Waveform::getTotalCharge(){

        float sum=0.0;
        for( auto _adc : m_waveform ){
          sum += _adc;
        }

        return sum * adc_to_pC;

      }

    //--------------------------------------------------------------------------

    TH1D *Waveform::makeWaveformHist(){

      double tstart = m_startbin*m_sampling_period;
      double tend = (m_startbin+m_nbins)*m_sampling_period;
      TH1D *m_wave = new TH1D( "","", tstart, tend, m_nbins );

      for( int bin=m_startbin; bin<m_startbin+m_nbins; bin++ ){
        m_wave->Fill( bin*m_sampling_period, m_waveform[bin] );
      }

      return m_wave;

    }

    //--------------------------------------------------------------------------


    void Waveform::clean()
    {

      m_raw_waveform.clear();
      m_waveform.clear();

      m_baseline_mean=0;
      m_baseline_width=0;

    }


} // end namespace

#endif
