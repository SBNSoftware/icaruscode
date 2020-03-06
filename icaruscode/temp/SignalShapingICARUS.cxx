//////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingICARUS.cxx
///
/// \brief  Generic signal shaping class.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "cetlib_except/exception.h"
#include "icaruscode/TPC/Utilities/SignalShapingICARUS.h"
#include "icaruscode/TPC/Utilities/ICARUSFFT.h"

namespace icarusutil
{

//----------------------------------------------------------------------
// Constructor.
//
SignalShapingICARUS::SignalShapingICARUS()
  : fResponseLocked(false)
  , fFilterLocked  (false)
  , fNorm (true)
{}

//----------------------------------------------------------------------
// Destructor.
//
SignalShapingICARUS::~SignalShapingICARUS()
{}

//----------------------------------------------------------------------
// Reset this class to its default-constructed state.
void SignalShapingICARUS::Reset()
{
    fResponseLocked = false;
    fFilterLocked = false;
    fResponse.clear();
    fConvKernel.clear();
    fFilter.clear();
    fDeconvKernel.clear();
    //Set deconvolution polarity to + as default
    fDeconvKernelPolarity = +1;
}


//----------------------------------------------------------------------
// Add a time domain response function.
void SignalShapingICARUS::AddResponseFunction(const std::vector<double>& resp, bool ResetResponse )
{
    // Make sure configuration is not locked.

    if(fResponseLocked)
        throw cet::exception("SignalShapingICARUS") << "Configuration locked.\n";

    // Copy new response function into fResponse attribute
    fResponse = resp;

    icarusutil::ICARUSFFT<double> fft(resp.size());
  
    // Is this the first response function?
    if ( fConvKernel.size() == 0 || ResetResponse ) 
    {
        // This is the first response function.
        // Just calculate the fourier transform.
        fConvKernel.resize(resp.size());
        fft.forwardFFT(fResponse,fConvKernel);
    }
    else 
    {
        if (resp.size() != fConvKernel.size())
            throw cet::exception("SignalShapingICARUS") << __func__ << ": inconsistent kernel size, "
              << resp.size() << " vs. " << fConvKernel.size();

        // Not the first response function.
        // Calculate the fourier transform of new response function.
        FrequencyVec kern(fResponse.size());

        fft.forwardFFT(fResponse, kern);

        for(unsigned int i=0; i<kern.size(); ++i) fConvKernel[i] *= kern[i];

        // Recalculate overall response function.
        fft.inverseFFT(fConvKernel, fResponse);
    }

    return;
}


//----------------------------------------------------------------------
// Add a frequency domain filter function to cumulative filter function.
void SignalShapingICARUS::AddFilterFunction(const FrequencyVec& filt)
{
    // Make sure configuration is not locked.
    if(fFilterLocked)
        throw cet::exception("SignalShapingICARUS") << "Configuration locked.\n";

    // If this is the first filter function, just copy the filter function.
    // Otherwise, update the overall filter function.
    if(fFilter.size() == 0) 
    {
        fFilter = filt;
    }
    else 
    {
        unsigned int n = std::min(fFilter.size(), filt.size());

        for(unsigned int i=0; i<n; ++i)              fFilter[i] *= filt[i];
        for(unsigned int i=n; i<fFilter.size(); ++i) fFilter[i] = 0.;
    }

    return;
}

//----------------------------------------------------------------------
// Add a DeconvKernel Polarity Flag to decide how to normalize
void SignalShapingICARUS::SetDeconvKernelPolarity(int pol)
{
    if ( (pol != 1) and (pol != -1) ) 
    {
        throw cet::exception("SignalShapingICARUS") << __func__
            << ": DeconvKernelPolarity should be +1 or -1 (got " << pol << "). Setting to +1\n";
        fDeconvKernelPolarity = +1;
    }
    else
        fDeconvKernelPolarity = pol;

    return;
}


//----------------------------------------------------------------------
// Test and lock the response and convolution kernel.
void SignalShapingICARUS::LockResponse() const
{
    // Do nothing if the response is already locked.
    if(!fResponseLocked) 
    {
        // Make sure response has been configured.
        if(fResponse.size() == 0)
            throw cet::exception("SignalShapingICARUS")
	                << "Response has not been configured.\n";

        // Set the lock flag.
        fResponseLocked = true;
    }

    return;
}


//----------------------------------------------------------------------
// Calculate the deconvolution kernel as the ratio
// of the filter function and convolution kernel.
void SignalShapingICARUS::CalculateDeconvKernel() const
{
    // Make sure configuration is not locked.
    if(fFilterLocked)
        throw cet::exception("SignalShapingICARUS") << "Configuration locked.\n";

    // Lock response configuration.
    LockResponse();

    // Get the FFT function
    icarusutil::ICARUSFFT<double> fft(fResponse.size());

    // Make sure filter function has been configured.
    if (fFilter.size() == 0)
        throw cet::exception("SignalShapingICARUS")
            << "Filter function has not been configured.\n";

    if (fFilter.size() != fResponse.size())
        throw cet::exception("SignalShapingICARUS") << "Filter function does not match. \n";

    // Calculate deconvolution kernel as the ratio of the
    // filter function and the convolution kernel.
    fDeconvKernel = fFilter;
   
    for(unsigned int i=0; i<fDeconvKernel.size(); ++i) 
    {
        if (std::abs(fConvKernel[i]) <= 0.0001) fDeconvKernel[i] = 0.;
        else                                    fDeconvKernel[i] /= fConvKernel[i];
    }

    // Normalize the deconvolution kernel.
    // Calculate the unnormalized deconvoluted response
    // (inverse FFT of filter function).
    icarusutil::TimeVec      deconv(fResponse.size(), 0.);
    icarusutil::FrequencyVec filter = fFilter;

    fft.inverseFFT(filter, deconv);

    if (fNorm)
    {
        // Find the peak value of the response
        // Should normally be at zero, but don't assume that.
        // Use DeconvKernelPolairty to find what to normalize to
        double peak_response = 0;
    
        if ( fDeconvKernelPolarity == -1 ) peak_response = 4096;
    
        for(unsigned int i = 0; i < fResponse.size(); ++i) 
        {
            if( (fResponse[i] > peak_response) and (fDeconvKernelPolarity == 1))          peak_response = fResponse[i];
            else if ( (fResponse[i] < peak_response) and ( fDeconvKernelPolarity == -1) ) peak_response = fResponse[i];
        }

        if ( fDeconvKernelPolarity == -1 ) peak_response *= -1;
    
        if (peak_response <= 0.)
            throw cet::exception("SignalShapingICARUS") << __func__
					    << ": peak should always be positive (got " << peak_response << ")\n";

        // Find the peak value of the deconvoluted response
        // Should normally be at zero, but don't assume that.
        double peak_deconv = 0.;

        for(unsigned int i = 0; i < deconv.size(); ++i) 
        {
            if(deconv[i] > peak_deconv) peak_deconv = deconv[i];
        }
    
        if (peak_deconv <= 0.)
            throw cet::exception("SignalShapingICARUS") << __func__
					    << ": deconvolution peak should always be positive (got " << peak_deconv << ")\n";

        // Multiply the deconvolution kernel by a factor such that
        // (Peak of response) = (Peak of deconvoluted response).
        double ratio = peak_response / peak_deconv;

        for(unsigned int i = 0; i < fDeconvKernel.size(); ++i) fDeconvKernel[i] *= ratio;
    }
  
    // Set the lock flag.
    fFilterLocked = true;

    return;
}

} // end of namespace
