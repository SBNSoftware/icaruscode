////////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingICARUS.h
///
/// \brief  Generic class for shaping signals on wires.
///
/// \author H. Greenlee
///
/// This is a generic class for shaping signals on wires during simulation
/// (convolution) and reconstruction (deconvolution).
///
/// This class acts as a repository for a consistent set of convolution
/// and deconvolution kernels.  It also supplies an interface for
/// convoluting either type of kernel with a time series signal.  All
/// FFT type calculations are done using LArFFT service.
///
/// This class has only a default constructor.  Configuration must be done
/// externally by calling configuration methods.  The proper method for
/// configuring this class is as follows.
///
/// 1.  Add one or more response functions using method AddReponseFunction.
/// 2.  Optionally call methods SetPeakResponseTime or ShiftResponseTime.
/// 3.  Add one or more filter functions using method AddFilterFunction.
/// 4.  Call method CalculateDeconvKernel once.
///
/// After the deconvolution kernel is calculated, the configuration is locked.
///
/// Notes on time and frequency series functions
/// ---------------------------------------------
///
/// Times and frequencies are measured in units of ticks and cycles/tick.
///
/// Time series are represented as vector<double> of length N, representing
/// sampled times on interval [0,N) ticks.   (N = LArFFT::FFTSize().)
///
/// Frequency series are represented as vector<TComplex> of length (N/2+1),
/// representing sampled frequencies on interval [0, 1/2] cycles/tick.
/// Negative frequencies (not stored) are complex conjugate of
/// corresponding positive frequency.
///
/// Update notes
/// -------------
///
/// * Yun-Tse Tsai (yuntse@slac.stanford.edu), July 17th, 2014<br/>
///     Modify
///     `void AddResponseFunction(const std::vector<double>& resp);`
///     to
///     `void AddResponseFunction(const std::vector<double>& resp, bool ResetResponse = false );`
///     If you want to reset your response, `fResponse` in this object, you can
///     do
///     `AddResponseFunction( yourResponse, true )`
///     The other part involving `AddResponseFunction` shouldn't be affected.
/// * X. Qian 2015/01/06 <br/>
///     Add the time offset variable<br/>
///     Need to add the set and extraction code
////////////////////////////////////////////////////////////////////////

#ifndef SignalShapingICARUS_H
#define SignalShapingICARUS_H

#include <vector>
#include <complex>
#include <cmath>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "icaruscode/Utilities/tools/SignalProcessingDefs.h"
#include "icaruscode/Utilities/ICARUSFFT.h"

namespace icarusutil {

class SignalShapingICARUS 
{
public:

    // Constructor, destructor.
    SignalShapingICARUS();

    virtual ~SignalShapingICARUS();

    // Accessors.
    const TimeVec&      Response()      const {return fResponse;}
    const TimeVec&      Response_save() const {return fResponse_save;}
    const FrequencyVec& ConvKernel()    const {return fConvKernel;}
    const FrequencyVec& Filter()        const {return fFilter;}
    const FrequencyVec& DeconvKernel()  const {return fDeconvKernel;}

    // Signal shaping methods.
    // Convolute a time series with convolution kernel.
    template <class T> void Convolute(std::vector<T>& func) const;

    // Convolute a time series with deconvolution kernel.
    template <class T> void Deconvolute(std::vector<T>& func) const;

    // Configuration methods.
    // Reset this class to default-constructed state.
    void Reset();

    void save_response()         {fResponse_save.clear(); fResponse_save=fResponse;}
    void set_normflag(bool flag) {fNorm = flag;}

    // Add a time domain response function.
    // Updates overall response function and convolution kernel.
    void AddResponseFunction(const TimeVec& resp, bool ResetResponse = false );

    // Add a filter function.
    void AddFilterFunction(const FrequencyVec& filt);

    //Add DeconvKernel Polarity switch to decide how to normalize
    //deconvoluted signal w.r.t. RawDigits. If +1 then normalize
    //to Max ADC, if -1 to Min ADC
    void SetDeconvKernelPolarity(int pol);

    // Test and lock the current response function.
    // Does not lock filter configuration.
    void LockResponse() const;

    // Calculate deconvolution kernel using current convolution kernel
    // and filter function.
    // Fully locks configuration.
    void CalculateDeconvKernel() const;

  private:

    // Lock flags.
    mutable bool         fResponseLocked;
    mutable bool         fFilterLocked;

    // Overall response.
    TimeVec              fResponse;
    TimeVec              fResponse_save;

    // Convolution kernel (fourier transform of response function).
    FrequencyVec         fConvKernel;

    // Overall filter function.
    FrequencyVec         fFilter;

    // Deconvolution kernel (= fFilter / fConvKernel).
    mutable FrequencyVec fDeconvKernel;

    // Deconvolution Kernel Polarity Flag
    // Set to +1 if deconv signal should be deconv to + ADC count
    // Set to -1 if one wants to normalize to - ADC count
    int                  fDeconvKernelPolarity;

    // Xin added */
    bool                 fNorm;
};

//----------------------------------------------------------------------
// Convolute a time series with current response.
template <class T> inline void SignalShapingICARUS::Convolute(std::vector<T>& func) const
{
    // Make sure that time series has the correct size.
    if (func.size() != fResponse.size())
    {   
        std::cout << "Convolute, input func size: " << func.size() << ", response size: " << fResponse.size() << std::endl;
        throw cet::exception("SignalShapingICARUS") << "Bad time series size = " << func.size() << "\n";
    }

    // Make sure response configuration is locked.
    if (!fResponseLocked) LockResponse();

    // Get an instance of the FFT
    icarusutil::ICARUSFFT<icarusutil::SigProcPrecision> fft(func.size());

    // Temporary work space
    FrequencyVec funcFFT(func.size());

    // Get FFT of input function
    fft.forwardFFT(funcFFT, func);

    // Convolute with the convolution kernel
    for(size_t idx = 0; idx < fConvKernel.size(); idx++)
        funcFFT[idx] *= fConvKernel[idx];

    // transform back to time domain
    fft.inverseFFT(func, funcFFT);

    return;
}

//----------------------------------------------------------------------
// Convolute a time series with deconvolution kernel.
template <class T> inline void SignalShapingICARUS::Deconvolute(std::vector<T>& func) const
{
    // Make sure that time series has the correct size.
    if(func.size() != fDeconvKernel.size())
    {   
        std::cout << "Convolute, input func size: " << func.size() << ", response size: " << fResponse.size() << std::endl;
        throw cet::exception("SignalShapingICARUS") << "Bad time series size = " << func.size() << "\n";
    }

    // Make sure deconvolution kernel is configured.
    if(!fFilterLocked) CalculateDeconvKernel();

    // Get an instance of the FFT
    icarusutil::ICARUSFFT<icarusutil::SigProcPrecision> fft(func.size());

    // Temporary work space
    FrequencyVec funcFFT(func.size());

    // Get FFT of input function
    fft.forwardFFT(func, funcFFT);

    // Convolute with the deconvolution kernel
    for(size_t idx = 0; idx < fDeconvKernel.size(); idx++)
        funcFFT[idx] *= fDeconvKernel[idx];

    // transform back to time domain
    fft.inverseFFT(funcFFT, func);

    return;
}

}

#endif
