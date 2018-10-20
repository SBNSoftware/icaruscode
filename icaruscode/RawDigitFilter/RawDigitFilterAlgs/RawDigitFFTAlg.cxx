
#include "RawDigitFFTAlg.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "icaruscode/Utilities/tools/IWaveformTool.h"
#include "icaruscode/Utilities/tools/IFilter.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include <cmath>
#include <algorithm>
#include <complex>

#include "TComplex.h"

namespace caldata
{

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitFFTAlg::RawDigitFFTAlg(fhicl::ParameterSet const & pset)
{
    reconfigure(pset);

    // Report.
    mf::LogInfo("RawDigitFFTAlg") << "RawDigitFFTAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitFFTAlg::~RawDigitFFTAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitFFTAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fTransformViewVec = pset.get<std::vector<bool>>  ("TransformViewVec", std::vector<bool>() = {true,false,false});
    fFillHistograms   = pset.get<bool             >  ("FillHistograms",                                      false);
    fHistDirName      = pset.get<std::string      >  ("HistDirName",                                   "FFT_hists");
    fLoWireByPlane    = pset.get<std::vector<size_t>>("LoWireByPlane",               std::vector<size_t>()={0,0,0});
    fHiWireByPlane    = pset.get<std::vector<size_t>>("HiWireByPlane",         std::vector<size_t>()={100,100,100});

    const fhicl::ParameterSet& waveformParamSet = pset.get<fhicl::ParameterSet>("WaveformTool");
    
    fWaveformTool     = art::make_tool<icarus_tool::IWaveformTool>(waveformParamSet);
    
    // Implement the tools for handling the responses
    const fhicl::ParameterSet& filterTools = pset.get<fhicl::ParameterSet>("FilterTools");
    
    for(const std::string& filterTool : filterTools.get_pset_names())
    {
        const fhicl::ParameterSet& filterToolParamSet = filterTools.get<fhicl::ParameterSet>(filterTool);
        size_t                     planeIdx           = filterToolParamSet.get<size_t>("Plane");
        
        fFilterToolMap.insert(std::pair<size_t,std::unique_ptr<icarus_tool::IFilter>>(planeIdx,art::make_tool<icarus_tool::IFilter>(filterToolParamSet)));
    }
}
    
//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitFFTAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fFillHistograms)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        
        // hijack hists here
        double sampleRate  = fDetectorProperties->SamplingRate();
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        double maxFreq     = 1.e6 / (2. * sampleRate);
        double minFreq     = 1.e6 / (2. * sampleRate * readOutSize);
        int    numSamples  = readOutSize / 2;

        // Make a directory for these histograms
        art::TFileDirectory dir = tfs->mkdir(fHistDirName.c_str());

        fFFTPowerVec.resize(3);
        fAveFFTPowerVec.resize(3);
        fConvFFTPowerVec.resize(3);
        fFilterFuncVec.resize(3);
        
        for(size_t plane = 0; plane < 3; plane++)
        {
            size_t numHists = fHiWireByPlane[plane] - fLoWireByPlane[plane];
            
            fFFTPowerVec[plane].resize(numHists);
           
            for(size_t idx = 0; idx < fHiWireByPlane[plane] - fLoWireByPlane[plane]; idx++)
            {
                std::string histName = "FFTPower_" + std::to_string(plane) + "-" + std::to_string(idx);
                
                fFFTPowerVec[plane][idx] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, minFreq, maxFreq, 0., 10000.);
            }
            
            std::string histName = "AveFFTPower_" + std::to_string(plane);
            
            fAveFFTPowerVec[plane] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, minFreq, maxFreq, 0., 1000.);
            
            histName = "ConvFFTPower_" + std::to_string(plane);
            
            fConvFFTPowerVec[plane] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, minFreq, maxFreq, 0., 1000.);
            
            histName = "FilterFunc_" + std::to_string(plane);
            
            fFilterFuncVec[plane] = dir.make<TProfile>(histName.c_str(),  "Filter Function;kHz;Power", numSamples, minFreq, maxFreq, 0., 1000.);
        }
    }
    
    return;
}
    
template <class T> void RawDigitFFTAlg::getFFTCorrection(std::vector<T>& corValVec, double minPowerThreshold) const
{
    // This version will take FFT of input waveform and then remove bins in the time domain with a power less
    // than the threshold input above.
    int fftDataSize = corValVec.size();
    
    Eigen::FFT<T> eigenFFT;
    
    std::vector<std::complex<T>> fftOutputVec(corValVec.size());
    
    eigenFFT.fwd(fftOutputVec, corValVec);
    
    size_t halfFFTDataSize(fftDataSize/2);
    
    std::vector<double> powerVec(halfFFTDataSize);
    
    std::transform(fftOutputVec.begin(), fftOutputVec.begin() + halfFFTDataSize, powerVec.begin(), [](const auto& val){return std::abs(val);});
    
    // Third step is to zap those bins under threshold
    for(size_t idx = 0; idx < halfFFTDataSize; idx++)
    {
        if (powerVec[idx] < minPowerThreshold)
        {
            fftOutputVec.at(idx)                   = std::complex<T>(0.,0.);
            fftOutputVec.at(fftDataSize - idx - 1) = std::complex<T>(0.,0.);
        }
    }
    
    eigenFFT.inv(corValVec, fftOutputVec);
    
    return;
}
    
template void RawDigitFFTAlg::getFFTCorrection<float>(std::vector<float>&, double) const;

template<class T> void RawDigitFFTAlg::getFFTCorrection(std::vector<T>& corValVec, size_t maxBin) const
{
    // This version will take FFT of input waveform and then remove bins in the time domain above the
    // cutoff frequency defined by maxBin passed in above
    int fftDataSize = corValVec.size();
    
    Eigen::FFT<T> eigenFFT;
    
    std::vector<std::complex<T>> fftOutputVec(corValVec.size());
    
    eigenFFT.fwd(fftOutputVec, corValVec);
    
    size_t halfFFTDataSize(fftDataSize/2);

    std::fill(fftOutputVec.begin() + maxBin,     fftOutputVec.begin() + halfFFTDataSize, std::complex<T>(0.,0.));
    std::fill(fftOutputVec.end()   - maxBin - 1, fftOutputVec.end(),                     std::complex<T>(0.,0.));

    eigenFFT.inv(corValVec, fftOutputVec);

    return;
}
    
void RawDigitFFTAlg::filterFFT(std::vector<short>& rawadc, size_t plane, size_t wire, float pedestal) const
{
    // Check there is something to do
    if (!fTransformViewVec.at(plane)) return;
    
    // Step one is to setup and then get the FFT transform of the input waveform
    int    fftDataSize = rawadc.size();
    double sampleRate  = fDetectorProperties->SamplingRate();
    double readOutSize = fDetectorProperties->ReadOutWindowSize();
    
    std::vector<float> fftInputVec;
    
    fftInputVec.resize(fftDataSize, 0.);
    
    std::transform(rawadc.begin(),rawadc.end(),fftInputVec.begin(),[pedestal](const auto& val){return float(float(val) - pedestal);});

    std::vector<std::complex<float>> fftOutputVec;
    
    fftOutputVec.resize(fftInputVec.size());
    
    Eigen::FFT<float> eigenFFT;

    eigenFFT.fwd(fftOutputVec,fftInputVec);
    
    size_t halfFFTDataSize(fftDataSize/2 + 1);

    std::vector<float> powerVec;
    powerVec.resize(halfFFTDataSize, 0.);
    
    std::transform(fftOutputVec.begin(), fftOutputVec.begin() + halfFFTDataSize, powerVec.begin(), [](const auto& complex){return std::abs(complex);});

    // Recover the filter function we are using...
    const std::vector<TComplex>& filter = fFilterToolMap.at(plane)->getResponseVec();
    
    // Make sure the filter has been correctly initialized
    if (filter.size() != halfFFTDataSize) fFilterToolMap.at(plane)->setResponse(fftDataSize,1.,1.);
    
    // Filter the FFT output
    // Ok, this is a pain...
    std::vector<std::complex<float>> filterVec;
    filterVec.reserve(filterVec.size());
    for(auto& rootComplex : filter) filterVec.emplace_back(rootComplex.Re(),rootComplex.Im());
    
    std::transform(fftOutputVec.begin(), fftOutputVec.begin() + fftOutputVec.size()/2, filterVec.begin(), fftOutputVec.begin(), std::multiplies<std::complex<float>>());
    
    for(size_t idx = 0; idx < fftOutputVec.size()/2; idx++) fftOutputVec.at(fftOutputVec.size() - idx - 1) = fftOutputVec.at(idx);

    eigenFFT.inv(fftInputVec, fftOutputVec);

    // Fill hists
    if (fFillHistograms)
    {
        // Fill any individual wire histograms we want to look at
        if (wire >= fLoWireByPlane[plane] && wire < fHiWireByPlane[plane])
        {
            // Fill the power spectrum histogram
            for(size_t idx = 0; idx < halfFFTDataSize; idx++)
            {
                float freq = 1.e6 * float(idx + 1)/ (sampleRate * readOutSize);
                fFFTPowerVec[plane][wire-fLoWireByPlane[plane]]->Fill(freq, std::min(powerVec[idx],float(999.)), 1.);
            }
        }
        
        for(size_t idx = 0; idx < halfFFTDataSize; idx++)
        {
            float freq = 1.e6 * float(idx + 1)/ (sampleRate * readOutSize);
            fAveFFTPowerVec[plane]->Fill(freq, std::min(powerVec[idx],float(999.)), 1.);
        }
        
        // Get the filter power vec
        std::transform(fftOutputVec.begin(), fftOutputVec.begin() + halfFFTDataSize, powerVec.begin(), [](const auto& val){return std::abs(val);});

        for(size_t idx = 0; idx < halfFFTDataSize; idx++)
        {
            float freq = 1.e6 * float(idx)/ (sampleRate * readOutSize);
            fConvFFTPowerVec[plane]->Fill(freq, std::min(powerVec[idx],float(999.)), 1.);
            fFilterFuncVec[plane]->Fill(freq, filter[idx], 1.);
        }
    }
    
    // Finally, we invert the resulting time domain values to recover the new waveform
    std::transform(fftInputVec.begin(), fftInputVec.end(), rawadc.begin(), [pedestal](const float& adc){return std::round(adc + pedestal);});

    return;
}

template void RawDigitFFTAlg::getFFTCorrection<float>(std::vector<float>& corValVec, size_t maxBin) const;
}
