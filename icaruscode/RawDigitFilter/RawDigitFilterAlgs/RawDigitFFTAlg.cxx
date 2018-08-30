
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

#include <cmath>
#include <algorithm>

#include "TVirtualFFT.h"
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
    
    TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
    
    // In the first step we copy the input array into a container of doubles to pass to the FFT
    std::vector<double> fftInputVec(fftDataSize);
    
    std::copy(corValVec.begin(),corValVec.end(),fftInputVec.begin());
    
    fftr2c->SetPoints(fftInputVec.data());
    fftr2c->Transform();
    
    // In the second step we recover the power spectrum
    std::vector<double> realVals(fftDataSize);
    std::vector<double> imaginaryVals(fftDataSize);
    
    fftr2c->GetPointsComplex(realVals.data(), imaginaryVals.data());
    
    size_t halfFFTDataSize(fftDataSize/2);
    
    std::vector<double> powerVec(halfFFTDataSize);
    
    std::transform(realVals.begin(), realVals.begin() + halfFFTDataSize, imaginaryVals.begin(), powerVec.begin(), [](const double& real, const double& imaginary){return std::sqrt(real*real + imaginary*imaginary);});
    
    // Third step is to zap those bins under threshold
    for(size_t idx = 0; idx < halfFFTDataSize; idx++)
    {
        if (powerVec[idx] < minPowerThreshold)
        {
            realVals[idx]                  = 0.;
            realVals[fftDataSize-idx]      = 0.;
            imaginaryVals[idx]             = 0.;
            imaginaryVals[fftDataSize-idx] = 0.;
        }
    }
    
    // Finally, we invert the resulting time domain values to recover the new waveform
    TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &fftDataSize, "C2R M K");
    
    fftc2r->SetPointsComplex(realVals.data(),imaginaryVals.data());
    fftc2r->Transform();
    
    double* fftOutputArray = fftc2r->GetPointsReal();
    
    double normFctr = 1. / double(fftDataSize);
    
    std::transform(fftOutputArray, fftOutputArray + fftDataSize, corValVec.begin(), [normFctr](const double& real){return real * normFctr;});
    
    delete fftc2r;
    delete fftr2c;
    
    return;
}
    
template void RawDigitFFTAlg::getFFTCorrection<float>(std::vector<float>&, double) const;

template<class T> void RawDigitFFTAlg::getFFTCorrection(std::vector<T>& corValVec, size_t maxBin) const
{
    // This version will take FFT of input waveform and then remove bins in the time domain above the
    // cutoff frequency defined by maxBin passed in above
    int fftDataSize = corValVec.size();
    
    TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
    
    // In the first step we copy the input array into a container of doubles to pass to the FFT
    std::vector<double> fftInputVec(fftDataSize);
    
    std::copy(corValVec.begin(),corValVec.end(),fftInputVec.begin());
    
    fftr2c->SetPoints(fftInputVec.data());
    fftr2c->Transform();
    
    // In the second step we recover the power spectrum
    std::vector<double> realVals(fftDataSize);
    std::vector<double> imaginaryVals(fftDataSize);
    
    fftr2c->GetPointsComplex(realVals.data(), imaginaryVals.data());
    
    size_t halfFFTDataSize(fftDataSize/2);
    
    std::vector<double> powerVec(halfFFTDataSize);
    
    std::transform(realVals.begin(), realVals.begin() + halfFFTDataSize, imaginaryVals.begin(), powerVec.begin(), [](const double& real, const double& imaginary){return std::sqrt(real*real + imaginary*imaginary);});
    
    // Zero all bins above selected frequency
    std::fill(realVals.begin()      + maxBin, realVals.begin()      + fftDataSize - maxBin, 0.);
    std::fill(imaginaryVals.begin() + maxBin, imaginaryVals.begin() + fftDataSize - maxBin, 0.);
    
    // Finally, we invert the resulting time domain values to recover the new waveform
    TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &fftDataSize, "C2R M K");
    
    fftc2r->SetPointsComplex(realVals.data(),imaginaryVals.data());
    fftc2r->Transform();
    
    double* fftOutputArray = fftc2r->GetPointsReal();
    
    double normFctr = 1. / double(fftDataSize);
    
    std::transform(fftOutputArray, fftOutputArray + fftDataSize, corValVec.begin(), [normFctr](const double& real){return real * normFctr;});
    
    delete fftc2r;
    delete fftr2c;
    
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

    TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C M");
    
    std::vector<double> fftInputVec;
    
    fftInputVec.resize(fftDataSize, 0.);
    
    std::transform(rawadc.begin(),rawadc.end(),fftInputVec.begin(),[pedestal](const auto& val){return double(float(val) - pedestal);});
    
    fftr2c->SetPoints(fftInputVec.data());
    fftr2c->Transform();
    
    // Now we set up and recover the FFT power spectrum
    std::vector<double>   realVals;
    std::vector<double>   imaginaryVals;
    std::vector<TComplex> complexVals;
    
    size_t halfFFTDataSize(fftDataSize/2 + 1);
    
    realVals.resize(halfFFTDataSize,0.);
    imaginaryVals.resize(halfFFTDataSize,0.);
    
    fftr2c->GetPointsComplex(realVals.data(), imaginaryVals.data());
    
    std::vector<double> powerVec;
    powerVec.resize(halfFFTDataSize, 0.);
    
    std::transform(realVals.begin(), realVals.begin() + halfFFTDataSize, imaginaryVals.begin(), powerVec.begin(), [](const double& real, const double& imaginary){return std::sqrt(real*real + imaginary*imaginary);});
    
    // Not sure the better way to do this...
    for(size_t complexIdx = 0; complexIdx < halfFFTDataSize; complexIdx++) complexVals.emplace_back(realVals.at(complexIdx),imaginaryVals.at(complexIdx));

    // Recover the filter function we are using...
    const std::vector<TComplex>& filter = fFilterToolMap.at(plane)->getResponseVec();
    
    // Make sure the filter has been correctly initialized
    if (filter.size() != halfFFTDataSize) fFilterToolMap.at(plane)->setResponse(fftDataSize,1.,1.);
    
    // Convolve this with the FFT of the input waveform
    std::transform(complexVals.begin(), complexVals.end(), filter.begin(), complexVals.begin(), std::multiplies<TComplex>());
    std::transform(complexVals.begin(), complexVals.end(), realVals.begin(),      [](const auto& val){return val.Re();});
    std::transform(complexVals.begin(), complexVals.end(), imaginaryVals.begin(), [](const auto& val){return val.Im();});

    // Fill hists
    if (fFillHistograms)
    {
        // Fill any individual wire histograms we want to look at
        if (wire >= fLoWireByPlane[plane] && wire < fHiWireByPlane[plane])
        {
            // Fill the power spectrum histogram
            for(size_t idx = 0; idx < halfFFTDataSize; idx++)
            {
                double freq = 1.e6 * double(idx + 1)/ (sampleRate * readOutSize);
                fFFTPowerVec[plane][wire-fLoWireByPlane[plane]]->Fill(freq, std::min(powerVec[idx],999.), 1.);
            }
        }
        
        for(size_t idx = 0; idx < halfFFTDataSize; idx++)
        {
            double freq = 1.e6 * double(idx + 1)/ (sampleRate * readOutSize);
            fAveFFTPowerVec[plane]->Fill(freq, std::min(powerVec[idx],999.), 1.);
        }
        
        // Get the filter power vec
        std::transform(complexVals.begin(), complexVals.end(), powerVec.begin(), [](const auto& val){return val.Rho();});

        for(size_t idx = 0; idx < halfFFTDataSize; idx++)
        {
            double freq = 1.e6 * double(idx)/ (sampleRate * readOutSize);
            fConvFFTPowerVec[plane]->Fill(freq, std::min(powerVec[idx],999.), 1.);
            fFilterFuncVec[plane]->Fill(freq, filter[idx], 1.);
        }
    }
    
    // Finally, we invert the resulting time domain values to recover the new waveform
    TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &fftDataSize, "C2R M");
    
    fftc2r->SetPointsComplex(realVals.data(),imaginaryVals.data());
    fftc2r->Transform();
    
    double* fftOutputArray = fftc2r->GetPointsReal();
    
    double normFctr = 1. / double(fftDataSize);
    
    std::transform(fftOutputArray, fftOutputArray + fftDataSize, rawadc.begin(), [normFctr,pedestal](const double& real){return std::round(real * normFctr + pedestal);});
    
    return;
}

template void RawDigitFFTAlg::getFFTCorrection<float>(std::vector<float>& corValVec, size_t maxBin) const;
/*
void RawDigitFFTAlg::triangleSmooth(std::vector<double>& smoothVec, size_t lowestBin) const
{
    std::vector<double>::iterator curItr  = smoothVec.begin() + 2 + lowestBin;
    std::vector<double>::iterator stopItr = smoothVec.end()   - 2;
    
    while(curItr != stopItr)
    {
        double newVal = (*(curItr - 2) + 2. * *(curItr - 1) + 3. * *curItr + 2. * *(curItr + 1) + *(curItr + 2)) / 9.;
        
        *curItr++ = newVal;
    }
    return;
}
    
void RawDigitFFTAlg::firstDerivative(std::vector<double>& inputVec, std::vector<double>& derivVec) const
{
    derivVec.resize(inputVec.size(), 0.);
    
    for(size_t idx = 1; idx < derivVec.size() - 1; idx++)
        derivVec.at(idx) = 0.5 * (inputVec.at(idx + 1) - inputVec.at(idx - 1));
    
    return;
}
    
void RawDigitFFTAlg::findPeaks(std::vector<double>& derivVec, PeakTupleVec& peakTupleVec, double threshold, size_t startBin) const
{
    peakTupleVec.clear();
    
    for(size_t currentBin = startBin; currentBin < derivVec.size(); currentBin++)
    {
        // Look for bin over threshold
        if (derivVec.at(currentBin) > threshold)
        {
            // search backward to find zero point
            size_t startBin = currentBin;
            
            while(startBin--)
                if (derivVec.at(startBin) < 0.) break;
            
            // Find the peak bin
            size_t negThreshBin = currentBin;
            
            while(++negThreshBin < derivVec.size())
                if (derivVec.at(negThreshBin) < -threshold) break;
            
            if (negThreshBin >= derivVec.size()) negThreshBin--;
            
            size_t peakBin = currentBin + (negThreshBin - currentBin) / 2;
            
            currentBin = negThreshBin;
            
            // Now find the end of the candidate peak
            while(++currentBin < derivVec.size())
                if (derivVec.at(currentBin) > 0.) break;
            
            if (currentBin >= derivVec.size()) currentBin--;
            
            peakTupleVec.push_back(PeakTuple(startBin,peakBin,currentBin));
        }
    }
    
    return;
}
*/    
}
