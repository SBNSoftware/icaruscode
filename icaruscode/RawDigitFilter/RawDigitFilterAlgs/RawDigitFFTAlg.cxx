
#include "RawDigitFFTAlg.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "icaruscode/Utilities/tools/IWaveformTool.h"

#include <cmath>
#include <algorithm>

#include "TVirtualFFT.h"

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
    fTransformViewVec = pset.get<std::vector<bool>>("TransformViewVec",     std::vector<bool>() = {true,false,false});
    fFillHistograms   = pset.get<bool             >("FillHistograms",                                          false);
    fHistDirName      = pset.get<std::string      >("HistDirName",                                       "FFT_hists");
    
    const fhicl::ParameterSet& waveformParamSet = pset.get<fhicl::ParameterSet>("WaveformTool");
    
    fWaveformTool     = art::make_tool<icarus_tool::IWaveformTool>(waveformParamSet);
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
        //    double sampleRate  = fDetectorProperties->SamplingRate();
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        //    double maxFreq     = 1000000. / (2. * sampleRate);
        //    double minFreq     = 1000000. / (2. * sampleRate * readOutSize);
        //    int    numSamples  = (readOutSize / 2 + 1) / 4;
        int numSamples     = readOutSize / 2;
        
        fCorValHistVec.resize(20);
        fFFTPowerVec.resize(20);
        fFFTPowerDerivVec.resize(20);
        fFFTRealVec.resize(20);
        fFFTImaginaryVec.resize(20);
        fFFTCorValHistVec.resize(20);
        fSmoothPowerVec.resize(20);
        
        // Make a directory for these histograms
        art::TFileDirectory dir = tfs->mkdir(fHistDirName.c_str());
        
        for(size_t idx = 0; idx < 20; idx++)
        {
            std::string histName = "RawWaveform_" + std::to_string(idx);
            
            fCorValHistVec[idx] = dir.make<TProfile>(histName.c_str(), "Raw Waveform;Tick", readOutSize, 0., readOutSize, -100., 100.);
            
            histName = "FFTPower_" + std::to_string(idx);
            
            fFFTPowerVec[idx] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, 0, numSamples, 0., 10000.);
            
            histName = "FFTPowerDeriv_" + std::to_string(idx);
            
            fFFTPowerDerivVec[idx] = dir.make<TProfile>(histName.c_str(),  "Power Deriv;kHz;Power", numSamples, 0, numSamples, -500., 500.);
            
            histName = "FFTReal_" + std::to_string(idx);
            
            fFFTRealVec[idx] = dir.make<TProfile>(histName.c_str(),  "Real values;kHz;Power", numSamples, 0, numSamples, -10000., 10000.);
            
            histName = "FFTImaginary_" + std::to_string(idx);
            
            fFFTImaginaryVec[idx] = dir.make<TProfile>(histName.c_str(),  "Imaginary values;kHz;Power", numSamples, 0, numSamples, -10000., 10000.);
            
            histName = "SmoothPWR_" + std::to_string(idx);
            
            fSmoothPowerVec[idx] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, 0, numSamples, 0., 10000.);
            
            histName = "FFTCorrected_" + std::to_string(idx);
            
            fFFTCorValHistVec[idx] = dir.make<TProfile>(histName.c_str(),  "Corrected Waveform;Tick", readOutSize, 0., readOutSize, -100., 100.);
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
    if (fTransformViewVec[plane])
    {
        size_t lowWire(350);
        size_t hiWire(370);
    
        //                double sampleFreq  = 1000000. / fDetectorProperties->SamplingRate();
        //                double readOutSize = fDetectorProperties->ReadOutWindowSize();
        //                double binSize     = sampleFreq / readOutSize;
        int    fftDataSize = rawadc.size();
    
        TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C M");
    
        std::vector<double> fftInputVec;
        
        fftInputVec.resize(fftDataSize, 0.);
    
        for(size_t tick = 0; tick < rawadc.size(); tick++)
        {
            fftInputVec[tick] = rawadc[tick] - pedestal;
        
            if (fFillHistograms && plane == 0 && wire >= lowWire && wire < hiWire)
                fCorValHistVec[wire-lowWire]->Fill(tick, fftInputVec[tick], 1.);
        }
    
        fftr2c->SetPoints(fftInputVec.data());
        fftr2c->Transform();
    
        // Recover the power spectrum...
        std::vector<double> realVals;
        std::vector<double> imaginaryVals;
        
        size_t halfFFTDataSize(fftDataSize/2 + 1);
        
        realVals.resize(halfFFTDataSize,0.);
        imaginaryVals.resize(halfFFTDataSize,0.);
    
        fftr2c->GetPointsComplex(realVals.data(), imaginaryVals.data());
    
        std::vector<double> powerVec;
        powerVec.resize(halfFFTDataSize, 0.);
    
        std::transform(realVals.begin(), realVals.begin() + halfFFTDataSize, imaginaryVals.begin(), powerVec.begin(), [](const double& real, const double& imaginary){return std::sqrt(real*real + imaginary*imaginary);});
    
        if (fFillHistograms && plane == 0 && wire >= lowWire && wire < hiWire)
        {
            // Fill the power spectrum histogram
            for(size_t idx = 0; idx < halfFFTDataSize; idx++)
                fFFTPowerVec[wire-lowWire]->Fill(idx, std::min(powerVec[idx],9999.), 1.);
        }
    
        if (fTransformViewVec[plane])
        {
            // Idea here is to run through the power spectrum and keep a running average of the n bins around the current bin
            size_t numBinsToAve(9);  // number bins either size of current bin
            size_t lowestBin(3);    //(275);   // Go no lower than this?

            size_t currentBin(halfFFTDataSize - numBinsToAve - 1);
//            size_t firstBin(currentBin - numBinsToAve - 1);
//            size_t lastBin(halfFFTDataSize - 1);
            
            // Initialization of running sum to start one bin past the first, include the "current" and "last"
//            double powerRunSum      = std::accumulate(powerVec.begin()      + firstBin + 1, powerVec.begin()      + lastBin + 1, 0.);
//            double realRunSum       = std::accumulate(realVals.begin()      + firstBin + 1, realVals.begin()      + lastBin + 1, 0.);
//            double imaginaryRunSum  = std::accumulate(imaginaryVals.begin() + firstBin + 1, imaginaryVals.begin() + lastBin + 1, 0.);
            
            fWaveformTool->triangleSmooth(powerVec, powerVec);
//            triangleSmooth(powerVec);
//            triangleSmooth(realVals);
//            triangleSmooth(imaginaryVals);
            
            std::vector<double> powerDerivVec;
            
            fWaveformTool->firstDerivative(powerVec, powerDerivVec);
            
            // Find the peaks...
            icarus_tool::IWaveformTool::PeakTupleVec peakTupleVec;
            
            fWaveformTool->findPeaks(powerDerivVec.begin() + 300, powerDerivVec.end(), peakTupleVec, 10., 0);
            
            // Try smoothing the peak regions
            for(const auto& peakTuple : peakTupleVec)
            {
                if (std::get<0>(peakTuple) >= powerVec.size() || std::get<2>(peakTuple) >= powerVec.size())
                {
                    std::cout << "indexing problem - first: " << std::get<0>(peakTuple) << ", last: " << std::get<2>(peakTuple) << std::endl;
                    continue;
                }
                
                size_t firstBin    = std::get<0>(peakTuple);
                size_t lastBin     = std::get<2>(peakTuple);
                double firstBinVal = powerVec.at(firstBin);
                double lastBinVal  = powerVec.at(lastBin);
                double stepVal     = (lastBinVal - firstBinVal) / double(lastBin - firstBin);
                double newBinVal   = firstBinVal + stepVal;
                
                while(++firstBin < lastBin)
                {
                    // Update the power first
                    powerVec.at(firstBin)  = newBinVal;
                    newBinVal             += stepVal;
                    
                    // Now scale the real and imaginary values...
                    double scaleFactor = 1. / sqrt(realVals.at(firstBin)*realVals.at(firstBin) + imaginaryVals.at(firstBin)*imaginaryVals.at(firstBin));
                    
                    realVals.at(firstBin)      *= powerVec.at(firstBin) * scaleFactor;
                    imaginaryVals.at(firstBin) *= powerVec.at(firstBin) * scaleFactor;
                }
            }
            
            while(currentBin > lowestBin)
            {
/*
                // update running sum
                powerRunSum     += powerVec.at(firstBin);
                powerRunSum     -= powerVec.at(lastBin);
                
                realRunSum      += realVals.at(firstBin);
                realRunSum      -= realVals.at(lastBin);
                
                imaginaryRunSum += imaginaryVals.at(firstBin);
                imaginaryRunSum -= imaginaryVals.at(lastBin);
                
                double runningCnt      = lastBin - firstBin - 1;
                double thisBinValue    = powerVec.at(currentBin);
                double avePowerThisBin = (powerRunSum - thisBinValue) / runningCnt;
                
                if (thisBinValue - avePowerThisBin > 150.)
                {
                    double aveRealThisBin      = (realRunSum      - realVals.at(currentBin))      / runningCnt;
                    double aveImaginaryThisBin = (imaginaryRunSum - imaginaryVals.at(currentBin)) / runningCnt;
                    
                    realRunSum                   -= realVals.at(currentBin) - aveRealThisBin;
                    realVals.at(currentBin)       = aveRealThisBin;
                    
                    imaginaryRunSum              -= imaginaryVals.at(currentBin) - aveImaginaryThisBin;
                    imaginaryVals.at(currentBin)  = aveImaginaryThisBin;
                    
                    powerRunSum                  -= thisBinValue - avePowerThisBin;
                    powerVec.at(currentBin)       = avePowerThisBin;
                }
*/
                double avePowerThisBin(powerVec.at(currentBin));
                
                if (fFillHistograms && plane == 0 && wire >= lowWire && wire < hiWire)
                {
                    fSmoothPowerVec[wire-lowWire]->Fill(currentBin, avePowerThisBin, 1.);
                    fFFTPowerDerivVec[wire-lowWire]->Fill(currentBin, powerDerivVec.at(currentBin), 1.);
                    fFFTRealVec[wire-lowWire]->Fill(currentBin, realVals.at(currentBin), 1.);
                    fFFTImaginaryVec[wire-lowWire]->Fill(currentBin, imaginaryVals.at(currentBin), 1.);
                }
                
                currentBin--;
//                firstBin--;
//                lastBin--;
            }
            
            // Remove the bins at the bottom
            for(size_t idx = 0; idx < numBinsToAve; idx++)
            {
                realVals[idx]      = 0.;
                imaginaryVals[idx] = 0.;
            }
        }
    
        // Finally, we invert the resulting time domain values to recover the new waveform
        TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &fftDataSize, "C2R M");
    
        fftc2r->SetPointsComplex(realVals.data(),imaginaryVals.data());
        fftc2r->Transform();
    
        double* fftOutputArray = fftc2r->GetPointsReal();
    
        double normFctr = 1. / double(fftDataSize);
    
        std::transform(fftOutputArray, fftOutputArray + fftDataSize, rawadc.begin(), [normFctr,pedestal](const double& real){return std::round(real * normFctr + pedestal);});
        
        if (fFillHistograms && plane == 0 && wire >= lowWire && wire < hiWire)
            for(int idx = 0; idx < fftDataSize; idx++) fFFTCorValHistVec[wire-lowWire]->Fill(idx, rawadc[idx] - pedestal, 1.);
    }
    
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
