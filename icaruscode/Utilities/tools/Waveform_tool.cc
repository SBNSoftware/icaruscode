////////////////////////////////////////////////////////////////////////
/// \file   Waveform.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IWaveformTool.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/SignalShaping.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "art/Utilities/make_tool.h"

#include "TVirtualFFT.h"

#include <fstream>
#include <iomanip>

namespace icarus_tool
{

class Waveform : IWaveformTool
{
public:
    explicit Waveform(const fhicl::ParameterSet& pset);
    
    ~Waveform() {}
    
    void configure(const fhicl::ParameterSet& pset)  override;
    
    using PeakTuple    = std::tuple<size_t,size_t,size_t>; // first bin, peak bin, last bin
    using PeakTupleVec = std::vector<PeakTuple>;
    
    void triangleSmooth(const std::vector<float>&,  std::vector<float>&,  size_t = 0)                           const override;
    void triangleSmooth(const std::vector<double>&, std::vector<double>&, size_t = 0)                           const override;
    void medianSmooth(  const std::vector<float>&,  std::vector<float>&,  size_t = 3)                           const override;
    void medianSmooth(  const std::vector<double>&, std::vector<double>&, size_t = 3)                           const override;
    void getTruncatedMeanRMS(const std::vector<double>&, double&, double&)                                      const override;
    void getTruncatedMeanRMS(const std::vector<float>&, float&, float&)                                         const override;
    void firstDerivative(const std::vector<float>&,  std::vector<float>&)                                       const override;
    void firstDerivative(const std::vector<double>&, std::vector<double>&)                                      const override;
    void findPeaks(std::vector<float>::iterator,  std::vector<float>::iterator,  PeakTupleVec&, float, size_t)  const override;
    void findPeaks(std::vector<double>::iterator, std::vector<double>::iterator, PeakTupleVec&, double, size_t) const override;
    void getFFTPower(const std::vector<float>& inputVec, std::vector<float>& outputPowerVec)                    const override;
    void getFFTPower(const std::vector<double>& inputVec, std::vector<double>& outputPowerVec)                  const override;
    
private:
    template <typename T> void triangleSmooth(const std::vector<T>&, std::vector<T>&, size_t = 0)                                        const;
    template <typename T> void medianSmooth(  const std::vector<T>&, std::vector<T>&, size_t = 3)                                        const;
    template <typename T> void getTruncatedMeanRMS(const std::vector<T>&, T&, T&)                                                        const;
    template <typename T> void firstDerivative(const std::vector<T>&,  std::vector<T>&)                                                  const;
    template <typename T> void findPeaks(typename std::vector<T>::iterator, typename std::vector<T>::iterator, PeakTupleVec&, T, size_t) const;
};
    
//----------------------------------------------------------------------
// Constructor.
Waveform::Waveform(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
void Waveform::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
//    fThisPlane       = pset.get<size_t>("Plane");
    
    return;
}
    
void Waveform::triangleSmooth(const std::vector<double>& inputVec, std::vector<double>& smoothVec, size_t lowestBin) const
{
    triangleSmooth<double>(inputVec, smoothVec, lowestBin);
    
    return;
}
    
void Waveform::triangleSmooth(const std::vector<float>& inputVec, std::vector<float>& smoothVec, size_t lowestBin) const
{
    triangleSmooth<float>(inputVec, smoothVec, lowestBin);
    
    return;
}

template <typename T> void Waveform::triangleSmooth(const std::vector<T>& inputVec, std::vector<T>& smoothVec, size_t lowestBin) const
{
    if (inputVec.size() != smoothVec.size()) smoothVec.resize(inputVec.size());
    
    std::copy(inputVec.begin(), inputVec.begin() + 2 + lowestBin, smoothVec.begin());
    std::copy(inputVec.end() - 2, inputVec.end(), smoothVec.end() - 2);
    
    typename std::vector<T>::iterator       curItr    = smoothVec.begin() + 2 + lowestBin;
    typename std::vector<T>::const_iterator curInItr  = inputVec.begin() + 1 + lowestBin;
    typename std::vector<T>::const_iterator stopInItr = inputVec.end()   - 2;
    
    while(curInItr++ != stopInItr)
    {
        // Take the weighted average of five consecutive points centered on current point
        T newVal = (*(curInItr - 2) + 2. * *(curInItr - 1) + 3. * *curInItr + 2. * *(curInItr + 1) + *(curInItr + 2)) / 9.;
        
        *curItr++ = newVal;
    }
    return;
}
    
void Waveform::medianSmooth(const std::vector<float>&inputVec, std::vector<float>& smoothVec, size_t nBins) const
{
    medianSmooth<float>(inputVec, smoothVec, nBins);
    
    return;
}

void Waveform::medianSmooth(const std::vector<double>& inputVec, std::vector<double>& smoothVec, size_t nBins) const
{
    medianSmooth<double>(inputVec, smoothVec, nBins);
    
    return;
}

template <typename T> void Waveform::medianSmooth(const std::vector<T>& inputVec, std::vector<T>& smoothVec, size_t nBins) const
{
    // For our purposes, nBins must be odd
    if (nBins % 2 == 0) nBins++;
    
    // Make sure the input vector is right sized
    if (inputVec.size() != smoothVec.size()) smoothVec.resize(inputVec.size());
    
    // Basic set up
    typename std::vector<T> medianVec(nBins);
    typename std::vector<T>::const_iterator startItr = inputVec.begin();
    typename std::vector<T>::const_iterator stopItr  = startItr;
    
    std::advance(stopItr, inputVec.size() - nBins);
    
    size_t medianBin = nBins/2;
    size_t smoothBin = medianBin;
    
    // First bins are not smoothed
    std::copy(startItr, startItr + medianBin, smoothVec.begin());
    
    while(std::distance(startItr,stopItr) > 0)
    {
        std::copy(startItr,startItr+nBins,medianVec.begin());
        std::sort(medianVec.begin(),medianVec.end());
        
        T medianVal = medianVec[medianBin];
        
        smoothVec[smoothBin++] = medianVal;
        
        startItr++;
    }
    
    // Last bins are not smoothed
    std::copy(startItr + medianBin, inputVec.end(), smoothVec.begin() + smoothBin);
    
    return;
}
    
void Waveform::getTruncatedMeanRMS(const std::vector<double>& waveform, double& mean, double& rms) const
{
    getTruncatedMeanRMS<double>(waveform, mean, rms);
}
    
void Waveform::getTruncatedMeanRMS(const std::vector<float>& waveform, float& mean, float& rms) const
{
    getTruncatedMeanRMS<float>(waveform, mean, rms);
}

template <typename T> void Waveform::getTruncatedMeanRMS(const std::vector<T>& waveform, T& mean, T& rms) const
{
    // We need to get a reliable estimate of the mean and can't assume the input waveform will be ~zero mean...
    // Basic idea is to find the most probable value in the ROI presented to us
    // From that we can develop an average of the true baseline of the ROI.
    // To do that we employ a map based scheme
    std::map<int,int> frequencyMap;
    int               mpCount(0);
    int               mpVal(0);
    
    for(const auto& val : waveform)
    {
        int intVal = std::round(4.*val);
        
        frequencyMap[intVal]++;
        
        if (frequencyMap.at(intVal) > mpCount)
        {
            mpCount = frequencyMap.at(intVal);
            mpVal   = intVal;
        }
    }
    
    // take a weighted average of two neighbor bins
    int meanCnt  = 0;
    int meanSum  = 0;
    int binRange = std::min(16, int(frequencyMap.size()/2 + 1));
    
    for(int idx = -binRange; idx <= binRange; idx++)
    {
        std::map<int,int>::iterator neighborItr = frequencyMap.find(mpVal+idx);
        
        if (neighborItr != frequencyMap.end() && 5 * neighborItr->second > mpCount)
        {
            meanSum += neighborItr->first * neighborItr->second;
            meanCnt += neighborItr->second;
        }
    }
    
    mean = 0.25 * float(meanSum) / float(meanCnt);  // Note that bins were expanded by a factor of 4 above
    
    // do rms calculation - the old fashioned way and over all adc values
    typename std::vector<T> locWaveform = waveform;
    
    std::transform(locWaveform.begin(), locWaveform.end(), locWaveform.begin(),std::bind(std::minus<T>(),std::placeholders::_1,mean));
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});

    // recalculate the rms
    T localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + meanCnt, locWaveform.begin(), 0.);
    rms = std::sqrt(std::max(T(0.),localRMS / T(meanCnt)));
    
    return;
}

void Waveform::firstDerivative(const std::vector<double>& inputVec, std::vector<double>& derivVec) const
{
    firstDerivative<double>(inputVec, derivVec);
    
    return;
}
    
void Waveform::firstDerivative(const std::vector<float>& inputVec, std::vector<float>& derivVec) const
{
    firstDerivative<float>(inputVec, derivVec);
    
    return;
}
    
template <typename T> void Waveform::firstDerivative(const std::vector<T>& inputVec, std::vector<T>& derivVec) const
{
    derivVec.resize(inputVec.size(), 0.);
    
    for(size_t idx = 1; idx < derivVec.size() - 1; idx++)
        derivVec.at(idx) = 0.5 * (inputVec.at(idx + 1) - inputVec.at(idx - 1));
    
    return;
}
    
void Waveform::findPeaks(std::vector<double>::iterator startItr, std::vector<double>::iterator stopItr, PeakTupleVec& peakTupleVec, double threshold, size_t firstTick) const
{
    findPeaks<double>(startItr, stopItr, peakTupleVec, threshold, firstTick);
    
    return;
}
    
void Waveform::findPeaks(std::vector<float>::iterator startItr, std::vector<float>::iterator stopItr, PeakTupleVec& peakTupleVec, float threshold, size_t firstTick) const
{
    findPeaks<float>(startItr, stopItr, peakTupleVec, threshold, firstTick);
    
    return;
}

template <typename T> void Waveform::findPeaks(typename std::vector<T>::iterator startItr,
                                               typename std::vector<T>::iterator stopItr,
                                               PeakTupleVec&                     peakTupleVec,
                                               T                                 threshold,
                                               size_t                            firstTick) const
{
    // Need a minimum distance or else nothing to do
    if (std::distance(startItr,stopItr) > 4)
    {
        // This is a divide and conquer algorithm, start by finding the maximum element.
        typename std::vector<T>::iterator firstItr = std::max_element(startItr,stopItr,[](float left, float right){return std::fabs(left) < std::fabs(right);});

        // Are we over threshold?
        if (std::fabs(*firstItr) > threshold)
        {
            // What am I thinking?
            // First task is to find the "other" lobe max point
            // Set one to the "first", the other to the "second"
            // Search backward from first to find start point, forward from second to find end point
            // Set mid point between first and second as "peak"?
            typename std::vector<T>::iterator secondItr = firstItr;
        
            // Assume if max bin is positive then second lobe is later
            if (*firstItr > 0)
            {
                typename std::vector<T>::iterator tempItr = secondItr;
            
                while(tempItr != stopItr)
                {
                    if (*++tempItr < -threshold)
                    {
                        if (*tempItr < *secondItr) secondItr = tempItr;
                    }
                    else if (secondItr != firstItr) break;
                }
            }
            // Otherwise it goes the other way
            else
            {
                typename std::vector<T>::iterator tempItr = secondItr;
            
                while(tempItr != startItr)
                {
                    if (*--tempItr > threshold)
                    {
                        if (*tempItr > *secondItr) secondItr = tempItr;
                    }
                    else if (secondItr != firstItr) break;
                }
            
                std::swap(firstItr,secondItr);
            }
        
            // It might that no real pulse was found
            if (firstItr != secondItr)
            {
                // Get the "peak" position
                size_t peakBin = std::distance(startItr,firstItr) + std::distance(firstItr,secondItr) / 2;
        
                // Advance (forward or backward) the first and second iterators to get back to zero crossing
                while(firstItr  != startItr) if (*--firstItr  < 0.) break;
                while(secondItr != stopItr)  if (*++secondItr > 0.) break;
        
                size_t firstBin = std::distance(startItr,firstItr);
                size_t lastBin  = std::distance(startItr,secondItr);
        
                // Find leading peaks
                findPeaks(startItr, firstItr, peakTupleVec, threshold, firstTick);
        
                // Save this peak
                peakTupleVec.push_back(PeakTuple(firstBin+firstTick,peakBin+firstTick,lastBin+firstTick));
        
                // Find downstream peaks
                findPeaks(secondItr, stopItr, peakTupleVec, threshold, firstTick + std::distance(startItr,secondItr));
            }
        }
    }

    return;
}
    
void Waveform::getFFTPower(const std::vector<float>& inputVec, std::vector<float>& outputPowerVec) const
{
    std::vector<double> inputDoubleVec(inputVec.size());
    std::vector<double> outputDoubleVec(inputVec.size()/2);
    
    std::copy(inputVec.begin(),inputVec.end(),inputDoubleVec.begin());
    
    getFFTPower(inputDoubleVec, outputDoubleVec);
    
    if (outputDoubleVec.size() != outputPowerVec.size()) outputPowerVec.resize(outputDoubleVec.size());
    
    std::copy(outputDoubleVec.begin(),outputDoubleVec.end(),outputPowerVec.begin());
    
    return;
}

void Waveform::getFFTPower(const std::vector<double>& inputVec, std::vector<double>& outputPowerVec) const
{
    // Get the FFT of the response
    int fftDataSize = inputVec.size();
    
    TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
    
    fftr2c->SetPoints(inputVec.data());
    fftr2c->Transform();
    
    // Recover the results so we can compute the power spectrum
    size_t halfFFTDataSize(fftDataSize/2 + 1);
    
    std::vector<double> realVals(halfFFTDataSize);
    std::vector<double> imaginaryVals(halfFFTDataSize);
    
    fftr2c->GetPointsComplex(realVals.data(), imaginaryVals.data());
    
    if (outputPowerVec.size() != halfFFTDataSize) outputPowerVec.resize(halfFFTDataSize,0.);
    
    std::transform(realVals.begin(), realVals.begin() + halfFFTDataSize, imaginaryVals.begin(), outputPowerVec.begin(), [](const double& real, const double& imaginary){return std::sqrt(real*real + imaginary*imaginary);});
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(Waveform)
}
