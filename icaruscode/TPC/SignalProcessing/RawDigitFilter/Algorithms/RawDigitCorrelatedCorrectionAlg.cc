
#include <cmath>
#include <algorithm>
#include <vector>

#include "RawDigitCorrelatedCorrectionAlg.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/LArFFTWPlan.h"
#include "lardata/Utilities/LArFFTW.h"

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
RawDigitCorrelatedCorrectionAlg::RawDigitCorrelatedCorrectionAlg(fhicl::ParameterSet const & pset)
{
    reconfigure(pset);

    // Report.
    mf::LogInfo("RawDigitCorrelatedCorrectionAlg") << "RawDigitCorrelatedCorrectionAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitCorrelatedCorrectionAlg::~RawDigitCorrelatedCorrectionAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitCorrelatedCorrectionAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                        0.15);
    fApplyCorSmoothing     = pset.get<bool>               ("ApplyCorSmoothing",                                        true);
    fApplyFFTCorrection    = pset.get<bool>               ("ApplyFFTCorrection",                                       true);
    fFillFFTHistograms     = pset.get<bool>               ("FillFFTHistograms",                                       false);
    fFFTHistsWireGroup     = pset.get<std::vector<size_t>>("FFTHistsWireGroup",         std::vector<size_t>() = {1, 33, 34});
    fFFTNumHists           = pset.get<std::vector<size_t>>("FFTNumWaveHistograms",       std::vector<size_t>() = {10,48,48});
    fFFTHistsStartTick     = pset.get<std::vector<double>>("FFTWaveHistsStartTick", std::vector<double>() = {96.,96.,7670.});
    fFFTMinPowerThreshold  = pset.get<std::vector<double>>("FFTPowerThreshold",     std::vector<double>() = {100.,75.,500.});
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",          std::vector<size_t>() = {48, 48, 96});
    fFillHistograms        = pset.get<bool>               ("FillHistograms",                                          false);
    fRunFFTCorrected       = pset.get<bool>               ("RunFFTCorrectedWires",                                    false);
    fNumRmsToSmoothVec     = pset.get<std::vector<float>> ("NumRmsToSmooth",          std::vector<float>() = {3.6, 3.6, 4.});
}

//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitCorrelatedCorrectionAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs)
{
}

void RawDigitCorrelatedCorrectionAlg::smoothCorrectionVec(std::vector<float>& corValVec, unsigned int& viewIdx) const
{
    // First get the truncated mean and rms for the input vector (noting that it is not in same format as raw data)
    // We need a local copy so we can sort it
    std::vector<float> localCorValVec = corValVec;

    std::sort(localCorValVec.begin(),localCorValVec.end());

    int   nTruncVal  = (1. - fTruncMeanFraction) * localCorValVec.size();
    float corValSum  = std::accumulate(localCorValVec.begin(),localCorValVec.begin() + nTruncVal,0.);
    float meanCorVal = corValSum / float(nTruncVal);

    std::vector<float> diffVec(nTruncVal);
    std::transform(localCorValVec.begin(),localCorValVec.begin() + nTruncVal, diffVec.begin(), std::bind(std::minus<float>(),std::placeholders::_1,meanCorVal));

    float rmsValSq   = std::inner_product(diffVec.begin(),diffVec.end(),diffVec.begin(),0.);
    float rmsVal     = std::sqrt(rmsValSq / float(nTruncVal));

    // Now set up to run through and do a "simple" interpolation over outliers
    std::vector<float>::iterator lastGoodItr = corValVec.begin();

    bool wasOutlier(false);

    for(std::vector<float>::iterator corValItr = lastGoodItr+1; corValItr != corValVec.end(); corValItr++)
    {
        if (fabs(*corValItr - meanCorVal) < fNumRmsToSmoothVec.at(viewIdx)*rmsVal)
        {
            if (wasOutlier)
            {
                float lastVal  = *lastGoodItr;
                float curVal   = *corValItr;
                float numTicks = std::distance(lastGoodItr,corValItr);
                float slope    = (curVal - lastVal) / numTicks;

                while(lastGoodItr != corValItr)
                {
                    *lastGoodItr++ = (numTicks - std::distance(lastGoodItr,corValItr)) * slope + lastVal;
                }
            }

            wasOutlier  = false;
            lastGoodItr = corValItr;
        }
        else wasOutlier = true;
    }

    return;
}

void RawDigitCorrelatedCorrectionAlg::removeCorrelatedNoise(RawDigitAdcIdxPair& digitIdxPair,
                                                            unsigned int        planeIdx,
                                                            std::vector<float>& truncMeanWireVec,
                                                            std::vector<float>& truncRmsWireVec,
                                                            std::vector<short>& minMaxWireVec,
                                                            std::vector<short>& meanWireVec,
                                                            std::vector<float>& skewnessWireVec,
                                                            std::vector<float>& neighborRatioWireVec,
                                                            std::vector<float>& pedCorWireVec,
                                                            unsigned int& fftSize, unsigned int& halfFFTSize,
							    void* fplan, void* rplan) const
{
    // This method represents and enhanced implementation of "Corey's Algorithm" for correcting the
    // correlated noise across a group of wires. The primary enhancement involves using a FFT to
    // "fit" for the underlying noise as a way to reduce the impact on the signal.
    WireToRawDigitVecMap& wireToRawDigitVecMap = digitIdxPair.first;
    WireToAdcIdxMap&      wireToAdcIdxMap      = digitIdxPair.second;

    size_t maxTimeSamples(wireToRawDigitVecMap.begin()->second.size());
    size_t baseWireIdx(wireToRawDigitVecMap.begin()->first - wireToRawDigitVecMap.begin()->first % fNumWiresToGroup[planeIdx]);

    std::vector<float> corValVec(maxTimeSamples);

    // First step is to get the correction values to apply to this set of input waveforms
    // Don't try to do correction if too few wires unless they have gaps
    if (wireToAdcIdxMap.size() > 2) // || largestGapSize > 2)
    {
        // Zero? This is probably not necessary
        std::fill(corValVec.begin(),corValVec.end(),0.);

        // Build the vector of corrections for each time bin
        for(size_t sampleIdx = 0; sampleIdx < maxTimeSamples; sampleIdx++)
        {
            // Define a vector for accumulating values...
            // Loop over the wires at this time bin and get their pedestal corrected ADC values
            // We'll use a simple stl vector for this
            std::vector<float> adcValuesVec;

            for(const auto& wireAdcItr : wireToAdcIdxMap)
            {
                // Check that we should be doing something in this range
                // Note that if the wire is not to be considered then the "start" bin will be after the last bin
                if (sampleIdx < wireAdcItr.second.first || sampleIdx >= wireAdcItr.second.second) continue;

                int wireIdx(wireAdcItr.first - baseWireIdx);

                // Accumulate
                adcValuesVec.push_back(float(wireToRawDigitVecMap.at(wireAdcItr.first)[sampleIdx]) - truncMeanWireVec[wireIdx]);
            }

            float medianValue = getMedian(adcValuesVec, float(-10000.));

            corValVec[sampleIdx] = medianValue;
        }

        // Try to eliminate any real outliers
        if (fApplyCorSmoothing) smoothCorrectionVec(corValVec, planeIdx);

        // Get the FFT correction
        if (fApplyFFTCorrection) {
          std::vector<std::complex<double>> fftOutputVec(halfFFTSize);
          util::LArFFTW lfftw(fftSize, fplan, rplan, 0);
          lfftw.DoFFT(corValVec, fftOutputVec);

          std::vector<double> powerVec(halfFFTSize);
          std::transform(fftOutputVec.begin(), fftOutputVec.begin() + halfFFTSize, powerVec.begin(), [](const auto& val){return std::abs(val);});

          // Want the first derivative
          std::vector<double> firstDerivVec(powerVec.size(), 0.);
    
          //fWaveformTool->firstDerivative(powerVec, firstDerivVec);
          for(size_t idx = 1; idx < firstDerivVec.size() - 1; idx++)
              firstDerivVec.at(idx) = 0.5 * (powerVec.at(idx + 1) - powerVec.at(idx - 1));

          // Find the peaks
          std::vector<std::tuple<size_t,size_t,size_t>> peakTupleVec;
    
          findPeaks(firstDerivVec.begin(),firstDerivVec.end(),peakTupleVec,fFFTMinPowerThreshold[planeIdx],0);
    
          if (!peakTupleVec.empty())
          {
              for(const auto& peakTuple : peakTupleVec)
              {
                  size_t startTick = std::get<0>(peakTuple);
                  size_t stopTick  = std::get<2>(peakTuple);
            
                  if (stopTick > startTick)
                  {
                      std::complex<double> slope = (fftOutputVec[stopTick] - fftOutputVec[startTick]) / double(stopTick - startTick);
                
                      for(size_t tick = startTick; tick < stopTick; tick++)
                      {
                          std::complex<double> interpVal = fftOutputVec[startTick] + double(tick - startTick) * slope;
                    
                          fftOutputVec[tick]                   = interpVal;
                          //fftOutputVec[fftDataSize - tick - 1] = interpVal;
                      }
                  }
              }
        
              std::vector<double> tmpVec(corValVec.size());
        
              lfftw.DoInvFFT(fftOutputVec, tmpVec);
        
              std::transform(corValVec.begin(),corValVec.end(),tmpVec.begin(),corValVec.begin(),std::minus<double>());
          }
        } // fApplyFFTCorrection

        // Now go through and apply the correction
        for(size_t sampleIdx = 0; sampleIdx < maxTimeSamples; sampleIdx++)
        {
            // Now run through and apply correction
            for (const auto& wireAdcItr : wireToAdcIdxMap)
            {
                float corVal(corValVec[sampleIdx]);
                int   wireIdx(wireAdcItr.first - baseWireIdx);

                // If the "start" bin is after the "stop" bin then we are meant to skip this wire in the averaging process
                // Or if the sample index is in a chirping section then no correction is applied.
                // Both cases are handled by looking at the sampleIdx
                if (sampleIdx < wireAdcItr.second.first || sampleIdx >= wireAdcItr.second.second)
                    corVal = 0.;

                //RawDigitVector& rawDataTimeVec = wireToRawDigitVecMap.at(wireIdx);
                short& rawDataTimeVal = wireToRawDigitVecMap.at(wireAdcItr.first)[sampleIdx];

                // Probably doesn't matter, but try to get slightly more accuracy by doing float math and rounding
                float newAdcValueFloat = float(rawDataTimeVal) - corVal - pedCorWireVec[wireIdx];
                rawDataTimeVal = std::round(newAdcValueFloat);
            }
        }
    }
    return;
}

template<class T> T RawDigitCorrelatedCorrectionAlg::getMedian(std::vector<T>& valuesVec, T defaultValue) const
{
    T medianValue(defaultValue);

    if (!valuesVec.empty())
    {
        std::sort(valuesVec.begin(),valuesVec.end());

        size_t medianIdx = valuesVec.size() / 2;

        medianValue = valuesVec[medianIdx];

        if (valuesVec.size() > 1 && medianIdx % 2) medianValue = (medianValue + valuesVec[medianIdx+1]) / 2;
    }

    return std::max(medianValue,defaultValue);
}

template <typename T> void RawDigitCorrelatedCorrectionAlg::findPeaks(typename std::vector<T>::iterator startItr,
                                                                        typename std::vector<T>::iterator stopItr,
                                                                        std::vector<std::tuple<size_t,size_t,size_t>>& peakTupleVec,
                                                                        T threshold,
                                                                        size_t firstTick) const
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
                peakTupleVec.push_back(std::tuple<size_t,size_t,size_t>(firstBin+firstTick,peakBin+firstTick,lastBin+firstTick));
        
                // Find downstream peaks
                findPeaks(secondItr, stopItr, peakTupleVec, threshold, firstTick + std::distance(startItr,secondItr));
            }
        }
    }

    return;
}



}
