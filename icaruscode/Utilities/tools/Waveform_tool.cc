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
    
    void triangleSmooth(std::vector<float>&,  std::vector<float>&,  size_t = 0)                                 const;
    void triangleSmooth(std::vector<double>&, std::vector<double>&, size_t = 0)                                 const;
    void firstDerivative(std::vector<float>&,  std::vector<float>&)                                             const;
    void firstDerivative(std::vector<double>&, std::vector<double>&)                                            const;
    void findPeaks(std::vector<float>::iterator,  std::vector<float>::iterator,  PeakTupleVec&, float, size_t)  const;
    void findPeaks(std::vector<double>::iterator, std::vector<double>::iterator, PeakTupleVec&, double, size_t) const;
    
private:
    template <typename T> void triangleSmooth(std::vector<T>&, std::vector<T>&, size_t = 0)                              const;
    template <typename T> void firstDerivative(std::vector<T>&,  std::vector<T>&)                                        const;
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
    
void Waveform::triangleSmooth(std::vector<double>& inputVec, std::vector<double>& smoothVec, size_t lowestBin) const
{
    triangleSmooth<double>(inputVec, smoothVec, lowestBin);
    
    return;
}
    
void Waveform::triangleSmooth(std::vector<float>& inputVec, std::vector<float>& smoothVec, size_t lowestBin) const
{
    triangleSmooth<float>(inputVec, smoothVec, lowestBin);
    
    return;
}

template <typename T> void Waveform::triangleSmooth(std::vector<T>& inputVec, std::vector<T>& smoothVec, size_t lowestBin) const
{
    if (inputVec.size() != smoothVec.size()) smoothVec.resize(inputVec.size());
    
    std::copy(inputVec.begin(), inputVec.begin() + 2 + lowestBin, smoothVec.begin());
    std::copy(inputVec.end() - 2, inputVec.end(), smoothVec.end() - 2);
    
    typename std::vector<T>::iterator curItr    = smoothVec.begin() + 2 + lowestBin;
    typename std::vector<T>::iterator curInItr  = inputVec.begin() + 1 + lowestBin;
    typename std::vector<T>::iterator stopInItr = inputVec.end()   - 2;
    
    while(curInItr++ != stopInItr)
    {
        // Take the weighted average of five consecutive points centered on current point
        T newVal = (*(curInItr - 2) + 2. * *(curInItr - 1) + 3. * *curInItr + 2. * *(curInItr + 1) + *(curInItr + 2)) / 9.;
        
        *curItr++ = newVal;
    }
    return;
}

void Waveform::firstDerivative(std::vector<double>& inputVec, std::vector<double>& derivVec) const
{
    firstDerivative<double>(inputVec, derivVec);
    
    return;
}
    
void Waveform::firstDerivative(std::vector<float>& inputVec, std::vector<float>& derivVec) const
{
    firstDerivative<float>(inputVec, derivVec);
    
    return;
}
    
template <typename T> void Waveform::firstDerivative(std::vector<T>& inputVec, std::vector<T>& derivVec) const
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

    
DEFINE_ART_CLASS_TOOL(Waveform)
}
