
#include "RawDigitBinAverageAlg.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <cmath>
#include <algorithm>

namespace caldata
{

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitBinAverageAlg::RawDigitBinAverageAlg(fhicl::ParameterSet const & pset)
{
    // Report.
    mf::LogInfo("RawDigitBinAverageAlg") << "RawDigitBinAverageAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitBinAverageAlg::~RawDigitBinAverageAlg()
{}

void RawDigitBinAverageAlg::doBinAverage(RawDigitVector& dataVec,
                                        size_t          binsToAverage) const
{
    size_t halfBinsToAverage(binsToAverage/2);
    
    float runningSum(0.);
    
    for(size_t idx = 0; idx < halfBinsToAverage; idx++) runningSum += dataVec[idx];
    
    // Declare erosion vector
    std::vector<float> meanVec;
    
    meanVec.resize(dataVec.size(), 0);
    std::vector<float>::iterator meanVecItr = meanVec.begin();
    
    // First pass through to build the erosion vector
    for(RawDigitVector::iterator dataItr = dataVec.begin(); dataItr != dataVec.end(); dataItr++)
    {
        size_t startOffset = std::distance(dataVec.begin(),dataItr);
        size_t stopOffset  = std::distance(dataItr,dataVec.end());
        size_t count       = std::min(2 * halfBinsToAverage, std::min(startOffset + halfBinsToAverage + 1, halfBinsToAverage + stopOffset - 1));
        
        if (startOffset >= halfBinsToAverage) runningSum -= *(dataItr - halfBinsToAverage);
        if (stopOffset  >  halfBinsToAverage) runningSum += *(dataItr + halfBinsToAverage);
        
        *meanVecItr++ = runningSum / float(count);
    }
    
    std::transform(meanVec.begin(),meanVec.end(),dataVec.begin(),[](const float& val){return std::round(val);});
    
    return;
}

void RawDigitBinAverageAlg::doTwoBinAverage(RawDigitVector& dataVec) const
{
    // This should be a straightforward transform
    std::transform(dataVec.begin(),dataVec.end()-1,dataVec.begin()+1,dataVec.begin(),[](const short& first, const short& last){return std::round(float(first+last)/2.);});
    
    return;
}
    
}
