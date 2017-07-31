#include "RawDigitFilterAlg.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
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
RawDigitFilterAlg::RawDigitFilterAlg(fhicl::ParameterSet const & pset) :
    fFirstEvent(false),
    fHistsInitialized(false),
    fCharacterizationAlg(pset),
    fBinAverageAlg(pset),
    fPedestalRetrievalAlg(art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider())
{
    reconfigure(pset);

    // Report.
    mf::LogInfo("RawDigitFilterAlg") << "RawDigitFilterAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitFilterAlg::~RawDigitFilterAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitFilterAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fTruncMeanFraction     = pset.get<float> ("TruncMeanFraction",  0.15);
    fStructuringElement    = pset.get<size_t>("StructuringElement",   30);
    fFillHistograms        = pset.get<bool>  ("FillHistograms",    false);
    
    fCharacterizationAlg.reconfigure(pset);
}

void RawDigitFilterAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fFillHistograms)
    {
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        
        fErosionVecHists.resize(4);
        fDilationVecHists.resize(4);
        
        for(size_t wireIdx = 0; wireIdx < 4; wireIdx++)
        {
            
            std::string vecName = "ErosionVec_" + std::to_string(wireIdx);
            
            fErosionVecHists[wireIdx] = tfs->make<TProfile>(vecName.c_str(), "Erosion Vector;tick", readOutSize, 0., readOutSize, -500., 500.);
            
            vecName = "DilationVec_" + std::to_string(wireIdx);
            
            fDilationVecHists[wireIdx] = tfs->make<TProfile>(vecName.c_str(), "Dilation Vector;tick", readOutSize, 0., readOutSize, -500., 500.);
        }
        
        fHistsInitialized = true;
    }
    
    return;
}
    
void RawDigitFilterAlg::doTopHatFilter(RawDigitVector& dataVec, size_t wire) const
{
    // We need to start by collecting mean/rms
    float truncMean;
    float rmsVal;
    float fracBins(2.*fTruncMeanFraction);
    
    fCharacterizationAlg.getMeanAndRms(dataVec, truncMean, rmsVal, fracBins);
    
    // round the mean to the nearest integer
    truncMean = std::round(truncMean);
    
    // Don't run baseline correction of rms is "good"
    if (rmsVal < 0.) return;
    
    RawDigitVector tempDataVec(dataVec);
    
    fBinAverageAlg.doBinAverage(tempDataVec, fStructuringElement/2);
    
    // Move the short int data into a vector of floats
    typedef std::vector<float> RawDigitFloatVec;
    
    RawDigitFloatVec floatDataVec(dataVec.size());
    
    // while moving, do a pedestal subtraction to handle this just once
    std::transform(tempDataVec.begin(),tempDataVec.end(),floatDataVec.begin(),[truncMean](const short& dataVal){return float(dataVal) - truncMean;});
    
    // Declare erosion and dilation vectors on the original data
    RawDigitFloatVec erosionVec(dataVec.size());
    RawDigitFloatVec dilationVec(dataVec.size());
    
    RawDigitFloatVec::iterator erosionItr  = erosionVec.begin();
    RawDigitFloatVec::iterator dilationItr = dilationVec.begin();
    
    size_t halfStructure(fStructuringElement/2);
    
    // Define sort function so we don't have to keep repeating it
    auto CompareToMean = [](const float& left, const float& right) {return abs(left) < abs(right);};
    //auto CompareToMean = [](const float& left, const float& right) {return left < right;};
    
    // Try to keep track of smallest element to reduce look ups
    RawDigitFloatVec::iterator smallestItr = std::min_element(floatDataVec.begin(), floatDataVec.begin() + halfStructure, CompareToMean);
    RawDigitFloatVec::iterator largestItr  = std::max_element(floatDataVec.begin(), floatDataVec.begin() + halfStructure, CompareToMean);
    
    // First pass through to build the erosion and dilation vectors for the first pass
    for(RawDigitFloatVec::iterator dataItr = floatDataVec.begin(); dataItr != floatDataVec.end(); dataItr++)
    {
        size_t startOffset = std::distance(floatDataVec.begin(), dataItr);
        size_t stopOffset  = std::min(halfStructure,size_t(std::distance(dataItr,floatDataVec.end())));
        
        // Erosion vector handled first
        // Sliding window - check case where the element being dropped from the window was the smallest element,
        // if so then we need to find a new smallest element
        if (startOffset > halfStructure && smallestItr == dataItr - halfStructure - 1)
            smallestItr = std::min_element(dataItr - halfStructure, dataItr + stopOffset, CompareToMean);
        
        // Check if the element being added is the smallest element
        if (stopOffset == halfStructure && abs(*smallestItr) >= abs(*(dataItr + halfStructure - 1)))
        //if (stopOffset == halfStructure && *smallestItr > *(dataItr + halfStructure - 1))
            smallestItr = dataItr + halfStructure - 1;
        
        *erosionItr++ = *smallestItr;
        
        // Now handle the dilation vector
        // Sliding window - check case where the element being dropped from the window was the largest element,
        // if so then we need to find a new largest element
        if (startOffset > halfStructure && largestItr == dataItr - halfStructure - 1)
            largestItr = std::max_element(dataItr - halfStructure, dataItr + stopOffset, CompareToMean);
        
        // Check if the element being added is the largest element
        if (stopOffset == halfStructure && abs(*largestItr) <= abs(*(dataItr + halfStructure - 1)))
        //if (stopOffset == halfStructure && *largestItr < *(dataItr + halfStructure - 1))
            largestItr = dataItr + halfStructure - 1;
        
        *dilationItr++ = *largestItr;
    }
    
    if (fHistsInitialized && (wire == 1596 || wire == 1613 || wire == 1619 || wire == 1623))
    {
        RawDigitFloatVec testAveVec;
        testAveVec.resize(dataVec.size());
        std::transform(erosionVec.begin(),erosionVec.end(),dilationVec.begin(),testAveVec.begin(),[](const float& left, const float& right){return 0.5 * (left + right);});
        
        size_t histIdx = 0;
        if      (wire == 1613) histIdx = 1;
        else if (wire == 1619) histIdx = 2;
        else if (wire == 1623) histIdx = 3;
        
        for(size_t tickIdx = 0; tickIdx < testAveVec.size(); tickIdx++) fErosionVecHists[histIdx]->Fill(tickIdx, testAveVec[tickIdx], 1.);
    }
    
    // Now declare the dilation vector for the erosion vector above
    RawDigitFloatVec openingVec(dataVec.size());
    
    RawDigitFloatVec::iterator openingItr = openingVec.begin();
    
    // Try to keep track of smallest element to reduce look ups
    largestItr = std::max_element(erosionVec.begin(), erosionVec.begin() + halfStructure, CompareToMean);
    
    // Go through the erosion vector to build the dilation vector
    for(RawDigitFloatVec::iterator dataItr = erosionVec.begin(); dataItr != erosionVec.end(); dataItr++)
    {
        size_t startOffset = std::distance(erosionVec.begin(), dataItr);
        size_t stopOffset  = std::min(halfStructure,size_t(std::distance(dataItr,erosionVec.end())));
        
        if (startOffset > halfStructure && largestItr == dataItr - halfStructure - 1)
            largestItr = std::max_element(dataItr - halfStructure, dataItr + stopOffset, CompareToMean);
        
        if (stopOffset == halfStructure && abs(*largestItr) <= abs(*(dataItr + halfStructure - 1)))
        //if (stopOffset == halfStructure && *largestItr < *(dataItr + halfStructure - 1))
            largestItr = dataItr + halfStructure - 1;
        
        *openingItr++ = *largestItr;
    }
    
    // Now declare the erosion vector for the dilation vector above
    RawDigitFloatVec closingVec(dataVec.size());
    
    RawDigitFloatVec::iterator closingItr = closingVec.begin();
    
    // Try to keep track of smallest element to reduce look ups
    smallestItr = std::min_element(dilationVec.begin(), dilationVec.begin() + halfStructure, CompareToMean);
    
    // Go through the erosion vector to build the dilation vector
    for(RawDigitFloatVec::iterator dataItr = dilationVec.begin(); dataItr != dilationVec.end(); dataItr++)
    {
        size_t startOffset = std::distance(dilationVec.begin(), dataItr);
        size_t stopOffset  = std::min(halfStructure,size_t(std::distance(dataItr,dilationVec.end())));
        
        if (startOffset > halfStructure && smallestItr == dataItr - halfStructure - 1)
            smallestItr = std::min_element(dataItr - halfStructure, dataItr + stopOffset, CompareToMean);
        
        if (stopOffset == halfStructure && abs(*smallestItr) >= abs(*(dataItr + halfStructure - 1)))
        //if (stopOffset == halfStructure && *smallestItr > *(dataItr + halfStructure - 1))
            smallestItr = dataItr + halfStructure - 1;
        
        *closingItr++ = *smallestItr;
    }
    
    // Now get the average between the two
    RawDigitFloatVec aveBaselineVec(dataVec.size());
    
    //std::transform(openingVec.begin(),openingVec.end(),closingVec.begin(),aveBaselineVec.begin(),[](const float& openVal, const float& closeVal){return 0.5 * (openVal + closeVal);});
    
    std::transform(openingVec.begin(),openingVec.end(),closingVec.begin(),aveBaselineVec.begin(),[](const float& openVal, const float& closeVal){return abs(openVal) < abs(closeVal) ? openVal : closeVal;});
    
    if (fHistsInitialized && (wire == 1596 || wire == 1613 || wire == 1619 || wire == 1623 || wire == 1108))
    {
        size_t histIdx = 0;
        if      (wire == 1613) histIdx = 1;
        else if (wire == 1619) histIdx = 2;
        else if (wire == 1623) histIdx = 3;
        
        for(size_t tickIdx = 0; tickIdx < aveBaselineVec.size(); tickIdx++) fDilationVecHists[histIdx]->Fill(tickIdx, aveBaselineVec[tickIdx], 1.);
        
        if (wire == 1108)
        {
            for(size_t tickIdx = 0; tickIdx < aveBaselineVec.size(); tickIdx++)
                std::cout << "TickIdx: " << tickIdx << ", erosion/dilation: " << erosionVec[tickIdx] << "/" << dilationVec[tickIdx] << ", open/close: " << openingVec[tickIdx] << "/" << closingVec[tickIdx] << ", correction: " << aveBaselineVec[tickIdx] << ", data val: " << floatDataVec[tickIdx] << std::endl;
        }
    }
    
    // Now do baseline correction
    std::transform(dataVec.begin(),dataVec.end(),aveBaselineVec.begin(),dataVec.begin(),[](const short& left, const short& right){return std::round(float(left) - right);});
    
    return;
}

void RawDigitFilterAlg::doAdaptiveFilter(RawDigitVector& dataVec) const
{
    // We need to start by collecting mean/rms
    float truncMean;
    float rmsVal;
    
    fCharacterizationAlg.getMeanAndRms(dataVec, truncMean, rmsVal, fTruncMeanFraction);
    
    // Declare erosion vector
    std::vector<float> meanVec;
    
    meanVec.resize(dataVec.size(), 0);
    std::vector<float>::iterator meanVecItr = meanVec.begin();
    
    size_t halfStructure(fStructuringElement/2);
    
    float runningSum(0.);
    
    for(size_t idx = 0; idx < halfStructure - 1; idx++) runningSum += dataVec[idx];
    
    // First pass through to build the erosion vector
    for(RawDigitVector::iterator dataItr = dataVec.begin(); dataItr != dataVec.end(); dataItr++)
    {
        size_t startOffset = std::distance(dataVec.begin(),dataItr);
        size_t stopOffset  = std::distance(dataItr,dataVec.end());
        size_t count       = std::min(2 * halfStructure, std::min(startOffset + halfStructure, halfStructure + stopOffset));
        
        if (startOffset >  halfStructure) runningSum -= *(dataItr - halfStructure - 1);
        if (stopOffset  >= halfStructure) runningSum += *(dataItr + halfStructure - 1);
        
        float adaptVal = runningSum / float(count) - truncMean;
        
        *meanVecItr++ = adaptVal;
    }
    
    // Now do baseline correction
    for(size_t idx = 0; idx < dataVec.size(); idx++)
        dataVec[idx] = dataVec[idx] - std::round(meanVec[idx]);

    return;
}
    
}
