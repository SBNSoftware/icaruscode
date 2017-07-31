
#include "RawDigitCharacterizationAlg.h"

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
RawDigitCharacterizationAlg::RawDigitCharacterizationAlg(fhicl::ParameterSet const & pset) :
                      fHistsInitialized(false),
                      fFirstEvent(true),
                      fChannelGroups(pset),
                      fPedestalRetrievalAlg(art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider())

{
    reconfigure(pset);

    // Report.
    mf::LogInfo("RawDigitCharacterizationAlg") << "RawDigitCharacterizationAlg configured\n";
}
    
//----------------------------------------------------------------------------
/// Destructor.
RawDigitCharacterizationAlg::~RawDigitCharacterizationAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitCharacterizationAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                        0.15);
    fRmsRejectionCutHi     = pset.get<std::vector<float>> ("RMSRejectionCutHi",     std::vector<float>() = {25.0,25.0,25.0});
    fRmsRejectionCutLow    = pset.get<std::vector<float>> ("RMSRejectionCutLow",    std::vector<float>() = {0.70,0.70,0.70});
    fRmsSelectionCut       = pset.get<std::vector<float>> ("RMSSelectionCut",       std::vector<float>() = {1.40,1.40,1.00});
    fMinMaxSelectionCut    = pset.get<std::vector<short>> ("MinMaxSelectionCut",        std::vector<short>() = {13, 13, 11});
    fTheChosenWire         = pset.get<unsigned int>       ("TheChosenWire",                                            1200);
    fMaxPedestalDiff       = pset.get<double>             ("MaxPedestalDiff",                                           10.);
    fHistsWireGroup        = pset.get<std::vector<size_t>>("FFTHistsWireGroup",         std::vector<size_t>() = {1, 33, 34});
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",          std::vector<size_t>() = {48, 48, 96});
    fFillHistograms        = pset.get<bool>               ("FillHistograms",                                           true);
}

//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitCharacterizationAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fFillHistograms)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        fAdcCntHist[0]    = tfs->make<TH1D>("CntUPlane", ";#adc",  200, 9000., 10000.);
        fAdcCntHist[1]    = tfs->make<TH1D>("CntVPlane", ";#adc",  200, 9000., 10000.);
        fAdcCntHist[2]    = tfs->make<TH1D>("CntWPlane", ";#adc",  200, 9000., 10000.);
        fAveValHist[0]    = tfs->make<TH1D>("AveUPlane", ";Ave",   120,  -30.,    30.);
        fAveValHist[1]    = tfs->make<TH1D>("AveVPlane", ";Ave",   120,  -30.,    30.);
        fAveValHist[2]    = tfs->make<TH1D>("AveWPlane", ";Ave",   120,  -30.,    30.);
        fRmsValHist[0]    = tfs->make<TH1D>("RmsUPlane", ";RMS",   200,    0.,    50.);
        fRmsValHist[1]    = tfs->make<TH1D>("RmsVPlane", ";RMS",   200,    0.,    50.);
        fRmsValHist[2]    = tfs->make<TH1D>("RmsWPlane", ";RMS",   200,    0.,    50.);
        fPedValHist[0]    = tfs->make<TH1D>("PedUPlane", ";Ped",   200,  1950,  2150.);
        fPedValHist[1]    = tfs->make<TH1D>("PedVPlane", ";Ped",   200,  1950,  2150.);
        fPedValHist[2]    = tfs->make<TH1D>("PedWPlane", ";Ped",   200,   350,   550.);
    
        fRmsValProf[0]    = tfs->make<TProfile>("RmsPlane0Prof",    ";Wire #",  1200, 0., 1200., 0., 100.);
        fRmsValProf[1]    = tfs->make<TProfile>("RmsPlane1Prof",    ";Wire #",  5000, 0., 5000., 0., 100.);
        fRmsValProf[2]    = tfs->make<TProfile>("RmsPlane2Prof",    ";Wire #",  5000, 0., 5000., 0., 100.);
    
        fMinMaxValProf[0] = tfs->make<TProfile>("MinMaxPlane0Prof", ";Wire #",  1200, 0., 1200., 0., 200.);
        fMinMaxValProf[1] = tfs->make<TProfile>("MinMaxPlane1Prof", ";Wire #",  5000, 0., 5000., 0., 200.);
        fMinMaxValProf[2] = tfs->make<TProfile>("MinMaxPlane2Prof", ";Wire #",  5000, 0., 5000., 0., 200.);

        fPedValProf[0]    = tfs->make<TProfile>("PedPlane0Prof",    ";Wire #",  1200, 0., 1200., 1500., 2500.);
        fPedValProf[1]    = tfs->make<TProfile>("PedPlane1Prof",    ";Wire #",  5000, 0., 5000., 1500., 2500.);
        fPedValProf[2]    = tfs->make<TProfile>("PedPlane2Prof",    ";Wire #",  5000, 0., 5000.,    0., 1000.);
    
        fAverageHist[0]   = tfs->make<TH1D>("Average0", ";Bin", 1000, 1500., 2500.);
        fAverageHist[1]   = tfs->make<TH1D>("Average1", ";Bin", 1000, 1500., 2500.);
        fAverageHist[2]   = tfs->make<TH1D>("Average2", ";Bin", 1000,    0., 1000.);
    
        fMinMaxProfiles.resize(3);
        fSkewnessProfiles.resize(3);
        fModeRatioProfiles.resize(3);
        
        for(size_t viewIdx = 0; viewIdx < 3; viewIdx++)
        {
            std::string minMaxName = "MinMax_" + std::to_string(viewIdx);
        
            fMinMaxProfiles[viewIdx] = tfs->make<TProfile>(minMaxName.c_str(), "Min/Max Profiles;Wire", fNumWiresToGroup[viewIdx], 0., fNumWiresToGroup[viewIdx], 0., 200.);
        
            minMaxName = "Skewness_" + std::to_string(viewIdx);
        
            fSkewnessProfiles[viewIdx] = tfs->make<TProfile>(minMaxName.c_str(), "Skewness Profiles;Wire", fNumWiresToGroup[viewIdx], 0., fNumWiresToGroup[viewIdx], -4., 4.);
        
            minMaxName = "ModeRatio_" + std::to_string(viewIdx);
        
            fModeRatioProfiles[viewIdx] = tfs->make<TProfile>(minMaxName.c_str(), "Mode Ratio;Wire", fNumWiresToGroup[viewIdx], 0., fNumWiresToGroup[viewIdx], 0., 1.2);
        }
    
        fHistsInitialized = true;
    }
    
    return;
}

// Basic waveform mean, rms and pedestal offset
void RawDigitCharacterizationAlg::getWaveformParams(const RawDigitVector& rawWaveform,
                                                    unsigned int          channel,
                                                    unsigned int          view,
                                                    unsigned int          wire,
                                                    float&                truncMean,
                                                    float&                truncRms,
                                                    short&                mean,
                                                    short&                median,
                                                    short&                mode,
                                                    float&                skewness,
                                                    float&                rms,
                                                    short&                minMax,
                                                    float&                neighborRatio,
                                                    float&                pedCorVal) const
{
    // We start by finding the most likely baseline which is most easily done by
    // finding the most populated bin and the average using the neighboring bins
    // To do this we'll use a map with key the bin number and data the count in that bin
    // Define the map first
    std::map<short,short> binAdcMap;
    
    // Populate the map
    for(const auto& adcVal : rawWaveform) binAdcMap[adcVal]++;
    
    // Find the max bin and the count
    std::map<short,short>::iterator maxBinItr = std::max_element(binAdcMap.begin(),binAdcMap.end(),[](const auto& lhs,const auto& rhs){return lhs.second < rhs.second;});
    short                           binMax(maxBinItr->first);
    short                           binMaxCnt(maxBinItr->second);
    
    // fill example hists - throw away code
    if (fHistsInitialized && fFirstEvent && wire == fTheChosenWire)
    {
        for(const auto& binAdcItr : binAdcMap)
        {
            fAverageHist[view]->Fill(binAdcItr.first, binAdcItr.second);
        }
    }
    
    // Armed with the max bin and its count, now set up to get an average
    // about this bin. We'll want to cut off at some user defined fraction
    // of the total bins on the wire
    size_t maxTimeSamples(rawWaveform.size());
    int    minNumBins = (1. - fTruncMeanFraction) * maxTimeSamples;
    int    curBinCnt(binMaxCnt);
    short  binOffset(1);
    
    truncMean = curBinCnt * binMax;
    
    // This loop to develop the average
    // In theory, we could also keep the sum of the squares for the rms but I had problems doing
    // it that way so will loop twice... (potential time savings goes here!)
    while(curBinCnt < minNumBins)
    {
        if (binAdcMap[binMax-binOffset])
        {
            curBinCnt += binAdcMap[binMax-binOffset];
            truncMean += double(binAdcMap[binMax-binOffset] * (binMax - binOffset));
        }
        
        if (binAdcMap[binMax+binOffset])
        {
            curBinCnt += binAdcMap[binMax+binOffset];
            truncMean += double(binAdcMap[binMax+binOffset] * (binMax + binOffset));
        }
        
        binOffset++;
    }
    
    truncMean /= float(curBinCnt);
    
    getTruncatedRMS(rawWaveform, truncMean, truncRms);
    
    // Recover the database version of the pedestal
    float pedestal = fPedestalRetrievalAlg.PedMean(channel);
    
    // The pedCorVal will transform from the average of waveform calculated here to the expected value of the pedestal.
    pedCorVal = truncMean - pedestal;
    
    // Determine the range of ADC values on this wire
    short minVal = *std::min_element(rawWaveform.begin(), rawWaveform.end());
    short maxVal = *std::max_element(rawWaveform.begin(), rawWaveform.end());
    
    minMax = std::min(maxVal - minVal + 1, 199);  // for the purposes of histogramming
    
    // We also want mean, median, rms, etc., for all ticks on the waveform
    std::vector<short> localTimeVec = rawWaveform;
    
    std::sort(localTimeVec.begin(),localTimeVec.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float realMean(float(std::accumulate(localTimeVec.begin(),localTimeVec.end(),0))/float(localTimeVec.size()));
    
    median = localTimeVec[localTimeVec.size()/2];
    mean   = std::round(realMean);
    
    std::vector<float> adcLessPedVec;
    
    adcLessPedVec.resize(localTimeVec.size());
    
    std::transform(localTimeVec.begin(),localTimeVec.end(),adcLessPedVec.begin(),std::bind2nd(std::minus<short>(),mean));
    
    rms      = std::sqrt(std::inner_product(adcLessPedVec.begin(), adcLessPedVec.end(), adcLessPedVec.begin(), 0.) / float(adcLessPedVec.size()));
    skewness = 3. * float(mean - median) / rms;
    
    // Final task is to get the mode and neighbor ratio
    // Fortunately, the map was already set up
    short neighborSum(0);
    short leftNeighbor(maxBinItr->second);
    short rightNeighbor(maxBinItr->second);
    short cnt(0);
    
    mode = maxBinItr->first;
    
    if (binAdcMap.find(mode-1) != binAdcMap.end())
    {
        leftNeighbor  = binAdcMap.find(mode-1)->second;
        neighborSum  += leftNeighbor;
        cnt++;
    }
    
    if (binAdcMap.find(mode+1) != binAdcMap.end())
    {
        rightNeighbor  = binAdcMap.find(mode+1)->second;
        neighborSum   += rightNeighbor;
        cnt++;
    }
    
    neighborRatio = float(neighborSum) / float(2*maxBinItr->second);

    neighborRatio = float(std::min(leftNeighbor,rightNeighbor)) / float(maxBinItr->second);
//    float leastNeighborRatio = float(std::min(leftNeighbor,rightNeighbor)) / float(maxBinItr->second);
    
    // Fill some histograms here
    if (fHistsInitialized)
    {
        short maxVal = *std::max_element(rawWaveform.begin(),rawWaveform.end());
        short minVal = *std::min_element(rawWaveform.begin(),rawWaveform.end());
        short minMax = std::min(maxVal - minVal,199);
        
        fAdcCntHist[view]->Fill(curBinCnt, 1.);
        fAveValHist[view]->Fill(std::max(-29.9, std::min(29.9,double(truncMean - pedestal))), 1.);
        fRmsValHist[view]->Fill(std::min(49.9, double(truncRms)), 1.);
        fRmsValProf[view]->Fill(wire, double(truncRms), 1.);
        fMinMaxValProf[view]->Fill(wire, double(minMax), 1.);
        fPedValProf[view]->Fill(wire, truncMean, 1.);
        fPedValHist[view]->Fill(truncMean, 1.);
    }
    
    
    if (wire / fNumWiresToGroup[view] == fHistsWireGroup[view])
    {
        float  leastNeighborRatio = float(std::min(leftNeighbor,rightNeighbor)) / float(maxBinItr->second);
        size_t wireIdx            = wire % fNumWiresToGroup[view];
        
        if (skewness > 0. && leastNeighborRatio < 0.7)
        {
            short threshold(6);
            
            RawDigitVector::const_iterator stopChirpItr = std::find_if(rawWaveform.begin(),rawWaveform.end(),[mean,threshold](const short& elem){return abs(elem - mean) > threshold;});
        }
        
        if (fHistsInitialized)
        {
            fMinMaxProfiles[view]->Fill(double(wireIdx+0.5), double(minMax), 1.);
            fSkewnessProfiles[view]->Fill(double(wireIdx+0.5), double(skewness), 1.);
            fModeRatioProfiles[view]->Fill(double(wireIdx+0.5), double(leastNeighborRatio), 1.);
        }
    }
    
    return;
}
void RawDigitCharacterizationAlg::getTruncatedRMS(const RawDigitVector& rawWaveform,
                                                  float&                pedestal,
                                                  float&                truncRms) const
{
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> adcLessPedVec;
    
    adcLessPedVec.resize(rawWaveform.size());
    
    // Fill the vector
    std::transform(rawWaveform.begin(),rawWaveform.end(),adcLessPedVec.begin(),std::bind2nd(std::minus<short>(),pedestal));
    
    // sort in ascending order so we can truncate the sume
    std::sort(adcLessPedVec.begin(), adcLessPedVec.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    int minNumBins = (1. - fTruncMeanFraction) * rawWaveform.size();
    
    // Get the truncated sum
    truncRms = std::inner_product(adcLessPedVec.begin(), adcLessPedVec.begin() + minNumBins, adcLessPedVec.begin(), 0.);
    truncRms = std::sqrt(std::max(0.,truncRms / double(minNumBins)));
    
    return;
}

void RawDigitCharacterizationAlg::getMeanRmsAndPedCor(const RawDigitVector& rawWaveform,
                                                      unsigned int          channel,
                                                      unsigned int          view,
                                                      unsigned int          wire,
                                                      float&                truncMean,
                                                      float&                rmsVal,
                                                      float&                pedCorVal) const
{
    // The strategy for finding the average for a given wire will be to
    // find the most populated bin and the average using the neighboring bins
    // To do this we'll use a map with key the bin number and data the count in that bin
    // Define the map first
    std::map<short,short> binAdcMap;
    
    // Populate the map
    for(const auto& adcVal : rawWaveform)
    {
        binAdcMap[adcVal]++;
    }
    
    // Find the max bin and the count
    std::map<short,short>::iterator maxBinItr = std::max_element(binAdcMap.begin(),binAdcMap.end(),[](const auto& lhs,const auto& rhs){return lhs.second < rhs.second;});
    short                           binMax(maxBinItr->first);
    short                           binMaxCnt(maxBinItr->second);
    
    // fill example hists - throw away code
    if (fHistsInitialized && fFirstEvent && wire == fTheChosenWire)
    {
        for(const auto& binAdcItr : binAdcMap)
        {
            fAverageHist[view]->Fill(binAdcItr.first, binAdcItr.second);
        }
    }
    
    // Armed with the max bin and its count, now set up to get an average
    // about this bin. We'll want to cut off at some user defined fraction
    // of the total bins on the wire
    size_t maxTimeSamples(rawWaveform.size());
    int    minNumBins = (1. - fTruncMeanFraction) * maxTimeSamples;
    int    curBinCnt(binMaxCnt);
    
    float peakValue(curBinCnt * binMax);
    short binOffset(1);
    
    truncMean = peakValue;
    
    // This loop to develop the average
    // In theory, we could also keep the sum of the squares for the rms but I had problems doing
    // it that way so will loop twice... (potential time savings goes here!)
    while(curBinCnt < minNumBins)
    {
        if (binAdcMap[binMax-binOffset])
        {
            curBinCnt += binAdcMap[binMax-binOffset];
            truncMean += double(binAdcMap[binMax-binOffset] * (binMax - binOffset));
        }
        
        if (binAdcMap[binMax+binOffset])
        {
            curBinCnt += binAdcMap[binMax+binOffset];
            truncMean += double(binAdcMap[binMax+binOffset] * (binMax + binOffset));
        }
        
        binOffset++;
    }
    
    truncMean /= double(curBinCnt);
    
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> adcLessPedVec;
    
    adcLessPedVec.resize(rawWaveform.size());
    
    size_t adcLessPedIdx(0);
    
    for(const auto& adcVal : rawWaveform)
    {
        double adcLessPed = adcVal - truncMean;
        adcLessPedVec[adcLessPedIdx++] = adcLessPed;
    }
    
    // sort in ascending order so we can truncate the sume
    std::sort(adcLessPedVec.begin(), adcLessPedVec.end());
    
    // Get the truncated sum
    rmsVal = std::inner_product(adcLessPedVec.begin(), adcLessPedVec.begin() + minNumBins, adcLessPedVec.begin(), 0.);
    rmsVal = std::sqrt(std::max(0.,rmsVal / double(minNumBins)));
    
    // Recover the database version of the pedestal
    float pedestal = fPedestalRetrievalAlg.PedMean(channel);
    
    pedCorVal = truncMean - pedestal;
    
    // Fill some histograms here
    if (fHistsInitialized)
    {
        short maxVal = *std::max_element(rawWaveform.begin(),rawWaveform.end());
        short minVal = *std::min_element(rawWaveform.begin(),rawWaveform.end());
        short minMax = std::min(maxVal - minVal,199);
        
        fAdcCntHist[view]->Fill(curBinCnt, 1.);
        fAveValHist[view]->Fill(std::max(-29.9, std::min(29.9,double(truncMean - pedestal))), 1.);
        fRmsValHist[view]->Fill(std::min(49.9, double(rmsVal)), 1.);
        fRmsValProf[view]->Fill(wire, double(rmsVal), 1.);
        fMinMaxValProf[view]->Fill(wire, double(minMax), 1.);
        fPedValProf[view]->Fill(wire, truncMean, 1.);
        fPedValHist[view]->Fill(truncMean, 1.);
    }
    
    // Output a message is there is significant different to the pedestal
    if (abs(truncMean - pedestal) > fMaxPedestalDiff)
    {
        mf::LogInfo("RawDigitCharacterizationAlg") << ">>> Pedestal mismatch, channel: " << channel << ", new value: " << truncMean << ", original: " << pedestal << ", rms: " << rmsVal << std::endl;
    }
    
    return;
}

void RawDigitCharacterizationAlg::getMeanAndRms(const RawDigitVector& rawWaveform,
                                                float&                aveVal,
                                                float&                rmsVal,
                                                float                 fracBins) const
{
    // The strategy for finding the average for a given wire will be to
    // find the most populated bin and the average using the neighboring bins
    // To do this we'll use a map with key the bin number and data the count in that bin
    // Define the map first
    std::map<short,short> binAdcMap;
    
    // Populate the map
    for(const auto& adcVal : rawWaveform)
    {
        binAdcMap[adcVal]++;
    }
    
    // Find the max bin
    short binMax(-1);
    short binMaxCnt(0);
    
    for(const auto& binAdcItr : binAdcMap)
    {
        if (binAdcItr.second > binMaxCnt)
        {
            binMax    = binAdcItr.first;
            binMaxCnt = binAdcItr.second;
        }
    }
    
    // Armed with the max bin and its count, now set up to get an average
    // about this bin. We'll want to cut off at some user defined fraction
    // of the total bins on the wire
    size_t maxTimeSamples(rawWaveform.size());
    int    minNumBins = std::round((1. - fracBins) * maxTimeSamples);
    int    curBinCnt(binMaxCnt);
    
    double peakValue(curBinCnt * double(binMax+0.5));
    double truncMean(peakValue);
    
    short binOffset(1);
    
    // This loop to develop the average
    // In theory, we could also keep the sum of the squares for the rms but I had problems doing
    // it that way so will loop twice... (potential time savings goes here!)
    while(curBinCnt < minNumBins)
    {
        if (binAdcMap[binMax-binOffset])
        {
            curBinCnt += binAdcMap[binMax-binOffset];
            truncMean += double(binAdcMap[binMax-binOffset] * double((binMax - binOffset)+0.5));
        }
        
        if (binAdcMap[binMax+binOffset])
        {
            curBinCnt += binAdcMap[binMax+binOffset];
            truncMean += double(binAdcMap[binMax+binOffset] * double((binMax + binOffset)+0.5));
        }
        
        binOffset++;
    }
    
    truncMean /= double(curBinCnt);
    
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> adcLessPedVec;
    
    adcLessPedVec.resize(rawWaveform.size());
    
    size_t adcLessPedIdx(0);
    
    for(const auto& adcVal : rawWaveform)
    {
        double adcLessPed = adcVal - truncMean;
        adcLessPedVec[adcLessPedIdx++] = adcLessPed;
    }
    
    // sort in ascending order so we can truncate the sum
    std::sort(adcLessPedVec.begin(), adcLessPedVec.end());
    
    // Get the truncated sum
    rmsVal = std::inner_product(adcLessPedVec.begin(), adcLessPedVec.begin() + minNumBins, adcLessPedVec.begin(), 0.);
    rmsVal = std::sqrt(std::max(0.,rmsVal / double(minNumBins)));
    aveVal = truncMean;
    
    return;
}

bool RawDigitCharacterizationAlg::classifyRawDigitVec(RawDigitVector&  rawWaveform,
                                                      unsigned int            viewIdx,
                                                      unsigned int            wire,
                                                      float                   truncRms,
                                                      short                   minMax,
                                                      short                   mean,
                                                      float                   skewness,
                                                      float                   neighborRatio,
                                                      GroupToDigitIdxPairMap& groupToDigitIdxPairMap) const
{
    // This simply classifies the input waveform:
    // a) determines if it should be added to the list of waveforms to process
    // b) if to be analyzed, places in the group of wires to process
    bool classified(false);
    
    // Dereference the selection/rejection cut
    float selectionCut = fMinMaxSelectionCut[viewIdx];
    float rejectionCut = fRmsRejectionCutHi[viewIdx];
    
    // Selection to process
    if (minMax > selectionCut && truncRms < rejectionCut)
    {
        size_t group = fChannelGroups.channelGroup(viewIdx,wire);
        
        if (groupToDigitIdxPairMap.find(group) == groupToDigitIdxPairMap.end())
            groupToDigitIdxPairMap.insert(std::pair<size_t,RawDigitAdcIdxPair>(group,RawDigitAdcIdxPair()));
        
        groupToDigitIdxPairMap.at(group).first.insert(WireToRawDigitVecPair(wire,rawWaveform));
        groupToDigitIdxPairMap.at(group).second.insert(std::pair<size_t,RawDigitVectorIdxPair>(wire,RawDigitVectorIdxPair(0,rawWaveform.size())));
        
        // Look for chirping wire sections. Confine this to only the V plane
        if (viewIdx == 1)
        {
            // Do wire shape corrections to look for chirping wires and other oddities to avoid
            // Recover our objects...
            WireToAdcIdxMap& wireToAdcIdxMap = groupToDigitIdxPairMap.at(group).second;

            // Set a threshold
            short threshold(6);
                    
            // If going from quiescent to on again, then the min/max will be large
            //if (skewnessWireVec[wireIdx] > 0. && minMaxWireVec[wireIdx] > 50 && truncRmsWireVec[wireIdx] > 2.)
            if (skewness > 0. && neighborRatio < 0.7 && minMax > 50)
            {
                RawDigitVector::iterator stopChirpItr = std::find_if(rawWaveform.begin(),rawWaveform.end(),[mean,threshold](const short& elem){return abs(elem - mean) > threshold;});
                
                size_t threshIndex = std::distance(rawWaveform.begin(),stopChirpItr);
                
                if (threshIndex > 60) wireToAdcIdxMap[wire].first = threshIndex;
            }
            // Check in the reverse direction?
            else if (minMax > 20 && neighborRatio < 0.7)
            {
                threshold = 3;
                
                RawDigitVector::reverse_iterator startChirpItr = std::find_if(rawWaveform.rbegin(),rawWaveform.rend(),[mean,threshold](const short& elem){return abs(elem - mean) > threshold;});
                
                size_t threshIndex = std::distance(rawWaveform.rbegin(),startChirpItr);
                
                if (threshIndex > 60) wireToAdcIdxMap[wire].second = rawWaveform.size() - threshIndex;
            }
        }
        
        classified = true;
    }
    
    return classified;
}

template<class T> T RawDigitCharacterizationAlg::getMedian(std::vector<T>& valuesVec, T defaultValue) const
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
    
}
