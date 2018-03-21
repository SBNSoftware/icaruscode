////////////////////////////////////////////////////////////////////////
/// \file   OpHitFinder.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

#include "icaruscode/Light/OpticalTools/IOpHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"

#include <cmath>
#include <fstream>

namespace light
{

class OpHitFinder : public IOpHitFinder
{
public:
    explicit OpHitFinder(const fhicl::ParameterSet& pset);
    
    ~OpHitFinder();
    
    void configure(const fhicl::ParameterSet& pset)                 override;
    void outputHistograms(art::TFileDirectory&)               const override;
    
    void FindOpHits(const raw::OpDetWaveform&, recob::OpHit&) const override;
    
private:
    // fhicl parameters
    int    fNumBinsToAverage;
    
    float getBaseline(const raw::OpDetWaveform&) const;
    
    std::unique_ptr<reco_tool::ICandidateHitFinder> fHitFinderTool;  ///< For finding candidate hits
};
    
//----------------------------------------------------------------------
// Constructor.
OpHitFinder::OpHitFinder(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
OpHitFinder::~OpHitFinder()
{
}
    
void OpHitFinder::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fNumBinsToAverage = pset.get<int>("NumBinsToAverage", 20);
    
    fHitFinderTool  = art::make_tool<reco_tool::ICandidateHitFinder>(pset.get<fhicl::ParameterSet>("CandidateHits"));

    return;
}

    
void OpHitFinder::FindOpHits(const raw::OpDetWaveform& opDetWaveform,
                             recob::OpHit&             opHit) const
{
    // The plan here:
    // 1) copy to a local vector
    // 2) Find the mean and rms
    // 3) Fill a map to find the most probable value
    float baseline = getBaseline(opDetWaveform);
    
    std::vector<float> locWaveform;
    
    locWaveform.resize(opDetWaveform.size());
    
    // The aim here is to baseline correct AND invert the waveform
    std::transform(opDetWaveform.begin(),opDetWaveform.end(),locWaveform.begin(),[baseline](const auto& val){return baseline - val;});
    
    reco_tool::ICandidateHitFinder::HitCandidateVec      hitCandidateVec;
    reco_tool::ICandidateHitFinder::MergeHitCandidateVec mergedCandidateHitVec;
    
    fHitFinderTool->findHitCandidates(locWaveform, 0, 0, hitCandidateVec);
    fHitFinderTool->MergeHitCandidates(locWaveform, hitCandidateVec, mergedCandidateHitVec);

    return;
}

float OpHitFinder::getBaseline(const raw::OpDetWaveform& locWaveform) const
{
    float meanVal = std::accumulate(locWaveform.begin(),locWaveform.end(),0.) / float(locWaveform.size());
    
    // now fill a map to determine the most probable value
    std::map<raw::ADC_Count_t,int> adcFrequencyMap;
    
    raw::ADC_Count_t maxBin(0);
    int              maxCount(0);
    
    for(const auto& adc : locWaveform)
    {
        int& adcFrequency = adcFrequencyMap[adc];
        
        if (++adcFrequency > maxCount)
        {
            maxBin   = adc;
            maxCount = adcFrequency;
        }
    }
    
    if (std::abs(float(maxBin) - meanVal) > 5.) std::cout << "screw up" << std::endl;
    
    float mostProbableBaseline(0.);
    int   mostProbableCount(0);
    
    for(raw::ADC_Count_t adcBin = maxBin - 3; adcBin <= maxBin + 3; adcBin++)
    {
        try{
            mostProbableBaseline += adcFrequencyMap.at(adcBin);
            mostProbableCount++;
        }
        catch(...) {}
    }
    
    mostProbableBaseline /= mostProbableCount;
   
    return mostProbableBaseline;
}

    
void OpHitFinder::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "OpHitFinderPlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fOpHitFinderVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "OpHitFinderPlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "OpHitFinder;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fOpHitFinderVec.at(bin).Re());
    }
*/
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(OpHitFinder)
}
