////////////////////////////////////////////////////////////////////////
/// \file   Baseline.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IBaseline.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "icaruscode/Utilities/SignalShapingICARUSService_service.h"

#include "TH1D.h"

#include <fstream>

namespace icarus_tool
{

class BaselineStandard : public IBaseline
{
public:
    explicit BaselineStandard(const fhicl::ParameterSet& pset);
    
    ~BaselineStandard();
    
    void configure(const fhicl::ParameterSet& pset)                                       override;
    void outputHistograms(art::TFileDirectory&)                                     const override;
    
    float GetBaseline(icarusutil::TimeVec const&, raw::ChannelID_t, size_t, size_t) const override;
    
private:
    // fhicl parameters
    int    fNumBinsToAverage;
    
    art::ServiceHandle<icarusutil::SignalShapingICARUSService> fSignalShaping;
};
    
//----------------------------------------------------------------------
// Constructor.
BaselineStandard::BaselineStandard(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
BaselineStandard::~BaselineStandard()
{
}
    
void BaselineStandard::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fNumBinsToAverage = pset.get<int>("NumBinsToAverage", 20);
    
    // Get signal shaping service.
    fSignalShaping = art::ServiceHandle<icarusutil::SignalShapingICARUSService>();
    
    return;
}

    
float BaselineStandard::GetBaseline(icarusutil::TimeVec const& holder,
                                    raw::ChannelID_t           channel,
                                    size_t                     roiStart,
                                    size_t                     roiLen) const
{
    float base=0;
    //1. Check Baseline match?
    // If not, include next ROI(if none, go to the end of signal)
    // If yes, proceed
    size_t nBinsToAve(fNumBinsToAverage);
    size_t roiStop(roiStart + roiLen);
    size_t newRoiStart(roiStart);
    size_t newRoiStop(roiStop);
    size_t nTries(4);
    
    // Calculate baslines from the very front of the deconvolution buffer and from the end
    float  basePre  = std::accumulate(holder.begin() + roiStart, holder.begin() + roiStart + nBinsToAve, 0.) / float(nBinsToAve);
    float  basePost = std::accumulate(holder.begin() + roiStop - nBinsToAve, holder.begin() + roiStop,   0.) / float(nBinsToAve);
    
    // emulate method for refining baseline from original version of CalWireROI
    float deconNoise = 1.26491 * fSignalShaping->GetDeconNoise(channel);    // 4./sqrt(10) * noise
    
    // If the estimated baseline from the front of the roi does not agree well with that from the end
    // of the roi then we'll extend the roi hoping for good agreement
    while(!(fabs(basePre - basePost) < deconNoise) && nTries++ < 3)
    {
        size_t nBinsToAdd(10);
        
        if (newRoiStart < nBinsToAdd)                newRoiStart  = 0;
        else                                         newRoiStart -= nBinsToAdd;
        if (newRoiStop + nBinsToAdd > holder.size()) newRoiStop   = holder.size();
        else                                         newRoiStop  += nBinsToAdd;
        
        basePre  = std::accumulate(holder.begin() + newRoiStart, holder.begin() + newRoiStart + nBinsToAve, 0.) / float(nBinsToAve);
        basePost = std::accumulate(holder.begin() + newRoiStop - nBinsToAve, holder.begin() + newRoiStop,   0.) / float(nBinsToAve);
    }
    
    // get spread in "pre" baseline
    float maxPre = *std::max_element(holder.begin() + roiStart, holder.begin() + roiStart + nBinsToAve);
    float minPre = *std::min_element(holder.begin() + roiStart, holder.begin() + roiStart + nBinsToAve);
    
    // Basically, we are hoping for a relatively smooth "front porch" for the waveform
    if ((maxPre - minPre) < 4.) base = basePre;
    else
    {
        // No success... see if the "back porch" is smooth
        size_t roiStop = roiStart + roiLen;
        
        float maxPost = *std::max_element(holder.begin() + roiStop - nBinsToAve, holder.begin() + roiStop);
        float minPost = *std::min_element(holder.begin() + roiStop - nBinsToAve, holder.begin() + roiStop);
        
        if ((maxPost - minPost) < 4.) base = basePost;
        // Starting to get desparate...
        else if ((maxPre  - minPre)  < 8. && std::fabs(basePre)  < std::fabs(basePost)) base = basePre;
        else if ((maxPost - minPost) < 8. && std::fabs(basePost) < std::fabs(basePre))  base = basePost;
        // Ok, apply brute force
        else
        {
            float min = *std::min_element(holder.begin()+roiStart,holder.begin()+roiStart+roiLen);
            float max = *std::max_element(holder.begin()+roiStart,holder.begin()+roiStart+roiLen);
            int nbin = std::ceil(max - min);
            if (nbin > 0){
                TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
                for (unsigned int bin = roiStart; bin < roiStart+roiLen; bin++){
                    h1->Fill(holder[bin]);
                }
                int   pedBin = h1->GetMaximumBin();
                float ped    = h1->GetBinCenter(pedBin);
                float ave=0,ncount = 0;
                for (unsigned int bin = roiStart; bin < roiStart+roiLen; bin++){
                    if (fabs(holder[bin]-ped)<2){
                        ave +=holder[bin];
                        ncount ++;
                    }
                }
                if (ncount==0) ncount=1;
                ave = ave/ncount;
                h1->Delete();
                base = ave;
            }
        }
    }
    
    return base;
}
    
void BaselineStandard::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "BaselinePlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fBaselineVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "BaselinePlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "Baseline;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fBaselineVec.at(bin).Re());
    }
*/
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(BaselineStandard)
}
