////////////////////////////////////////////////////////////////////////
/// \file   OpHitFinder.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/Light/OpticalTools/IOpHitFinder.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "icaruscode/Utilities/SignalShapingServiceICARUS.h"

#include "TH1D.h"

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
    
    art::ServiceHandle<util::SignalShapingServiceICARUS> fSignalShaping;
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
    
    // Get signal shaping service.
    fSignalShaping = art::ServiceHandle<util::SignalShapingServiceICARUS>();
    
    return;
}

    
void OpHitFinder::FindOpHits(const raw::OpDetWaveform& opDetWaveform,
                             recob::OpHit&             opHit) const
{
    
    return;
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
