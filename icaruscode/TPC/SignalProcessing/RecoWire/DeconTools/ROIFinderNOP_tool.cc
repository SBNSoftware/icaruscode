////////////////////////////////////////////////////////////////////////
/// \file   ROIFinder.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IROIFinder.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <fstream>

namespace icarus_tool
{

class ROIFinderNOP : public IROIFinder
{
public:
    explicit ROIFinderNOP(const fhicl::ParameterSet& pset);
    
    ~ROIFinderNOP();
    
    void   configure(const fhicl::ParameterSet& pset)                              override;
    void   initializeHistograms(art::TFileDirectory&)                        const override;
    size_t plane()                                                           const override {return fPlane;}

    void FindROIs(const Waveform&, size_t, size_t, double, CandidateROIVec&) const override;
    
private:
    // Member variables from the fhicl file
    size_t fPlane;                      ///< Tool can be plane dependent
};
    
//----------------------------------------------------------------------
// Constructor.
ROIFinderNOP::ROIFinderNOP(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ROIFinderNOP::~ROIFinderNOP()
{
}
    
void ROIFinderNOP::configure(const fhicl::ParameterSet& pset)
{    
    fPlane = pset.get< size_t >("Plane");

    return;
}
    
void ROIFinderNOP::FindROIs(const Waveform& waveform, size_t channel, size_t cnt, double rmsNoise, CandidateROIVec& roiVec) const
{
    // We just make the entire waveform an "ROI"
    roiVec.clear();
    roiVec.emplace_back(size_t(0),waveform.size()-1);
    
    return;
}

void ROIFinderNOP::initializeHistograms(art::TFileDirectory& histDir) const
{    
    return;
}
    
DEFINE_ART_CLASS_TOOL(ROIFinderNOP)
}
