////////////////////////////////////////////////////////////////////////
/// \file   NoBaseline.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/TPC/SignalProcessing/RecoWire/DeconTools/IBaseline.h"
#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include <fstream>

namespace icarus_tool
{

class NoBaseline : public IBaseline
{
public:
    explicit NoBaseline(const fhicl::ParameterSet& pset);
    
    ~NoBaseline();
    
    void configure(const fhicl::ParameterSet& pset)                                       override;
    void outputHistograms(art::TFileDirectory&)                                     const override;
    
    float GetBaseline(const icarusutil::TimeVec&, raw::ChannelID_t, size_t, size_t) const override;
    
private:
};
    
//----------------------------------------------------------------------
// Constructor.
NoBaseline::NoBaseline(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
NoBaseline::~NoBaseline()
{
}
    
void NoBaseline::configure(const fhicl::ParameterSet& pset)
{
    return;
}

    
float NoBaseline::GetBaseline(const icarusutil::TimeVec& holder,
                              raw::ChannelID_t           channel,
                              size_t                     roiStart,
                              size_t                     roiLen) const
{
    return 0.;
}
    
void NoBaseline::outputHistograms(art::TFileDirectory& histDir) const
{
    return;
}
    
DEFINE_ART_CLASS_TOOL(NoBaseline)
}
