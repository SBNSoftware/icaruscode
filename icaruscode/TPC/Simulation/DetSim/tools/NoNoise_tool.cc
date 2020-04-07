////////////////////////////////////////////////////////////////////////
/// \file   NoNoise.cc
/// \author F. Varanini
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IGenNoise.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

namespace icarus_tool
{

class NoNoise : IGenNoise
{
public:
    explicit NoNoise(const fhicl::ParameterSet& pset);
    
    ~NoNoise();
    
    void configure(const fhicl::ParameterSet& pset)               override;
    
    void nextEvent() override  {return;};

    void generateNoise(CLHEP::HepRandomEngine&,
                       CLHEP::HepRandomEngine&,
                       icarusutil::TimeVec&,
                       detinfo::DetectorPropertiesData const&,
                       double, unsigned int) override;
    
private:

};
    
//----------------------------------------------------------------------
// Constructor.
NoNoise::NoNoise(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
NoNoise::~NoNoise()
{
}
    
void NoNoise::configure(const fhicl::ParameterSet& pset)
{
    // Nothing to do here
    return;
}

void NoNoise::generateNoise(CLHEP::HepRandomEngine&,
                            CLHEP::HepRandomEngine&,
                            icarusutil::TimeVec& noise,
                            detinfo::DetectorPropertiesData const&,
                            double noise_factor,
                            unsigned int channel)
{
    // Set all values to 0
    std::fill(noise.begin(), noise.end(), 0.);

    return;
}
    
DEFINE_ART_CLASS_TOOL(NoNoise)
}
