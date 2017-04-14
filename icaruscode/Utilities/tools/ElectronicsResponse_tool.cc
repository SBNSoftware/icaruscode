////////////////////////////////////////////////////////////////////////
/// \file   ElectronicsResponse.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IElectronicsResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <fstream>

namespace icarus_tool
{

class ElectronicsResponse : IElectronicsResponse
{
public:
    explicit ElectronicsResponse(const fhicl::ParameterSet& pset);
    
    ~ElectronicsResponse() {}
    
    void configure(const fhicl::ParameterSet& pset)   override;
    void setResponse(size_t numBins, double binWidth) override;
    
    size_t                     getPlane()           const override {return fPlane;}
    double                     getASICGain()        const override {return fASICGain;}
    double                     getASICShapingTime() const override {return fASICShapingTime;}
    const std::vector<double>& getResponseVec()     const override {return fElectronicsResponseVec;}
    
private:
    // Member variables from the fhicl file
    size_t              fPlane;
    double              fASICGain;
    double              fASICShapingTime;
    double              fADCPerPCAtLowestASICGain;
    
    // Container for the field response "function"
    std::vector<double> fElectronicsResponseVec;
};
    
//----------------------------------------------------------------------
// Constructor.
ElectronicsResponse::ElectronicsResponse(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
void ElectronicsResponse::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fPlane                    = pset.get<size_t>("Plane");
    fASICGain                 = pset.get<double>("ASICGainInMVPerFC");
    fASICShapingTime          = pset.get<double>("ASICShapingTime");
    fADCPerPCAtLowestASICGain = pset.get<double>("ADCPerPCAtLowestASICGain");
    
    return;
}
    
void ElectronicsResponse::setResponse(size_t numBins, double binWidth)
{
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    
    // what is this number supposed to be?
    samplingRate = 1000.;
    
    fElectronicsResponseVec.resize(numBins, 0.);
    
    // This note from Filippo:
    // The following sets the ICARUS electronics response function in
    // time-space. Function comes from BNL SPICE simulation of ICARUS
    // electronics. SPICE gives the electronics transfer function in
    // frequency-space. The inverse laplace transform of that function
    // (in time-space) was calculated in Mathematica and is what is being
    // used below. Parameters Ao and To are cumulative gain/timing parameters
    // from the full (ASIC->Intermediate amp->Receiver->ADC) electronics chain.
    // They have been adjusted to make the SPICE simulation to match the
    // actual electronics response. Default params are Ao=1.4, To=0.5us.
    
    for(size_t timeIdx = 0; timeIdx < numBins; timeIdx++)
    {
        double time = double(timeIdx) * binWidth / samplingRate;
        
        fElectronicsResponseVec.at(timeIdx) = time / fASICShapingTime * exp(-time / fASICShapingTime);
    }
    
    double maxValue = *std::max_element(fElectronicsResponseVec.begin(),fElectronicsResponseVec.end());
    
    // normalize fElectResponse[i], before the convolution
    // Put in overall normalization in a pedantic way:
    // first put in the pulse area per eleectron at the lowest gain setting,
    // then normalize by the actual ASIC gain setting used.
    // This code is executed only during initialization of service,
    // so don't worry about code inefficiencies here.
    //    double last_integral=0;
    //    double last_max=0;
    
    //Normalization are the following
    // Peak is firstly normalized to 1
    // thus we expect peak to be 1 * 9390 (fADCPerPCtAtLowestAsicGain) * 1.602e-7 * (1 fC) = 9.39 ADC
    // At 4.7 mV/fC, the ADC value should be 4.7 (mV/fC) * 2 (ADC/mV) ~ 9.4 ADC/fC
    // so the normalization are consistent
    
    for (auto& element : fElectronicsResponseVec)
        element *= (fASICGain / 6.5) * (fADCPerPCAtLowestASICGain * 1.60217657e-7) / maxValue;
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(ElectronicsResponse)
}
