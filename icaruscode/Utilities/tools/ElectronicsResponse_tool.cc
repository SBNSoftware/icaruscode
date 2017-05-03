////////////////////////////////////////////////////////////////////////
/// \file   ElectronicsResponse.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IElectronicsResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TH1D.h"

#include <fstream>

namespace icarus_tool
{

class ElectronicsResponse : IElectronicsResponse
{
public:
    explicit ElectronicsResponse(const fhicl::ParameterSet& pset);
    
    ~ElectronicsResponse() {}
    
    void configure(const fhicl::ParameterSet& pset)         override;
    void setResponse(size_t numBins, double binWidth)       override;
    void outputHistograms(art::TFileDirectory&)       const override;
    
    size_t                     getPlane()           const override {return fPlane;}
    double                     getFCperADCMicroS()  const override {return fFCperADCMicroS;}
    double                     getASICShapingTime() const override {return fASICShapingTime;}
    const std::vector<double>& getResponseVec()     const override {return fElectronicsResponseVec;}
    
private:
    // Member variables from the fhicl file
    size_t              fPlane;
    double              fFCperADCMicroS;
    double              fASICShapingTime;
    double              fADCPerPCAtLowestASICGain;
    
    // Container for the electronics response "function"
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
    fFCperADCMicroS           = pset.get<double>("FCperADCMicroS");
    fASICShapingTime          = pset.get<double>("ASICShapingTime");
    fADCPerPCAtLowestASICGain = pset.get<double>("ADCPerPCAtLowestASICGain");
    
    return;
}
    
void ElectronicsResponse::setResponse(size_t numBins, double binWidth)
{
//    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
//    double      samplingRate = detprop->SamplingRate();
    double      timeCorrect  = 1.e-3;
    
    binWidth *= timeCorrect;
    
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
        double time = double(timeIdx) * binWidth;
        
        fElectronicsResponseVec.at(timeIdx) = time / fASICShapingTime * exp(-time / fASICShapingTime);
    }
    
//    double maxValue = *std::max_element(fElectronicsResponseVec.begin(),fElectronicsResponseVec.end());
    
    // normalize fElectResponse[i], before the convolution
    // Put in overall normalization in a pedantic way:
    // first put in the pulse area per eleectron at the lowest gain setting,
    // then normalize by the actual ASIC gain setting used.
    // This code is executed only during initialization of service,
    // so don't worry about code inefficiencies here.
    //    double last_integral=0;
    //    double last_max=0;
    
    // ICARUS Normalization are the following
    // Field response is normalized to 1 electron. Shaping function written as (t/tau)*exp(-t/tau) is normalized to 1
    // From test pulse measurement with FLIC@CERN we have 0.027 fC/(ADC*us)
    // Therefore 0.027*6242 electrons/(ADC*us)
    
    for (auto& element : fElectronicsResponseVec)
        element /= (fFCperADCMicroS * 6242.);
    
    return;
}
    
void ElectronicsResponse::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
    std::string dirName = "ElectronicsPlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    std::string histName = "ElectronicsResponse_" + std::to_string(fPlane);
    
    TH1D* hist = dir.make<TH1D>(histName.c_str(), "Response", fElectronicsResponseVec.size(), 0., fElectronicsResponseVec.size());
    
    for(size_t idx = 0; idx < fElectronicsResponseVec.size(); idx++) hist->Fill(idx, fElectronicsResponseVec.at(idx));
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(ElectronicsResponse)
}
