////////////////////////////////////////////////////////////////////////
/// \file   Filter.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IFilter.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TF1.h"
#include "TComplex.h"

#include <fstream>

namespace icarus_tool
{

class Filter : IFilter
{
public:
    explicit Filter(const fhicl::ParameterSet& pset);
    
    ~Filter();
    
    void configure(const fhicl::ParameterSet& pset)                          override;
    void setResponse(size_t numBins, double correct3D, double timeScaleFctr) override;
    
    size_t                       getPlane()       const override {return fPlane;}
    const std::vector<TComplex>& getResponseVec() const override {return fFilterVec;}
    
private:
    // Member variables from the fhicl file
    size_t                fPlane;
    std::vector<double>   fParameters;
    std::string           fFunctionString;
    double                fFilterWidthCorrectionFactor;
    
    // The root version of the function
    TF1*                  fFunction;
    
    // Container for the field response "function"
    std::vector<TComplex> fFilterVec;
};
    
//----------------------------------------------------------------------
// Constructor.
Filter::Filter(const fhicl::ParameterSet& pset) :
    fFunction(0)
{
    configure(pset);
}
    
Filter::~Filter()
{
    if (fFunction) delete fFunction;
}
    
void Filter::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fPlane                       = pset.get<size_t>("Plane");
    fParameters                  = pset.get<std::vector<double>>("FilterParametersVec");
    fFunctionString              = pset.get<std::string>("FilterFunction");
    fFilterWidthCorrectionFactor = pset.get<double>("FilterWidthCorrectionFactor");

    std::string functionName = "Filter_plane" + std::to_string(fPlane) + "wr00";
    
    // Get rid of any existing function
    if (fFunction) delete fFunction;
    
    fFunction = new TF1(functionName.c_str(), fFunctionString.c_str());
    
    size_t idx(0);
    for(const auto& parameter : fParameters) fFunction->SetParameter(idx++, parameter);
    
    return;
}
    
void Filter::setResponse(size_t numBins, double correct3D, double timeScaleFctr)
{
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins2     = numBins / 2;
    
    // Set the range on the function
    fFunction->SetRange(0., double(numBins2));
    
    // now to scale the filter function!
    // only scale params 1,2 &3
    double timeFactor = 1. / (timeScaleFctr * correct3D * fFilterWidthCorrectionFactor);
    size_t paramIdx(0);
    
    while(paramIdx++ < 3)
        fFunction->SetParameter(paramIdx, timeFactor * fParameters[paramIdx]);

    // Now ready to set the response vector
    for(size_t bin = 0; bin <= size_t(numBins2); bin++)
    {
        // This takes a sampling rate in ns -> gives a frequency in cycles/us
        double freq = 500. * bin / (samplingRate * numBins2);
        double f = fFunction->Eval(freq);
        
        fFilterVec.push_back(TComplex(f, 0.));
    }
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(Filter)
}
