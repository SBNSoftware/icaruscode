////////////////////////////////////////////////////////////////////////
/// \file   Filter.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IFilter.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TF1.h"
#include "TProfile.h"
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
    void outputHistograms(art::TFileDirectory&)                        const override;
    
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
    // Note that here we are working in frequency space, not in the time domain...
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = 1.e-3 * detprop->SamplingRate(); // Note sampling rate is in ns, convert to us
    double      maxFreq      = 1. / (2. * samplingRate);        // highest frequency in cycles/us
    double      freqRes      = maxFreq / double(numBins/2);     // frequency resolution in cycles/us

    // Set the range on the function
    fFunction->SetRange(0, double(numBins/2));

    // now to scale the filter function!
    // only scale params 1,2 &3
    double timeFactor = 1. / (timeScaleFctr * correct3D * fFilterWidthCorrectionFactor);
    size_t paramIdx(0);

    for(const auto& parameter : fParameters) fFunction->SetParameter(paramIdx++, timeFactor * parameter);
    
    // Don't assume that the filter vec has not already been initialized...
    fFilterVec.clear();
    
    // Keep track of the peak value
    double peakVal(std::numeric_limits<double>::min());

    // Now ready to set the response vector
    for(size_t bin = 0; bin <= numBins/2 + 1; bin++)
    {
        // This takes a sampling rate in ns -> gives a frequency in cycles/us
        double freq = bin * freqRes;

        double f = fFunction->Eval(freq);
        
        peakVal = std::max(peakVal, f);
        
        fFilterVec.push_back(TComplex(f, 0.));
    }
    
    // "Normalize" to peak value
    for(auto& filterValue : fFilterVec) filterValue = filterValue / peakVal;

    return;
}
    
void Filter::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
    std::string dirName = "FilterPlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      numBins      = fFilterVec.size();
    double      samplingRate = detprop->SamplingRate(); // Sampling time in us
    double      maxFreq      = 1.e6 / (2. * samplingRate);
    double      minFreq      = maxFreq / numBins;
    std::string histName     = "FilterPlane_" + std::to_string(fPlane);
    
    TProfile*   hist         = dir.make<TProfile>(histName.c_str(), "Filter;Frequency(MHz)", numBins, minFreq, maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = bin * minFreq;
        
        hist->Fill(freq, fFilterVec.at(bin).Re(), 1.);
    }
    
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(Filter)
}
