////////////////////////////////////////////////////////////////////////
/// \file   Filter.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IFilter.h"
#include "art/Utilities/ToolMacros.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TF1.h"
#include "TProfile.h"

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
    
    size_t                          getPlane()                         const override {return fPlane;}
    const icarusutil::FrequencyVec& getResponseVec()                   const override {return fFilterVec;}
    
private:
    // Member variables from the fhicl file
    size_t                   fPlane;
    std::vector<double>      fParameters;
    std::string              fFunctionString;
    double                   fFilterWidthCorrectionFactor;
    
    // Container for the field response "function"
    icarusutil::FrequencyVec fFilterVec;
};
    
//----------------------------------------------------------------------
// Constructor.
Filter::Filter(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
Filter::~Filter()
{
    return;
}
    
void Filter::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fPlane                       = pset.get<size_t>("Plane");
    fParameters                  = pset.get<std::vector<double>>("FilterParametersVec");
    fFunctionString              = pset.get<std::string>("FilterFunction");
    fFilterWidthCorrectionFactor = pset.get<double>("FilterWidthCorrectionFactor");
    
    return;
}
    
void Filter::setResponse(size_t numBins, double correct3D, double timeScaleFctr)
{
    // Note that here we are working in frequency space, not in the time domain...
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = 1.e-3 * detprop->SamplingRate(); // Note sampling rate is in ns, convert to us
    double      maxFreq      = 1.e3 / (2. * samplingRate);      // highest frequency in cycles/us (MHz)
    double      freqRes      = maxFreq / double(numBins/2);     // frequency resolution in cycles/us
    
    std::string funcName = "tempFilter";
    TF1 function(funcName.c_str(),fFunctionString.c_str());

    // Set the range on the function, probably not really necessary
    function.SetRange(0, maxFreq);

    // now to scale the filter function!
    // only scale params 1,2 &3
    double timeFactor = 1. / (timeScaleFctr * correct3D * fFilterWidthCorrectionFactor);
    size_t paramIdx(0);

    for(const auto& parameter : fParameters) function.SetParameter(paramIdx++, timeFactor * parameter);
    
    // Don't assume that the filter vec has not already been initialized...
    // Note we are setting up the FFT to be the "full" so add a bin for folding over
    fFilterVec.resize(numBins+1,icarusutil::ComplexVal(0.,0.));
    
    // Keep track of the peak value
    double peakVal(std::numeric_limits<double>::min());

    // Now ready to set the response vector
    for(size_t bin = 0; bin < numBins/2 + 1; bin++)
    {
        // This takes a sampling rate in ns -> gives a frequency in cycles/us
        double freq = bin * freqRes;
        double f    = function.Eval(freq);

        peakVal = std::max(peakVal, f);
        
        fFilterVec[bin] = icarusutil::ComplexVal(f, 0.);
    }

    // Since we are returning the "full" FFT... folder over the first half
    for(size_t bin = 0; bin < numBins/2; bin++)
        fFilterVec[numBins-bin] = fFilterVec[bin];
    
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
    double      samplingRate = 1.e-3 * detprop->SamplingRate(); // Sampling time in us
    double      maxFreq      = 1.e3 / (2. * samplingRate);      // Max frequency in MHz
    double      minFreq      = maxFreq / numBins;
    std::string histName     = "FilterPlane_" + std::to_string(fPlane);
    TProfile*   hist         = dir.make<TProfile>(histName.c_str(), "Filter;Frequency(MHz)", numBins, minFreq, maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = bin * minFreq;
        
        hist->Fill(freq, std::abs(fFilterVec.at(bin)), 1.);
    }
    
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(Filter)
}
