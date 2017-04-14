////////////////////////////////////////////////////////////////////////
/// \file   FieldResponse.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IFieldResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "TFile.h"
#include "TH1D.h"

#include <fstream>
#include <iomanip>

namespace icarus_tool
{

class FieldResponse : IFieldResponse
{
public:
    explicit FieldResponse(const fhicl::ParameterSet& pset);
    
    ~FieldResponse() {}
    
    void configure(const fhicl::ParameterSet& pset)                         override;
    void setResponse(double weight, double correct3D, double timeScaleFctr) override;
    
    size_t                     getPlane()             const override;
    size_t                     getNumBins()           const override;
    double                     getBinCenter(int bin)  const override;
    double                     getBinContent(int bin) const override;
    double                     getLowEdge()           const override;
    double                     getHighEdge()          const override;
    double                     getBinWidth()          const override;
    double                     getTOffset()           const override;
    double                     getIntegral()          const override;
    double                     interpolate(double x)  const override;
    
    const std::vector<double>& getResponseVec()       const override {return fFieldResponseVec;}
    
private:
    // Utility routine for converting numbers to strings
    std::string         numberToString(int number);
    
    // Make sure we have been initialized
    bool                fIsValid;
    
    // Member variables from the fhicl file
    size_t              fThisPlane;
    geo::SigType_t      fSignalType;
    std::string         fFieldResponseFileName;
    std::string         fFieldResponseFileVersion;
    std::string         fFieldResponseHistName;
    double              fFieldResponseAmplitude;
    double              fTimeCorrectionFactor;
    
    // Pointer to the input histogram
    TH1D*               fFieldResponseHist;
    
    // Container for the field response "function"
    std::vector<double> fFieldResponseVec;
    
    // Derived variables
    double              fT0Offset;
};
    
//----------------------------------------------------------------------
// Constructor.
FieldResponse::FieldResponse(const fhicl::ParameterSet& pset) :
    fIsValid(false), fSignalType(geo::kMysteryType)
{
    configure(pset);
}
    
void FieldResponse::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fThisPlane                = pset.get<size_t>("Plane");
    fSignalType               = pset.get<size_t>("SignalType") == 0 ? geo::kInduction : geo::kCollection;
    fFieldResponseFileName    = pset.get<std::string>("FieldResponseFileName");
    fFieldResponseFileVersion = pset.get<std::string>("FieldResponseFileVersion");
    fFieldResponseHistName    = pset.get<std::string>("FieldResponseHistName");
    fFieldResponseAmplitude   = pset.get<double>("FieldResponseAmplitude");
    fTimeCorrectionFactor     = pset.get<double>("TimeCorrectionFactor");
    
    // Recover the input field response histogram
    std::string fileName = fFieldResponseFileName + "_vw" + numberToString(fThisPlane) + "_" + fFieldResponseFileVersion + ".root";
    
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file(fileName, fullFileName);
    
    TFile inputFile(fullFileName.c_str(), "READ");
    
    if (!inputFile.IsOpen())
        throw cet::exception("FieldResponse::configure") << "Unable to open input file: " << fileName << std::endl;
    
    std::string histName = fFieldResponseHistName + "_vw" + numberToString(fThisPlane) + "_" + fFieldResponseFileVersion + ".root";
    
    fFieldResponseHist = (TH1D*)inputFile.Get(histName.c_str());
    
    // Calculation of the T0 offset depends on the signal type
    // For collection planes it is as simple as finding the maximum bin
    if (fSignalType == geo::kCollection)
    {
        int binWithMaxValue = fFieldResponseHist->GetMaximumBin();
        
        fT0Offset = fFieldResponseHist->GetXaxis()->GetBinCenter(binWithMaxValue) - fFieldResponseHist->GetXaxis()->GetBinCenter(1);
    }
    else
    {
        // The technique below from Leon Rochester
        int binWithMinimumValue = fFieldResponseHist->GetMinimumBin();
        
        // Do a search backward to find the zero point cross over
        while(1)
        {
            // Did we go too far?
            if (binWithMinimumValue < 0)
                throw cet::exception("FieldResponse::configure") << "Cannot find zero-point crossover for induction response!" << std::endl;
                
            double content = fFieldResponseHist->GetBinContent(binWithMinimumValue);
            
            if (content > 0.) break;
            
            binWithMinimumValue--;
        }
        
        fT0Offset = fFieldResponseHist->GetXaxis()->GetBinCenter(binWithMinimumValue) - fFieldResponseHist->GetXaxis()->GetBinCenter(1);
    }
    
    fIsValid = true;
    
    return;
}
    
void FieldResponse::setResponse(double weight, double correction3D, double timeScaleFctr)
{
    double timeFactor    = correction3D * timeScaleFctr;
    size_t numBins       = getNumBins();
    size_t nResponseBins = numBins * timeFactor;
    
    fFieldResponseVec.resize(nResponseBins);

    double x0     = getBinCenter(1);
    double xf     = getBinCenter(numBins);
    double deltaX = (xf - x0)/(numBins-1);
    
    for(size_t bin = 1; bin <= nResponseBins; bin++)
    {
        double xVal = x0 + deltaX * (bin-1) / timeFactor;
        
        fFieldResponseVec.at(bin-1) = interpolate(xVal) * fFieldResponseAmplitude * weight;
    }
    
    return;
}
    
size_t FieldResponse::getPlane() const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fThisPlane;
}
    
size_t FieldResponse::getNumBins() const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fFieldResponseHist->GetXaxis()->GetNbins();
}
    
double FieldResponse::getLowEdge() const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fFieldResponseHist->GetBinCenter(1) - 0.5 * fFieldResponseHist->GetBinWidth(1);
}
    
double FieldResponse::getBinCenter(int bin) const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fFieldResponseHist->GetBinCenter(bin);
}
    
double FieldResponse::getBinContent(int bin) const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fFieldResponseHist->GetBinContent(bin);
}
    
double FieldResponse::getHighEdge() const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    size_t nBins = getNumBins();
    
    return fFieldResponseHist->GetBinCenter(nBins) + 0.5 * fFieldResponseHist->GetBinWidth(nBins);
}
    
double FieldResponse::getBinWidth() const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fFieldResponseHist->GetBinWidth(1) * fTimeCorrectionFactor;
}
    
double FieldResponse::getTOffset() const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fT0Offset;
}
    
double FieldResponse::getIntegral() const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fFieldResponseHist->Integral();
}
    
double FieldResponse::interpolate(double x) const
{
    if (!fIsValid)
        throw cet::exception("FieldResponse::getPlane") << "Attempting to access plane info when tool is invalid state" << std::endl;
    
    return fFieldResponseHist->Interpolate(x);
}

std::string FieldResponse::numberToString(int number)
{
    std::ostringstream string;
    
    string << std::setfill('0') << std::setw(2) << number;
    
    return string.str();
}

    
DEFINE_ART_CLASS_TOOL(FieldResponse)
}
