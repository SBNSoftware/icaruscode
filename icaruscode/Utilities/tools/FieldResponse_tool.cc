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
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
//#include "lardata/DetectorInfoServices/DetectorClocksService.h"
//#include "lardata/Utilities/LArFFT.h"
#include "TFile.h"
#include "TH1D.h"

#include <fstream>

namespace icarus_tool
{

class FieldResponse : IFieldResponse
{
public:
    explicit FieldResponse(const fhicl::ParameterSet& pset);
    
    ~FieldResponse() {}
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    size_t getPlane()             const;
    size_t getNumBins()           const;
    double getBinCenter(int bin)  const;
    double getBinContent(int bin) const;
    double getLowEdge()           const;
    double getHighEdge()          const;
    double getBinWidth()          const;
    double getTOffset()           const;
    double getIntegral()          const;
    double interpolate(double x)  const;
    
private:
    // Make sure we have been initialized
    bool           fIsValid;
    
    // Member variables from the fhicl file
    size_t         fThisPlane;
    geo::SigType_t fSignalType;
    std::string    fFieldResponseFileName;
    std::string    fFieldResponseFileVersion;
    std::string    fFieldResponseHistName;
    
    // Pointer to the input histogram
    TH1D*          fFieldResponseHist;
    
    // Derived variables
    double         fT0Offset;
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
    
    // Recover the input field response histogram
    std::string fileName = fFieldResponseFileName + Form("_vw%02i_", int(fThisPlane)) + fFieldResponseFileVersion + ".root";
    
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file(fileName, fullFileName);
    
    TFile inputFile(fullFileName.c_str(), "READ");
    
    if (!inputFile.IsOpen())
        throw cet::exception("FieldResponse::configure") << "Unable to open input file: " << fileName << std::endl;
    
    std::string histName = fFieldResponseHistName + Form("_vw%02i_", int(fThisPlane)) + fFieldResponseFileVersion + ".root";
    
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
    
    return fFieldResponseHist->GetBinWidth(1) * 1000.;
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
    

    
DEFINE_ART_CLASS_TOOL(FieldResponse)
}
