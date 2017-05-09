////////////////////////////////////////////////////////////////////////
/// \file   FieldResponse.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IFieldResponse.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "icaruscode/Utilities/tools/IWaveformTool.h"
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
    
    void configure(const fhicl::ParameterSet&)        override;
    void setResponse(double, double, double)          override;
    void outputHistograms(art::TFileDirectory&) const override;
    
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
    int binOfInterest = fFieldResponseHist->GetMinimumBin();
    
    // For collection planes it is as simple as finding the maximum bin
    if (fSignalType == geo::kCollection) binOfInterest = fFieldResponseHist->GetMaximumBin();
    
    // Do a backwards search to find the first positive bin
    while(1)
    {
        // Did we go too far?
        if (binOfInterest < 0)
            throw cet::exception("FieldResponse::configure") << "Cannot find zero-point crossover for induction response!" << std::endl;
            
        double content = fFieldResponseHist->GetBinContent(binOfInterest);
        
        if (content > 0.) break;
        
        binOfInterest--;
    }
    
    fT0Offset = -(fFieldResponseHist->GetXaxis()->GetBinCenter(binOfInterest) - fFieldResponseHist->GetXaxis()->GetBinCenter(1)) * fTimeCorrectionFactor;
    
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
    
void FieldResponse::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
    art::TFileDirectory dir = histDir.mkdir(fFieldResponseHistName.c_str());
    
    TAxis* xAxis = fFieldResponseHist->GetXaxis();
    
    TH1D* hist = dir.make<TH1D>(fFieldResponseHistName.c_str(), "Field Response; Time(us)", xAxis->GetNbins(), xAxis->GetXmin(), xAxis->GetXmax());
    
    double binWidth = (xAxis->GetXmax() - xAxis->GetXmin()) / double(xAxis->GetNbins());
    
    std::vector<double> histResponseVec(xAxis->GetNbins());
    
    for(int idx = 0; idx < xAxis->GetNbins(); idx++)
    {
        double xBin   = xAxis->GetXmin() + idx * binWidth;
        double binVal = fFieldResponseHist->GetBinContent(idx);
        
        hist->Fill(xBin, binVal);
        histResponseVec.at(idx) = binVal;
    }
    
    // Let's apply some smoothing as an experiment... first let's get the tool we need
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    std::unique_ptr<icarus_tool::IWaveformTool> waveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);

    // Make a copy of the response vec
    std::vector<double> smoothedResponseVec;
    
    // Run the triangulation smoothing
    waveformTool->triangleSmooth(histResponseVec, smoothedResponseVec);

    // Now make histogram of this
    std::string histName = "Smooth_" + fFieldResponseHistName;
    
    TH1D* smoothHist = dir.make<TH1D>(histName.c_str(), "Field Response; Time(us)", xAxis->GetNbins(), xAxis->GetXmin(), xAxis->GetXmax());
    
    for(size_t idx = 0; idx < smoothedResponseVec.size(); idx++)
    {
        double xBin = xAxis->GetXmin() + idx * binWidth;
        
        smoothHist->Fill(xBin,smoothedResponseVec.at(idx));
    }
    
    // Get the FFT of the response
    std::vector<double> powerVec;
    
    waveformTool->getFFTPower(histResponseVec, powerVec);
    
    // Now we can plot this...
    double maxFreq   = 0.5 / binWidth;   // binWidth will be in us, maxFreq will be units of MHz
    double freqWidth = maxFreq / powerVec.size();
    
    histName = "FFT_" + fFieldResponseHistName;
    
    TH1D* fftHist = dir.make<TH1D>(histName.c_str(), "Field Response FFT; Frequency(MHz)", powerVec.size(), 0., maxFreq);
    
    for(size_t idx = 0; idx < powerVec.size(); idx++)
    {
        double bin = (idx + 0.5) * freqWidth;
        
        fftHist->Fill(bin, powerVec.at(idx));
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
