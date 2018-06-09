////////////////////////////////////////////////////////////////////////
/// \file   PeakFitterICARUS.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "icaruscode/HitFinder/HitFinderTools/IPeakFitter.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include <cmath>
#include <fstream>
#include "TH1F.h"
#include "TF1.h"

namespace reco_tool
{

class PeakFitterICARUS : IPeakFitter
{
public:
    explicit PeakFitterICARUS(const fhicl::ParameterSet& pset);
    
    ~PeakFitterICARUS();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    void findPeakParameters(const std::vector<float>&,
                            const ICandidateHitFinder::HitCandidateVec&,
                            PeakParamsVec&,
                            double&,
                            int&) const override;
    
     static Double_t fitf(Double_t *x, Double_t *par);
    
private:
    // Member variables from the fhicl file
    double                   fMinWidth;     ///< minimum initial width for ICARUS fit
    double                   fMaxWidthMult; ///< multiplier for max width for ICARUS fit
    double                   fPeakRange;    ///< set range limits for peak center
    double                   fAmpRange;     ///< set range limit for peak amplitude
    
    mutable TH1F             fHistogram;
    
    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
};
    
//----------------------------------------------------------------------
// Constructor.
PeakFitterICARUS::PeakFitterICARUS(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
PeakFitterICARUS::~PeakFitterICARUS()
{
}
    
void PeakFitterICARUS::configure(const fhicl::ParameterSet& pset)
{
    // Start by recovering the parameters
    fMinWidth     = pset.get<double>("MinWidth",      0.5);
    fMaxWidthMult = pset.get<double>("MaxWidthMult",  3.);
    fPeakRange    = pset.get<double>("PeakRangeFact", 2.);
    fAmpRange     = pset.get<double>("PeakAmpRange",  2.);
    
    fHistogram    = TH1F("PeakFitterHitSignal","",500,0.,500.);
    
    fHistogram.Sumw2();
    
    std::string function = "Gaus(0)";
    
    return;
}
    
// --------------------------------------------------------------------------------------------
void PeakFitterICARUS::findPeakParameters(const std::vector<float>&                   roiSignalVec,
                                            const ICandidateHitFinder::HitCandidateVec& hitCandidateVec,
                                            PeakParamsVec&                              peakParamsVec,
                                            double&                                     chi2PerNDF,
                                            int&                                        NDF) const
{
    // The following is a translation of the original FitICARUSs function in the original
    // GausHitFinder module originally authored by Jonathan Asaadi
    //
    // *** NOTE: this algorithm assumes the reference time for input hit candidates is to
    //           the first tick of the input waveform (ie 0)
    //
    if (hitCandidateVec.empty()) return;
    
    // in case of a fit failure, set the chi-square to infinity
    chi2PerNDF = std::numeric_limits<double>::infinity();
    
    int startTime = hitCandidateVec.front().startTick;
    int endTime   = hitCandidateVec.back().stopTick;
    int roiSize   = endTime - startTime;
    
    // Check to see if we need a bigger histogram for fitting
    if (roiSize > fHistogram.GetNbinsX())
    {
        std::string histName = "PeakFitterHitSignal_" + std::to_string(roiSize);
        fHistogram = TH1F(histName.c_str(),"",roiSize,0.,roiSize);
        fHistogram.Sumw2();
    }
    
    for(int idx = 0; idx < roiSize; idx++) fHistogram.SetBinContent(idx+1,roiSignalVec.at(startTime+idx));
    
   // for(int idx = 0; idx < roiSize; idx++)
     //   std::cout << " bin " << idx << " fHistogram " << fHistogram.GetBinContent(idx+1) << std::endl;

    
    // Build the string to describe the fit formula
    std::string equation = "gaus(0)";
    

    // Now define the complete function to fit
    // Now define the complete function to fit
    TF1 Func("ICARUSfunc",fitf,roiSize,0,roiSize);
    
    // ### Setting the parameters for the ICARUS Fit ###
    //int parIdx(0);
    for(auto& candidateHit : hitCandidateVec)
    {
        double peakMean   = candidateHit.hitCenter - float(startTime);
        double peakWidth  = candidateHit.hitSigma;
       // std::cout << " peakWidth " << peakWidth << std::endl;
       // std::cout << " hitcenter " << candidateHit.hitCenter << " starttime " << startTime <<std::endl;
        

        double amplitude  = candidateHit.hitHeight;
        double meanLowLim = std::max(peakMean - fPeakRange * peakWidth,              0.);
        double meanHiLim  = std::min(peakMean + fPeakRange * peakWidth, double(roiSize));
      //  std::cout << " amplitude " << amplitude << std::endl;
      //  std::cout << " peakMean " << peakMean << std::endl;

        Func.SetParameter(0,0);
        Func.SetParameter(1, amplitude);
        Func.SetParameter(2, peakMean);
        Func.SetParameter(3,peakWidth/2);
        Func.SetParameter(4,peakWidth/2);
        
        Func.SetParLimits(0, -5, 5);
        Func.SetParLimits(1, 0.1 * amplitude,  10. * amplitude);
        Func.SetParLimits(2, meanLowLim,meanHiLim);
        Func.SetParLimits(3, std::max(fMinWidth, 0.1 * peakWidth), fMaxWidthMult * peakWidth);
        Func.SetParLimits(4, std::max(fMinWidth, 0.1 * peakWidth), fMaxWidthMult * peakWidth);
        
  //      std::cout << " peakWidth limits " << std::max(fMinWidth, 0.1 * peakWidth) <<" " <<  fMaxWidthMult * peakWidth << std::endl;
    //      std::cout << " peakMean limits " << meanLowLim <<" " <<  meanHiLim << std::endl;
        
    }
    
    int fitResult(-1);
    
    try
    { fitResult = fHistogram.Fit(&Func,"QNRWB","", 0., roiSize);}
    catch(...)
    {mf::LogWarning("GausHitFinder") << "Fitter failed finding a hit";}
    
   if(fitResult!=0)
       std::cout << " fit cannot converge " << std::endl;
        // ##################################################
        // ### Getting the fitted parameters from the fit ###
        // ##################################################
     NDF        = roiSize-5;
        chi2PerNDF = (Func.GetChisquare() / NDF);
    
 //   std::cout << " chi2 " << Func.GetChisquare() << std::endl;
  //  std::cout << " ndf " << NDF << std::endl;

    //      std::cout << " chi2ndf " << chi2PerNDF<< std::endl;
        //parIdx = 0;
        for(size_t idx = 0; idx < hitCandidateVec.size(); idx++)
        {
            PeakFitParams_t peakParams;
            
            peakParams.peakAmplitude      = Func.GetParameter(1);
            peakParams.peakAmplitudeError = Func.GetParError(1);
            peakParams.peakCenter         = Func.GetParameter(2) + float(startTime);
            peakParams.peakCenterError    = Func.GetParError(2);
    //std::cout << " rising time " << Func.GetParameter(3) << " falling time " <<Func.GetParameter(4) << std::endl;
            peakParams.peakTauLeft        = Func.GetParameter(3);
            peakParams.peakTauLeftError   = Func.GetParError(3);
            peakParams.peakTauRight       = Func.GetParameter(4);
            peakParams.peakTauRightError  = Func.GetParError(4);
            peakParams.peakBaseline       = Func.GetParameter(0);
            peakParams.peakBaselineError  = Func.GetParError(0);
            
            peakParamsVec.emplace_back(peakParams);
            
    }
    
    //Gaus.Delete();
    Func.Delete();
    return;
}
    Double_t PeakFitterICARUS::fitf(Double_t *x, Double_t *par)
    {
        // Double_t arg = 0;
        
        Double_t fitval = par[0]+par[1]*TMath::Exp(-(x[0]-par[2])/par[3])/(1+TMath::Exp(-(x[0]-par[3])/par[4]));
        return fitval;
    }
DEFINE_ART_CLASS_TOOL(PeakFitterICARUS)
}
