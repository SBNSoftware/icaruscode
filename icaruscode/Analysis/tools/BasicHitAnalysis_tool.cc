
#include "icaruscode/Analysis/tools/IHitHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <cmath>
#include <algorithm>

namespace BasicHitAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       BasicHitAnalysis
    // Module Type: producer
    // File:        BasicHitAnalysis.h
    //
    //              The intent of this module is to provide methods for
    //              "analyzing" hits on waveforms
    //
    // Configuration parameters:
    //
    // TruncMeanFraction     - the fraction of waveform bins to discard when
    //
    // Created by Tracy Usher (usher@slac.stanford.edu) on February 19, 2016
    //
    ////////////////////////////////////////////////////////////////////////
    
// The following typedefs will, obviously, be useful
using HitPtrVec       = std::vector<art::Ptr<recob::Hit>>;
using ViewHitMap      = std::map<size_t,HitPtrVec>;
using TrackViewHitMap = std::map<int,ViewHitMap>;

class BasicHitAnalysis : virtual public IHitHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit BasicHitAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~BasicHitAnalysis();
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset) override;

    /**
     *  @brief Interface for initializing the histograms to be filled
     *
     *  @param TFileService   handle to the TFile service
     *  @param string         subdirectory to store the hists in
     */
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&) override;
    
    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    void endJob(int numEvents) override;
    
    /**
     *  @brief Interface for filling histograms
     */
    void fillHistograms(const HitPtrVec&)  const override;
    
private:
    
    // Fcl parameters.
    std::string fLocalDirName;     ///< Fraction for truncated mean
    
    // Pointers to the histograms we'll create.
    TH1D*     fHitsByWire[3];
    TH1D*     fDriftTimes[3];
    TH1D*     fHitsByTime[3];
    TH1D*     fPulseHeight[3];
    TH1D*     fPulseHeightSingle[3];
    TH1D*     fPulseHeightMulti[3];
    TH1D*     fChi2DOF[3];
    TH1D*     fNumDegFree[3];
    TH1D*     fChi2DOFSingle[3];
    TH1D*     fHitMult[3];
    TH1D*     fHitCharge[3];
    TH1D*     fFitWidth[3];
    TH1D*     fHitSumADC[3];
    TH2D*     fNDFVsChi2[3];
    TH2D*     fPulseHVsWidth[3];
    TH2D*     fPulseHVsCharge[3];
    TH1D*     fBadWPulseHeight;
    TH2D*     fBadWPulseHVsWidth;
    TH1D*     fBadWHitsByWire;
    TProfile* fPulseHVsHitNo[3];
    TProfile* fChargeVsHitNo[3];
    TProfile* fChargeVsHitNoS[3];
    
    TH2D*     fSPHvsIdx[3];
    TH2D*     fSWidVsIdx[3];
    TH2D*     f1PPHvsWid[3];
    TH2D*     fSPPHvsWid[3];
    TH2D*     fSOPHvsWid[3];
    TH2D*     fPHRatVsIdx[3];
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
BasicHitAnalysis::BasicHitAnalysis(fhicl::ParameterSet const & pset)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    configure(pset);
    
    // Report.
    mf::LogInfo("BasicHitAnalysis") << "BasicHitAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
BasicHitAnalysis::~BasicHitAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void BasicHitAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fLocalDirName = pset.get<std::string>("LocalDirName", std::string("wow"));
}

//----------------------------------------------------------------------------
/// Begin job method.
void BasicHitAnalysis::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());

    fHitsByWire[0]            = dir.make<TH1D>("HitsByWire0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0));
    fHitsByWire[1]            = dir.make<TH1D>("HitsByWire1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1));
    fHitsByWire[2]            = dir.make<TH1D>("HitsByWire2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    
    fDriftTimes[0]            = dir.make<TH1D>("DriftTime0",  ";time(ticks)", 3200, 0., 3200.);
    fDriftTimes[1]            = dir.make<TH1D>("DriftTime1",  ";time(ticks)", 3200, 0., 3200.);
    fDriftTimes[2]            = dir.make<TH1D>("DriftTime2",  ";time(ticks)", 3200, 0., 3200.);
    
    fHitsByTime[0]            = dir.make<TH1D>("HitsByTime0", ";Tick",   1600, 0., 3200.);
    fHitsByTime[1]            = dir.make<TH1D>("HitsByTime1", ";Tick",   1600, 0., 3200.);
    fHitsByTime[2]            = dir.make<TH1D>("HitsByTime2", ";Tick",   1600, 0., 3200.);
    
    fPulseHeight[0]           = dir.make<TH1D>("PulseHeight0",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeight[1]           = dir.make<TH1D>("PulseHeight1",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeight[2]           = dir.make<TH1D>("PulseHeight2",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[0]     = dir.make<TH1D>("PulseHeightS0", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[1]     = dir.make<TH1D>("PulseHeightS1", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[2]     = dir.make<TH1D>("PulseHeightS2", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[0]      = dir.make<TH1D>("PulseHeightM0", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[1]      = dir.make<TH1D>("PulseHeightM1", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[2]      = dir.make<TH1D>("PulseHeightM2", "PH (ADC)",  300,  0.,  150.);
    fChi2DOF[0]               = dir.make<TH1D>("Chi2DOF0",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOF[1]               = dir.make<TH1D>("Chi2DOF1",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOF[2]               = dir.make<TH1D>("Chi2DOF2",      "Chi2DOF",   502, -1.,  250.);
    fNumDegFree[0]            = dir.make<TH1D>("NumDegFree0",   "NDF",       200,  0.,  200.);
    fNumDegFree[1]            = dir.make<TH1D>("NumDegFree1",   "NDF",       200,  0.,  200.);
    fNumDegFree[2]            = dir.make<TH1D>("NumDegFree2",   "NDF",       200,  0.,  200.);
    fChi2DOFSingle[0]         = dir.make<TH1D>("Chi2DOFS0",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[1]         = dir.make<TH1D>("Chi2DOFS1",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[2]         = dir.make<TH1D>("Chi2DOFS2",     "Chi2DOF",   502, -1.,  250.);
    fHitMult[0]               = dir.make<TH1D>("HitMult0",      "# hits",     15,  0.,   15.);
    fHitMult[1]               = dir.make<TH1D>("HitMult1",      "# hits",     15,  0.,   15.);
    fHitMult[2]               = dir.make<TH1D>("HitMult2",      "# hits",     15,  0.,   15.);
    fHitCharge[0]             = dir.make<TH1D>("HitCharge0",    "Charge",   1000,  0., 2000.);
    fHitCharge[1]             = dir.make<TH1D>("HitCharge1",    "Charge",   1000,  0., 2000.);
    fHitCharge[2]             = dir.make<TH1D>("HitCharge2",    "Charge",   1000,  0., 2000.);
    fFitWidth[0]              = dir.make<TH1D>("FitWidth0",     "Width",     100,  0.,   20.);
    fFitWidth[1]              = dir.make<TH1D>("FitWidth1",     "Width",     100,  0.,   20.);
    fFitWidth[2]              = dir.make<TH1D>("FitWidth2",     "Width",     100,  0.,   20.);
    fHitSumADC[0]             = dir.make<TH1D>("SumADC0",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADC[1]             = dir.make<TH1D>("SumADC1",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADC[2]             = dir.make<TH1D>("SumADC2",       "Sum ADC",  1000,  0., 2000.);
    
    fNDFVsChi2[0]             = dir.make<TH2D>("NDFVsChi20",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    fNDFVsChi2[1]             = dir.make<TH2D>("NDFVsChi21",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    fNDFVsChi2[2]             = dir.make<TH2D>("NDFVsChi22",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    
    fPulseHVsWidth[0]         = dir.make<TH2D>("PHVsWidth0",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fPulseHVsWidth[1]         = dir.make<TH2D>("PHVsWidth1",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fPulseHVsWidth[2]         = dir.make<TH2D>("PHVsWidth2",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    
    fPulseHVsCharge[0]        = dir.make<TH2D>("PHVsChrg0",     ";PH;Q",     100,  0.,  100., 100,  0., 2000.);
    fPulseHVsCharge[1]        = dir.make<TH2D>("PHVsChrg1",     ";PH;Q",     100,  0.,  100., 100,  0., 2000.);
    fPulseHVsCharge[2]        = dir.make<TH2D>("PHVsChrg2",     ";PH;Q",     100,  0.,  100., 100,  0., 2000.);
    
    fPulseHVsHitNo[0]         = dir.make<TProfile>("PHVsNo0",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);
    fPulseHVsHitNo[1]         = dir.make<TProfile>("PHVsNo1",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);
    fPulseHVsHitNo[2]         = dir.make<TProfile>("PHVsNo2",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);
    
    fChargeVsHitNo[0]         = dir.make<TProfile>("QVsNo0",    ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    fChargeVsHitNo[1]         = dir.make<TProfile>("QVsNo1",    ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    fChargeVsHitNo[2]         = dir.make<TProfile>("QVsNo2",    ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    
    fChargeVsHitNoS[0]        = dir.make<TProfile>("QVsNoS0",   ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    fChargeVsHitNoS[1]        = dir.make<TProfile>("QVsNoS1",   ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    fChargeVsHitNoS[2]        = dir.make<TProfile>("QVsNoS2",   ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    
    fBadWPulseHeight          = dir.make<TH1D>("BWPulseHeight", "PH (ADC)",  300,  0.,  150.);
    fBadWPulseHVsWidth        = dir.make<TH2D>("BWPHVsWidth",   ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    fBadWHitsByWire           = dir.make<TH1D>("BWHitsByWire",  ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    
    fSPHvsIdx[0]              = dir.make<TH2D>("SPHVsIdx0",     ";PH;Idx", 30,  0.,  30., 100,  0., 100.);
    fSPHvsIdx[1]              = dir.make<TH2D>("SPHVsIdx1",     ";PH;Idx", 30,  0.,  30., 100,  0., 100.);
    fSPHvsIdx[2]              = dir.make<TH2D>("SPHVsIdx2",     ";PH;Idx", 30,  0.,  30., 100,  0., 100.);
    
    fSWidVsIdx[0]             = dir.make<TH2D>("SWidsIdx0",     ";Width;Idx", 30,  0.,  30., 100,  0., 10.);
    fSWidVsIdx[1]             = dir.make<TH2D>("SWidsIdx1",     ";Width;Idx", 30,  0.,  30., 100,  0., 10.);
    fSWidVsIdx[2]             = dir.make<TH2D>("SWidsIdx2",     ";Width;Idx", 30,  0.,  30., 100,  0., 10.);
    
    f1PPHvsWid[0]             = dir.make<TH2D>("1PPHVsWid0",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    f1PPHvsWid[1]             = dir.make<TH2D>("1PPHVsWid1",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    f1PPHvsWid[2]             = dir.make<TH2D>("1PPHVsWid2",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    
    fSPPHvsWid[0]             = dir.make<TH2D>("SPPHVsWid0",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    fSPPHvsWid[1]             = dir.make<TH2D>("SPPHVsWid1",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    fSPPHvsWid[2]             = dir.make<TH2D>("SPPHVsWid2",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    
    fSOPHvsWid[0]             = dir.make<TH2D>("SOPHVsWid0",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    fSOPHvsWid[1]             = dir.make<TH2D>("SOPHVsWid1",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    fSOPHvsWid[2]             = dir.make<TH2D>("SOPHVsWid2",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    
    fPHRatVsIdx[0]            = dir.make<TH2D>("PHRatVsIdx0",   ";PHRat;Idx", 30,  0.,  30., 51,  0., 1.02);
    fPHRatVsIdx[1]            = dir.make<TH2D>("PHRatVsIdx1",   ";PHRat;Idx", 30,  0.,  30., 51,  0., 1.02);
    fPHRatVsIdx[2]            = dir.make<TH2D>("PHRatVsIdx2",   ";PHRat;Idx", 30,  0.,  30., 51,  0., 1.02);

    return;
}
    
void BasicHitAnalysis::fillHistograms(const HitPtrVec& hitPtrVec) const
{
    // Keep track of number of hits per view
    size_t nHitsPerView[] = {0,0,0};
    size_t negCount(0);
    
    std::vector<const recob::Hit*> hitSnippetVec;
    
    // Loop the hits and make some plots
    for(const auto& hitPtr : hitPtrVec)
    {
        // Extract interesting hit parameters
        const geo::WireID& wireID   = hitPtr->WireID();
        float              chi2DOF  = std::min(hitPtr->GoodnessOfFit(),float(249.8));
        int                numDOF   = hitPtr->DegreesOfFreedom();
        int                hitMult  = hitPtr->Multiplicity();
        float              peakTime = std::min(float(3199.),hitPtr->PeakTime());
        float              charge   = hitPtr->Integral();
        float              sumADC   = hitPtr->SummedADC();
        float              hitPH    = std::min(hitPtr->PeakAmplitude(),float(149.8));
        float              hitSigma = hitPtr->RMS();
        
        size_t             plane    = wireID.Plane;
        size_t             wire     = wireID.Wire;
        
        if (hitPH < 0.)
        {
            negCount++;
            std::cout << "Hit Plane: " << plane << ", wire: " << wire << ", T: " << peakTime << ", PH: " << hitPH << ", charge: " << charge << ", sumADC: " << sumADC << std::endl;
        }
        
        nHitsPerView[plane]++;
        
        fHitsByWire[plane]->Fill(wire,1.);
        fHitsByTime[plane]->Fill(peakTime, 1.);
        fPulseHeight[plane]->Fill(hitPH, 1.);
        fChi2DOF[plane]->Fill(chi2DOF, 1.);
        fNumDegFree[plane]->Fill(numDOF, 1.);
        fHitMult[plane]->Fill(hitMult, 1.);
        fHitCharge[plane]->Fill(charge, 1.);
        fFitWidth[plane]->Fill(std::min(float(19.99),hitSigma), 1.);
        fHitSumADC[plane]->Fill(sumADC, 1.);
        fNDFVsChi2[plane]->Fill(numDOF, chi2DOF, 1.);
        fDriftTimes[plane]->Fill(peakTime, 1.);
        
        if (hitMult == 1)
        {
            fPulseHeightSingle[plane]->Fill(hitPH, 1.);
            fChi2DOFSingle[plane]->Fill(chi2DOF, 1.);
            fPulseHVsWidth[plane]->Fill(std::min(float(99.9),hitPH), std::min(float(19.99),hitSigma), 1.);
            fPulseHVsCharge[plane]->Fill(std::min(float(99.9),hitPH), std::min(float(1999.),charge), 1.);
            
            if (plane == 2 && hitPH < 5 && hitSigma < 2.2)
            {
                std::cout << "++> Plane: " << plane << ", wire: " << wire << ", peakTime: " << peakTime << ", ph: " << hitPH << ", w: " << hitSigma << std::endl;
                
                fBadWPulseHeight->Fill(hitPH,1.);
                fBadWPulseHVsWidth->Fill(std::min(float(99.9),hitPH), std::min(float(9.99),hitSigma), 1.);
                fBadWHitsByWire->Fill(wire,1.);
            }
        }
        else
            fPulseHeightMulti[plane]->Fill(hitPH, 1.);
        
        // Look at hits on snippets
        if (!hitSnippetVec.empty() && hitSnippetVec.back()->LocalIndex() >= hitPtr->LocalIndex())
        {
            // Only worried about multi hit snippets
            if (hitSnippetVec.size() > 1)
            {
                // Sort in order of largest to smallest pulse height
                std::sort(hitSnippetVec.begin(),hitSnippetVec.end(),[](const auto* left, const auto* right){return left->PeakAmplitude() > right->PeakAmplitude();});
                
                float maxPulseHeight = hitSnippetVec.front()->PeakAmplitude();
                
                for(size_t idx = 0; idx < hitSnippetVec.size(); idx++)
                {
                    float pulseHeight      = hitSnippetVec.at(idx)->PeakAmplitude();
                    float pulseWid         = hitSnippetVec.at(idx)->RMS();
                    float pulseHeightRatio = pulseHeight / maxPulseHeight;
                    
                    size_t plane = hitSnippetVec.at(idx)->WireID().Plane;
                    
                    fSPHvsIdx[plane]->Fill(idx, std::min(float(99.9),pulseHeight), 1.);
                    fSWidVsIdx[plane]->Fill(idx, std::min(float(9.99),pulseWid), 1.);
                    fPHRatVsIdx[plane]->Fill(idx, pulseHeightRatio, 1.);
                    
                    if (idx == 0) fSPPHvsWid[plane]->Fill(std::min(float(99.9),pulseHeight), std::min(float(9.99),pulseWid), 1.);
                    else          fSOPHvsWid[plane]->Fill(std::min(float(99.9),pulseHeight), std::min(float(9.99),pulseWid), 1.);
                }
            }
            else
            {
                float  pulseHeight = hitSnippetVec.front()->PeakAmplitude();
                float  pulseWid    = hitSnippetVec.front()->RMS();
                size_t plane       = hitSnippetVec.front()->WireID().Plane;
                
                f1PPHvsWid[plane]->Fill(std::min(float(99.9),pulseHeight), std::min(float(9.99),pulseWid), 1.);
            }
            
            hitSnippetVec.clear();
        }
        
        hitSnippetVec.push_back(hitPtr.get());
    }
    
    return;
}
    
// Useful for normalizing histograms
void BasicHitAnalysis::endJob(int numEvents)
{
    // Normalize wire profiles to be hits/event
    double normFactor(1./numEvents);
    
    for(size_t idx = 0; idx < 3; idx++) fHitsByWire[idx]->Scale(normFactor);
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(BasicHitAnalysis)
}
