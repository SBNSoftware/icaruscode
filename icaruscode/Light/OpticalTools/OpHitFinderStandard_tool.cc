////////////////////////////////////////////////////////////////////////
/// \file   OpHitFinderStandard.cc
/// \author A. Falcone
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/Light/OpticalTools/IOpHitFinder.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"

#include <fstream>
#include <algorithm> // std::minmax_element()

namespace light
{

class OpHitFinderStandard : IOpHitFinder
{
public:
    explicit OpHitFinderStandard(const fhicl::ParameterSet& pset);
    
    ~OpHitFinderStandard();
    
    void configure(const fhicl::ParameterSet& pset)                  override;
    void outputHistograms(art::TFileDirectory&)                const override;
    
    void FindOpHits(const raw::OpDetWaveform&, OpHitVec&)      const override;
    
private:
    double fSPEArea;         //conversion between phe and Adc*ns
    double fHitThreshold;
    int    fBaselineSample;
};
    
//----------------------------------------------------------------------
// Constructor.
OpHitFinderStandard::OpHitFinderStandard(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
OpHitFinderStandard::~OpHitFinderStandard()
{
}
    
void OpHitFinderStandard::configure(const fhicl::ParameterSet& pset)
{
    fHitThreshold    = pset.get< double >("HitThreshold");
    fSPEArea         = pset.get< double >("SPEArea");
    fBaselineSample  = pset.get< int    >("BaselineSample");
    
    return;
}

    
void OpHitFinderStandard::FindOpHits(const raw::OpDetWaveform& opDetWaveform,
                                     OpHitVec&                 opHitVec) const
{
    int chNumber = opDetWaveform.ChannelNumber();
    std::cout << "Photon channel: " << chNumber << std::endl;
    // Load pulses into WaveformVector
    //std::vector < short_t >  WaveformVector = wvf.Waveform();
    
    int grsize = opDetWaveform.size();
    
    //std::vector<Double_t> TimeVector;
    //std::vector<Double_t> ADCVector;
    
    //TimeVector.size()= grsize;
    //ADCVector.size()= grsize;

    double TimeVector[10000];
    double ADCVector[10000];

    unsigned short frame;
    double         Area;
    double         fasttotal;
    double         time_abs;
    double         FWHM;
    double         phelec;

    double baseline=0;
    
    for (int btime =0; btime< fBaselineSample; btime++)
    {
        baseline = baseline+ opDetWaveform[btime];
    }
    
    baseline = baseline/fBaselineSample;
    
    std::cout << "Baseline " << baseline << std::endl;
    
    
    for (int wtime=0; wtime< grsize; wtime++)
    {
        TimeVector[wtime] = wtime;
        ADCVector[wtime]  = -(opDetWaveform[wtime]-baseline);
    }
    
    TGraph *gr = new TGraph(grsize,TimeVector,ADCVector);
    
    int n_graph = gr->GetN();
    double *y_graph = gr->GetY();
    
    int    min_time        = TMath::LocMax(n_graph,y_graph);
    double min_time_to_put = TMath::LocMax(n_graph,y_graph);
    double min             = y_graph[min_time];
    double min_to_put      = min;
    
    std::cout << "Min " << min << std::endl;
    
    if (min>fHitThreshold)
    {
        TF1 *funz= new TF1("funz", "pol1(0)", min_time-2.0, min_time);
        
        gr->Fit("funz","R");
        
        double start_moment = funz->GetX(baseline, 0, min_time);
        
        std::cout << "Start " << start_moment << std::endl;
        
        TF1 *gauss_start= new TF1("gauss_start", "gaus", min_time-5.0, min_time);
        TF1 *gauss_end  = new TF1("gauss_end", "gaus", min_time, min_time+10);
        
        gauss_start->SetParameter(1,min_time);
        gauss_end->SetParameter(1,min_time);
        
        gr->Fit("gauss_start","R");
        
        double Constant1 = gauss_start->GetParameter(0);
        //double Mean1 = gauss_start->GetParameter(1);
        double Sigma1 = gauss_start->GetParameter(2);
        
        std::cout << "GaussParam 00 " << gauss_start->GetParameter(0) << "GaussParam 01 " << gauss_start->GetParameter(1) << "GaussParam 02 " << gauss_start->GetParameter(2) << std::endl;
        
        gr->Fit("gauss_end","R");
        
        double Constant2 = gauss_end->GetParameter(0);
        //double Mean2 = gauss_end->GetParameter(1);
        double Sigma2 = gauss_end->GetParameter(2);
        
        std::cout << "GaussParam 10 " << gauss_end->GetParameter(0) << "GaussParam 11 " << gauss_end->GetParameter(1) << "GaussParam 12 " << gauss_end->GetParameter(2) << std::endl;
        
        frame     = 1;
        Area      =  ((Constant1*Sigma1)/2 + (Constant2*Sigma2)/2)*sqrt(2*3.14159);
        fasttotal = 3/4;
        time_abs  = sqrt(min_time_to_put);
        FWHM      = 2.35*((Sigma1+Sigma2)/2);
        phelec    = Area/fSPEArea;
        
        //    recob::OpHit adcVec(chNumber, min_time_to_put, time_abs, frame, FWHM, Area, min_to_put, phelec, fasttotal);//including hit info
        //    pulseVecPtr->emplace_back(std::move(adcVec));
        
        funz->~TF1();
        gauss_start->~TF1();
        gauss_end->~TF1();
        gr->~TGraph();
    }
    else
    {
        //std::cout << "No OpHit in channel " << chNumber << std::endl;
        min_time_to_put=0;
        min_to_put=0;
        Area=0;
        frame=0;
        fasttotal=0;
        time_abs=0;
        FWHM=0;
        phelec=0;
    }
    
    recob::OpHit opHit(chNumber, min_time_to_put, time_abs, frame, FWHM, Area, min_to_put, phelec, fasttotal);
    
    opHitVec.push_back(recob::OpHit());
    opHitVec.back() = opHit;

    return;
}
    
void OpHitFinderStandard::outputHistograms(art::TFileDirectory& histDir) const
{
    // It is assumed that the input TFileDirectory has been set up to group histograms into a common
    // folder at the calling routine's level. Here we create one more level of indirection to keep
    // histograms made by this tool separate.
/*
    std::string dirName = "OpHitFinderPlane_" + std::to_string(fPlane);
    
    art::TFileDirectory dir = histDir.mkdir(dirName.c_str());
    
    auto const* detprop      = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double      samplingRate = detprop->SamplingRate();
    double      numBins      = fOpHitFinderVec.size();
    double      maxFreq      = 500. / samplingRate;
    std::string histName     = "OpHitFinderPlane_" + std::to_string(fPlane);
    
    TH1D*       hist         = dir.make<TH1D>(histName.c_str(), "OpHitFinder;Frequency(MHz)", numBins, 0., maxFreq);
    
    for(int bin = 0; bin < numBins; bin++)
    {
        double freq = maxFreq * double(bin + 0.5) / double(numBins);
        
        hist->Fill(freq, fOpHitFinderVec.at(bin).Re());
    }
*/
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(OpHitFinderStandard)
}
