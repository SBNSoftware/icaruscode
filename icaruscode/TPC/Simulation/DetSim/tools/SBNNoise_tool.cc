////////////////////////////////////////////////////////////////////////
/// \file   SBNNoise.cc
/// \author F. Varanini
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IGenNoise.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art_root_io/TFileService.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "icaruscode/TPC/Utilities/tools/IWaveformTool.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TComplex.h"

#include <complex.h>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include <fstream>

namespace icarus_tool
{

class SBNNoise : IGenNoise
{
public:
    explicit SBNNoise(const fhicl::ParameterSet& pset);
    
    ~SBNNoise();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    void nextEvent() override;

    void generateNoise(CLHEP::HepRandomEngine& noise_engine,
                       CLHEP::HepRandomEngine& cornoise_engine,
                       icarusutil::TimeVec& noise,
                       double noise_factor,
                       unsigned int wire) override;
    
private:
    void GenerateCorrelatedNoise(CLHEP::HepRandomEngine&, icarusutil::TimeVec&, double, unsigned int);
    void GenerateUncorrelatedNoise(CLHEP::HepRandomEngine&, icarusutil::TimeVec&, double, unsigned int);
    void GenNoise(std::function<void (double[])>&, const icarusutil::TimeVec&, icarusutil::TimeVec&, float);
    void ComputeRMSs();
    void makeHistograms();
    
    // Member variables from the fhicl file
    size_t                                      fPlane;
    int                                         fMedianNumBins;
    float                                       fNoiseRand;
    long                                        fCorrelatedSeed;
    long                                        fUncorrelatedSeed;
    float                                       fIncoherentNoiseFrac;
    bool                                        fStoreHistograms;
    std::string                                 fInputNoiseHistFileName;
    std::string                                 fHistogramName;
    std::string                                 fCorrelatedHistogramName;
    std::string                                 fUncorrelatedHistogramName;
    std::string                                 fCorrelatedRMSHistoName;
    std::string                                 fUncorrelatedRMSHistoName;
    std::string                                 fTotalRMSHistoName;
    std::unique_ptr<icarus_tool::IWaveformTool> fWaveformTool;

    // We'll recover the bin contents and store in a vector
    // with the likely false hope this will be faster...
    icarusutil::TimeVec                          fCoherentNoiseVec;       //< Input full noise frequency distribution
    icarusutil::TimeVec                          fIncoherentNoiseVec;       //< Input full noise frequency distribution

    std::vector<double>                         fCorrAmpDistVec;     //< Keep track of motherboard contributions
    
    double                                      fIncoherentNoiseRMS; //< RMS of full noise waveform
    double                                      fCoherentNoiseRMS;   //< RMS of full noise waveform

    // Container for doing the work
    icarusutil::FrequencyVec                    fNoiseFrequencyVec;
    
    // Keep track of seed initialization for uncorrelated noise
    bool                                        fNeedFirstSeed=true;
    
    // Histograms
    TProfile*                                   fInputNoiseHist;
    TProfile*                                   fMedianNoiseHist;
    TProfile*                                   fPeakNoiseHist;
    TProfile*                                   fCorAmpDistHist;
    TH1D*                                       corrRMSHistPtr;
    TH1D*                                       uncorrRMSHistPtr;
    TH1D*                                       totalRMSHistPtr;

float totalRMS;

    
    // Keep instance of the eigen FFT
    Eigen::FFT<double>                          fEigenFFT;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();   ///< Detector properties service
};
    
//----------------------------------------------------------------------
// Constructor.
SBNNoise::SBNNoise(const fhicl::ParameterSet& pset)
{
    // Recover the configuration of the tool from the input fhicl file and set up
    configure(pset);
ComputeRMSs();
    
    // Output some histograms to catalogue what's been done
    makeHistograms();
}
    
SBNNoise::~SBNNoise()
{
}
    
void SBNNoise::configure(const fhicl::ParameterSet& pset)
{
    // Recover the histogram used for noise generation
    fPlane                     = pset.get< size_t      >("Plane");
    fMedianNumBins             = pset.get< int         >("MedianNumBins");
    fNoiseRand                 = pset.get< float       >("NoiseRand");
    fCorrelatedSeed            = pset.get< long        >("CorrelatedSeed",1000);
    fUncorrelatedSeed          = pset.get< long        >("UncorrelatedSeed",5000);
    fIncoherentNoiseFrac       = pset.get< float       >("IncoherentNoiseFraction",0.5);
    fStoreHistograms           = pset.get< bool        >("StoreHistograms");
    fInputNoiseHistFileName    = pset.get< std::string >("NoiseHistFileName");
    fCorrelatedHistogramName   = pset.get< std::string >("CorrelatedHistogramName");
    fUncorrelatedHistogramName = pset.get< std::string >("UncorrelatedHistogramName");
    fCorrelatedRMSHistoName    = pset.get< std::string >("CorrelatedRMSHistoName");
    fUncorrelatedRMSHistoName  = pset.get< std::string >("UncorrelatedRMSHistoName");
    fTotalRMSHistoName         = pset.get< std::string >("TotalRMSHistoName");
    // Initialize the work vector
    fNoiseFrequencyVec.resize(fDetectorProperties->NumberTimeSamples(),std::complex<float>(0.,0.));

    // Set up to input the histogram with the overall noise spectrum
    std::string fullFileName;
    std::string corrAmpFileName;

    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file(fInputNoiseHistFileName, fullFileName);
    
    TFile inputFile(fullFileName.c_str(), "READ");
    
    if (!inputFile.IsOpen())
        throw cet::exception("NoiseFromHist::configure") << "Unable to open input file: " << fInputNoiseHistFileName << std::endl;
    
    TH1D* corrHistPtr = (TH1D*)inputFile.Get(fCorrelatedHistogramName.c_str());    
    if (!corrHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fCorrelatedHistogramName << std::endl;

  TH1D* uncorrHistPtr = (TH1D*)inputFile.Get(fUncorrelatedHistogramName.c_str());    
    if (!uncorrHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fUncorrelatedHistogramName << std::endl;

    corrRMSHistPtr = (TH1D*)inputFile.Get(fCorrelatedRMSHistoName.c_str());  
    if (!corrRMSHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fCorrelatedRMSHistoName << std::endl;

  uncorrRMSHistPtr = (TH1D*)inputFile.Get(fUncorrelatedRMSHistoName.c_str());  
    if (!uncorrRMSHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fUncorrelatedRMSHistoName << std::endl;

 totalRMSHistPtr = (TH1D*)inputFile.Get(fTotalRMSHistoName.c_str());  
    if (!totalRMSHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fTotalRMSHistoName << std::endl;

    // Close the input file
    inputFile.Close();
std::cout << " corr nbins " << corrHistPtr->GetNbinsX() << std::endl;
std::cout << " uncorr nbins " << uncorrHistPtr->GetNbinsX() << std::endl;
    fCoherentNoiseVec.resize(corrHistPtr->GetNbinsX(), 0.);
    
    // Recover the bin contents into local vectors
    for(size_t histIdx = 0; histIdx < size_t(corrHistPtr->GetNbinsX()); histIdx++)
        fCoherentNoiseVec[histIdx] = corrHistPtr->GetBinContent(histIdx+1);
fIncoherentNoiseVec.resize(uncorrHistPtr->GetNbinsX(), 0.);
    
    // Recover the bin contents into local vectors
    for(size_t histIdx = 0; histIdx < size_t(uncorrHistPtr->GetNbinsX()); histIdx++)
        fIncoherentNoiseVec[histIdx] = uncorrHistPtr->GetBinContent(histIdx+1);
    // Should we store hists?

std::cout << " after filling vectors " << std::endl;
    if (fStoreHistograms)
    {
        // Define histograms
        art::ServiceHandle<art::TFileService> tfs;
    
        art::TFileDirectory* histDirectory = tfs.get();
        
        // Make a directory for these histograms
        art::TFileDirectory dir = histDirectory->mkdir(Form("CorNoisePlane%1zu",fPlane));
        
        float sampleRate  = fDetectorProperties->SamplingRate();
        float readOutSize = fDetectorProperties->ReadOutWindowSize();
        float maxFreq     = 1.e6 / (2. * sampleRate);
        float minFreq     = 1.e6 / (2. * sampleRate * readOutSize);
        int   numSamples  = readOutSize / 2;
        
        fInputNoiseHist   = dir.make<TProfile>("InNoise",   ";freq(kHz)", numSamples, minFreq, maxFreq);
        fMedianNoiseHist  = dir.make<TProfile>("MedNoise",  ";freq(kHz)", numSamples, minFreq, maxFreq);;
        fPeakNoiseHist    = dir.make<TProfile>("PeakNoise", ";freq(kHz)", numSamples, minFreq, maxFreq);;
        
        fCorAmpDistHist   = dir.make<TProfile>("CorAmp",    ";Motherboard", fCorrAmpDistVec.size(),0.,fCorrAmpDistVec.size());
    }

    // Recover an instance of the waveform tool
    // Here we just make a parameterset to pass to it...
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    fWaveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);
    
    return;
}
    
void SBNNoise::nextEvent()
{
    // We update the correlated seed because we want to see different noise event-by-event
    fCorrelatedSeed   = (333 * fCorrelatedSeed) % 900000000;
    
    return;
}

void SBNNoise::generateNoise(CLHEP::HepRandomEngine& engine_unc,
                                    CLHEP::HepRandomEngine& engine_corr,
                                    icarusutil::TimeVec&     noise,
                                    double                  noise_factor,
                                    unsigned int            channel)
{
noise_factor=totalRMS;
//std::cout <<  " generating noise " << std::endl;
    // Define a couple of vectors to hold intermediate work
    icarusutil::TimeVec noise_unc(noise.size(),0.);
    icarusutil::TimeVec noise_corr(noise.size(),0.);
    
    // Make sure the work vector is size right with the output
    if (fNoiseFrequencyVec.size() != noise.size()) fNoiseFrequencyVec.resize(noise.size(),std::complex<float>(0.,0.));
    //std::cout <<  " generating uncorrelated noise " << std::endl;
    // If applying incoherent noise call the generator
    if (fIncoherentNoiseFrac > 0.) GenerateUncorrelatedNoise(engine_unc,noise_unc,noise_factor,channel);
    
    int board=channel/64;
//std::cout <<  " generating correlated noise " << std::endl;
    // If applying coherent noise call the generator
    if (fIncoherentNoiseFrac < 1.) GenerateCorrelatedNoise(engine_corr, noise_corr, noise_factor, board);
   // std::cout <<  " summing noise " << std::endl;
    // Take the noise as the simple sum of the two contributions
    std::transform(noise_unc.begin(),noise_unc.end(),noise_corr.begin(),noise.begin(),std::plus<float>());
    
    return;
}
    
void SBNNoise::GenerateUncorrelatedNoise(CLHEP::HepRandomEngine& engine, icarusutil::TimeVec &noise, double noise_factor, unsigned int channel)
{
    // Here we aim to produce a waveform consisting of incoherent noise
    // Note that this is expected to be the dominate noise contribution
    // Check for seed initialization
    if (fNeedFirstSeed)
    {
        engine.setSeed(fUncorrelatedSeed,0);
        fNeedFirstSeed = false;
    }
    
    // Get the generator
    CLHEP::RandFlat noiseGen(engine,0,1);
    
    std::function<void (double[])> randGenFunc = [&noiseGen](double randArray[]){noiseGen.fireArray(2,randArray);};

    float  scaleFactor = fIncoherentNoiseFrac * noise_factor / fIncoherentNoiseRMS;
    
    GenNoise(randGenFunc, fIncoherentNoiseVec, noise, scaleFactor);

    return;
}
    
void SBNNoise::GenerateCorrelatedNoise(CLHEP::HepRandomEngine& engine, icarusutil::TimeVec &noise, double noise_factor, unsigned int board)
{    
   
        // Set the engine seed to the board being considered
        engine.setSeed(fCorrelatedSeed+board,0);
        
        CLHEP::RandFlat noiseGen(engine,0,1);
        
        std::function<void (double[])> randGenFunc = [&noiseGen](double randArray[]){noiseGen.fireArray(2,randArray);};
        
        float fraction    = std::sqrt(1. - fIncoherentNoiseFrac * fIncoherentNoiseFrac);
        float scaleFactor = fraction * noise_factor / fCoherentNoiseRMS;
        
        GenNoise(randGenFunc, fCoherentNoiseVec, noise, scaleFactor);
    
    
    return;
}
    
void SBNNoise::GenNoise(std::function<void (double[])>& gen,const icarusutil::TimeVec& freqDist, icarusutil::TimeVec& noise, float scaleFactor)
{
    double rnd_corr[2] = {0.,0.};
    
    // Build out the frequency vector
    for(size_t i=0; i< noise.size()/2; ++i)
    {
        // exponential noise spectrum
        gen(rnd_corr);
        
        float pval  = freqDist[i] * ((1-fNoiseRand) + 2 * fNoiseRand*rnd_corr[0]) * scaleFactor;
        float phase = rnd_corr[1] * 2. * TMath::Pi();
        
        std::complex<float> tc(pval*cos(phase),pval*sin(phase));
        
        fNoiseFrequencyVec[i] = tc;
    }
    
    // inverse FFT MCSignal
    fEigenFFT.inv(noise, fNoiseFrequencyVec);
    
    return;
}

void SBNNoise::makeHistograms()
{
    
    return;
}
   
void SBNNoise::ComputeRMSs()
{
    // Let's get the rms we expect from the incoherent noise contribution to the waveform
    // A couple of ways to do this, let's basically invert the frequency spectrum to
    // produce a waveform and then get the rms from that
    std::function<void (double[])> randGenFunc = [](double randArray[]){randArray[0]=0.5; randArray[1]=0.5;};
    
    icarusutil::TimeVec waveNoise(fIncoherentNoiseVec.size());
    float              scaleFactor = 1.;
    
    GenNoise(randGenFunc, fIncoherentNoiseVec, waveNoise, scaleFactor);
    
    // Now get the details...
    double nSig(3.);
    double mean,rmsTrunc;
    int    nTrunc;
    
    // Use the waveform tool to recover the full rms
    fWaveformTool->getTruncatedMeanRMS(waveNoise, nSig, mean, fIncoherentNoiseRMS, rmsTrunc, nTrunc);
    
    // Do the same for the coherent term
    GenNoise(randGenFunc, fCoherentNoiseVec, waveNoise, scaleFactor);
    
    fWaveformTool->getTruncatedMeanRMS(waveNoise, nSig, mean, fCoherentNoiseRMS, rmsTrunc, nTrunc);

    if (fStoreHistograms)
    {
        for(size_t idx = 0; idx < fCorrAmpDistVec.size(); idx++)
            fCorAmpDistHist->Fill(idx, fCorrAmpDistVec[idx], 1.);
    }

 totalRMS=totalRMSHistPtr->GetMean();
double rmsUnc=uncorrRMSHistPtr->GetMean();
//double rmsCorr=corrRMSHistPtr->GetMean();
fIncoherentNoiseFrac=rmsUnc/totalRMS;
    return;
}
 
DEFINE_ART_CLASS_TOOL(SBNNoise)
}
