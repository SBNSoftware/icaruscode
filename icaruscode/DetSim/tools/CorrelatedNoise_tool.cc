////////////////////////////////////////////////////////////////////////
/// \file   UncorrelatedNoise.cc
/// \author F. Varanini
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IGenNoise.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Optional/TFileService.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

#include "icaruscode/Utilities/tools/IWaveformTool.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TComplex.h"

#include <fstream>

namespace icarus_tool
{

class CorrelatedNoise : IGenNoise
{
public:
    explicit CorrelatedNoise(const fhicl::ParameterSet& pset);
    
    ~CorrelatedNoise();
    
    void configure(const fhicl::ParameterSet& pset)                          override;

    void GenerateNoise(CLHEP::HepRandomEngine&, std::vector<float> &noise, double noise_factor, unsigned int wire) const override;
    void GenerateUncorrelatedNoise(std::vector<float> &noise, double noise_factor, unsigned int wire) const ;
    void GenerateCorrelatedNoise(std::vector<float> &noise, double noise_factor, unsigned int wire) const ;
    
    void ExtractCorrelatedAmplitude(double &, int) const;


    void SelectContinuousSpectrum() ;
    void FindPeaks() ;
    
private:
    void makeHistograms();
    
    // Member variables from the fhicl file
    size_t              fPlane;
    int                 fMedianNumBins;
    double              fNoiseRand;
    long                fCorrelatedSeed;
    long                fUncorrelatedSeed;
    bool                fStoreHistograms;
    std::string         fInputNoiseHistFileName;
    std::string         fHistogramName;
    std::string         fCorrAmpHistFileName;
    std::string         fCorrAmpHistogramName;
    
    double              fHistNormFactor;

    std::vector<int>    fPeakVec;

    std::unique_ptr<icarus_tool::IWaveformTool> fWaveformTool;

    // We'll recover the bin contents and store in a vector
    // with the likely false hope this will be faster...
    std::vector<double> fNoiseHistVec;
    std::vector<double> fNoiseContVec;
    
    std::vector<double> fCorrAmpDistVec;
    
    // Histograms
    TProfile*           fInputNoiseHist;
    TProfile*           fMedianNoiseHist;
    TProfile*           fPeakNoiseHist;

    // Local random generators
    std::unique_ptr<CLHEP::RandFlat> fCorrelatedGen;
    std::unique_ptr<CLHEP::RandFlat> fUncorrelatedGen;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();   ///< Detector properties service
};
    
//----------------------------------------------------------------------
// Constructor.
CorrelatedNoise::CorrelatedNoise(const fhicl::ParameterSet& pset)
{
    std::cout << " Correlated Noise " << std::endl;
    configure(pset);
//    FindPeaks();
    SelectContinuousSpectrum();
    makeHistograms();
    
    // Set seeds for the two random engines

    // FIXME: When art 3.02 is released, the 'getEngine' call will be
    // deprecated.  One of the reasons for this is that it can be
    // difficult to determine which module is active whenever
    // 'getEngine' is called.  To solve this, the random-number
    // engines should be passed in to the 'CorrelatedNoise'
    // constructor from the modules that creates the engines.
    // Interface will be added to nutools so that the 'getEngine'
    // function never needs to be called.
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine_unc = rng->getEngine(art::ScheduleID::first(),pset.get<std::string>("module_label"),"noise");
    engine_unc.setSeed(fUncorrelatedSeed,0);
    fUncorrelatedGen = std::make_unique<CLHEP::RandFlat>(engine_unc,-1,1);
    
    CLHEP::HepRandomEngine &engine_corr = rng->getEngine(art::ScheduleID::first(),pset.get<std::string>("module_label"),"cornoise");
    engine_corr.setSeed(fCorrelatedSeed,0);
    fCorrelatedGen = std::make_unique<CLHEP::RandFlat>(engine_corr,-1,1);
}
    
CorrelatedNoise::~CorrelatedNoise()
{
}
    
void CorrelatedNoise::configure(const fhicl::ParameterSet& pset)
{
    // Recover the histogram used for noise generation
    fPlane                  = pset.get<size_t>("Plane");
    fMedianNumBins          = pset.get<int>("MedianNumBins");
    fNoiseRand              = pset.get< double>("NoiseRand");
    fCorrelatedSeed         = pset.get< long >("CorrelatedSeed",1000);
    fUncorrelatedSeed       = pset.get< long >("UncorrelatedSeed",5000);
    fStoreHistograms        = pset.get<bool>("StoreHistograms");
    fInputNoiseHistFileName = pset.get<std::string>("NoiseHistFileName");
    fHistogramName          = pset.get<std::string>("HistogramName");
    fHistNormFactor         = pset.get<double>("HistNormFactor");
    fCorrAmpHistFileName    = pset.get<std::string>("CorrAmpHistFileName");
    fCorrAmpHistogramName   = pset.get<std::string>("CorrAmpHistogramName");
    
    std::string fullFileName;
    std::string corrAmpFileName;

    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file(fInputNoiseHistFileName, fullFileName);
    
    TFile inputFile(fullFileName.c_str(), "READ");
    
    if (!inputFile.IsOpen())
        throw cet::exception("NoiseFromHist::configure") << "Unable to open input file: " << fInputNoiseHistFileName << std::endl;
    
    TH1D* histPtr = (TH1D*)inputFile.Get(fHistogramName.c_str());
    
    if (!histPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fHistogramName << std::endl;
    
    // Close the input file
    inputFile.Close();
    
    searchPath.find_file(fCorrAmpHistFileName, corrAmpFileName);
    
    TFile corrAmpInputFile(corrAmpFileName.c_str(), "READ");
    
    if (!corrAmpInputFile.IsOpen())
        throw cet::exception("CorrelatedNoise::configure") << "Unable to open input file: " << fCorrAmpHistFileName << std::endl;
    
    TH1D* corrAmpHistPtr = (TH1D*)corrAmpInputFile.Get(fCorrAmpHistogramName.c_str());
    
    if (!corrAmpHistPtr)
        throw cet::exception("CorrelatedNoise::configure") << "Unable to recover desired histogram: " << fCorrAmpHistogramName << std::endl;
    
    fNoiseHistVec.resize(histPtr->GetNbinsX(), 0.);
    fCorrAmpDistVec.resize(corrAmpHistPtr->GetNbinsX(),0.);
    
    
    for(size_t histIdx = 0; histIdx < size_t(histPtr->GetNbinsX()); histIdx++) {
        fNoiseHistVec[histIdx] = histPtr->GetBinContent(histIdx+1);
    }
    
    for(size_t histIdx = 0; histIdx < 200; histIdx++) {
        fCorrAmpDistVec[histIdx] = corrAmpHistPtr->GetBinContent(histIdx+1);
    }
    
    // Close the input file
    inputFile.Close();
    
    // Should we store hists?
    if (fStoreHistograms)
    {
        // Define histograms
        art::ServiceHandle<art::TFileService> tfs;
    
        art::TFileDirectory* histDirectory = tfs.get();
        
        // Make a directory for these histograms
        art::TFileDirectory dir = histDirectory->mkdir(Form("CorNoisePlane%1zu",fPlane));
        
        double sampleRate  = fDetectorProperties->SamplingRate();
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        double maxFreq     = 1.e6 / (2. * sampleRate);
        double minFreq     = 1.e6 / (2. * sampleRate * readOutSize);
        int    numSamples  = readOutSize / 2;
        
        fInputNoiseHist   = dir.make<TProfile>("InNoise",   ";freq(kHz)", numSamples, minFreq, maxFreq);
        fMedianNoiseHist  = dir.make<TProfile>("MedNoise",  ";freq(kHz)", numSamples, minFreq, maxFreq);;
        fPeakNoiseHist    = dir.make<TProfile>("PeakNoise", ";freq(kHz)", numSamples, minFreq, maxFreq);;
    }

    // Recover an instance of the waveform tool
    // Here we just make a parameterset to pass to it...
    fhicl::ParameterSet waveformToolParams;
    
    waveformToolParams.put<std::string>("tool_type","Waveform");
    
    fWaveformTool = art::make_tool<icarus_tool::IWaveformTool>(waveformToolParams);
    
    return;
}

void CorrelatedNoise::GenerateNoise(CLHEP::HepRandomEngine&, std::vector<float> &noise, double noise_factor, unsigned int channel) const
{
    //std::cout << " generating noise " << std::endl;
    std::vector<float> noise_unc;
    std::vector<float> noise_corr;
    
    //ExtractCorrelatedAmplitude(corr_factor);
    
    noise_unc.resize(noise.size(),0.);
    
    GenerateUncorrelatedNoise(noise_unc,noise_factor,channel);
    
    int board=channel/32;
    
    noise_corr.resize(noise.size(), 0.);
    
    GenerateCorrelatedNoise(noise_corr,noise_factor,board);
    
    std::transform(noise_unc.begin(),noise_unc.end(),noise_corr.begin(),noise.begin(),std::plus<float>());
    
    return;
}
    
void CorrelatedNoise::GenerateUncorrelatedNoise(std::vector<float> &noise, double noise_factor, unsigned int channel) const
{
    //std::cout << " generating uncorrelated noise " << std::endl;

    art::ServiceHandle<util::LArFFT> fFFT;
 
    size_t nFFTTicks = fFFT->FFTSize();
        //std::cout << " ticks " <<nFFTTicks << std::endl;
       // std::cout << " noise size " <<nFFTTicks << std::endl;

    if(noise.size() != nFFTTicks)
        throw cet::exception("SimWireICARUS")
        << "\033[93m"
        << "Frequency noise vector length must match FFT Ticks (FFT size)"
        << " ... " << noise.size() << " != " << nFFTTicks
        << "\033[00m"
        << std::endl;
    
    // noise in frequency space
    std::vector<TComplex> noiseFrequency(nFFTTicks/2+1, 0.);
    
    double pval       = 0.;
    double phase      = 0.;
    double rnd_unc[2] = {0.,0.};

    double scaleFactor = fHistNormFactor * noise_factor;
    
    // uncorrelated, smooth frequency spectrum
     // std::cout << " before loop " <<nFFTTicks << std::endl;
    for(size_t i=0; i< nFFTTicks/2 + 1; ++i)
    {
        // exponential noise spectrum
        fUncorrelatedGen->fireArray(2,rnd_unc,0,1);
      //  if(i<10) std::cout << " channel " << channel << " uncorrelated bin " << i << " rnd0 " << rnd_unc[0] << std::endl;
        pval = fNoiseContVec[i] * ((1-fNoiseRand) + 2 * fNoiseRand*rnd_unc[0]) * scaleFactor;

        phase = rnd_unc[1] * 2. * TMath::Pi();

        TComplex tc(pval*cos(phase),pval*sin(phase));

        noiseFrequency.at(i) += tc;
    }
    
    // inverse FFT MCSignal
    fFFT->DoInvFFT(noiseFrequency, noise);
    
    noiseFrequency.clear();

    return;
}
    
void CorrelatedNoise::GenerateCorrelatedNoise(std::vector<float> &noise, double noise_factor, unsigned int board) const
{
    art::ServiceHandle<util::LArFFT> fFFT;
    
    double rnd_corr[2] = {0.,0.};
    size_t nFFTTicks   = fFFT->FFTSize();
    
    double cf=1;
    ExtractCorrelatedAmplitude(cf,board);
    
   // std::cout << " noise size " << noise.size() << " ticks " << nFFTTicks << std::endl;
    
    if(noise.size() != nFFTTicks)
        throw cet::exception("SimWireICARUS")
        << "\033[93m"
        << "Frequency noise vector length must match FFT Ticks (FFT size)"
        << " ... " << noise.size() << " != " << nFFTTicks
        << "\033[00m"
        << std::endl;
    
   // std::cout << " after ticks check " << std::endl;
    
    // noise in frequency space
    std::vector<TComplex> noiseFrequency(nFFTTicks/2+1, 0.);
    
    double pval        = 0.;
    double phase       = 0.;
    double scaleFactor = (1. - fHistNormFactor) * noise_factor*cf;

    for(size_t i=0; i< nFFTTicks/2 + 1; ++i)
    {
        if(!fPeakVec[i]) continue;
        
        // exponential noise spectrum
        fCorrelatedGen->fireArray(2,rnd_corr);
          
           // std::cout << " board " << board << " random corr" << rnd_corr[0] << std::endl;
//        pval = fNoiseHistVec[i] * ((1-fNoiseRand) + 2 * fNoiseRand*rnd_corr[0]) * scaleFactor;
        pval = fPeakVec[i] * ((1-fNoiseRand) + 2 * fNoiseRand*rnd_corr[0]) * scaleFactor;

        phase = rnd_corr[1] * 2. * TMath::Pi();
        
        TComplex tc(pval*cos(phase),pval*sin(phase));
        
        noiseFrequency.at(i) += tc;
    }
    
    // inverse FFT MCSignal
    fFFT->DoInvFFT(noiseFrequency, noise);
    
    noiseFrequency.clear();
    
    return;
}

void CorrelatedNoise::FindPeaks()
{
    //std::cout << " FindPeaks " << std::endl;

    float dThr=50;
    std::vector<double> dHeight;
    dHeight.resize(fNoiseHistVec.size(), 0.);

    dHeight[0]=0;
    fPeakVec.resize(fNoiseHistVec.size(),0);
    
    //std::cout << " before loop " << fNoiseHistVec.size() << std::endl;
    
    for(size_t histIdx = 1; histIdx < fNoiseHistVec.size(); histIdx++) {
        dHeight[histIdx]=fNoiseHistVec[histIdx]-fNoiseHistVec[histIdx-1];
    //    std::cout << " histidx " << histIdx << " dheight " << dHeight[histIdx] << std::endl;
    }
    
    for(size_t histIdx = 1; histIdx < fNoiseHistVec.size(); histIdx++)
        if(dHeight[histIdx]>dThr&&dHeight[histIdx+1]<dThr&&histIdx>5) fPeakVec[histIdx]=1;
    
    return;
}

void CorrelatedNoise::SelectContinuousSpectrum()
{
    fNoiseContVec.resize(fNoiseHistVec.size(), 0.);
    
    fWaveformTool->medianSmooth(fNoiseHistVec, fNoiseContVec, fMedianNumBins);
    
    std::vector<double> peakVec;
    
    peakVec.resize(fNoiseHistVec.size());
    
    std::transform(fNoiseHistVec.begin(),fNoiseHistVec.end(),fNoiseContVec.begin(),peakVec.begin(),std::minus<double>());
    
    std::vector<double>::iterator maxItr = std::max_element(peakVec.begin(),peakVec.end());
    
    double maxPeak   = *maxItr;
    double threshold = std::min(0.1 * maxPeak, 20.);
    
    fPeakVec.clear();
    fPeakVec.resize(peakVec.size(),0.);
    
    for(size_t idx = 0; idx < peakVec.size(); idx++)
    {
        if (peakVec[idx] > threshold) fPeakVec[idx] = peakVec[idx];
        
        if (fStoreHistograms)
        {
            double freq = 1.e6 * double(idx)/ (2. * fDetectorProperties->SamplingRate() * fPeakVec.size());
            
            fInputNoiseHist->Fill(freq,fNoiseHistVec.at(idx),1.);
            fMedianNoiseHist->Fill(freq,fNoiseContVec.at(idx),1.);
            fPeakNoiseHist->Fill(freq,peakVec.at(idx),1.);
        }
    }
    
    return;
}
    
void CorrelatedNoise::ExtractCorrelatedAmplitude(double& corrFactor, int board) const
{
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    //CLHEP::HepRandomEngine &engine_corr = rng->getEngine("cornoise");
    //CLHEP::RandFlat flat_corr(engine_corr,0,1);
    CLHEP::RandGeneral amp_corr(fCorrAmpDistVec.data(),fCorrAmpDistVec.size(),0);
    amp_corr.setTheSeed(board);
    double rnd_corr[1]      = {0.};
    
    amp_corr.fireArray(1,rnd_corr);
    
    double cfmedio=0.2287;
    corrFactor=rnd_corr[0]/cfmedio;
    //corrFactor=10;
    //  if(corrFactor>3)
    //   std::cout << " corr noise factor " << corrFactor << std::endl;
   
    // corrFactor=1;
}
    
void CorrelatedNoise::makeHistograms()
{
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(CorrelatedNoise)
}
