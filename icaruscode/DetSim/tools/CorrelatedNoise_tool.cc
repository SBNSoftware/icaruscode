////////////////////////////////////////////////////////////////////////
/// \file   UncorrelatedNoise.cc
/// \author F. Varanini
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IGenNoise.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "lardata/Utilities/LArFFT.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH1D.h"
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

    void GenerateNoise(std::vector<float> &noise, double noise_factor, unsigned int wire) const override;
    void GenerateUncorrelatedNoise(std::vector<float> &noise, double noise_factor, unsigned int wire) const ;
    void GenerateCorrelatedNoise(std::vector<float> &noise, double noise_factor, unsigned int wire) const ;
    
    void ExtractCorrelatedAmplitude(double &) const;


    void SelectContinuousSpectrum() ;
    void FindPeaks() ;
    
private:
    // Member variables from the fhicl file
    double              fNoiseRand;
    std::string         fInputNoiseHistFileName;
    std::string         fHistogramName;
    std::string         fCorrAmpHistFileName;
    std::string         fCorrAmpHistogramName;
    
    double              fHistNormFactor;

    std::vector<int> peaks;

    // We'll recover the bin contents and store in a vector
    // with the likely false hope this will be faster...
    std::vector<double> fNoiseHistVec;
    std::vector<double> fNoiseContVec;
    
    double corrAmpDist[200];

};
    
//----------------------------------------------------------------------
// Constructor.
CorrelatedNoise::CorrelatedNoise(const fhicl::ParameterSet& pset)
{
    std::cout << " Correlated Noise " << std::endl;
    configure(pset);
    FindPeaks();
    SelectContinuousSpectrum();
}
    
CorrelatedNoise::~CorrelatedNoise()
{
}
    
void CorrelatedNoise::configure(const fhicl::ParameterSet& pset)
{
    // Recover the histogram used for noise generation
    fNoiseRand              = pset.get< double>("NoiseRand");
    fInputNoiseHistFileName = pset.get<std::string>("NoiseHistFileName");
    fHistogramName          = pset.get<std::string>("HistogramName");
    fHistNormFactor         = pset.get<double>("HistNormFactor");
    fCorrAmpHistFileName = pset.get<std::string>("CorrAmpHistFileName");
    fCorrAmpHistogramName          = pset.get<std::string>("CorrAmpHistogramName");
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
    
    fNoiseHistVec.resize(histPtr->GetNbinsX(), 0.);
    peaks.resize(histPtr->GetNbinsX(), 0.);

    
    for(size_t histIdx = 0; histIdx < size_t(histPtr->GetNbinsX()); histIdx++) {
        fNoiseHistVec[histIdx] = histPtr->GetBinContent(histIdx+1);
      //  std::cout << " histidx " << histIdx << " noisehist " << fNoiseHistVec[histIdx] << std::endl;
    }
    
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
    peaks.resize(histPtr->GetNbinsX(), 0.);
    
    
    for(size_t histIdx = 0; histIdx < size_t(histPtr->GetNbinsX()); histIdx++) {
        fNoiseHistVec[histIdx] = histPtr->GetBinContent(histIdx+1);
        //  std::cout << " histidx " << histIdx << " noisehist " << fNoiseHistVec[histIdx] << std::endl;
    }
    
    for(size_t histIdx = 0; histIdx < 200; histIdx++) {
        corrAmpDist[histIdx] = histPtr->GetBinContent(histIdx+1);
        //  std::cout << " histidx " << histIdx << " noisehist " << fNoiseHistVec[histIdx] << std::endl;
    }
    
    // Close the input file
    inputFile.Close();

    
    
    return;
}

void CorrelatedNoise::GenerateNoise(std::vector<float> &noise, double noise_factor, unsigned int channel) const
    {
        //std::cout << " generating noise " << std::endl;
        std::vector<float> noise_unc  = noise;
        std::vector<float> noise_corr = noise;
        
        GenerateUncorrelatedNoise(noise_unc,noise_factor,channel);
        GenerateCorrelatedNoise(noise_corr,noise_factor,channel);
        
        for(size_t j=0;j<noise_corr.size();j++) {
            noise[j]=noise_corr[j];
        if(j<10) std::cout << "channel " << channel << " j " << j << " noise_unc " << noise_unc[j] << std::endl;
        if(j<10) std::cout << "channel " << channel << " j " << j << " noise_corr " << noise_corr[j] << std::endl;
        if(j<10) std::cout << "channel " << channel << " j " << j << " noise " << noise[j] << std::endl;
        }
    }
void CorrelatedNoise::GenerateUncorrelatedNoise(std::vector<float> &noise, double noise_factor, unsigned int channel) const
    {
        //std::cout << " generating uncorrelated noise " << std::endl;

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    art::ServiceHandle<util::LArFFT>               fFFT;

    CLHEP::HepRandomEngine &engine_unc = rng->getEngine("noise");
        engine_unc.setSeed(channel,0);

    CLHEP::RandFlat flat_unc(engine_unc,-1,1);
 
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
    double rnd_unc[2] = {0.};

    double scaleFactor = fHistNormFactor * noise_factor;
    
    // uncorrelated, smooth frequency spectrum
     // std::cout << " before loop " <<nFFTTicks << std::endl;
    for(size_t i=0; i< nFFTTicks/2 + 1; ++i)
    {
        // exponential noise spectrum
        flat_unc.fireArray(2,rnd_unc,0,1);
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
    void CorrelatedNoise::GenerateCorrelatedNoise(std::vector<float> &noise, double noise_factor, unsigned int channel) const
    {
       // std::cout << " generate correlated noise " << std::endl;
        art::ServiceHandle<art::RandomNumberGenerator> rng;
        art::ServiceHandle<util::LArFFT>               fFFT;
        
        int board = channel/32;
        CLHEP::HepRandomEngine &engine_corr = rng->getEngine("cornoise");
        engine_corr.setSeed(board,0);
        CLHEP::RandGeneral amp_corr(corrAmpDist,200,1);
        double rnd_corr[2]      = {0.};
        size_t nFFTTicks = fFFT->FFTSize();
        
        if(noise.size() != nFFTTicks)
            throw cet::exception("SimWireICARUS")
            << "\033[93m"
            << "Frequency noise vector length must match FFT Ticks (FFT size)"
            << " ... " << noise.size() << " != " << nFFTTicks
            << "\033[00m"
            << std::endl;
        
        // noise in frequency space
        std::vector<TComplex> noiseFrequency(nFFTTicks/2+1, 0.);
        
        double pval        = 0.;
        double phase       = 0.;
        double scaleFactor = (1. - fHistNormFactor) * noise_factor;

        for(size_t i=0; i< nFFTTicks/2 + 1; ++i)
        {
            if(!peaks[i]) continue;
        
        // exponential noise spectrum
        amp_corr.fireArray(2,rnd_corr);
        
        pval = fNoiseHistVec[i] * ((1-fNoiseRand) + 2 * fNoiseRand*rnd_corr[0]) * scaleFactor;
        
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
    peaks[0]=0;

    //std::cout << " before loop " << fNoiseHistVec.size() << std::endl;
    
    for(size_t histIdx = 1; histIdx < fNoiseHistVec.size(); histIdx++) {
        dHeight[histIdx]=fNoiseHistVec[histIdx]-fNoiseHistVec[histIdx-1];
    //    std::cout << " histidx " << histIdx << " dheight " << dHeight[histIdx] << std::endl;
    }
    
    for(size_t histIdx = 1; histIdx < fNoiseHistVec.size(); histIdx++)
        if(dHeight[histIdx]>dThr&&dHeight[histIdx+1]<dThr&&histIdx>5) {
          peaks[histIdx]=1;
           // std::cout << " peak at bin " << histIdx << std::endl;
        }
    
}

void CorrelatedNoise::SelectContinuousSpectrum()
{
    fNoiseContVec.resize(fNoiseHistVec.size(), 0.);

    
for(size_t histIdx = 1; histIdx < fNoiseHistVec.size(); histIdx++)
    if(!peaks[histIdx])
        fNoiseContVec[histIdx]=fNoiseHistVec[histIdx];
    else
 fNoiseContVec[histIdx]=0.5*(fNoiseHistVec[histIdx-1]+fNoiseHistVec[histIdx-1]);
    
}
    
void CorrelatedNoise::ExtractCorrelatedAmplitude(double& corrFactor) const
    {
        art::ServiceHandle<art::RandomNumberGenerator> rng;
        CLHEP::HepRandomEngine &engine_corr = rng->getEngine("cornoise");
        CLHEP::RandFlat flat_corr(engine_corr,0,1);
        double rnd_corr[1]      = {0.};
        
        flat_corr.fireArray(1,rnd_corr,0,1);
        
        corrFactor=rnd_corr[0]*10;
        //std::cout << " corr noise factor " << corrFactor << std::endl;
    }
    
DEFINE_ART_CLASS_TOOL(CorrelatedNoise)
}
