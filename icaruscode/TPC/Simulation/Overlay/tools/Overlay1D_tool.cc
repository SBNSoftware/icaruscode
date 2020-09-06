////////////////////////////////////////////////////////////////////////
/// \file   Overlay1D.cc
/// \author F. Varanini
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IOverlay.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art_root_io/TFileService.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "icarus_signal_processing/WaveformTools.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH1F.h"
#include "TProfile.h"
#include "TFile.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include <fstream>

namespace icarus_tool
{

class Overlay1D : IOverlay
{
public:
    explicit Overlay1D(const fhicl::ParameterSet& pset);
    
    ~Overlay1D();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    void nextEvent() override;

    void generateNoise(CLHEP::HepRandomEngine& noise_engine,
                       CLHEP::HepRandomEngine& cornoise_engine,
                       icarusutil::TimeVec& noise,
                       double noise_factor,
                       unsigned int wire) override;
    
private:
    void GenerateOverlay1D(CLHEP::HepRandomEngine&, icarusutil::TimeVec&, double, unsigned int);
    void GenNoise(std::function<void (double[])>&, const icarusutil::TimeVec&, icarusutil::TimeVec&, double);
    void ExtractCorrelatedAmplitude(float&, int) const;
    void SelectContinuousSpectrum() ;
    void FindPeaks() ;
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
    std::string                                 fCorrAmpHistFileName;
    std::string                                 fCorrAmpHistogramName;

    using WaveformTools = icarus_signal_processing::WaveformTools<icarusutil::SigProcPrecision>;

    WaveformTools                               fWaveformTool;

    // We'll recover the bin contents and store in a vector
    // with the likely false hope this will be faster...
    icarusutil::TimeVec                         fNoiseHistVec;       //< Input full noise frequency distribution
    icarusutil::TimeVec                         fCoherentNoiseVec;   //< Peak distribution for coherent noise
    icarusutil::TimeVec                         fIncoherentNoiseVec; //< Incoherent frequency distribution
    icarusutil::TimeVec                         fCorrAmpDistVec;     //< Keep track of motherboard contributions
    
    icarusutil::SigProcPrecision                fIncoherentNoiseRMS; //< RMS of full noise waveform
    icarusutil::SigProcPrecision                fCoherentNoiseRMS;   //< RMS of full noise waveform

    // Container for doing the work
    icarusutil::FrequencyVec                    fNoiseFrequencyVec;
    
    // Keep track of seed initialization for uncorrelated noise
    bool                                        fNeedFirstSeed=true;
    
    // Histograms
    TProfile*                                   fInputNoiseHist;
    TProfile*                                   fMedianNoiseHist;
    TProfile*                                   fPeakNoiseHist;
    TProfile*                                   fCorAmpDistHist;
    
    // Keep instance of the eigen FFT
    Eigen::FFT<double>                          fEigenFFT;
    
};
    
//----------------------------------------------------------------------
// Constructor.
Overlay1D::Overlay1D(const fhicl::ParameterSet& pset)
{
    // Recover the configuration of the tool from the input fhicl file and set up
    configure(pset);
    
    // Now break out the coherent from the incoherent using the input overall spectrum
    SelectContinuousSpectrum();
    
    // Output some histograms to catalogue what's been done
    makeHistograms();
}
    
Overlay1D::~Overlay1D()
{
}
    
void Overlay1D::configure(const fhicl::ParameterSet& pset)
{
    // Recover the histogram used for noise generation
    fPlane                  = pset.get< size_t      >("Plane");
    fMedianNumBins          = pset.get< int         >("MedianNumBins");
    fNoiseRand              = pset.get< float       >("NoiseRand");
    fCorrelatedSeed         = pset.get< long        >("CorrelatedSeed",1000);
    fUncorrelatedSeed       = pset.get< long        >("UncorrelatedSeed",5000);
    fIncoherentNoiseFrac    = pset.get< float       >("IncoherentNoiseFraction",0.5);
    fStoreHistograms        = pset.get< bool        >("StoreHistograms");
    fInputNoiseHistFileName = pset.get< std::string >("NoiseHistFileName");
    fHistogramName          = pset.get< std::string >("HistogramName");
    fCorrAmpHistFileName    = pset.get< std::string >("CorrAmpHistFileName");
    fCorrAmpHistogramName   = pset.get< std::string >("CorrAmpHistogramName");
    
    // Initialize the work vector
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
    fNoiseFrequencyVec.resize(detProp.NumberTimeSamples(),std::complex<float>(0.,0.));

    // Set up to input the histogram with the overall noise spectrum
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
    
    // Now grab the histogram givng the observed "strenght" of the coherent noise
    searchPath.find_file(fCorrAmpHistFileName, corrAmpFileName);
    
    TFile corrAmpInputFile(corrAmpFileName.c_str(), "READ");
    
    if (!corrAmpInputFile.IsOpen())
        throw cet::exception("Overlay1D::configure") << "Unable to open input file: " << fCorrAmpHistFileName << std::endl;
    
    TH1D* corrAmpHistPtr = (TH1D*)corrAmpInputFile.Get(fCorrAmpHistogramName.c_str());
    
    if (!corrAmpHistPtr)
        throw cet::exception("Overlay1D::configure") << "Unable to recover desired histogram: " << fCorrAmpHistogramName << std::endl;
    
    // Close the input file
    corrAmpInputFile.Close();

    fNoiseHistVec.resize(histPtr->GetNbinsX(), 0.);
    fCorrAmpDistVec.resize(corrAmpHistPtr->GetNbinsX(),0.);
    
    // Recover the bin contents into local vectors
    for(size_t histIdx = 0; histIdx < size_t(histPtr->GetNbinsX()); histIdx++)
        fNoiseHistVec[histIdx] = histPtr->GetBinContent(histIdx+1);
    
    for(size_t histIdx = 0; histIdx < size_t(corrAmpHistPtr->GetNbinsX()); histIdx++)  // was 200
        fCorrAmpDistVec[histIdx] = corrAmpHistPtr->GetBinContent(histIdx+1);
    
    // Should we store hists?
    if (fStoreHistograms)
    {
        // Define histograms
        art::ServiceHandle<art::TFileService> tfs;
    
        art::TFileDirectory* histDirectory = tfs.get();
        
        // Make a directory for these histograms
        art::TFileDirectory dir = histDirectory->mkdir(Form("CorNoisePlane%1zu",fPlane));
        
        float sampleRate  = sampling_rate(clockData);
        float readOutSize = detProp.ReadOutWindowSize();
        float maxFreq     = 1.e6 / (2. * sampleRate);
        float minFreq     = 1.e6 / (2. * sampleRate * readOutSize);
        int   numSamples  = readOutSize / 2;
        
        fInputNoiseHist   = dir.make<TProfile>("InNoise",   ";freq(kHz)", numSamples, minFreq, maxFreq);
        fMedianNoiseHist  = dir.make<TProfile>("MedNoise",  ";freq(kHz)", numSamples, minFreq, maxFreq);;
        fPeakNoiseHist    = dir.make<TProfile>("PeakNoise", ";freq(kHz)", numSamples, minFreq, maxFreq);;
        
        fCorAmpDistHist   = dir.make<TProfile>("CorAmp",    ";Motherboard", fCorrAmpDistVec.size(),0.,fCorrAmpDistVec.size());
    }
    
    return;
}
    
void Overlay1D::nextEvent()
{
    // We update the correlated seed because we want to see different noise event-by-event
    fCorrelatedSeed   = (333 * fCorrelatedSeed) % 900000000;
    
    return;
}

void Overlay1D::generateNoise(CLHEP::HepRandomEngine& engine_unc,
                                    CLHEP::HepRandomEngine& engine_corr,
                                    icarusutil::TimeVec&    noise,
                                    double                  noise_factor,
                                    unsigned int            channel)
{
    // Define a couple of vectors to hold intermediate work
    icarusutil::TimeVec noise_unc(noise.size(),0.);
    icarusutil::TimeVec noise_corr(noise.size(),0.);
    
    // Make sure the work vector is size right with the output
    if (fNoiseFrequencyVec.size() != noise.size()) fNoiseFrequencyVec.resize(noise.size(),std::complex<float>(0.,0.));
    
    // If applying incoherent noise call the generator
    if (fIncoherentNoiseFrac > 0.) GenerateOverlay1D(engine_unc,noise_unc,noise_factor,channel);
    
    int board=channel/32;

    // If applying coherent noise call the generator
    if (fIncoherentNoiseFrac < 1.) GenerateOverlay1D(engine_corr, noise_corr, noise_factor, board);
    
    // Take the noise as the simple sum of the two contributions
    std::transform(noise_unc.begin(),noise_unc.end(),noise_corr.begin(),noise.begin(),std::plus<float>());
    
    return;
}
    
void Overlay1D::GenerateOverlay1D(CLHEP::HepRandomEngine& engine, icarusutil::TimeVec& noise, double noise_factor, unsigned int channel)
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

    double scaleFactor = fIncoherentNoiseFrac * noise_factor / fIncoherentNoiseRMS;
    
    GenNoise(randGenFunc, fIncoherentNoiseVec, noise, scaleFactor);

    return;
}
    
void Overlay1D::GenNoise(std::function<void (double[])>& gen,const icarusutil::TimeVec& freqDist, icarusutil::TimeVec& noise, double scaleFactor)
{
    double rnd_corr[2] = {0.,0.};
    
    // Build out the frequency vector
    for(size_t i=0; i< noise.size()/2; ++i)
    {
        // exponential noise spectrum
        gen(rnd_corr);
        
        double pval  = freqDist[i] * ((1-fNoiseRand) + 2 * fNoiseRand*rnd_corr[0]) * scaleFactor;
        double phase = rnd_corr[1] * 2. * M_PI;
        
        std::complex<double> tc(pval*cos(phase),pval*sin(phase));
        
        fNoiseFrequencyVec[i] = tc;
    }
    
    // inverse FFT MCSignal
    fEigenFFT.inv(noise, fNoiseFrequencyVec);
    
    return;
}

void Overlay1D::FindPeaks()
{
    float dThr=50;
    std::vector<float> dHeight;
    dHeight.resize(fNoiseHistVec.size(), 0.);

    dHeight[0]=0;
    fCoherentNoiseVec.resize(fNoiseHistVec.size(),0);
    
    for(size_t histIdx = 1; histIdx < fNoiseHistVec.size(); histIdx++) {
        dHeight[histIdx]=fNoiseHistVec[histIdx]-fNoiseHistVec[histIdx-1];
    }
    
    for(size_t histIdx = 1; histIdx < fNoiseHistVec.size(); histIdx++)
        if(dHeight[histIdx]>dThr&&dHeight[histIdx+1]<dThr&&histIdx>5) fCoherentNoiseVec[histIdx]=1;
    
    return;
}

void Overlay1D::SelectContinuousSpectrum()
{
    // In the below we take the overall noise spectrum and perform a long bin smoothing to get its
    // general shape. We can then subtract this from the input spectrum to return a spectrum containing
    // the "spikes" which are associated to the coherent noise.
    // Generally, the inchorent noise model is then that with the spikes subtracted,
    // the coherent model is derived from the spikes...
    fIncoherentNoiseVec.resize(fNoiseHistVec.size(), 0.);
    
    // This does a median smoothing based on the number of bins we input
    fWaveformTool.medianSmooth(fNoiseHistVec, fIncoherentNoiseVec, fMedianNumBins);
    
    icarusutil::TimeVec peakVec(fNoiseHistVec.size());
    
    // Subtract the smoothed from the input to get the peaks
    std::transform(fNoiseHistVec.begin(),fNoiseHistVec.end(),fIncoherentNoiseVec.begin(),peakVec.begin(),std::minus<float>());
    
    // Try to clean up the "peak" vector a bit
    icarusutil::TimeVec::iterator maxItr = std::max_element(peakVec.begin(),peakVec.end());
    
    float maxPeak   = *maxItr;
    float threshold = std::min(0.1 * maxPeak, 20.);

    fCoherentNoiseVec.clear();
    fCoherentNoiseVec.resize(peakVec.size(),0.);
    
    auto const samplingRate = sampling_rate(art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob());
    for(size_t idx = 0; idx < peakVec.size(); idx++)
    {
        if (peakVec[idx] > threshold) fCoherentNoiseVec[idx] = peakVec[idx];
        
        if (fStoreHistograms)
        {
            float freq = 1.e6 * float(idx)/ (2. * samplingRate * fCoherentNoiseVec.size());
            
            fInputNoiseHist->Fill(freq,fNoiseHistVec.at(idx),1.);
            fMedianNoiseHist->Fill(freq,fIncoherentNoiseVec.at(idx),1.);
            fPeakNoiseHist->Fill(freq,peakVec.at(idx),1.);
        }
    }
    
    // Let's get the rms we expect from the incoherent noise contribution to the waveform
    // A couple of ways to do this, let's basically invert the frequency spectrum to
    // produce a waveform and then get the rms from that
    std::function<void (double[])> randGenFunc = [](double randArray[]){randArray[0]=0.5; randArray[1]=0.5;};
    
    icarusutil::TimeVec waveNoise(fIncoherentNoiseVec.size());
    float               scaleFactor = 1.;
    
    GenNoise(randGenFunc, fIncoherentNoiseVec, waveNoise, scaleFactor);
    
    // Now get the details...
    icarusutil::SigProcPrecision nSig(3.);
    icarusutil::SigProcPrecision mean,rmsTrunc;
    int                          nTrunc;
    int                          range;
    
    // Use the waveform tool to recover the full rms
    fWaveformTool.getTruncatedMean(waveNoise, mean, nTrunc, range);
    fWaveformTool.getTruncatedRMS(waveNoise, nSig, fIncoherentNoiseRMS, rmsTrunc, nTrunc);
    
    // Do the same for the coherent term
    GenNoise(randGenFunc, fCoherentNoiseVec, waveNoise, scaleFactor);
    
    fWaveformTool.getTruncatedMean(waveNoise, mean, nTrunc, range);
    fWaveformTool.getTruncatedRMS(waveNoise, nSig, fCoherentNoiseRMS, rmsTrunc, nTrunc);

    if (fStoreHistograms)
    {
        for(size_t idx = 0; idx < fCorrAmpDistVec.size(); idx++)
            fCorAmpDistHist->Fill(idx, fCorrAmpDistVec[idx], 1.);
    }

    return;
}
    
void Overlay1D::ExtractCorrelatedAmplitude(float& corrFactor, int board) const
{
    CLHEP::RandGeneral amp_corr(fCorrAmpDistVec.data(),fCorrAmpDistVec.size(),0);
    amp_corr.setTheSeed(board+1);
    double rnd_corr[1] = {0.};
    
    amp_corr.fireArray(1,rnd_corr);
    
    float cfmedio=0.2287;
    corrFactor=rnd_corr[0]/cfmedio;
}
    
void Overlay1D::makeHistograms()
{
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(Overlay1D)
}
