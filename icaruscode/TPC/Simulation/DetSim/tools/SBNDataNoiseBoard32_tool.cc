////////////////////////////////////////////////////////////////////////
/// \file   SBNDataNoiseBoard32.cc
/// \author F. Varanini
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IGenNoise.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art_root_io/TFileService.h"

#include "larcore/Geometry/Geometry.h"

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

class SBNDataNoiseBoard32 : IGenNoise
{
public:
    explicit SBNDataNoiseBoard32(const fhicl::ParameterSet& pset);
    
    ~SBNDataNoiseBoard32();
    
    void configure(const fhicl::ParameterSet& pset) override ;
    
    void nextEvent() override ;
   
    void generateNoise(CLHEP::HepRandomEngine&,
                       CLHEP::HepRandomEngine& ,
                       CLHEP::HepRandomEngine& ,
                       icarusutil::TimeVec&,
                       detinfo::DetectorPropertiesData const&,
                       double,
                       const geo::PlaneID&,
                       unsigned int, int )  override;
    
private:
    void GenerateCorrelatedNoise(CLHEP::HepRandomEngine&, icarusutil::TimeVec&, double, unsigned int);
void Generate32Noise(CLHEP::HepRandomEngine&, icarusutil::TimeVec&, double, unsigned int, int);
    void GenerateUncorrelatedNoise(CLHEP::HepRandomEngine&, icarusutil::TimeVec&, double, unsigned int);
    void GenNoise(std::function<void (double[])>&, const icarusutil::TimeVec&, const icarusutil::TimeVec&, icarusutil::TimeVec&, float);
    void ComputeRMSs();
    void makeHistograms();
    void SampleCorrelatedRMSs() ;
    void ExtractUncorrelatedRMS(float&, int, int) const;   
    void makeBoardHistos(unsigned int) ; 

    // Member variables from the fhicl file
    size_t                                      fPlane;
    int                                         fMedianNumBins;
    float                                       fNoiseRand;
    long                                        fCorrelatedSeed;
    long                                        fUncorrelatedSeed;
    std::vector<float>                          fIncoherentNoiseFrac;
    bool                                        fStoreHistograms;
    std::string                                 fInputNoiseHistFileName;
    std::string                                 fHistogramName;
    std::string                                 fCorrelatedHistogramName;
    std::string                                 fUncorrelatedHistogramName;
std::string                                 f32HistogramName;
    std::string                                 fCorrelatedRMSHistoName;
    std::string                                 fUncorrelatedRMSHistoName;
    std::string                                 fTotalRMSHistoName;

    using WaveformTools = icarus_signal_processing::WaveformTools<icarusutil::SigProcPrecision>;

    WaveformTools                               fWaveformTool;

    // We'll recover the bin contents and store in a vector
    // with the likely false hope this will be faster...
    using BoardToNoiseVecMap = std::unordered_map<unsigned int, icarusutil::TimeVec>;

    BoardToNoiseVecMap                          fCoherentBoardToNoiseVecMap;
    BoardToNoiseVecMap                          fIncoherentBoardToNoiseVecMap;
    BoardToNoiseVecMap                          fCoherentBoardToNoiseWidthVecMap;
     BoardToNoiseVecMap                          fIncoherentBoardToNoiseWidthVecMap;
  BoardToNoiseVecMap                          f32BoardToNoiseVecMap;
    BoardToNoiseVecMap                          f32BoardToNoiseWidthVecMap;

    std::vector<icarusutil::TimeVec>            fCoherentNoiseVec;       //< Input full noise frequency distribution
    std::vector<icarusutil::TimeVec>            fIncoherentNoiseVec;       //< Input full noise frequency distribution

    using BoardToHistMap = std::unordered_map<unsigned int, TH1D*>;

    // Container for doing the work
    icarusutil::FrequencyVec                    fNoiseFrequencyVec;
    
    // Keep track of seed initialization for uncorrelated noise
    bool                                        fNeedFirstSeed=true;
    
    // Histograms
    TProfile*                                   fInputNoiseHist;
    TH1D*                                       fMediaNoiseHist;
    TProfile*                                   fPeakNoiseHist;

    TFile*                                      fHistogramFile;

    std::vector<float>                          totalRMS;
    std::vector<float>                          rmsUnc;
    std::vector<float>                          rmsCorr;
    
    // Keep instance of the eigen FFT
    Eigen::FFT<double>                          fEigenFFT;
    
};
    
//----------------------------------------------------------------------
// Constructor.
SBNDataNoiseBoard32::SBNDataNoiseBoard32(const fhicl::ParameterSet& pset)
{
    // Recover the configuration of the tool from the input fhicl file and set up
    configure(pset);
    ComputeRMSs();
    
    // Output some histograms to catalogue what's been done
    makeHistograms();
    std::cout << " after making histos " << std::endl;
}
    
SBNDataNoiseBoard32::~SBNDataNoiseBoard32()
{
    // Close the input file
    fHistogramFile->Close();
}
    
void SBNDataNoiseBoard32::configure(const fhicl::ParameterSet& pset)
{
std::cout << " configuring tool " << std::endl;
    // Recover the histogram used for noise generation
    fPlane                          = pset.get< size_t             >("Plane");
    fMedianNumBins                  = pset.get< int                >("MedianNumBins");
    fNoiseRand                      = pset.get< float              >("NoiseRand");
    fCorrelatedSeed                 = pset.get< long               >("CorrelatedSeed",1000);
    fUncorrelatedSeed               = pset.get< long               >("UncorrelatedSeed",5000);
    std::vector<float> noiseFracVec = pset.get< std::vector<float> >("IncoherentNoiseFraction", std::vector<float>());
    fStoreHistograms                = pset.get< bool               >("StoreHistograms");
    fInputNoiseHistFileName         = pset.get< std::string        >("NoiseHistFileName");
    fCorrelatedHistogramName        = pset.get< std::string        >("CorrelatedHistogramName");
    fUncorrelatedHistogramName      = pset.get< std::string        >("UncorrelatedHistogramName");
    f32HistogramName      = pset.get< std::string        >("Corr32HistogramName");
    fCorrelatedRMSHistoName         = pset.get< std::string        >("CorrelatedRMSHistoName");
    fUncorrelatedRMSHistoName       = pset.get< std::string        >("UncorrelatedRMSHistoName");
    fTotalRMSHistoName              = pset.get< std::string        >("TotalRMSHistoName");
 
    for(auto& noiseFrac : noiseFracVec) {
        fIncoherentNoiseFrac.push_back(float(noiseFrac));
    }

    std::string fullFileName;

    cet::search_path searchPath("FW_SEARCH_PATH");

    searchPath.find_file(fInputNoiseHistFileName, fullFileName);

    fHistogramFile = TFile::Open(fullFileName.c_str(), "READ");
    
    if (!fHistogramFile->IsOpen())
        throw cet::exception("NoiseFromHist::configure") << "Unable to open input file: " << fInputNoiseHistFileName << std::endl;



    return;
}

void SBNDataNoiseBoard32::makeBoardHistos(unsigned int board)
{
    if (fCoherentBoardToNoiseVecMap.find(board) == fCoherentBoardToNoiseVecMap.end())
    {
    //   std::cout << "opening input file: " << fCorrelatedHistogramName << std::endl;    

        TH1D* corrHistPtr = (TH1D*)fHistogramFile->Get((fCorrelatedHistogramName+std::to_string(board)).c_str());    
        if (!corrHistPtr)
            throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fCorrelatedHistogramName+std::to_string(board) << std::endl;

        icarusutil::TimeVec& noiseVec = fCoherentBoardToNoiseVecMap[board];
        icarusutil::TimeVec& widthVec = fCoherentBoardToNoiseWidthVecMap[board];

        noiseVec.resize(corrHistPtr->GetNbinsX(),0.);
        widthVec.resize(corrHistPtr->GetNbinsX(),0.);
  //for(size_t histIdx = 0; histIdx < size_t(corrHistPtr->GetNbinsX()); histIdx++) 
   // std::cout << " idx " << histIdx << " histo content " << corrHistPtr->GetBinContent(histIdx+1) << " histo width " << corrHistPtr->GetBinError(histIdx+1) << std::endl;
        for(size_t histIdx = 0; histIdx < size_t(corrHistPtr->GetNbinsX()); histIdx++) noiseVec[histIdx] = corrHistPtr->GetBinContent(histIdx+1);
        for(size_t histIdx = 0; histIdx < size_t(corrHistPtr->GetNbinsX()); histIdx++) widthVec[histIdx] = corrHistPtr->GetBinError(histIdx+1);
    }
    if (f32BoardToNoiseVecMap.find(board) == f32BoardToNoiseVecMap.end())
    {
    //   std::cout << "opening input file: " << fCorrelatedHistogramName << std::endl;    

        TH1D* corr32HistPtr = (TH1D*)fHistogramFile->Get((f32HistogramName+std::to_string(board)).c_str());    
        if (!corr32HistPtr)
            throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << f32HistogramName+std::to_string(board) << std::endl;

        icarusutil::TimeVec& noiseVec = f32BoardToNoiseVecMap[board];
        icarusutil::TimeVec& widthVec = f32BoardToNoiseWidthVecMap[board];

        noiseVec.resize(corr32HistPtr->GetNbinsX(),0.);
        widthVec.resize(corr32HistPtr->GetNbinsX(),0.);
  //for(size_t histIdx = 0; histIdx < size_t(corrHistPtr->GetNbinsX()); histIdx++) 
   // std::cout << " idx " << histIdx << " histo content " << corrHistPtr->GetBinContent(histIdx+1) << " histo width " << corrHistPtr->GetBinError(histIdx+1) << std::endl;
        for(size_t histIdx = 0; histIdx < size_t(corr32HistPtr->GetNbinsX()); histIdx++) noiseVec[histIdx] = corr32HistPtr->GetBinContent(histIdx+1);
        for(size_t histIdx = 0; histIdx < size_t(corr32HistPtr->GetNbinsX()); histIdx++) widthVec[histIdx] = corr32HistPtr->GetBinError(histIdx+1);
    }
    if(fIncoherentBoardToNoiseVecMap.find(board) == fIncoherentBoardToNoiseVecMap.end())
    {
        TH1D* uncorrHistPtr = (TH1D*)fHistogramFile->Get((fUncorrelatedHistogramName+std::to_string(board)).c_str());    
        if (!uncorrHistPtr)
            throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fUncorrelatedHistogramName << std::endl;

icarusutil::TimeVec& noiseVec = fIncoherentBoardToNoiseVecMap[board];
        icarusutil::TimeVec& widthVec = fIncoherentBoardToNoiseWidthVecMap[board];

        noiseVec.resize(uncorrHistPtr->GetNbinsX(),0.);
        widthVec.resize(uncorrHistPtr->GetNbinsX(),0.);
        for(size_t histIdx = 0; histIdx < size_t(uncorrHistPtr->GetNbinsX()); histIdx++) noiseVec[histIdx] = uncorrHistPtr->GetBinContent(histIdx+1);
        for(size_t histIdx = 0; histIdx < size_t(uncorrHistPtr->GetNbinsX()); histIdx++) widthVec[histIdx] = uncorrHistPtr->GetBinError(histIdx+1);
    }

    // Should we store hists?
    if (fStoreHistograms)
    {
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
        auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);

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
        fMediaNoiseHist   = dir.make<TH1D>("MedNoise",  ";ADC", 100, -10., -10.);;
        fPeakNoiseHist    = dir.make<TProfile>("PeakNoise", ";freq(kHz)", numSamples, minFreq, maxFreq);;
    }
    //SampleCorrelatedRMSs();
}

void SBNDataNoiseBoard32::nextEvent()
{
    // We update the correlated seed because we want to see different noise event-by-event
    fCorrelatedSeed   = (333 * fCorrelatedSeed) % 900000000;
    SampleCorrelatedRMSs();
std::cout << " end nextevent " << std::endl;
    return;
}

void SBNDataNoiseBoard32::generateNoise(CLHEP::HepRandomEngine&                engine_unc,
                                      CLHEP::HepRandomEngine&                engine_corr,
                                      CLHEP::HepRandomEngine&                engine_32,
                                      icarusutil::TimeVec&                   noise,
                                      detinfo::DetectorPropertiesData const&,
                                      double                                 noise_factor,
                                      const geo::PlaneID&                    planeID,
                                      unsigned int                           board, int hboard)
{
    makeBoardHistos(board);

    std::cout << " generating noise board " << board << std::endl;
    //GET THE GEOMETRY.
    art::ServiceHandle<geo::Geometry> geom;

    // We need to know the plane to look up parameters
//    size_t cryostat=planeID.Cryostat;
//    size_t tpc=planeID.TPC;

    //std::cout <<  " index " << index << " generating noise totalRMS " << totalRMS[index] << std::endl;
    // Define a couple of vectors to hold intermediate work
    icarusutil::TimeVec noise_unc(noise.size(),0.);
    icarusutil::TimeVec noise_corr(noise.size(),0.);
    icarusutil::TimeVec noise_32(noise.size(),0.);
    icarusutil::TimeVec noise_temp(noise.size(),0.);
    
    // Make sure the work vector is size right with the output
    if (fNoiseFrequencyVec.size() != noise.size()) fNoiseFrequencyVec.resize(noise.size(),std::complex<float>(0.,0.));

    // If applying incoherent noise call the generator
    GenerateUncorrelatedNoise(engine_unc,noise_unc,noise_factor, board); 

    GenerateCorrelatedNoise(engine_corr, noise_corr, noise_factor, board);
std::cout << " before generating 32 " << std::endl;

    Generate32Noise(engine_32, noise_32, noise_factor, board, hboard);

    std::cout <<  " summing noise " << std::endl;
    // Take the noise as the simple sum of the two contributions
    std::transform(noise_unc.begin(),noise_unc.end(),noise_corr.begin(),noise_temp.begin(),std::plus<float>());
    std::transform(noise_temp.begin(),noise_temp.end(),noise_32.begin(),noise.begin(),std::plus<float>());
    
//    float mediaNoise=0;

//    mediaNoise/=(noise.size());

    //std::cout << " media noise size " << noise.size() << std::endl;
    //fMediaNoiseHist->Fill(mediaNoise);
    //std::cout << " media noise " << mediaNoise << std::endl;
    return;
}
    
void SBNDataNoiseBoard32::GenerateUncorrelatedNoise(CLHEP::HepRandomEngine& engine, icarusutil::TimeVec &noise, double noise_factor, unsigned int board)
{
    // Here we aim to produce a waveform consisting of incoherent noise
    // Note that this is expected to be the dominate noise contribution
    // Check for seed initialization
    if (fNeedFirstSeed)
    {
        engine.setSeed(fUncorrelatedSeed,0);
        fNeedFirstSeed = false;
    }
    std::cout << " unc seed " << fUncorrelatedSeed << std::endl;
    // Get the generator
    CLHEP::RandFlat noiseGen(engine,0,1);
    
    std::function<void (double[])> randGenFunc = [&noiseGen](double randArray[]){noiseGen.fireArray(2,randArray);};

    GenNoise(randGenFunc, fIncoherentBoardToNoiseVecMap[board], fIncoherentBoardToNoiseWidthVecMap[board], noise, noise_factor);

    return;
}
    
void SBNDataNoiseBoard32::GenerateCorrelatedNoise(CLHEP::HepRandomEngine& engine, icarusutil::TimeVec &noise, double noise_factor, unsigned int board)
{    
   
    // Set the engine seed to the board being considered
    engine.setSeed(fCorrelatedSeed+board,0);
    
    CLHEP::RandFlat noiseGen(engine,0,1);
    
    std::function<void (double[])> randGenFunc = [&noiseGen](double randArray[]){noiseGen.fireArray(2,randArray);};

    GenNoise(randGenFunc, fCoherentBoardToNoiseVecMap[board], fCoherentBoardToNoiseWidthVecMap[board], noise, noise_factor);
    
    return;
}
void SBNDataNoiseBoard32::Generate32Noise(CLHEP::HepRandomEngine& engine, icarusutil::TimeVec &noise, double noise_factor, unsigned int board, int hboard)
{    
   
    // Set the engine seed to the board being considered
    engine.setSeed(fCorrelatedSeed+2*board+hboard,0);
   // std::cout << " after setting seed " << std::endl;
    CLHEP::RandFlat noiseGen(engine,0,1);
    
    std::function<void (double[])> randGenFunc = [&noiseGen](double randArray[]){noiseGen.fireArray(2,randArray);};

    GenNoise(randGenFunc, f32BoardToNoiseVecMap[board], f32BoardToNoiseWidthVecMap[board], noise, noise_factor);
  // std::cout << " noise " << noise.at(0) << std::endl;
    return;
}
    
void SBNDataNoiseBoard32::GenNoise(std::function<void (double[])>& gen,const icarusutil::TimeVec& freqDist,const icarusutil::TimeVec& width, icarusutil::TimeVec& noise, float scaleFactor)
{
    double rnd_corr[2] = {0.,0.};
    // std::cout << " noise size " << noise.size() << std::endl;
    // Build out the frequency vector
    for(size_t i=0; i< noise.size()/2; ++i)
    {
       // std::cout << " i " << i << " freqdist " << freqDist[i] << " width " << width[i] << std::endl;
        // exponential noise spectrum
        gen(rnd_corr);
        // if(i!=10) continue;
float relWidth=width[i]/freqDist[i];
//std::cout << " relWidth " << relWidth << std::endl;
       float pval  = freqDist[i] * ((1-relWidth) + 2 * relWidth*rnd_corr[0]) * scaleFactor;
//std::cout << " checking random scale " << ((1-relWidth) + 2 * relWidth*rnd_corr[0]) << std::endl; 
        float phase = rnd_corr[1] * 2. * M_PI;
        //   float pval  = freqDist[i]*scaleFactor ;
        //  float phase = 0;
        std::complex<float> tc(pval*cos(phase),pval*sin(phase));
        
        fNoiseFrequencyVec[i] = tc;
        //std::cout << " i " << i << " noise freqvec " << fNoiseFrequencyVec[i] << std::endl; 
    }
    
    // inverse FFT MCSignal
    fEigenFFT.inv(noise, fNoiseFrequencyVec);
    //    for(unsigned int jn=0;jn<noise.size();jn++) std::cout << " jn " << jn << " noise sum " << noise.at(jn) << std::endl; 
    //exit(22);
    //std::cout << " end gen noise " << std::endl;
    return;
}

void SBNDataNoiseBoard32::makeHistograms()
{
    
    return;
}
   
void SBNDataNoiseBoard32::ComputeRMSs()
{


}
void SBNDataNoiseBoard32::SampleCorrelatedRMSs() 
{

}
void SBNDataNoiseBoard32::ExtractUncorrelatedRMS(float& cf, int channel, int index) const
{
}
    
 
DEFINE_ART_CLASS_TOOL(SBNDataNoiseBoard32)
}
