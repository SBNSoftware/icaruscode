////////////////////////////////////////////////////////////////////////
/// \file   SBNDataNoiseBoard.cc
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

class SBNDataNoiseBoard : IGenNoise
{
public:
    explicit SBNDataNoiseBoard(const fhicl::ParameterSet& pset);
    
    ~SBNDataNoiseBoard();
    
    void configure(const fhicl::ParameterSet& pset) override;
    
    void nextEvent() override;
   
    void generateNoise(CLHEP::HepRandomEngine& noise_engine,
                       CLHEP::HepRandomEngine& cornoise_engine,
                       icarusutil::TimeVec& noise,
                       detinfo::DetectorPropertiesData const&,
                       double noise_factor,
                       unsigned int wire,
                       unsigned int board) override;
    
private:
    void GenerateCorrelatedNoise(CLHEP::HepRandomEngine&, icarusutil::TimeVec&, double, unsigned int);
    void GenerateUncorrelatedNoise(CLHEP::HepRandomEngine&, icarusutil::TimeVec&, double, unsigned int);
    void GenNoise(std::function<void (double[])>&, const icarusutil::TimeVec&, icarusutil::TimeVec&, float);
    void ComputeRMSs();
    void makeHistograms();
void SampleCorrelatedRMSs() ;
void ExtractUncorrelatedRMS(float&, int, int) const;   
 void makeBoardHistos(int) ; 

    // Member variables from the fhicl file
    size_t                                      fPlane;
    int                                         fMedianNumBins;
    float                                       fNoiseRand;
    long                                        fCorrelatedSeed;
    long                                        fUncorrelatedSeed;
    std::vector<float>                                      fIncoherentNoiseFrac;
    bool                                        fStoreHistograms;
   std::string                              fInputNoiseHistFileName;
    std::string                                 fHistogramName;
    std::string                                 fCorrelatedHistogramName;
    std::string                                 fUncorrelatedHistogramName;
    std::string                                 fCorrelatedRMSHistoName;
    std::string                                 fUncorrelatedRMSHistoName;
    std::string                                 fTotalRMSHistoName;

    float corrFactors[216][4];  // this will be sparse, could use a map here I bet

    using WaveformTools = icarus_signal_processing::WaveformTools<icarusutil::SigProcPrecision>;

    WaveformTools                               fWaveformTool;

    // We'll recover the bin contents and store in a vector
    // with the likely false hope this will be faster...
    std::vector<icarusutil::TimeVec>                         fCoherentNoiseVec;       //< Input full noise frequency distribution
    std::vector<icarusutil::TimeVec>                          fIncoherentNoiseVec;       //< Input full noise frequency distribution
    
    double                                      fIncoherentNoiseRMS; //< RMS of full noise waveform
    double                                      fCoherentNoiseRMS;   //< RMS of full noise waveform

    // Container for doing the work
    icarusutil::FrequencyVec                    fNoiseFrequencyVec;
    
    // Keep track of seed initialization for uncorrelated noise
    bool                                        fNeedFirstSeed=true;
    
    // Histograms
    TProfile*                                   fInputNoiseHist;
   TH1D*                                       fMediaNoiseHist;
   TProfile*                                   fPeakNoiseHist;
  
   TH1D*                                      corrRMSHistPtr;
 TH1D*                                    uncorrRMSHistPtr;
   TH1D*                                      totalRMSHistPtr;

std::vector<float> totalRMS;
std::vector<float> rmsUnc;
std::vector<float> rmsCorr;
    
    // Keep instance of the eigen FFT
    Eigen::FFT<double>                          fEigenFFT;
    
};
    
//----------------------------------------------------------------------
// Constructor.
SBNDataNoiseBoard::SBNDataNoiseBoard(const fhicl::ParameterSet& pset)
{
    // Recover the configuration of the tool from the input fhicl file and set up
    configure(pset);
ComputeRMSs();
    
    // Output some histograms to catalogue what's been done
    makeHistograms();
std::cout << " after making histos " << std::endl;
}
    
SBNDataNoiseBoard::~SBNDataNoiseBoard()
{
}
    
void SBNDataNoiseBoard::configure(const fhicl::ParameterSet& pset)
{
std::cout << " configuring tool " << std::endl;
    // Recover the histogram used for noise generation
    fPlane                     = pset.get< size_t      >("Plane");
    fMedianNumBins             = pset.get< int         >("MedianNumBins");
    fNoiseRand                 = pset.get< float       >("NoiseRand");
    fCorrelatedSeed            = pset.get< long        >("CorrelatedSeed",1000);
    fUncorrelatedSeed          = pset.get< long        >("UncorrelatedSeed",5000);
 std::vector<float> noiseFracVec = pset.get< std::vector<float> >("IncoherentNoiseFraction", std::vector<float>());
 for(auto& noiseFrac : noiseFracVec) {
        fIncoherentNoiseFrac.push_back(float(noiseFrac));
    }
   // fIncoherentNoiseFrac       = pset.get< float       >("IncoherentNoiseFraction",0.5);
    fStoreHistograms           = pset.get< bool        >("StoreHistograms");

 //   std::vector<std::string> noiseToolParamSetVec = pset.get< std::vector<std::string> >("NoiseHistFileName", std::vector<std::string>());
 /*  for(auto& noiseToolParams : noiseToolParamSetVec) {
        fInputNoiseHistFileName.push_back(std::string(noiseToolParams));
    }*/
   fInputNoiseHistFileName    = pset.get< std::string >("NoiseHistFileName");
    fCorrelatedHistogramName   = pset.get< std::string >("CorrelatedHistogramName");
    fUncorrelatedHistogramName = pset.get< std::string >("UncorrelatedHistogramName");
    fCorrelatedRMSHistoName    = pset.get< std::string >("CorrelatedRMSHistoName");
    fUncorrelatedRMSHistoName  = pset.get< std::string >("UncorrelatedRMSHistoName");
    fTotalRMSHistoName         = pset.get< std::string >("TotalRMSHistoName");

std::cout << " end configuring tool " << std::endl;
    return;
}
void SBNDataNoiseBoard::makeBoardHistos(int board)
{
    // Initialize the work vector
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
    fNoiseFrequencyVec.resize(detProp.NumberTimeSamples(),std::complex<float>(0.,0.));

//for(auto& filename : fInputNoiseHistFileName) {
    // Set up to input the histogram with the overall noise spectrum
    std::string fullFileName;

    cet::search_path searchPath("FW_SEARCH_PATH");

    searchPath.find_file(fInputNoiseHistFileName, fullFileName);

    TFile inputFile(fullFileName.c_str(), "READ");
    
    if (!inputFile.IsOpen())
        throw cet::exception("NoiseFromHist::configure") << "Unable to open input file: " << fInputNoiseHistFileName << std::endl;
    
    TH1D* corrHistPtr = (TH1D*)inputFile.Get((fCorrelatedHistogramName+std::to_string(board)).c_str());    
    if (!corrHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fCorrelatedHistogramName+std::to_string(board) << std::endl;

  TH1D* uncorrHistPtr = (TH1D*)inputFile.Get((fUncorrelatedHistogramName+std::to_string(board)).c_str());    
    if (!uncorrHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fUncorrelatedHistogramName << std::endl;

  TH1D*  corrRMSHistPtr=(TH1D*)inputFile.Get((fCorrelatedRMSHistoName+std::to_string(board)).c_str());    
    if (!corrRMSHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fCorrelatedRMSHistoName << std::endl;

 TH1D* uncorrRMSHistPtr=(TH1D*)inputFile.Get((fCorrelatedRMSHistoName+std::to_string(board)).c_str()); 
    if (!uncorrRMSHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fUncorrelatedRMSHistoName << std::endl;

 TH1D* totalRMSHistPtr=(TH1D*)inputFile.Get((fTotalRMSHistoName+std::to_string(board)).c_str()); 
    if (!totalRMSHistPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fTotalRMSHistoName << std::endl;

    // Close the input file
    inputFile.Close();

icarusutil::TimeVec corvec;
    corvec.resize(corrHistPtr->GetNbinsX(), 0.);
    
    // Recover the bin contents into local vectors
    for(size_t histIdx = 0; histIdx < size_t(corrHistPtr->GetNbinsX()); histIdx++)
        corvec[histIdx] = corrHistPtr->GetBinContent(histIdx+1);
fCoherentNoiseVec.push_back(corvec);
icarusutil::TimeVec incvec;
incvec.resize(uncorrHistPtr->GetNbinsX(), 0.);
    
    // Recover the bin contents into local vectors
    for(size_t histIdx = 0; histIdx < size_t(uncorrHistPtr->GetNbinsX()); histIdx++)
        incvec[histIdx] = uncorrHistPtr->GetBinContent(histIdx+1);
fIncoherentNoiseVec.push_back(incvec);
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
       fMediaNoiseHist  = dir.make<TH1D>("MedNoise",  ";ADC", 100, -10., -10.);;
        fPeakNoiseHist    = dir.make<TProfile>("PeakNoise", ";freq(kHz)", numSamples, minFreq, maxFreq);;
       
    }
   //SampleCorrelatedRMSs();
}    
void SBNDataNoiseBoard::nextEvent()
{
    // We update the correlated seed because we want to see different noise event-by-event
    fCorrelatedSeed   = (333 * fCorrelatedSeed) % 900000000;
       SampleCorrelatedRMSs();
    return;
}

void SBNDataNoiseBoard::generateNoise(CLHEP::HepRandomEngine& engine_unc,
                                    CLHEP::HepRandomEngine& engine_corr,
                                    icarusutil::TimeVec&     noise,
                             detinfo::DetectorPropertiesData const&,
                                    double                  noise_factor,
                                    unsigned int            channel,
                                    unsigned int            board)
{
 makeBoardHistos(board);
//std::cout << " generating noise channel " << channel << std::endl;
   //GET THE GEOMETRY.
    art::ServiceHandle<geo::Geometry> geom;
    // get the WireID for this hit
          std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
          // for now, just take the first option returned from ChannelToWire
          geo::WireID wid  = wids[0];
          // We need to know the plane to look up parameters
         
          size_t cryostat=wid.Cryostat;
          size_t tpc=wid.TPC;
          size_t iWire=wid.Wire;
//std::cout << " generating noise cryostat " << cryostat << " tpc " << tpc << " wire " << iWire << std::endl;
/*int index=-1;
if(cryostat==0&&tpc<2) index=0;
if(cryostat==0&&tpc>1) index=1;
if(cryostat==1&&tpc<2) index=2;
if(cryostat==1&&tpc>1) index=3;*/
//std::cout << " index " << index << std::endl;
//std::cout << " totalrms size " << totalRMS.size() << std::endl;

//std::cout <<  " index " << index << " generating noise totalRMS " << totalRMS[index] << std::endl;
    // Define a couple of vectors to hold intermediate work
    icarusutil::TimeVec noise_unc(noise.size(),0.);
    icarusutil::TimeVec noise_corr(noise.size(),0.);
    
    // Make sure the work vector is size right with the output
    if (fNoiseFrequencyVec.size() != noise.size()) fNoiseFrequencyVec.resize(noise.size(),std::complex<float>(0.,0.));
 //   std::cout <<  " generating uncorrelated noise " << std::endl;
    // If applying incoherent noise call the generator
   GenerateUncorrelatedNoise(engine_unc,noise_unc,noise_factor,channel); 


//std::cout << " before corrfactors " << std::endl;
float cf=1;
//std::cout << " board" << board << " cf " << cf << std::endl;

   GenerateCorrelatedNoise(engine_corr, noise_corr, noise_factor*cf, board);

   // std::cout <<  " summing noise " << std::endl;
    // Take the noise as the simple sum of the two contributions
    std::transform(noise_unc.begin(),noise_unc.end(),noise_corr.begin(),noise.begin(),std::plus<float>());
    
 float mediaNoise=0;
 for(unsigned int jn=0;jn<noise.size();jn++) {
if(!cryostat&&!tpc&&!fPlane&&iWire<2) 
{
//std::cout << " jn " << jn << " noise sum " << noise.at(jn) << std::endl; 
//std::cout << " jn " << jn << " noise unc " << noise_unc.at(jn) << std::endl; 
//std::cout << " jn " << jn << " noise corr " << noise_corr.at(jn) << std::endl; 
}
  mediaNoise+=noise.at(jn);
}

mediaNoise/=(noise.size());
//std::cout << " media noise size " << noise.size() << std::endl;
//fMediaNoiseHist->Fill(mediaNoise);
//std::cout << " media noise " << mediaNoise << std::endl;
    return;
}
    
void SBNDataNoiseBoard::GenerateUncorrelatedNoise(CLHEP::HepRandomEngine& engine, icarusutil::TimeVec &noise, double noise_factor, unsigned int channel)
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
float cf=1;
//std::cout << " before extracting cf " << std::endl;
//ExtractUncorrelatedRMS(cf,channel,index);
    float  scaleFactor = cf*noise_factor;
   //std::cout << " fraction " << cf <<" unc scale Factor " << scaleFactor << " index " << index << std::endl;
//std::cout << "incnoisevec size " << fIncoherentNoiseVec.size() << std::endl;
    GenNoise(randGenFunc, fIncoherentNoiseVec[0], noise, scaleFactor);

    return;
}
    
void SBNDataNoiseBoard::GenerateCorrelatedNoise(CLHEP::HepRandomEngine& engine, icarusutil::TimeVec &noise, double noise_factor, unsigned int board)
{    
   
        // Set the engine seed to the board being considered
        engine.setSeed(fCorrelatedSeed+board,0);
        
        CLHEP::RandFlat noiseGen(engine,0,1);
        
        std::function<void (double[])> randGenFunc = [&noiseGen](double randArray[]){noiseGen.fireArray(2,randArray);};
    

        float scaleFactor = noise_factor;
    //      std::cout << " fraction " << fraction << " corr scale Factor " << scaleFactor << std::endl;
        GenNoise(randGenFunc, fCoherentNoiseVec[0], noise, scaleFactor);
    
    
    return;
}
    
void SBNDataNoiseBoard::GenNoise(std::function<void (double[])>& gen,const icarusutil::TimeVec& freqDist, icarusutil::TimeVec& noise, float scaleFactor)
{
    double rnd_corr[2] = {0.,0.};
   // std::cout << " noise size " << noise.size() << std::endl;
    // Build out the frequency vector
    for(size_t i=0; i< noise.size()/2; ++i)
    {
//std::cout << " i " << i << " freqdist " << freqDist[i] << std::endl;
        // exponential noise spectrum
        gen(rnd_corr);
     // if(i!=10) continue;
        float pval  = freqDist[i] * ((1-fNoiseRand) + 2 * fNoiseRand*rnd_corr[0]) * scaleFactor;
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

void SBNDataNoiseBoard::makeHistograms()
{
    
    return;
}
   
void SBNDataNoiseBoard::ComputeRMSs()
{


}
void SBNDataNoiseBoard::SampleCorrelatedRMSs() 
{

}
void SBNDataNoiseBoard::ExtractUncorrelatedRMS(float& cf, int channel, int index) const
{
}
    
 
DEFINE_ART_CLASS_TOOL(SBNDataNoiseBoard)
}
