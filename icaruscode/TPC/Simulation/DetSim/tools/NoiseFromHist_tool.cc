////////////////////////////////////////////////////////////////////////
/// \file   NoiseFromHist.cc
/// \author F. Varanini
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "IGenNoise.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Provenance/ModuleContext.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// FFT
#include "icarus_signal_processing/ICARUSFFT.h"

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH1D.h"
#include "TFile.h"

#include <fstream>

namespace icarus_tool
{

class NoiseFromHist : IGenNoise
{
public:
    explicit NoiseFromHist(const fhicl::ParameterSet& pset);
    
    ~NoiseFromHist();
    
    void configure(const fhicl::ParameterSet& pset)               override;
    
    void nextEvent() override  {return;};

    void generateNoise(CLHEP::HepRandomEngine&,
                       CLHEP::HepRandomEngine&,
                       icarusutil::TimeVec&, double, unsigned int) override;
    
private:
    // Member variables from the fhicl file
    double              fNoiseRand;
    std::string         fInputNoiseHistFileName;
    std::string         fHistogramName;
    
    double              fHistNormFactor;

    // We'll recover the bin contents and store in a vector
    // with the likely false hope this will be faster...
    std::vector<double> fNoiseHistVec;

    const detinfo::DetectorProperties*                fDetector;              //< Pointer to the detector properties
    std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>> fFFT;
};
    
//----------------------------------------------------------------------
// Constructor.
NoiseFromHist::NoiseFromHist(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
NoiseFromHist::~NoiseFromHist()
{
}
    
void NoiseFromHist::configure(const fhicl::ParameterSet& pset)
{
    // Recover the histogram used for noise generation
    fNoiseRand              = pset.get< double>("NoiseRand");
    fInputNoiseHistFileName = pset.get<std::string>("NoiseHistFileName");
    fHistogramName          = pset.get<std::string>("HistogramName");
    fHistNormFactor         = pset.get<double>("HistNormFactor");
    
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file(fInputNoiseHistFileName, fullFileName);
    
    TFile inputFile(fullFileName.c_str(), "READ");
    
    if (!inputFile.IsOpen())
        throw cet::exception("NoiseFromHist::configure") << "Unable to open input file: " << fInputNoiseHistFileName << std::endl;
    
    TH1D* histPtr = (TH1D*)inputFile.Get(fHistogramName.c_str());
    
    if (!histPtr)
        throw cet::exception("NoiseFromHist::configure") << "Unable to recover desired histogram: " << fHistogramName << std::endl;
    
    fNoiseHistVec.resize(histPtr->GetNbinsX(), 0.);
    
    for(size_t histIdx = 0; histIdx < size_t(histPtr->GetNbinsX()); histIdx++)
        fNoiseHistVec[histIdx] = histPtr->GetBinContent(histIdx+1);
    
    // Close the input file
    inputFile.Close();

    fDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Now set up our plans for doing the convolution
    int numberTimeSamples = fDetector->NumberTimeSamples();

    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(numberTimeSamples);
   
    return;
}

void NoiseFromHist::generateNoise(CLHEP::HepRandomEngine& engine,
                                  CLHEP::HepRandomEngine&,
                                  icarusutil::TimeVec& noise,
                                  double noise_factor,
                                  unsigned int channel)
{    
    CLHEP::RandFlat flat(engine,-1,1);
    
    size_t nFFTTicks = fDetector->NumberTimeSamples();
    
    if(noise.size() != nFFTTicks)
        throw cet::exception("SimWireICARUS")
        << "\033[93m"
        << "Frequency noise vector length must match FFT Ticks (FFT size)"
        << " ... " << noise.size() << " != " << nFFTTicks
        << "\033[00m"
        << std::endl;
    
    // noise in frequency space
    std::vector<std::complex<double>> noiseFrequency(nFFTTicks/2+1, 0.);
    
    double pval        = 0.;
    double phase       = 0.;
    double rnd[2]      = {0.};
    double scaleFactor = fHistNormFactor * noise_factor;
    
    // width of frequencyBin in kHz
    
    for(size_t i=0; i< nFFTTicks/2 + 1; ++i)
    {
        // exponential noise spectrum
        flat.fireArray(2,rnd,0,1);
        
        pval = fNoiseHistVec[i] * ((1-fNoiseRand) + 2 * fNoiseRand*rnd[0]) * scaleFactor;

        phase = rnd[1] * 2. * M_PI;

        std::complex<double> tc(pval*cos(phase),pval*sin(phase));

        noiseFrequency.at(i) += tc;
    }
    
    // inverse FFT MCSignal
    fFFT->inverseFFT(noiseFrequency, noise);
    
    noiseFrequency.clear();

    return;
}
    
DEFINE_ART_CLASS_TOOL(NoiseFromHist)
}
