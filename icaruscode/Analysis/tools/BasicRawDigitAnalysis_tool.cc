#include "icaruscode/Analysis/tools/IRawDigitHistogramTool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"
#include "icaruscode/TPC/Utilities/SignalShapingICARUSService_service.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

#include "icaruscode/TPC/SignalProcessing/RawDigitFilter/Algorithms/RawDigitCharacterizationAlg.h"
#include "icaruscode/TPC/Utilities/tools/SignalProcessingDefs.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"

#include "icarus_signal_processing/Filters/ICARUSFFT.h"
#include "icarus_signal_processing/WaveformTools.h"

#include <cmath>
#include <algorithm>
#include <vector>

namespace BasicRawDigitAnalysis
{
    ////////////////////////////////////////////////////////////////////////
    //
    // Class:       BasicRawDigitAnalysis
    // Module Type: producer
    // File:        BasicRawDigitAnalysis.h
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

class BasicRawDigitAnalysis : virtual public IRawDigitHistogramTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit BasicRawDigitAnalysis(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~BasicRawDigitAnalysis();
    
    // provide for initialization
    void configure(fhicl::ParameterSet const & pset) override;

    /**
     *  @brief Interface for initializing the histograms to be filled
     *
     *  @param TFileService   handle to the TFile service
     *  @param string         subdirectory to store the hists in
     */
    void initializeHists(detinfo::DetectorClocksData const& clockData,
                         detinfo::DetectorPropertiesData const& detProp,
                         art::ServiceHandle<art::TFileService>&, const std::string&) override;
    
    /**
     *  @brief Interface for method to executve at the end of run processing
     *
     *  @param int            number of events to use for normalization
     */
    void endJob(int numEvents) override;
    
    /**
     *  @brief Interface for filling histograms
     */
    void fillHistograms(const detinfo::DetectorClocksData& clockData,
                        const detinfo::DetectorPropertiesData& detProp,
                        const IRawDigitHistogramTool::RawDigitPtrVec&,
                        const IRawDigitHistogramTool::SimChannelMap&) const override;
    
private:
    void filterFFT(const detinfo::DetectorClocksData& clockData,
                   const detinfo::DetectorPropertiesData& detProp,
                   std::vector<short>&, raw::ChannelID_t, size_t, size_t, float, bool) const;

    // Fcl parameters.
    std::vector<size_t>                  fLoWireByPlane;    ///< Low wire for individual wire histograms
    std::vector<size_t>                  fHiWireByPlane;    ///< Hi wire for individual wire histograms
    std::vector<std::string>             fFFTFitFuncVec;    ///< Function definitions for fitting the average FFT power spectra
    std::vector<std::vector<double>>     fParameterVec;     ///< Initial parameters for fit function

    // Pointers to the histograms we'll create.
    std::vector<TH1D*>                   fTruncMeanHist;
    std::vector<TH1D*>                   fTruncRmsHist;
    std::vector<TH1D*>                   fFullRmsHist;

    std::vector<std::vector<TProfile*>>  fFFTPowerVec;
    std::vector<std::vector<TProfile*>>  fFFTPowerDerivVec;
    std::vector<std::vector<TProfile*>>  fFFTRealVec;
    std::vector<std::vector<TProfile*>>  fFFTImaginaryVec;
    std::vector<std::vector<TProfile*>>  fSmoothPowerVec;
    
    std::vector<TProfile*>               fAveFFTPowerVec;
    std::vector<TProfile*>               fConvFFTPowerVec;
    std::vector<TProfile*>               fConvKernelVec;
    std::vector<TProfile*>               fFilterFuncVec;
    std::vector<TProfile*>               fAveFFTPowerDerivVec;
    std::vector<TProfile*>               fAveFFTRealVec;
    std::vector<TProfile*>               fAveFFTImaginaryVec;
    std::vector<TProfile*>               fAveSmoothPowerVec;

    caldata::RawDigitCharacterizationAlg     fCharacterizationAlg;

    using WaveformTools = icarus_signal_processing::WaveformTools<double>;
    
    WaveformTools                            fWaveformTool;

    using FFTPointer = std::unique_ptr<icarus_signal_processing::ICARUSFFT<double>>;

    FFTPointer                               fFFT;                   //< Object to handle thread safe FFT

    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::WireReadoutGeom&                fChannelMapAlg;        ///< pointer to ChannelMapAlg
    icarusutil::SignalShapingICARUSService&  fSignalServices;       ///< The signal shaping service
    const lariov::DetPedestalProvider&       fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
BasicRawDigitAnalysis::BasicRawDigitAnalysis(fhicl::ParameterSet const & pset) :
    fCharacterizationAlg(pset.get<fhicl::ParameterSet>("CharacterizationAlg")),
    fChannelMapAlg(art::ServiceHandle<geo::WireReadout const>()->Get()),
    fSignalServices(*art::ServiceHandle<icarusutil::SignalShapingICARUSService>()),
    fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())
{
    // Now set up our plans for doing the convolution
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
    int numberTimeSamples = detProp.NumberTimeSamples();

    fFFT = std::make_unique<icarus_signal_processing::ICARUSFFT<double>>(numberTimeSamples);
    
    configure(pset);
    
    // Report.
    mf::LogInfo("BasicRawDigitAnalysis") << "BasicRawDigitAnalysis configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
BasicRawDigitAnalysis::~BasicRawDigitAnalysis()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void BasicRawDigitAnalysis::configure(fhicl::ParameterSet const & pset)
{
    fLoWireByPlane  = pset.get<std::vector<size_t>>             ("LoWireByPlane",                         std::vector<size_t>()={0,0,0});
    fHiWireByPlane  = pset.get<std::vector<size_t>>             ("HiWireByPlane",                   std::vector<size_t>()={100,100,100});
    fFFTFitFuncVec  = pset.get<std::vector<std::string>>        ("FFTFunctionVec",             std::vector<std::string>()={"1","1","1"});
    fParameterVec   = pset.get<std::vector<std::vector<double>>>("FFTFuncParamsVec", std::vector<std::vector<double>>() = {{1},{1},{1}});

    //const fhicl::ParameterSet& waveformParamSet = pset.get<fhicl::ParameterSet>("WaveformTool");
}

//----------------------------------------------------------------------------
/// Begin job method.
void BasicRawDigitAnalysis::initializeHists(detinfo::DetectorClocksData const& clockData,
                                            detinfo::DetectorPropertiesData const& detProp,
                                            art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());
    
    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    
    // hijack hists here
    double sampleRate  = sampling_rate(clockData);
    double readOutSize = detProp.ReadOutWindowSize();
    double maxFreq     = 1.e6 / (2. * sampleRate);
    size_t numSamples  = readOutSize / 2;
    
    fTruncMeanHist.resize(3);
    fTruncRmsHist.resize(3);
    fFullRmsHist.resize(3);

    fFFTPowerVec.resize(3);
    fFFTPowerDerivVec.resize(3);
    fFFTRealVec.resize(3);
    fFFTImaginaryVec.resize(3);
    fSmoothPowerVec.resize(3);
    
    fAveFFTPowerVec.resize(3);
    fConvFFTPowerVec.resize(3);
    fConvKernelVec.resize(3);
    fFilterFuncVec.resize(3);
    fAveFFTPowerDerivVec.resize(3);
    fAveFFTRealVec.resize(3);
    fAveFFTImaginaryVec.resize(3);
    fAveSmoothPowerVec.resize(3);

    for(size_t plane = 0; plane < fChannelMapAlg.Nplanes(); plane++)
    {
        size_t numHists = fHiWireByPlane[plane] - fLoWireByPlane[plane];
        
        fFFTPowerVec[plane].resize(numHists);
        fFFTPowerDerivVec[plane].resize(numHists);
        fFFTRealVec[plane].resize(numHists);
        fFFTImaginaryVec[plane].resize(numHists);
        fSmoothPowerVec[plane].resize(numHists);

        for(size_t idx = 0; idx < 20; idx++)
        {
            std::string histName = "FFTPower_" + std::to_string(plane) + "-" + std::to_string(idx);
        
            fFFTPowerVec[plane][idx] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, 0., maxFreq, 0., 10000.);
        
            histName = "FFTPowerDeriv_" + std::to_string(plane) + "-" + std::to_string(idx);
        
            fFFTPowerDerivVec[plane][idx] = dir.make<TProfile>(histName.c_str(),  "Power Deriv;kHz;Power", numSamples, 0., maxFreq, -500., 500.);
        
            histName = "FFTReal_" + std::to_string(plane) + "-" + std::to_string(idx);
        
            fFFTRealVec[plane][idx] = dir.make<TProfile>(histName.c_str(),  "Real values;kHz;Power", numSamples, 0., maxFreq, -10000., 10000.);
        
            histName = "FFTImaginary_" + std::to_string(plane) + "-" + std::to_string(idx);
        
            fFFTImaginaryVec[plane][idx] = dir.make<TProfile>(histName.c_str(),  "Imaginary values;kHz;Power", numSamples, 0., maxFreq, -10000., 10000.);
        
            histName = "SmoothPWR_" + std::to_string(plane) + "-" + std::to_string(idx);
        
            fSmoothPowerVec[plane][idx] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, 0., maxFreq, 0., 10000.);
        }
        
        std::string histName = "AveFFTPower_" + std::to_string(plane);
        
        fAveFFTPowerVec[plane] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, 0., maxFreq, 0., 1000.);
        
        histName = "ConvFFTPower_" + std::to_string(plane);
        
        fConvFFTPowerVec[plane] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, 0., maxFreq, 0., 1000.);
        
        histName = "ConvKernel_" + std::to_string(plane);
        
        fConvKernelVec[plane] = dir.make<TProfile>(histName.c_str(),  "Convolution Kernel;kHz;Power", numSamples, 0., maxFreq, 0., 1000.);
        
        histName = "FilterFunc_" + std::to_string(plane);
        
        fFilterFuncVec[plane] = dir.make<TProfile>(histName.c_str(),  "Filter Function;kHz;Power", numSamples, 0., maxFreq, 0., 1000.);

        histName = "AveFFTPowerDeriv_" + std::to_string(plane);

        fAveFFTPowerDerivVec[plane] = dir.make<TProfile>(histName.c_str(),  "Power Deriv;kHz;Power", numSamples, 0., maxFreq, -500., 500.);
        
        histName = "AveFFTReal_" + std::to_string(plane);
        
        fAveFFTRealVec[plane] = dir.make<TProfile>(histName.c_str(),  "Real values;kHz;Power", numSamples, 0., maxFreq, -10000., 1000.);
        
        histName = "AveFFTImaginary_" + std::to_string(plane);
        
        fAveFFTImaginaryVec[plane] = dir.make<TProfile>(histName.c_str(),  "Imaginary values;kHz;Power", numSamples, 0., maxFreq, -1000., 1000.);
        
        histName = "AveSmoothPWR_" + std::to_string(plane);
        
        fAveSmoothPowerVec[plane] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, 0., maxFreq, 0., 1000.);
        
        histName = "TruncMean_" + std::to_string(plane);
        
        fTruncMeanHist[plane] = dir.make<TH1D>(histName.c_str(), ";ADC", 200, -50., 50.);
        
        histName = "TruncRMS_" + std::to_string(plane);
        
        fTruncRmsHist[plane] = dir.make<TH1D>(histName.c_str(), ";ADC", 100, 0., 20.);
        
        histName = "FullRMS_" + std::to_string(plane);
        
        fFullRmsHist[plane] = dir.make<TH1D>(histName.c_str(), ";ADC", 100, 0., 20.);
        
        // Need a channel...
        raw::ChannelID_t channel = fChannelMapAlg.PlaneWireToChannel(geo::WireID(0, 0, plane, 0));
        
        // Recover the filter from signal shaping services...
        const icarusutil::FrequencyVec& response = fSignalServices.GetResponse(channel).getConvKernel();
        const icarusutil::FrequencyVec& filter   = fSignalServices.GetResponse(channel).getFilter()->getResponseVec();
        
        for(size_t idx = 0; idx < numSamples; idx++)
        {
            double freq = 1.e6 * double(idx)/ (sampleRate * readOutSize);
            fConvKernelVec[plane]->Fill(freq, std::abs(response.at(idx)), 1.);
            fFilterFuncVec[plane]->Fill(freq, std::abs(filter.at(idx)), 1.);
        }
    }

    return;
}
    
void BasicRawDigitAnalysis::fillHistograms(const detinfo::DetectorClocksData& clockData,
                                           const detinfo::DetectorPropertiesData& detProp,
                                           const IRawDigitHistogramTool::RawDigitPtrVec& rawDigitPtrVec,
                                           const IRawDigitHistogramTool::SimChannelMap&  wireReadout) const
{
    // Sadly, the RawDigits come to us in an unsorted condition which is not optimal for
    // what we want to do here. So we make a vector of pointers to the input raw digits and sort them
    std::vector<const raw::RawDigit*> rawDigitVec;
    
    // Ugliness to fill the pointer vector...
    for(size_t idx = 0; idx < rawDigitPtrVec.size(); idx++) rawDigitVec.push_back(rawDigitPtrVec.at(idx).get());
    
    // Sort (use a lambda to sort by channel id)
    std::sort(rawDigitVec.begin(),rawDigitVec.end(),[](const raw::RawDigit* left, const raw::RawDigit* right) {return left->Channel() < right->Channel();});
    
    // Ok, to do the correlated noise removal we are going to need a rather impressive data structure...
    // Because we need to unpack each wire's data, we will need to "explode" it out into a data structure
    // here... with the good news that we'll release the memory at the end of the module so should not
    // impact downstream processing (I hope!).
    // What we are going to do is make a vector over planes of vectors over wires of vectors over time samples
    //std::vector<RawDigitVector> rawDataWireTimeVec;
    std::vector<caldata::RawDigitVector> rawDataWireTimeVec;
    std::vector<float>                   truncMeanWireVec;
    std::vector<float>                   truncRmsWireVec;
    std::vector<short>                   meanWireVec;
    std::vector<short>                   medianWireVec;
    std::vector<short>                   modeWireVec;
    std::vector<float>                   skewnessWireVec;
    std::vector<float>                   fullRmsWireVec;
    std::vector<short>                   minMaxWireVec;
    std::vector<float>                   neighborRatioWireVec;
    std::vector<float>                   pedCorWireVec;
    
    // We're stealing this from the raw digit filter, here we are going to set "number wires to group" to 1
    size_t numWiresToGroup(1);
    
    rawDataWireTimeVec.resize(numWiresToGroup);
    truncMeanWireVec.resize(numWiresToGroup);
    truncRmsWireVec.resize(numWiresToGroup);
    meanWireVec.resize(numWiresToGroup);
    medianWireVec.resize(numWiresToGroup);
    modeWireVec.resize(numWiresToGroup);
    skewnessWireVec.resize(numWiresToGroup);
    fullRmsWireVec.resize(numWiresToGroup);
    minMaxWireVec.resize(numWiresToGroup);
    neighborRatioWireVec.resize(numWiresToGroup);
    pedCorWireVec.resize(numWiresToGroup);

    // Commence looping over raw digits
    for(const auto& rawDigit : rawDigitVec)
    {
        raw::ChannelID_t channel = rawDigit->Channel();
        
        bool goodChan(true);
        
        // The below try-catch block may no longer be necessary
        // Decode the channel and make sure we have a valid one
        std::vector<geo::WireID> wids;
        try {
            wids = fChannelMapAlg.ChannelToWire(channel);
        }
        catch(...)
        {
            //std::cout << "===>> Found illegal channel with id: " << channel << std::endl;
            goodChan = false;
        }
        
        if (!goodChan) continue;
        
        // Recover plane and wire in the plane
        unsigned int plane = wids[0].Plane;
        unsigned int wire  = wids[0].Wire;
        
        unsigned int dataSize = rawDigit->Samples();
        
        if (dataSize < 1)
        {
            std::cout << "****>> Found zero length raw digit buffer, channel: " << channel << ", plane: " << plane << ", wire: " << wire << std::endl;
            continue;
        }
        
        // If MC, does this channel have signal?
        bool hasSignal = wireReadout.find(channel) != wireReadout.end();

        // vector holding uncompressed adc values
        std::vector<short>& rawadc = rawDataWireTimeVec[0];
        
        if (rawadc.size() != dataSize) rawadc.resize(dataSize);
            
        // And now uncompress
        raw::Uncompress(rawDigit->ADCs(), rawadc, rawDigit->Compression());
        
        // Recover the database version of the pedestal
        float pedestal = fPedestalRetrievalAlg.PedMean(channel);
        
        filterFFT(clockData, detProp, rawadc, channel, plane, wire, pedestal, hasSignal);
        
        // Only rest if no signal on wire
        if (!hasSignal)
        {
            // Get the kitchen sink
            fCharacterizationAlg.getWaveformParams(rawadc,
                                                   channel,
                                                   plane,
                                                   wire,
                                                   truncMeanWireVec[0],
                                                   truncRmsWireVec[0],
                                                   meanWireVec[0],
                                                   medianWireVec[0],
                                                   modeWireVec[0],
                                                   skewnessWireVec[0],
                                                   fullRmsWireVec[0],
                                                   minMaxWireVec[0],
                                                   neighborRatioWireVec[0],
                                                   pedCorWireVec[0]);
        
            // Now fill histograms...
            fTruncMeanHist[plane]->Fill(truncMeanWireVec[0] - pedestal, 1.);
            fTruncRmsHist[plane]->Fill(truncRmsWireVec[0], 1.);
            fFullRmsHist[plane]->Fill(fullRmsWireVec[0], 1.);
        }
    }
    
    return;
}
    
void BasicRawDigitAnalysis::filterFFT(detinfo::DetectorClocksData const& clockData,
                                      detinfo::DetectorPropertiesData const& detProp,
                                      std::vector<short>& rawadc,
                                      raw::ChannelID_t channel,
                                      size_t plane, size_t wire, float pedestal, bool hasSignal) const
{
    double sampleRate  = sampling_rate(clockData);
    double readOutSize = detProp.ReadOutWindowSize();
    //       double binSize     = sampleFreq / readOutSize;
    
    // Step one is to setup and then get the FFT transform of the input waveform
    size_t fftDataSize = rawadc.size();
    size_t halfFFTDataSize = fftDataSize / 2 + 1;
    
    icarusutil::TimeVec inputVec(fftDataSize, 0.);
    icarusutil::TimeVec powerVec(fftDataSize, 0.);

    std::transform(rawadc.begin(),rawadc.end(),inputVec.begin(),[pedestal](const auto& val){return double(val) - pedestal;});

    fFFT->getFFTPower(inputVec,powerVec);

    // Fill any individual wire histograms we want to look at
    if (wire >= fLoWireByPlane[plane] && wire < fHiWireByPlane[plane])
    {
        // Fill the power spectrum histogram
        for(size_t idx = 0; idx < halfFFTDataSize; idx++)
        {
            double freq = 1.e6 * double(idx + 1)/ (sampleRate * readOutSize);
            fFFTPowerVec[plane][wire-fLoWireByPlane[plane]]->Fill(freq, std::min(powerVec[idx],999.), 1.);
        }
    }
    
    for(size_t idx = 0; idx < halfFFTDataSize; idx++)
    {
        double freq = 1.e6 * double(idx + 1)/ (sampleRate * readOutSize);
        fAveFFTPowerVec[plane]->Fill(freq, std::min(powerVec[idx],999.), 1.);
    }

    // Idea here is to run through the power spectrum and keep a running average of the n bins around the current bin
    size_t numBinsToAve(9);  // number bins either size of current bin
    size_t lowestBin(3);    //(275);   // Go no lower than this?
    
    size_t currentBin(halfFFTDataSize - numBinsToAve - 1);
    
    fWaveformTool.triangleSmooth(powerVec, powerVec);
    
    icarusutil::TimeVec powerDerivVec;
    
    fWaveformTool.firstDerivative(powerVec, powerDerivVec);
    
    // Find the peaks...
    icarus_signal_processing::WaveformTools<float>::PeakTupleVec peakTupleVec;
    
    fWaveformTool.findPeaks(powerDerivVec.begin() + 300, powerDerivVec.end(), peakTupleVec, 10., 0);

    icarusutil::TimeVec smoothPowerVec = powerVec;
    
    // Try smoothing the peak regions
    for(const auto& peakTuple : peakTupleVec)
    {
        if (std::get<0>(peakTuple) >= powerVec.size() || std::get<2>(peakTuple) >= powerVec.size())
        {
            std::cout << "indexing problem - first: " << std::get<0>(peakTuple) << ", last: " << std::get<2>(peakTuple) << std::endl;
            continue;
        }
        
        size_t firstBin    = std::get<0>(peakTuple);
        size_t lastBin     = std::get<2>(peakTuple);
        double firstBinVal = powerVec.at(firstBin);
        double lastBinVal  = powerVec.at(lastBin);
        double stepVal     = (lastBinVal - firstBinVal) / double(lastBin - firstBin);
        double newBinVal   = firstBinVal + stepVal;
        
        while(++firstBin < lastBin)
        {
            // Update the power first
            smoothPowerVec[firstBin]  = newBinVal;
            newBinVal                += stepVal;
        }
    }
    
    while(currentBin > lowestBin)
    {
        double avePowerThisBin(smoothPowerVec.at(currentBin));
        double freq = 1.e6 * double(currentBin)/ (sampleRate * readOutSize);

        if (wire >= fLoWireByPlane[plane] && wire < fHiWireByPlane[plane])
        {
            fSmoothPowerVec[plane][wire-fLoWireByPlane[plane]]->Fill(freq, avePowerThisBin, 1.);
            fFFTPowerDerivVec[plane][wire-fLoWireByPlane[plane]]->Fill(freq, powerDerivVec.at(currentBin), 1.);
        }
        
        fAveSmoothPowerVec[plane]->Fill(freq, avePowerThisBin, 1.);
        fAveFFTPowerDerivVec[plane]->Fill(freq, powerDerivVec.at(currentBin), 1.);

        currentBin--;
    }
    
    // Recover the filter from signal shaping services...
    const icarusutil::FrequencyVec& filter = fSignalServices.GetResponse(channel).getFilter()->getResponseVec();
    
    // Convolve this with the FFT of the input waveform
    fFFT->convolute(inputVec,filter,0);
    fFFT->getFFTPower(inputVec,powerVec);

    for(size_t idx = 0; idx < halfFFTDataSize; idx++)
    {
        double freq = 1.e6 * double(idx)/ (sampleRate * readOutSize);
        fConvFFTPowerVec[plane]->Fill(freq, std::min(powerVec[idx],999.), 1.);
    }

    return;
}

// Useful for normalizing histograms
void BasicRawDigitAnalysis::endJob(int numEvents)
{
    // Nothing to do if nothing was initialized
    if (fAveFFTPowerVec.empty()) return;
    
    // A task to complete is to fit the average power displays with aim to develop a "good" filter function and
    // get the signal to noise ratio
    for(size_t planeIdx = 0; planeIdx < fChannelMapAlg.Nplanes(); planeIdx++)
    {
        TH1* avePowerHist = fAveFFTPowerVec[planeIdx];
        
        // Create the fitting function, use the histogram name to help
        std::string funcName = std::string(avePowerHist->GetName()) + "_func";
        
        // Create the function object
        TF1 fitFunc(funcName.c_str(),fFFTFitFuncVec.at(planeIdx).c_str(),avePowerHist->GetMinimum(),avePowerHist->GetMaximum());
        
        // Set initial parameters
        int paramIdx(0);
        
        for(const auto& param : fParameterVec.at(planeIdx)) fitFunc.SetParameter(paramIdx++, param);
        
        int fitResult(-1);
        
        try
        { fitResult = avePowerHist->Fit(&fitFunc,"QNRWB","", avePowerHist->GetMinimum(),avePowerHist->GetMaximum());}
        catch(...)
        {
            std::cout << "******* FFT power vec fit failure, skipping *******" << std::endl;
            continue;
        }
        
        if (!fitResult)
        {
            double chi2PerNDF = (fitFunc.GetChisquare() / fitFunc.GetNDF());
            int    NDF        = fitFunc.GetNDF();
            
            std::cout << "******************** Fit of " << avePowerHist->GetName() << " ********************" << std::endl;
            std::cout << "-- Function:   " << fFFTFitFuncVec.at(planeIdx) << ", chi2PerNDF: " << chi2PerNDF << ", NDF: " << NDF << std::endl;
            std::cout << "-- Parameters - 0: " << fitFunc.GetParameter(0);
            
            for(size_t idx = 1; idx < fParameterVec.at(planeIdx).size(); idx++) std::cout << ", " << idx << ": " << fitFunc.GetParameter(idx);
            
            std::cout << std::endl;
        }
    }
    
    return;
}
    
DEFINE_ART_CLASS_TOOL(BasicRawDigitAnalysis)
}
