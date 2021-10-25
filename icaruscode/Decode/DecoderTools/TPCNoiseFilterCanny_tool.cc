/**
 *  @file   TPCNoiseFilterCannyMC_tool.cc
 *
 *  @brief  This tool converts from daq to LArSoft format with noise filtering
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"

#include "sbndaq-artdaq-core/Overlays/ICARUS/PhysCrateFragment.hh"

#include "icaruscode/Decode/DecoderTools/INoiseFilter.h"

#include "icarus_signal_processing/ICARUSSigProcDefs.h"
#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Filters/FFTFilterFunctions.h"
#include "icarus_signal_processing/Filters/ImageFilters.h"
#include "icarus_signal_processing/Denoising.h"
#include "icarus_signal_processing/Detection/EdgeDetection.h"
#include "icarus_signal_processing/Filters/BilateralFilters.h"
#include "icarus_signal_processing/ROIFinder2D.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace daq {
/**
 *  @brief  TPCNoiseFilterCannyMC class definiton
 */
class TPCNoiseFilterCannyMC : virtual public INoiseFilter
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit TPCNoiseFilterCannyMC(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~TPCNoiseFilterCannyMC();

    /**
     *  @brief Interface for configuring the particular algorithm tool
     *
     *  @param ParameterSet  The input set of parameters for configuration
     */
    virtual void configure(const fhicl::ParameterSet&) override;

    /**
     *  @brief Given a set of recob hits, run DBscan to form 3D clusters
     *
     *  @param fragment            The artdaq fragment to process
     */
    virtual void process_fragment(detinfo::DetectorClocksData const&,
                                  const daq::INoiseFilter::ChannelVec&,
                                  const icarus_signal_processing::ArrayFloat&) override;

    /**
     *  @brief Recover the channels for the processed fragment
     */
    const icarus_signal_processing::VectorInt&  getChannelIDs()       const override {return fChannelIDVec;}

    /**
     *  @brief Recover the selection values
     */
    const icarus_signal_processing::ArrayBool&  getSelectionVals()    const override {return fSelectVals;};

    /**
     *  @brief Recover the ROI values
     */
    const icarus_signal_processing::ArrayBool&  getROIVals()          const override {return fROIVals;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarus_signal_processing::ArrayFloat& getRawWaveforms()     const override {return fRawWaveforms;};

    /**
     *  @brief Recover the pedestal subtracted waveforms
     */
    const icarus_signal_processing::ArrayFloat& getPedCorWaveforms()  const override {return fPedCorWaveforms;};

    /**
     *  @brief Recover the "intrinsic" RMS
     */
    const icarus_signal_processing::ArrayFloat& getIntrinsicRMS()     const override {return fIntrinsicRMS;};

    /**
     *  @brief Recover the correction median values
     */
    const icarus_signal_processing::ArrayFloat& getCorrectedMedians()  const override {return fCorrectedMedians;};

    /**
     *  @brief Recover the waveforms less coherent noise
     */
    const icarus_signal_processing::ArrayFloat& getWaveLessCoherent()  const override {return fWaveLessCoherent;};

    /**
     *  @brief Recover the morphological filter waveforms
     */
    const icarus_signal_processing::ArrayFloat& getMorphedWaveforms()  const override {return fMorphedWaveforms;};

    /**
     *  @brief Recover the pedestals for each channel
     */
    const icarus_signal_processing::VectorFloat& getPedestalVals()     const override {return fPedestalVals;};

    /**
     *  @brief Recover the full RMS before coherent noise
     */
    const icarus_signal_processing::VectorFloat& getFullRMSVals()      const override {return fFullRMSVals;};

    /**
     *  @brief Recover the truncated RMS noise
     */
    const icarus_signal_processing::VectorFloat& getTruncRMSVals()     const override {return fTruncRMSVals;};

    /**
     *  @brief Recover the number of bins after truncation
     */
    const icarus_signal_processing::VectorInt&   getNumTruncBins()     const override {return fNumTruncBins;};

private:

    using FloatPairVec = std::vector<std::pair<float,float>>;

    float                                          fSigmaForTruncation;         //< Selection cut for truncated rms calculation
    bool                                           fUseFFTFilter;               //< Turn on/off the use of the FFT filter
    bool                                           fDiagnosticOutput;           //< If true will spew endless messages to output
    FloatPairVec                                   fFFTSigmaValsVec;            //< Gives the sigmas for the filter function
    FloatPairVec                                   fFFTCutoffValsVec;           //< Gives the cutoffs for the filter function
      
    std::vector<std::string>                       fFilterModeVec;              //< Allowed modes for the filter

    // fhicl parameters
    std::vector<size_t>                            fStructuringElement;         ///< Structuring element for morphological filter
    std::vector<float>                             fThreshold;                  ///< Threshold to apply for saving signal

    // Start with parameters for Butterworth Filter
    unsigned int                                   fButterworthOrder;           ///< Order parameter for Butterworth filter
    unsigned int                                   fButterworthThreshold;       ///< Threshold for Butterworth filter

    // Parameters for the 2D morphological filter
    unsigned int                                   fMorph2DStructuringElementX; ///< Structuring element in X
    unsigned int                                   fMorph2DStructuringElementY; ///< Structuring element in Y

    // Parameters for the denoiser
    unsigned int                                   fCoherentNoiseGrouping;      ///< Number of consecutive channels in coherent noise subtraction
    unsigned int                                   fCoherentNoiseOffset;        ///< Offset for the midplane...
    unsigned int                                   fMorphologicalWindow;        ///< Window size for filter
//    bool                                           fOutputStats;                ///< Output of timiing statistics?
    float                                          fCoherentThresholdFactor;    ///< Threshold factor for coherent noise removal

    // Parameters for the ROI finding
    unsigned int                                   fADFilter_SX;                ///< 
    unsigned int                                   fADFilter_SY;                ///< 
    float                                          fSigma_x;                    ///<
    float                                          fSigma_y;                    ///<
    float                                          fSigma_r;                    ///<
    float                                          fLowThreshold;               ///<
    float                                          fHighThreshold;              ///<
    unsigned int                                   fBinaryDilation_SX;          ///<
    unsigned int                                   fBinaryDilation_SY;          ///<

    // Allocate containers for noise processing
    icarus_signal_processing::VectorInt            fChannelIDVec;
    icarus_signal_processing::ArrayBool            fSelectVals;
    icarus_signal_processing::ArrayBool            fROIVals;
    icarus_signal_processing::ArrayFloat           fRawWaveforms;
    icarus_signal_processing::ArrayFloat           fPedCorWaveforms;
    icarus_signal_processing::ArrayFloat           fIntrinsicRMS;
    icarus_signal_processing::ArrayFloat           fCorrectedMedians;
    icarus_signal_processing::ArrayFloat           fWaveLessCoherent;
    icarus_signal_processing::ArrayFloat           fMorphedWaveforms;
      
    icarus_signal_processing::VectorFloat          fPedestalVals;
    icarus_signal_processing::VectorFloat          fFullRMSVals;
    icarus_signal_processing::VectorFloat          fTruncRMSVals;
    icarus_signal_processing::VectorInt            fNumTruncBins;
    icarus_signal_processing::VectorInt            fRangeBins;

    icarus_signal_processing::VectorFloat          fThresholdVec;

    icarus_signal_processing::FilterFunctionVec    fFilterFunctionVec;
    
    const geo::Geometry*                           fGeometry;              //< pointer to the Geometry service

    // Keep track of the FFT 
    icarus_signal_processing::FFTFilterFunctionVec fFFTFilterFunctionVec;

    // Our denoising functions
    std::unique_ptr<icarus_signal_processing::IFFTFilterFunction>        fButterworthFilter;
    std::unique_ptr<icarus_signal_processing::IMorphologicalFunctions2D> fMorphologicalFilter;
    std::unique_ptr<icarus_signal_processing::IDenoiser2D>               fDenoiser2D;
    std::unique_ptr<icarus_signal_processing::BilateralFilters>          fBilateralFilters;
    std::unique_ptr<icarus_signal_processing::EdgeDetection>             fEdgeDetection;
    std::unique_ptr<icarus_signal_processing::IROIFinder2D>              fROIFinder2D;

};

TPCNoiseFilterCannyMC::TPCNoiseFilterCannyMC(fhicl::ParameterSet const &pset)
{
    this->configure(pset);

    fSelectVals.clear();
    fROIVals.clear();
    fRawWaveforms.clear();
    fPedCorWaveforms.clear();
    fIntrinsicRMS.clear();
    fCorrectedMedians.clear();
    fWaveLessCoherent.clear();
    fMorphedWaveforms.clear();

    fPedestalVals.clear();
    fFullRMSVals.clear();
    fTruncRMSVals.clear();
    fNumTruncBins.clear();
    fRangeBins.clear();

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TPCNoiseFilterCannyMC::~TPCNoiseFilterCannyMC()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void TPCNoiseFilterCannyMC::configure(fhicl::ParameterSet const &pset)
{
    fSigmaForTruncation         = pset.get<float                   >("NSigmaForTrucation",  3.5);
    fCoherentNoiseGrouping      = pset.get<size_t                  >("CoherentGrouping",     64);
    fUseFFTFilter               = pset.get<bool                    >("UseFFTFilter",        true);
    fDiagnosticOutput           = pset.get<bool                    >("DiagnosticOutput",    false);
    fFilterModeVec              = pset.get<std::vector<std::string>>("FilterModeVec",       std::vector<std::string>()={"e","g","d"});
     
    fFFTSigmaValsVec            = pset.get<FloatPairVec            >("FFTSigmaVals",        FloatPairVec()={{1.5,20.}, {1.5,20.}, {2.0,20.}});
    fFFTCutoffValsVec           = pset.get<FloatPairVec            >("FFTCutoffVals",       FloatPairVec()={{8.,800.}, {8.,800.}, {0.0,800.}});
     
    // Recover parameters for noise/ROI
    fStructuringElement         = pset.get<std::vector<size_t>     >("StructuringElement",                       std::vector<size_t>()={8,16});
    fThreshold                  = pset.get<std::vector<float>      >("Threshold",                       std::vector<float>()={2.75,2.75,2.75});
     
    fButterworthOrder           = pset.get<unsigned int            >("ButterworthOrder",     2);
    fButterworthThreshold       = pset.get<unsigned int            >("ButterworthThreshld", 30);

    //fButterworthFilter = std::make_unique<icarus_signal_processing::HighPassButterworthFilter>(fButterworthThreshold,fButterworthOrder,4096);
    fButterworthFilter = std::make_unique<icarus_signal_processing::NoFFTFilter>();

    fMorph2DStructuringElementX = pset.get<unsigned int            >("Morph2DStructuringElementX", 7);
    fMorph2DStructuringElementY = pset.get<unsigned int            >("Morph2DStructuringElementX", 28);

    fMorphologicalFilter = std::make_unique<icarus_signal_processing::Dilation2D>(fMorph2DStructuringElementX,fMorph2DStructuringElementY);

    fCoherentNoiseGrouping      = pset.get<unsigned int            >("CoherentNoiseGrouping",    32);
    fCoherentNoiseOffset        = pset.get<unsigned int            >("CoherentNoiseOffset",      24);
    fMorphologicalWindow        = pset.get<unsigned int            >("MorphologicalWindow",      10);
    fCoherentThresholdFactor    = pset.get<float                   >("CoherentThresholdFactor", 2.5);

    fThresholdVec.resize(6560/fCoherentNoiseGrouping,fCoherentThresholdFactor);

    //fDenoiser2D = std::make_unique<icarus_signal_processing::Denoiser2D_Hough>(fMorphologicalFilter.get(), fThresholdVec, fCoherentNoiseGrouping, fCoherentNoiseOffset, fMorphologicalWindow);
    fDenoiser2D = std::make_unique<icarus_signal_processing::Denoiser2D>(fMorphologicalFilter.get(), fThresholdVec, fCoherentNoiseGrouping, fMorphologicalWindow);

    fADFilter_SX                = pset.get<unsigned int            >("ADFilter_SX",         7);
    fADFilter_SY                = pset.get<unsigned int            >("ADFilter_SY",         7);
    fSigma_x                    = pset.get<float                   >("Sigma_x",          10.0);
    fSigma_y                    = pset.get<float                   >("Sigma_y",          10.0);
    fSigma_r                    = pset.get<float                   >("Sigma_r",          30.0);
                   
    fLowThreshold               = pset.get<float                   >("LowThreshold",     10.0);
    fHighThreshold              = pset.get<float                   >("HighThreshold",    20.0); 
                   
    fBinaryDilation_SX          = pset.get<unsigned int            >("BinaryDilation_SX",  31);
    fBinaryDilation_SY          = pset.get<unsigned int            >("BinaryDilation_SY",  31);

    fBilateralFilters = std::make_unique<icarus_signal_processing::BilateralFilters>();
    fEdgeDetection    = std::make_unique<icarus_signal_processing::EdgeDetection>();

    fROIFinder2D = std::make_unique<icarus_signal_processing::ROICannyFilter>(fButterworthFilter.get(), 
                                                                              fDenoiser2D.get(), 
                                                                              fBilateralFilters.get(),
                                                                              fEdgeDetection.get(), 
                                                                              fADFilter_SX,
                                                                              fADFilter_SY,
                                                                              fSigma_x,
                                                                              fSigma_y,
                                                                              fSigma_r,
                                                                              fLowThreshold,
                                                                              fHighThreshold,
                                                                              fBinaryDilation_SX,
                                                                              fBinaryDilation_SY);

    fGeometry   = art::ServiceHandle<geo::Geometry const>{}.get();

    fFFTFilterFunctionVec.clear();

    if (fUseFFTFilter)
    {
        for(int plane = 0; plane < 3; plane++)
        {
            if (plane < 2) fFFTFilterFunctionVec.emplace_back(std::make_unique<icarus_signal_processing::WindowFFTFilter>(fFFTSigmaValsVec[plane], fFFTCutoffValsVec[plane]));
            else           fFFTFilterFunctionVec.emplace_back(std::make_unique<icarus_signal_processing::NoFFTFilter>());
        }
    }

    return;
}

void TPCNoiseFilterCannyMC::process_fragment(detinfo::DetectorClocksData const&,
                                               const daq::INoiseFilter::ChannelVec&  channelVec,
                                               const icarus_signal_processing::ArrayFloat& dataArray)
{
    cet::cpu_timer theClockTotal;

    theClockTotal.start();

    // Recover the number of channels and ticks
    unsigned int numChannels = dataArray.size();
    unsigned int numTicks    = dataArray[0].size();

    if (fSelectVals.size()       < numChannels)  fSelectVals.resize(numChannels, icarus_signal_processing::VectorBool(numTicks));
    if (fROIVals.size()          < numChannels)  fROIVals.resize(numChannels,  icarus_signal_processing::VectorBool(numTicks));
    if (fRawWaveforms.size()     < numChannels)  fRawWaveforms.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));
    if (fPedCorWaveforms.size()  < numChannels)  fPedCorWaveforms.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));
    if (fIntrinsicRMS.size()     < numChannels)  fIntrinsicRMS.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));
    if (fCorrectedMedians.size() < numChannels)  fCorrectedMedians.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));
    if (fWaveLessCoherent.size() < numChannels)  fWaveLessCoherent.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));
    if (fMorphedWaveforms.size() < numChannels)  fMorphedWaveforms.resize(numChannels, icarus_signal_processing::VectorFloat(numTicks));

    if (fChannelIDVec.size()     < numChannels)  fChannelIDVec.resize(numChannels);
    if (fPedestalVals.size()     < numChannels)  fPedestalVals.resize(numChannels);
    if (fFullRMSVals.size()      < numChannels)  fFullRMSVals.resize(numChannels);
    if (fTruncRMSVals.size()     < numChannels)  fTruncRMSVals.resize(numChannels);
    if (fNumTruncBins.size()     < numChannels)  fNumTruncBins.resize(numChannels);
    if (fRangeBins.size()        < numChannels)  fRangeBins.resize(numChannels);

    if (fThresholdVec.size()     < numChannels)  fThresholdVec.resize(numChannels / fCoherentNoiseGrouping);

    if (fFilterFunctionVec.size() < numChannels) fFilterFunctionVec.resize(numChannels);

    std::cout <<"  -->process_fragment with " << numChannels << " channels and " << numTicks << " ticks, array sizes: " << fCorrectedMedians.size() << ", " << fCorrectedMedians[1].size() <<  std::endl;

    icarus_signal_processing::Denoiser1D           denoiser;
    icarus_signal_processing::WaveformTools<float> waveformTools;

    // Make a pass throught to do pedestal corrections and get raw waveform information
    for(size_t idx = 0; idx < numChannels; idx++)
    {
        icarus_signal_processing::VectorFloat& pedCorDataVec = fPedCorWaveforms[idx];

        // Keep track of the channel
        fChannelIDVec[idx] = channelVec[idx];

        // We need to recover info on which plane we have
        std::vector<geo::WireID> widVec = fGeometry->ChannelToWire(fChannelIDVec[idx]);

        // Handle the filter function to use for this channel
        unsigned int plane = widVec[0].Plane;

        // Set the threshold which toggles between planes
        fThresholdVec[idx / fCoherentNoiseGrouping] = fThreshold[plane];

        switch(fFilterModeVec[plane][0])
        {
            case 'd' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Dilation1D>(fStructuringElement[1]);
                break;
            case 'e' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Erosion1D>(fStructuringElement[1]);
                break;
            case 'g' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Gradient1D>(fStructuringElement[1]);
                break;
            case 'a' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Average1D>(fStructuringElement[1]);
                break;
            case 'm' :
                fFilterFunctionVec[idx] = std::make_unique<icarus_signal_processing::Median1D>(fStructuringElement[1]);
                break;
            default:
                std::cout << "***** FOUND NO MATCH FOR TYPE: " << fFilterModeVec[plane] << ", plane " << plane << " DURING INITIALIZATION OF FILTER FUNCTIONS IN TPCNoiseFilterCannyMC" << std::endl;
                break;
        }

        // Now determine the pedestal and correct for it
        waveformTools.getPedestalCorrectedWaveform(dataArray[idx],
                                                   pedCorDataVec,
                                                   fSigmaForTruncation,
                                                   fPedestalVals[idx],
                                                   fFullRMSVals[idx],
                                                   fTruncRMSVals[idx],
                                                   fNumTruncBins[idx],
                                                   fRangeBins[idx]);

        // Convolve with a filter function
        if (fUseFFTFilter) (*fFFTFilterFunctionVec[plane])(pedCorDataVec);
    }

    icarus_signal_processing::ArrayFloat finalErosion(numChannels,icarus_signal_processing::VectorFloat(numTicks,0.));

    std::cout << "  --> calling icarus_signal_processing code" << std::endl;

    // Now pass the entire data array to the denoisercoherent
    (*fROIFinder2D)(fPedCorWaveforms,fRawWaveforms,fROIVals,fWaveLessCoherent,fCorrectedMedians,fIntrinsicRMS,fMorphedWaveforms,finalErosion);

    std::cout << "  --> have returned from denoising" << std::endl;

    theClockTotal.stop();

    double totalTime = theClockTotal.accumulated_real_time();

    mf::LogInfo("TPCNoiseFilterCannyMC") << "    *totalTime: " << totalTime << std::endl;

    return;
}


DEFINE_ART_CLASS_TOOL(TPCNoiseFilterCannyMC)
} // namespace lar_cluster3d
