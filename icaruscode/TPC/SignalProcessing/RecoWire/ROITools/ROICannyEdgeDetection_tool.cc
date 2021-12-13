////////////////////////////////////////////////////////////////////////
/// \file   ROICannyEdgeDetection.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "icaruscode/TPC/SignalProcessing/RecoWire/ROITools/IROILocator.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "cetlib/cpu_timer.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Denoising.h"
#include "icarus_signal_processing/ROIFinder2D.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <fstream>
#include <chrono>

namespace icarus_tool
{

class ROICannyEdgeDetection : public IROILocator
{
public:
    explicit ROICannyEdgeDetection(const fhicl::ParameterSet& pset);
    
    ~ROICannyEdgeDetection();
    
    void configure(const fhicl::ParameterSet& pset) override;
    void initializeHistograms(art::TFileDirectory&) override {return;}
    
    void FindROIs(const art::Event&, const ArrayFloat&, const geo::PlaneID&, ArrayFloat&, ArrayBool&) override;
    
private:

    using FloatPairVec = std::vector<std::pair<float,float>>;

    bool                                           fDiagnosticOutput;           //< Enable diagnostic output if desired

    // Start with parameters for Butterworth Filter
    unsigned int                                   fButterworthOrder;           ///< Order parameter for Butterworth filter
    unsigned int                                   fButterworthThreshold;       ///< Threshold for Butterworth filter

//    // Parameters for the 2D morphological filter
//    unsigned int                                   fMorph2DStructuringElementX; ///< Structuring element in X
//    unsigned int                                   fMorph2DStructuringElementY; ///< Structuring element in Y
//
//    // Parameters for the denoiser
//    unsigned int                                   fCoherentNoiseGrouping;      ///< Number of consecutive channels in coherent noise subtraction
//    unsigned int                                   fCoherentNoiseOffset;        ///< Offset for the midplane...
//    unsigned int                                   fMorphologicalWindow;        ///< Window size for filter
////    bool                                           fOutputStats;                ///< Output of timiing statistics?
//    float                                          fCoherentThresholdFactor;    ///< Threshold factor for coherent noise removal

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

    icarus_signal_processing::VectorFloat          fThresholdVec;
    
    const geo::Geometry*                           fGeometry;              //< pointer to the Geometry service

    // Our denoising functions
    std::unique_ptr<icarus_signal_processing::IFFTFilterFunction>        fButterworthFilter;
//    std::unique_ptr<icarus_signal_processing::IMorphologicalFunctions2D> fMorphologicalFilter;
//    std::unique_ptr<icarus_signal_processing::IDenoiser2D>               fDenoiser2D;
    std::unique_ptr<icarus_signal_processing::BilateralFilters>          fBilateralFilters;
    std::unique_ptr<icarus_signal_processing::EdgeDetection>             fEdgeDetection;
    std::unique_ptr<icarus_signal_processing::IROIFinder2D>              fROIFinder2D;
};
    
//----------------------------------------------------------------------
// Constructor.
ROICannyEdgeDetection::ROICannyEdgeDetection(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
ROICannyEdgeDetection::~ROICannyEdgeDetection()
{
}
    
void ROICannyEdgeDetection::configure(const fhicl::ParameterSet& pset)
{
    fDiagnosticOutput           = pset.get<bool                    >("DiagnosticOutput",    false);
     
    fButterworthOrder           = pset.get<unsigned int            >("ButterworthOrder",     2);
    fButterworthThreshold       = pset.get<unsigned int            >("ButterworthThreshld", 30);

    fButterworthFilter = std::make_unique<icarus_signal_processing::HighPassButterworthFilter>(fButterworthThreshold,fButterworthOrder,4096);
    //fButterworthFilter = std::make_unique<icarus_signal_processing::NoFFTFilter>();

//    fMorph2DStructuringElementX = pset.get<unsigned int            >("Morph2DStructuringElementX", 7);
//    fMorph2DStructuringElementY = pset.get<unsigned int            >("Morph2DStructuringElementX", 28);
//
//    fMorphologicalFilter = std::make_unique<icarus_signal_processing::Dilation2D>(fMorph2DStructuringElementX,fMorph2DStructuringElementY);
//
//    fCoherentNoiseGrouping      = pset.get<unsigned int            >("CoherentNoiseGrouping",    32);
//    fCoherentNoiseOffset        = pset.get<unsigned int            >("CoherentNoiseOffset",      24);
//    fMorphologicalWindow        = pset.get<unsigned int            >("MorphologicalWindow",      10);
//    fCoherentThresholdFactor    = pset.get<float                   >("CoherentThresholdFactor", 2.5);
//
//    fThresholdVec.resize(6560/fCoherentNoiseGrouping,fCoherentThresholdFactor);
//
//    fDenoiser2D = std::make_unique<icarus_signal_processing::Denoiser2D_Hough>(fMorphologicalFilter.get(), fThresholdVec, fCoherentNoiseGrouping, fCoherentNoiseOffset, fMorphologicalWindow);
//    //fDenoiser2D = std::make_unique<icarus_signal_processing::Denoiser2D>(fMorphologicalFilter.get(), fThresholdVec, fCoherentNoiseGrouping, fMorphologicalWindow);

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
                                                                              nullptr, //fDenoiser2D.get(), 
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
    
    return;
}

void ROICannyEdgeDetection::FindROIs(const art::Event& event, const ArrayFloat& inputImage, const geo::PlaneID& planeID, ArrayFloat& output, ArrayBool& outputROIs)
{
    cet::cpu_timer theClockTotal;

    theClockTotal.start();

    std::cout << "  --> calling icarus_signal_processing canny edge finder" << std::endl;

    // Now pass the entire data array to the denoisercoherent
    (*fROIFinder2D)(inputImage,output,outputROIs); //,fWaveLessCoherent,fCorrectedMedians,fIntrinsicRMS,fMorphedWaveforms,finalErosion);

    std::cout << "  --> have returned from canny" << std::endl;

    theClockTotal.stop();

    double totalTime = theClockTotal.accumulated_real_time();

    mf::LogInfo("TPCNoiseFilterCannyMC") << "    *totalTime: " << totalTime << std::endl;

    std::cout << "--> ROICannyEdgeDetection finished!" << std::endl;
    std::cout << "    - Total time: " << totalTime << std::endl;
     
    return;
}

DEFINE_ART_CLASS_TOOL(ROICannyEdgeDetection)
}
