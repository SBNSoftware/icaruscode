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
    
    void FindROIs(const ArrayFloat&, const geo::PlaneID&, ArrayFloat&, ArrayBool&)    const override;
    
private:

    std::unique_ptr<icarus_signal_processing::BilateralFilters> fBilateralFilters;
    std::unique_ptr<icarus_signal_processing::EdgeDetection>    fEdgeDetection;
    std::unique_ptr<icarus_signal_processing::IROIFinder2D>     fROIFinder2D;

    // Parameters for the ROI finding
    unsigned int         fADFilter_SX;                ///< 
    unsigned int         fADFilter_SY;                ///< 
    float                fSigma_x;                    ///<
    float                fSigma_y;                    ///<
    float                fSigma_r;                    ///<
    float                fLowThreshold;               ///<
    float                fHighThreshold;              ///<
    unsigned int         fBinaryDilation_SX;          ///<
    unsigned int         fBinaryDilation_SY;          ///<
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
    fADFilter_SX           = pset.get<unsigned int>("ADFilter_SX",         7);
    fADFilter_SY           = pset.get<unsigned int>("ADFilter_SY",         7);
    fSigma_x               = pset.get<float       >("Sigma_x",          10.0);
    fSigma_y               = pset.get<float       >("Sigma_y",          10.0);
    fSigma_r               = pset.get<float       >("Sigma_r",          30.0);

    fLowThreshold          = pset.get<float       >("LowThreshold",     10.0);
    fHighThreshold         = pset.get<float       >("HighThreshold",    20.0); 

    fBinaryDilation_SX     = pset.get<unsigned int>("BinaryDilation_SX",  31);
    fBinaryDilation_SY     = pset.get<unsigned int>("BinaryDilation_SY",  31);

    fBilateralFilters = std::make_unique<icarus_signal_processing::BilateralFilters>();
    fEdgeDetection    = std::make_unique<icarus_signal_processing::EdgeDetection>();
    
    return;
}

void ROICannyEdgeDetection::FindROIs(const ArrayFloat& inputImage, const geo::PlaneID& planeID, ArrayFloat& output, ArrayBool& outputROIs) const
{
    unsigned int numChannels = inputImage.size();
    unsigned int numTicks    = inputImage[0].size();
  
    std::chrono::high_resolution_clock::time_point funcStartTime = std::chrono::high_resolution_clock::now();

    // 5. Directional Smoothing
    std::cout << "++> Step 5: Directional smoothing" << std::endl;
    icarus_signal_processing::ArrayFloat buffer   (numChannels, icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat sobelX   (numChannels, icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat sobelY   (numChannels, icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat gradient (numChannels, icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayFloat direction(numChannels, icarus_signal_processing::VectorFloat(numTicks,0.));
    icarus_signal_processing::ArrayBool  rois     (numChannels, icarus_signal_processing::VectorBool( numTicks,false));

    std::chrono::high_resolution_clock::time_point sobelStartTime = std::chrono::high_resolution_clock::now();

    fEdgeDetection->SepSobel(inputImage, sobelX, sobelY, gradient, direction);

    std::chrono::high_resolution_clock::time_point sobelStopTime      = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point bilateralStartTime = sobelStopTime;

    std::cout << "==> Step 6: Apply bilateral filter" << std::endl;

    fBilateralFilters->directional(inputImage, direction, buffer, fADFilter_SX, fADFilter_SY, fSigma_x, fSigma_y, fSigma_r, 360);

    std::chrono::high_resolution_clock::time_point bilateralStopTime = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point dilationStartTime = bilateralStopTime;

    std::cout << "==> Step 7: Apply Second Morphological Enhancing" << std::endl;

    icarus_signal_processing::Dilation2D(fADFilter_SX,fADFilter_SY)(buffer.begin(), numChannels, output.begin());

    std::chrono::high_resolution_clock::time_point dilationStopTime    = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point edgeInterpStartTime = dilationStopTime;

    std::cout << "==> Step 8: Perform Canny Edge Detection" << std::endl;

    // 6. Apply Canny Edge Detection
    for(auto& gradVec : gradient) std::fill(gradVec.begin(),gradVec.end(),0.);

    fEdgeDetection->SepSobel(output, sobelX, sobelY, gradient, direction);
    fEdgeDetection->EdgeNMSInterpolation(gradient, sobelX, sobelY, direction, output);

    std::chrono::high_resolution_clock::time_point edgeInterpStopTime  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point hysterisisStartTime = edgeInterpStopTime;

    std::vector<int> strongEdgeRows;
    std::vector<int> strongEdgeCols;
    std::vector<int> weakEdgeRows;
    std::vector<int> weakEdgeCols;

    fEdgeDetection->SparseHysteresisThresholding(output, fLowThreshold, fHighThreshold, rois);

    std::cout << "==> Final Step: get dilation, numChannels: " << numChannels << ", rois: " << rois.size() << ", output: " << outputROIs.size() << std::endl;

    std::chrono::high_resolution_clock::time_point hysterisisStopTime      = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point binaryDilationStartTime = hysterisisStopTime;

    icarus_signal_processing::Dilation2D(fBinaryDilation_SX,fBinaryDilation_SY)(rois.begin(), numChannels, outputROIs.begin());

    std::chrono::high_resolution_clock::time_point binaryDilationStopTime  = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point funcStopTime            = hysterisisStopTime;
  
    std::chrono::duration<double> funcTime       = std::chrono::duration_cast<std::chrono::duration<double>>(funcStopTime - funcStartTime);
    std::chrono::duration<double> sobelTime      = std::chrono::duration_cast<std::chrono::duration<double>>(sobelStopTime - sobelStartTime);
    std::chrono::duration<double> bilateralTime  = std::chrono::duration_cast<std::chrono::duration<double>>(bilateralStopTime - bilateralStartTime);
    std::chrono::duration<double> dilationTime   = std::chrono::duration_cast<std::chrono::duration<double>>(dilationStopTime - dilationStartTime);
    std::chrono::duration<double> edgeInterpTime = std::chrono::duration_cast<std::chrono::duration<double>>(edgeInterpStopTime - edgeInterpStartTime);
    std::chrono::duration<double> hysterisisTime = std::chrono::duration_cast<std::chrono::duration<double>>(hysterisisStopTime - hysterisisStartTime);
    std::chrono::duration<double> binaryDilTime  = std::chrono::duration_cast<std::chrono::duration<double>>(binaryDilationStopTime - binaryDilationStartTime);

    std::cout << "--> ROICannyEdgeDetection finished!" << std::endl;
    std::cout << "    - Total time: " << funcTime.count() << ", sobel: " << sobelTime.count() << ", bilateral: " << bilateralTime.count() << ", dilation: " << dilationTime.count() << ", Edge: " << edgeInterpTime.count() << ", hysterisis: " << hysterisisTime.count() << ", dilate: " << binaryDilTime.count() << std::endl;
     
    return;
}

DEFINE_ART_CLASS_TOOL(ROICannyEdgeDetection)
}
