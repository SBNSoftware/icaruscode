///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingICARUSService.h
///
/// \brief  Service to provide ICARUS-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution)./Users/Lsrea/newSim/SignalShapingICARUSService.h
///
/// \author H. Greenlee, major mods by L. Rochester
///
/// This service inherits from SignalShaping and supplies
/// ICARUS-specific configuration.  It is intended that SimWire and
/// CalWire modules will access this service.
///
/// FCL parameters:
///
/// FieldBins       - Number of bins of field response (generated from the histogram).
/// Col3DCorrection - 3D path length correction for collection plane. (not invoked)
/// Ind3DCorrection - 3D path length correction for induction plane.  (not invoked)
/// FieldRespAmpVec - vector of response amplitudes, one for each view
/// ShapeTimeConst  - Time constants for exponential shaping.
/// FilterVec       - vector of filter function parameters, one for each view
/// FilterParamsVec - Vector of filter function parameters.
///
/// \update notes: Leon Rochester (lsrea@slac.stanford.edu, Jan 12, 2015
///                many changes, need to be documented better
///                 1. the three (or n) views are now represented by a vector of views
///                 2. There are separate SignalShaping objects for convolution and
///                    deconvolution
///
///                Yun-Tse Tsai (yuntse@slac.stanford.edu), July 17th, 2014
///                 1. Read in field responses from input histograms
///                 2. Allow different sampling rates in the input
///                    field response
///                    NOTE: The InputFieldRespSamplingRate parameter has
///                    NOT implemented for the field response input
///                    as a function (UseFunctionFieldShape)
///                 3. Allow different electron drift velocities from
///                    which the input field responses are obtained
///                 4. Convolute the field and electronic responses,
///                    and then sample the convoluted function with
///                    the nominal sampling rate (detinfo::DetectorPropertiesService).
///                    NOTE: Currently this doesn't include the filter
///                    function and the deconvolution kernel.
///                    We may want to include them later?
///                 5. Disable fColSignalShaping.SetPeakResponseTime(0.),
///                    so that the peak time in the input field response
///                    is preserved.
///                 6. Somebody needs to unify the units of time (microsec
///                    or nanosec); I'm fainting!
///
/// New function:   void SetResponseSampling();
///
/// Modified functions: void init();
///                     void SetFieldResponse();
///
/// New FCL parameters:
/// DefaultDriftVelocity       - The electron drift velocity used to obtain
///                              the input field response waveforms
/// InputFieldRespSamplingRate - The sampling rate in the input field response
/// UseHistogramFieldShape     - Use the field response from an input histogram,
///                              if both UseFunctionFieldShape and
///                              UseHistogramFieldShape are false, we will
///                              use the toy field responses (a bipolar square
///                              function for induction planes, a ramp function
///                              for collection planes.)
/// FieldResponseFname         - Name of the file containing the input field
///                              response histograms
/// FieldResponseHistoName     - Name of the field response histograms,
///                              the format in the code will be
///                              FieldResponseHistoName_U(V,Y)
///update notes: Jyoti Joshi (jjoshi@bnl.gov), Jan 13, 2015
//               1. Modification to GetShapingTime function to read in different
//                  shaping time for different planes
//               2. Modification to GetASICGain fucntion to read in different gain
//                  settings for different planes
////////////////////////////////////////////////////////////////////////

#ifndef SignalShapingICARUSService_H
#define SignalShapingICARUSService_H

#include <vector>
#include <map>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "icaruscode/TPC/Utilities/tools/IResponse.h"
#include "TH1D.h"

using DoubleVec  = std::vector<double>;
using DoubleVec2 = std::vector<DoubleVec>;

namespace icarusutil {
    
class SignalShapingICARUSService //: private lar::EnsureOnlyOneSchedule
{
public:
    
    // Constructor, destructor.
    SignalShapingICARUSService(const fhicl::ParameterSet& pset,
                                   art::ActivityRegistry& reg);
    ~SignalShapingICARUSService();
    
    // Update configuration parameters.
    void                          reconfigure(const fhicl::ParameterSet& pset);
    
    // Accessors.
    DoubleVec2                    GetNoiseFactVec() { return fNoiseFactVec; }
    
    double                        GetASICGain(unsigned int const channel)            const;
    double                        GetShapingTime(unsigned int const channel)         const;
    
    double                        GetRawNoise(unsigned int const channel)            const;
    double                        GetDeconNoise(unsigned int const channel)          const;
    
    const icarus_tool::IResponse& GetResponse(size_t channel)                        const;
     
    int                           ResponseTOffset(unsigned int const channel)        const;
    
    void                          SetDecon(size_t fftsize, size_t channel);
    double                        GetDeconNorm() {return fDeconNorm;};
    
    
private:
    
    // Private configuration methods.
    using IResponsePtr             = std::unique_ptr<icarus_tool::IResponse>;
    using ResponseVec              = std::vector<IResponsePtr>;
    using PlaneToResponseMap       = std::map<size_t, ResponseVec>;
    
    // Post-constructor initialization.
    void init() const{const_cast<SignalShapingICARUSService*>(this)->init();}
    void init();
    
    // Attributes.
    bool fInit;                                                 ///< Initialization flag
    
    // Fcl parameters.
    size_t             fPlaneForNormalization; ///< Normalize responses to this plane
    double             fDeconNorm;             ///< Final normalization to apply
    DoubleVec2         fNoiseFactVec;          ///< RMS noise in ADCs for lowest gain setting
    size_t             fInitialFFTSize;        ///< Size we initially initalize the responses
    bool               fPrintResponses;
    bool               fStoreHistograms;
    
    // Field response tools
    PlaneToResponseMap fPlaneToResponseMap;
};

} // end of namespace

DECLARE_ART_SERVICE(icarusutil::SignalShapingICARUSService, SHARED)
#endif


