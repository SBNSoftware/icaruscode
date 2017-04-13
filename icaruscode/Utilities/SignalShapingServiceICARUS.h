///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingServiceICARUS.h
///
/// \brief  Service to provide ICARUS-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution)./Users/Lsrea/newSim/SignalShapingServiceICARUS.h
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

#ifndef SIGNALSHAPINGSERVICEICARUS_H
#define SIGNALSHAPINGSERVICEICARUS_H

#include <vector>
#include <map>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "lardata/Utilities/SignalShaping.h"
#include "TF1.h"
#include "TH1D.h"

// LArSoft include
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"

using DoubleVec  = std::vector<double>;
using DoubleVec2 = std::vector<DoubleVec>;
using TH1FVec2   = std::vector<std::vector<TH1F*>>;

namespace icarus_tool
{
    class IFieldResponse;
}

namespace util {
    class SignalShapingServiceICARUS {
    public:
        
        // Constructor, destructor.
        
        SignalShapingServiceICARUS(const fhicl::ParameterSet& pset,
                                       art::ActivityRegistry& reg);
        ~SignalShapingServiceICARUS();
        
        // Update configuration parameters.
        
        void reconfigure(const fhicl::ParameterSet& pset);
        
        // Accessors.
        DoubleVec2 GetNoiseFactVec()              { return fNoiseFactVec; }
        
        double GetASICGain(unsigned int const channel) const;
        double GetShapingTime(unsigned int const channel) const;
        
        double GetRawNoise(unsigned int const channel) const ;
        double GetDeconNoise(unsigned int const channel) const;
        
        const util::SignalShaping& SignalShaping(size_t channel, size_t ktype = 0) const;
        
        int FieldResponseTOffset(unsigned int const channel, size_t ktype) const;
        
        // Do convolution calcution (for simulation).
        
        template <class T> void Convolute(size_t channel, std::vector<T>& func) const;
        
        // Do deconvolution calcution (for reconstruction).
        
        template <class T> void Deconvolute(size_t channel, std::vector<T>& func) const;
        
        void SetDecon(size_t fftsize, size_t channel);
        double GetDeconNorm(){return fDeconNorm;};
        
        
    private:
        
        // Private configuration methods.
        using IFieldResponsePtr       = std::unique_ptr<icarus_tool::IFieldResponse>;
        using FieldResponseVec        = std::vector<IFieldResponsePtr>;
        using PlaneToFieldResponseMap = std::map<size_t, FieldResponseVec>;
        
        // Post-constructor initialization.
        
        void init() const{const_cast<SignalShapingServiceICARUS*>(this)->init();}
        void init();
        
        // Calculate response functions.
        // Copied from SimWireICARUS.
        void SetTimeScaleFactor();
        
        void SetFieldResponse(size_t ktype);
        
        void SetElectResponse(size_t ktype, size_t plane, double shapingtime, double gain);  //changed to read different peaking time for different planes
        
        // Calculate filter functions.
        void SetFilters();
        
        // Attributes.
        bool fInit;               ///< Initialization flag
        
        // Sample the response function, including a configurable
        // drift velocity of electrons
        void SetResponseSampling(size_t ktype, int mode=0, size_t channel=0);
        
        // Fcl parameters.
        size_t                              fViewForNormalization;
        
        double                              fDeconNorm;
        double                              fADCPerPCAtLowestASICGain; ///< Pulse amplitude gain for a 1 pc charge impulse after convoluting it the with field and electronics response with the lowest ASIC gain setting of 4.7 mV/fC
        
        DoubleVec2                          fNoiseFactVec;       ///< RMS noise in ADCs for lowest gain setting
        DoubleVec                           fASICGainInMVPerFC;       ///< Cold electronics ASIC gain setting in mV/fC
        DoubleVec                           fDefaultDriftVelocity;  ///< Default drift velocity of electrons in cm/usec
        bool                                fUseCalibratedResponses;         //Flag to use Calibrated Responses for 70kV
        
        DoubleVec                           fCalibResponseTOffset; // calibrated time offset to align U/V/Y Signals
        
        // Field response tools
        PlaneToFieldResponseMap             fPlaneToFieldResponseVec;
        
        // test
        
        DoubleVec                           f3DCorrectionVec;  ///< correction factor to account for 3D path of electrons, 1 for each plane (default = 1.0)
        
        double                              fTimeScaleFactor;
        bool                                fStretchFullResponse;
        
        DoubleVec                           fFieldRespAmpVec;
        DoubleVec                           fShapeTimeConst; ///< time constants for exponential shaping
        std::vector<int>                    fDeconvPol;         ///< switch for DeconvKernel normalization sign (+ -> max pos ADC, - -> max neg ADC). Entry 0,1,2 = U,V,Y plane settings
        std::vector<TF1*>                   fFilterTF1Vec;     ///< Vector of Parameterized filter functions
        std::vector<std::string>            fFilterFuncVec;
        std::vector<std::vector<TComplex> > fFilterVec;
        DoubleVec2                          fFilterParamsVec;
        DoubleVec                           fFilterWidthCorrectionFactor;  // a knob

        bool                                fGetFilterFromHisto;   		///< Flag that allows to use a filter function from a histogram instead of the functional dependency
        
        double                              fDefaultEField;
        double                              fDefaultTemperature;
        
        DoubleVec                           fTimeScaleParams;
        
        std::vector<TH1D*>                  fFilterHistVec;
        
        // Following attributes hold the convolution and deconvolution kernels
        
        std::vector<std::vector<util::SignalShaping> > fSignalShapingVec;
        // Field response.
        
        std::vector<DoubleVec2 >            fFieldResponseVec;
        
        // Electronics response.
        
        std::vector<DoubleVec>              fElectResponse;
        
        bool                                fPrintResponses;
        bool                                fManualInterpolation;
        
        // some diagnostic histograms
        
        bool fHistDoneF[3];
    };
}




//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceICARUS::Convolute(size_t channel, std::vector<T>& func) const
{
    SignalShaping(channel, 0).Convolute(func);
    
    //negative number
    int time_offset = FieldResponseTOffset(channel,0);
    
    std::vector<T> temp;
    if (time_offset <=0){
        temp.assign(func.begin(),func.begin()-time_offset);
        func.erase(func.begin(),func.begin()-time_offset);
        func.insert(func.end(),temp.begin(),temp.end());
    }else{
        temp.assign(func.end()-time_offset,func.end());
        func.erase(func.end()-time_offset,func.end());
        func.insert(func.begin(),temp.begin(),temp.end());
    }
    
}

//----------------------------------------------------------------------
// Do deconvolution.
template <class T> inline void util::SignalShapingServiceICARUS::Deconvolute(size_t channel, std::vector<T>& func) const
{
    size_t ktype = 1;
    SignalShaping(channel, ktype).Deconvolute(func);
    
    int time_offset = FieldResponseTOffset(channel,ktype);
    
    std::vector<T> temp;
    if (time_offset <=0){
        temp.assign(func.end()+time_offset,func.end());
        func.erase(func.end()+time_offset,func.end());
        func.insert(func.begin(),temp.begin(),temp.end());
    }else{
        temp.assign(func.begin(),func.begin()+time_offset);
        func.erase(func.begin(),func.begin()+time_offset);
        func.insert(func.end(),temp.begin(),temp.end());
        
        
    }
}

DECLARE_ART_SERVICE(util::SignalShapingServiceICARUS, LEGACY)
#endif
