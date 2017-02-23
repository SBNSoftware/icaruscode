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
using DoubleVec2 = std::vector< DoubleVec >;
using DoubleVec3 = std::vector< DoubleVec2 >;
using DoubleVec4 = std::vector< DoubleVec3 >;
using TH1FVec4   = std::vector<std::vector<std::vector<std::vector<TH1F*> > > >;


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
        
        DoubleVec2 GetNoiseFactVec()                { return fNoiseFactVec; }
        
        std::vector<std::vector<size_t> > GetNResponses()       { return fNResponses; }
        std::vector<std::vector<size_t> > GetNYZResponses()     { return fNYZResponses; }
        std::vector<std::vector<size_t> > GetNdatadrivenResponses()     { return fNdatadrivenResponses; }
        std::vector<std::vector<size_t> > GetNActiveResponses() { return fNActiveResponses; }
        std::vector<std::vector<size_t> > GetNYZActiveResponses() { return fNYZActiveResponses; }
        std::vector<std::vector<size_t> > GetNdatadrivenActiveResponses() { return fNdatadrivenActiveResponses; }
        std::vector<std::vector<double> > GetYZchargeScaling()  { return fYZchargeScaling; }
        //std::vector<std::vector<std::vector<int> > > GetYZwireOverlap() { return fYZwireOverlap; }
        std::vector<std::vector<int> > GetMisconfiguredU()      { return fMisconfiguredU; }
        
        
        std::vector<size_t> GetViewIndex()       { return fViewIndex; }
        
        bool IsResponseYZDependent()    { return fYZdependentResponse; }
        bool IsdatadrivenResponse()     { return fdatadrivenResponse; }
        bool IsMisconfiguredUIncluded() { return fIncludeMisconfiguredU; }
        
        double GetASICGain(unsigned int const channel) const;
        double GetShapingTime(unsigned int const channel) const;
        
        double GetRawNoise(unsigned int const channel) const ;
        double GetDeconNoise(unsigned int const channel) const;
        
        const std::vector<TComplex>& GetConvKernel(unsigned int channel, unsigned int wire) const;  // M. Mooney
        double Get2DFilterVal(size_t planeNum, size_t freqDimension, double binFrac) const;  // M. Mooney
        double Get2DFilterNorm(size_t planeNum) const;  // M. Mooney
        
        const util::SignalShaping& SignalShaping(size_t channel, size_t wire = 0, size_t ktype = 0) const;
        
        int FieldResponseTOffset(unsigned int const channel, size_t ktype) const;
        
        // Do convolution calcution (for simulation).
        
        template <class T> void Convolute(size_t channel, std::vector<T>& func) const;
        template <class T> void Convolute(size_t channel, size_t wire, std::vector<T>& func) const;
        
        // Do deconvolution calcution (for reconstruction).
        
        template <class T> void Deconvolute(size_t channel, std::vector<T>& func) const;
        template <class T> void Deconvolute(size_t channel, size_t wire, std::vector<T>& func) const;
        
        void SetDecon(size_t fftsize, size_t channel);
        double GetDeconNorm(){return fDeconNorm;};
        
        
    private:
        
        // Private configuration methods.
        
        // Post-constructor initialization.
        
        void init() const{const_cast<SignalShapingServiceICARUS*>(this)->init();}
        void init();
        
        // Calculate response functions.
        // Copied from SimWireICARUS.
        
        void SetTimeScaleFactor();
        
        void SetFieldResponse(size_t ktype);
        // void SetElectResponse(size_t ktype);
        
        void SetElectResponse(size_t ktype, double shapingtime, double gain);  //changed to read different peaking time for different planes
        
        // Calculate filter functions.
        
        void SetFilters();
        
        // Pick the electronics configuration
        size_t GetConfig(size_t channel) const;
        
        // Attributes.
        
        bool fInit;               ///< Initialization flag
        bool fInitConfigMap;
        
        mutable std::map<size_t, size_t> fConfigMap;
        mutable size_t fConfigMapFirstChannel;
        mutable size_t fConfigMapLastChannel;
        
        // Sample the response function, including a configurable
        // drift velocity of electrons
        
        void SetResponseSampling(size_t ktype, size_t config, int mode=0, size_t channel=0);
        
        void SetFieldResponseTOffsets( const TH1F* resp, const size_t ktype);
        
        DoubleVec3 fTestParams;
        
        size_t fNConfigs;
        size_t fNPlanes;
        size_t fNViews;
        
        // Fcl parameters.
        std::vector<size_t>      fViewIndex;
        std::map<size_t, size_t> fViewMap;
        size_t                   fViewForNormalization;
        
        double fDeconNorm;
        double fADCPerPCAtLowestASICGain; ///< Pulse amplitude gain for a 1 pc charge impulse after convoluting it the with field and electronics response with the lowest ASIC gain setting of 4.7 mV/fC
        
        DoubleVec2 fNoiseFactVec;       ///< RMS noise in ADCs for lowest gain setting
        
        std::vector<std::vector<size_t> > fNResponses;
        std::vector<std::vector<size_t> > fNYZResponses;
        std::vector<std::vector<size_t> > fNdatadrivenResponses;
        std::vector<std::vector<size_t> > fNActiveResponses;
        std::vector<std::vector<size_t> > fNYZActiveResponses;
        std::vector<std::vector<size_t> > fNdatadrivenActiveResponses;
        
        std::vector<std::vector<double> > fYZchargeScaling;
        //std::vector<std::vector<std::vector<int> > > fYZwireOverlap;
        std::vector<std::vector<int> > fMisconfiguredU;
        
        bool fYZdependentResponse;
        bool fdatadrivenResponse;
        bool fIncludeMisconfiguredU;
        
        DoubleVec2 fASICGainInMVPerFC;       ///< Cold electronics ASIC gain setting in mV/fC
        
        DoubleVec fDefaultDriftVelocity;  ///< Default drift velocity of electrons in cm/usec
        DoubleVec2  fFieldResponseTOffset;  ///< Time offset for field response in ns
        
        bool fUseCalibratedResponses;         //Flag to use Calibrated Responses for 70kV
        
        DoubleVec fCalibResponseTOffset; // calibrated time offset to align U/V/Y Signals
        
        // test
        
        int fNFieldBins[2]; // BR
        //size_t fNFieldBins[2];         		///< number of bins for field response
        double fFieldLowEdge[2];           ///< low edge of the field response histo (for test output)
        double fFieldBin1Center[2];
        double fFieldBinWidth[2];       ///<  Bin with of the input field response.
        
        DoubleVec f3DCorrectionVec;  ///< correction factor to account for 3D path of electrons, 1 for each plane (default = 1.0)
        
        double fTimeScaleFactor;
        bool   fStretchFullResponse;
        
        DoubleVec fFieldRespAmpVec;
        DoubleVec2 fShapeTimeConst; ///< time constants for exponential shaping
        std::vector<int> fDeconvPol;         ///< switch for DeconvKernel normalization sign (+ -> max pos ADC, - -> max neg ADC). Entry 0,1,2 = U,V,Y plane settings
        std::vector<TF1*> fFilterTF1Vec;     ///< Vector of Parameterized filter functions
        std::vector<std::string> fFilterFuncVec;
        std::vector<std::vector<TComplex> > fFilterVec;
        DoubleVec2 fFilterParamsVec;
        DoubleVec fFilterWidthCorrectionFactor;  // a knob
        
        // Induced charge deconvolution additions (M. Mooney)
        std::vector<TF1*> fFilterTF1VecICTime;
        std::vector<std::string> fFilterFuncVecICTime;
        std::vector<TF1*> fFilterTF1VecICWire;
        std::vector<std::string> fFilterFuncVecICWire;
        DoubleVec fFilterScaleVecICTime;
        DoubleVec fFilterScaleVecICWire;
        DoubleVec fFilterNormVecIC;
        
        std::vector<double> fFilterICTimeMaxFreq;
        DoubleVec fFilterICTimeMaxVal;
        
        DoubleVec fFilterICWireMaxFreq;
        DoubleVec fFilterICWireMaxVal;
        
        
        
        bool fGetFilterFromHisto;   		///< Flag that allows to use a filter function from a histogram instead of the functional dependency
        
        TH1FVec4 fFieldResponseHistVec;
        
        double fDefaultEField;
        double fDefaultTemperature;
        
        DoubleVec fTimeScaleParams;
        
        std::vector<TH1D*> fFilterHistVec;
        
        // Following attributes hold the convolution and deconvolution kernels
        
        std::vector<std::vector<std::vector<std::vector<util::SignalShaping> > > > fSignalShapingVec;
        // Field response.
        
        std::vector<DoubleVec4 >fFieldResponseVec;
        
        // Electronics response.
        
        std::vector<DoubleVec> fElectResponse;
        
        bool fPrintResponses;
        bool fManualInterpolation;
        
        // some diagnostic histograms
        
        TH1D* fHRawResponse[3];
        TH1D* fHStretchedResponse[3];
        TH1D* fHFullResponse[3];
        TH1D* fHSampledResponse[3];
        
        bool fHistDone[3];
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

// Do convolution.
template <class T> inline void util::SignalShapingServiceICARUS::Convolute(size_t channel, size_t wire, std::vector<T>& func) const
{
    std::cout << " before convoluting " << channel << " wire " << wire << " funcsize " << func.size() << std::endl;
    for (int jj=0;jj<4096;jj++)
       if(func[jj])
        std::cout << " printing func " << jj << " " << func[jj] << std::endl;
    SignalShaping(channel, wire).Convolute(func);
    std::cout << " after convoluting " << channel << " wire " << std::endl;

    //negative number
    int time_offset = FieldResponseTOffset(channel,0);
    std::cout << " after offset " << channel << " offset " << time_offset << " size " << func.size() << std::endl;
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
    
     std::cout << " end convoluting " << channel << " wire " << std::endl;
}

//----------------------------------------------------------------------
// Do deconvolution.
template <class T> inline void util::SignalShapingServiceICARUS::Deconvolute(size_t channel, std::vector<T>& func) const
{
    size_t ktype = 1;
    SignalShaping(channel, 0, ktype).Deconvolute(func);
    
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

//----------------------------------------------------------------------
// Do deconvolution.

template <class T> inline void util::SignalShapingServiceICARUS::Deconvolute(size_t channel, size_t wire, std::vector<T>& func) const
{
    size_t ktype = 1;
    SignalShaping(channel, wire, ktype).Deconvolute(func);
    
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
