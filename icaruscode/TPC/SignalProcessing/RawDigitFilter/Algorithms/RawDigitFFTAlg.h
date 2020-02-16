#ifndef RAWDIGITFFTALG_H
#define RAWDIGITFFTALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       RawDigitFFTAlg
// Module Type: producer
// File:        RawDigitFFTAlg.h
//
//              This module provides some basic Fast Fourier Transform
//              algorithms for operating on RawDigit waveforms
//
// Configuration parameters:
//
// FillHistograms        - Turn on histogram filling for diagnostics
// RunFFTInputWires      - FFT analyze the input RawDigits if true - diagnostics
// RunFFTCorrectedWires  - FFT analyze the output RawDigits if true - diagnostics
//
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
//
////////////////////////////////////////////////////////////////////////

#include "RawDigitNoiseFilterDefs.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include "TProfile.h"

namespace icarus_tool
{
    class IWaveformTool;
    class IFilter;
}

namespace caldata
{
    
class RawDigitFFTAlg
{
public:

    // Copnstructors, destructor.
    RawDigitFFTAlg(fhicl::ParameterSet const & pset);
    ~RawDigitFFTAlg();

    // Provide for reinitialization if necessary
    void reconfigure(fhicl::ParameterSet const & pset);
    void initializeHists(art::ServiceHandle<art::TFileService>&);
    
    template <class T> void getFFTCorrection(std::vector<T>&, double) const;
    
    template <class T> void getFFTCorrection(std::vector<T>&, size_t) const;
    
    void filterFFT(std::vector<short>&, size_t, size_t, float pedestal=0.);
    
private:
    
    std::vector<bool>                                      fTransformViewVec;      ///< apply FFT transform to this view
    bool                                                   fFillHistograms;        ///< if true then will fill diagnostic hists
    std::string                                            fHistDirName;           ///< If writing histograms, the folder name
    std::vector<size_t>                                    fLoWireByPlane;         ///< Low wire for individual wire histograms
    std::vector<size_t>                                    fHiWireByPlane;         ///< Hi wire for individual wire histograms
    
    // Try to optimize the filter FFT function with static memory...
    std::map<size_t,std::vector<std::complex<float>>>      fFilterVecMap;
    std::vector<float>                                     fFFTInputVec;
    std::vector<std::complex<float>>                       fFFTOutputVec;
    std::vector<float>                                     fPowerVec;

    std::vector<std::vector<TProfile*>>                    fFFTPowerVec;
    std::vector<TProfile*>                                 fAveFFTPowerVec;
    std::vector<TProfile*>                                 fConvFFTPowerVec;
    std::vector<TProfile*>                                 fFilterFuncVec;

    std::unique_ptr<icarus_tool::IWaveformTool>            fWaveformTool;
    std::map<size_t,std::unique_ptr<icarus_tool::IFilter>> fFilterToolMap;
    
    std::unique_ptr<Eigen::FFT<float>>                     fEigenFFT;

    // Useful services, keep copies for now (we can update during begin run periods)
    detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();   ///< Detector properties service
};
    
} // end caldata namespace

#endif
