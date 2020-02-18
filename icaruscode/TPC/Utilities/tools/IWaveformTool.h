///////////////////////////////////////////////////////////////////////
///
/// \file   IWaveformTool.h
///
/// \brief  This is the interface class for a tool to aid in the analysis
///         of waveforms.
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IWaveformTool_H
#define IWaveformTool_H

#include "fhiclcpp/ParameterSet.h"
#include "TProfile.h"

#include "icaruscode/TPC/Utilities/tools/SignalProcessingDefs.h"

namespace icarus_tool
{
    template <class T> using Waveform = std::vector<T>;
    
    enum HistogramType : int
    {
        WAVEFORM,
        WAVELESSAVE,
        EROSION,
        DILATION,
        AVERAGE,
        DIFFERENCE,
        OPENING,
        CLOSING,
        DOPENCLOSING,
        LASTELEMENT
    };
    
    using HistogramMap = std::map<int, TProfile*>;
    
    class IWaveformTool
    {
    public:
        virtual ~IWaveformTool() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)  = 0;
        
        using PeakTuple    = std::tuple<size_t,size_t,size_t>;   // first bin, peak bin, last bin
        using PeakTupleVec = std::vector<PeakTuple>;             // The collection of candidate peaks
        
        virtual void triangleSmooth(const std::vector<float>&,  std::vector<float>&,  size_t = 0)                           const = 0;
        virtual void triangleSmooth(const std::vector<double>&, std::vector<double>&, size_t = 0)                           const = 0;
        virtual void medianSmooth(  const std::vector<float>&,  std::vector<float>&,  size_t = 3)                           const = 0;
        virtual void medianSmooth(  const std::vector<double>&, std::vector<double>&, size_t = 3)                           const = 0;
        virtual void getTruncatedMean(const std::vector<float>&, float&, int&)                                              const = 0;
        virtual void getTruncatedMean(const std::vector<double>&, double&, int&)                                            const = 0;
        virtual void getTruncatedMeanRMS(const std::vector<float>&, float, float&, float&, float&, int&)                    const = 0;
        virtual void getTruncatedMeanRMS(const std::vector<double>&, double, double&, double&, double&, int&)               const = 0;
        virtual void firstDerivative(const std::vector<float>&,  std::vector<float>&)                                       const = 0;
        virtual void firstDerivative(const std::vector<double>&, std::vector<double>&)                                      const = 0;
        virtual void findPeaks(std::vector<float>::iterator,  std::vector<float>::iterator,  PeakTupleVec&, float,  size_t) const = 0;
        virtual void findPeaks(std::vector<double>::iterator, std::vector<double>::iterator, PeakTupleVec&, double, size_t) const = 0;
        virtual void getFFTPower(const std::vector<float>& inputVec, std::vector<float>& outputPowerVec)                    const = 0;
        virtual void getFFTPower(const std::vector<double>& inputVec, std::vector<double>& outputPowerVec)                  const = 0;
        
        virtual void getErosionDilationAverageDifference(const Waveform<short>&,                //< Input waveform
                                                         int,                                   //< Structuring element
                                                         HistogramMap&,                         //< Map of histograms to fill
                                                         Waveform<short>&,                      //< Output erosion vector
                                                         Waveform<short>&,                      //< Output dilation vector
                                                         Waveform<short>&,                      //< Output ave erosion & dilation
                                                         Waveform<short>&)         const = 0;   //< Output diff erosion and dilation
       virtual void getErosionDilationAverageDifference(const Waveform<float>&,
                                                         int,
                                                         HistogramMap&,
                                                         Waveform<float>&,
                                                         Waveform<float>&,
                                                         Waveform<float>&,
                                                         Waveform<float>&)         const = 0;
        virtual void getErosionDilationAverageDifference(const Waveform<double>&,
                                                         int,
                                                         HistogramMap&,
                                                         Waveform<double>&,
                                                         Waveform<double>&,
                                                         Waveform<double>&,
                                                         Waveform<double>&)         const = 0;
        
        virtual void getOpeningAndClosing(const Waveform<short>&,                                //< Input erosions vector
                                          const Waveform<short>&,                                //< Input dilation vector
                                          int,                                                   //< Structuring element
                                          HistogramMap&,                                         //< Map of histograms to fill
                                          Waveform<short>&,                                      //< Output closing vector
                                          Waveform<short>&)         const = 0;                   //< Output opening vector
        virtual void getOpeningAndClosing(const Waveform<float>&,                                //< Input erosions vector
                                          const Waveform<float>&,                                //< Input dilation vector
                                          int,                                                   //< Structuring element
                                          HistogramMap&,                                         //< Map of histograms to fill
                                          Waveform<float>&,                                      //< Output closing vector
                                          Waveform<float>&)         const = 0;                   //< Output opening vector
        virtual void getOpeningAndClosing(const Waveform<double>&,                               //< Input erosions vector
                                          const Waveform<double>&,                               //< Input dilation vector
                                          int,                                                   //< Structuring element
                                          HistogramMap&,                                         //< Map of histograms to fill
                                          Waveform<double>&,                                     //< Output closing vector
                                          Waveform<double>&)         const = 0;                  //< Output opening vector
    };
}

#endif
