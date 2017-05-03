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

namespace icarus_tool
{
    class IWaveformTool
    {
    public:
        virtual ~IWaveformTool() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)  = 0;
        
        using PeakTuple    = std::tuple<size_t,size_t,size_t>;   // first bin, peak bin, last bin
        using PeakTupleVec = std::vector<PeakTuple>;             // The collection of candidate peaks
        
        virtual void triangleSmooth(std::vector<float>&,  std::vector<float>&,  size_t = 0)                                 const = 0;
        virtual void triangleSmooth(std::vector<double>&, std::vector<double>&, size_t = 0)                                 const = 0;
        virtual void firstDerivative(std::vector<float>&,  std::vector<float>&)                                             const = 0;
        virtual void firstDerivative(std::vector<double>&, std::vector<double>&)                                            const = 0;
        virtual void findPeaks(std::vector<float>::iterator,  std::vector<float>::iterator,  PeakTupleVec&, float,  size_t) const = 0;
        virtual void findPeaks(std::vector<double>::iterator, std::vector<double>::iterator, PeakTupleVec&, double, size_t) const = 0;
    };
}

#endif
