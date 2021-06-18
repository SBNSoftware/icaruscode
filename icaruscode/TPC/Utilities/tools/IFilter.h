///////////////////////////////////////////////////////////////////////
///
/// \file   IFilter.h
///
/// \brief  This is the interface class for a tool to handle a filter
///         for the overall response
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IFilter_H
#define IFilter_H

#include "fhiclcpp/ParameterSet.h"
#include "icaruscode/TPC/Utilities/tools/SignalProcessingDefs.h"

namespace art
{
    class TFileDirectory;
}

namespace icarus_tool
{
    class IFilter
    {
    public:
        virtual ~IFilter() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)                          = 0;
        virtual void setResponse(size_t numBins, double correct3D, double timeScaleFctr) = 0;
        virtual void outputHistograms(art::TFileDirectory&)                        const = 0;
        
        virtual size_t                          getPlane()                         const = 0;
        virtual const icarusutil::FrequencyVec& getResponseVec()                   const = 0;
    };
}

#endif
