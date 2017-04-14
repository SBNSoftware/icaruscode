///////////////////////////////////////////////////////////////////////
///
/// \file   IElectronicsResponse.h
///
/// \brief  This is the interface class for a tool to handle the
///         electronics response.
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IElectronicsResponse_H
#define IElectronicsResponse_H

#include "fhiclcpp/ParameterSet.h"

namespace icarus_tool
{
    class IElectronicsResponse
    {
    public:
        virtual ~IElectronicsResponse() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)   = 0;
        virtual void setResponse(size_t numBins, double binWidth) = 0;
        
        virtual size_t                     getPlane()           const = 0;
        virtual double                     getASICGain()        const = 0;
        virtual double                     getASICShapingTime() const = 0;
        virtual const std::vector<double>& getResponseVec()     const = 0;
    };
}

#endif
