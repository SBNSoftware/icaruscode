///////////////////////////////////////////////////////////////////////
///
/// \file   IResponse.h
///
/// \brief  This is the interface class for a tool to handle the field response
///         It is assumed that the field response is described in a to be
///         input histogram, the member
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef shprIResponse_H
#define shprIResponse_H

#include "fhiclcpp/ParameterSet.h"

namespace util
{
    class SignalShaper;
}

namespace art
{
    class TFileDirectory;
}

namespace icarus_tool
{
    class IFieldResponse;
    class IElectronicsResponse;
    class shprIFilter;
    
    class shprIResponse
    {
    public:
        virtual ~shprIResponse() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)   = 0;
        virtual void setResponse(double weight)                   = 0;
        virtual void outputHistograms(art::TFileDirectory&) const = 0;
        
        virtual size_t                      getPlane()               const = 0;
        
        virtual const IFieldResponse*       getFieldResponse()       const = 0;
        virtual const IElectronicsResponse* getElectronicsResponse() const = 0;
        virtual const shprIFilter*            getFilter()              const = 0;
        
        virtual const util::SignalShaper&   getSignalShaper()        const = 0;
    };
}

#endif
