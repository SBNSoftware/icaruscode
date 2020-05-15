///////////////////////////////////////////////////////////////////////
///
/// \file   IFieldResponse.h
///
/// \brief  This is the interface class for a tool to handle the field response
///         It is assumed that the field response is described in a to be
///         input histogram, the member
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IFieldResponse_H
#define IFieldResponse_H

#include "fhiclcpp/ParameterSet.h"
#include "icaruscode/TPC/Utilities/tools/SignalProcessingDefs.h"

namespace art
{
    class TFileDirectory;
}

namespace icarus_tool
{
    class IFieldResponse
    {
    public:
        virtual ~IFieldResponse() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)                         = 0;
        virtual void setResponse(double weight, double correct3D, double timeScaleFctr) = 0;
        virtual void outputHistograms(art::TFileDirectory&)                       const = 0;
        
        virtual size_t                          getPlane()                        const = 0;
        virtual size_t                          getNumBins()                      const = 0;
        virtual double                          getBinCenter(int bin)             const = 0;
        virtual double                          getBinContent(int bin)            const = 0;
        virtual double                          getLowEdge()                      const = 0;
        virtual double                          getHighEdge()                     const = 0;
        virtual double                          getBinWidth()                     const = 0;
        virtual double                          getTOffset()                      const = 0;
        virtual double                          getIntegral()                     const = 0;
        virtual double                          interpolate(double x)             const = 0;
        
        virtual const icarusutil::TimeVec&      getResponseVec()                  const = 0;
        virtual const icarusutil::FrequencyVec& getResponseFFTVec()               const = 0;
    };
}

#endif
