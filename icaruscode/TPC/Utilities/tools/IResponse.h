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

#ifndef IResponse_H
#define IResponse_H

#include "fhiclcpp/ParameterSet.h"
#include "icaruscode/TPC/Utilities/tools/SignalProcessingDefs.h"
#include "icaruscode/TPC/Utilities/tools/IFieldResponse.h"
#include "icaruscode/TPC/Utilities/tools/IElectronicsResponse.h"
#include "icaruscode/TPC/Utilities/tools/IFilter.h"

namespace icarusutil
{
    class SignalShapingICARUS;
}

namespace art
{
    class TFileDirectory;
}

namespace icarus_tool
{

class IResponse
{
public:
    virtual ~IResponse() noexcept = default;

    /**
     *  @brief Setup routines
     *
     *  @param output                the object containting the art output
     *  @param clusHitPairVector     List of 3D hits to output as "extreme" space points
     */
    virtual void configure(const fhicl::ParameterSet& pset)   = 0;
    virtual void setResponse(double sampling_rate, double weight) = 0; // rate in ns
    virtual void outputHistograms(double sampling_rate, art::TFileDirectory&) const = 0;

    /**
     *  @brief Return the plane these functions represent
     */
    virtual size_t                                  getPlane()                             const = 0;

    /**
     *  @brief Recover the individual response elements
     */
    virtual const IFieldResponse*                   getFieldResponse()                     const = 0;
    virtual const IElectronicsResponse*             getElectronicsResponse()               const = 0;
    virtual const IFilter*                          getFilter()                            const = 0;
    
    /**
     *  @brief here recover the combined response elements
     *
     *  @param output                the object containting the art output
     *  @param clusHitPairVector     List of 3D hits to output as "extreme" space points
     */
    virtual size_t                                   getNumberTimeSamples()                 const = 0;
    virtual const icarusutil::TimeVec&               getResponse()                          const = 0;
    virtual const icarusutil::FrequencyVec&          getConvKernel()                        const = 0;
    virtual const icarusutil::FrequencyVec&          getDeconvKernel()                      const = 0;
    virtual double                                   getTOffset()                           const = 0;
};
}

#endif
