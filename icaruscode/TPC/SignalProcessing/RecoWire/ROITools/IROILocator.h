///////////////////////////////////////////////////////////////////////
///
/// \file   IROILocator.h
///
/// \brief  This provides an interface for tools which are tasked with
///         finding ROI's in input waveforms. This allows different
///         approaches to be tried interchangeably
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IROILocator_H
#define IROILocator_H

#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"

namespace icarus_tool
{
    class IROILocator
    {
    public:
        virtual ~IROILocator() noexcept = default;
        
        virtual void   configure(const fhicl::ParameterSet& pset) = 0;
        
        // Define the waveform container
        using VectorBool  = std::vector<bool>;
        using VectorFloat = std::vector<float>;
        using ArrayBool   = std::vector<VectorBool>;
        using ArrayFloat  = std::vector<VectorFloat>;

        using PlaneIDVec  = std::vector<geo::PlaneID>;
        
        // Find the ROI's
        virtual void FindROIs(const ArrayFloat&, const geo::PlaneID&, ArrayFloat&, ArrayBool&) const = 0;
    };
}

#endif
