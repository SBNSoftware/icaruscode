////////////////////////////////////////////////////////////////////////////////
/// \file IcarusGeometryHelper.h
/// \brief Geometry helper service for Icarus geometries.
///
/// Handles Icarus-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class
///
/// \verion $Id
/// \author rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef Icarus_ExptGeoHelperInterface_h
#define Icarus_ExptGeoHelperInterface_h

#include "larcore/Geometry/ExptGeoHelperInterface.h"

#include <memory>

// Forward declarations
//
namespace geo
{
  class ChannelMapAlg;
}

namespace Icarus
{
  class IcarusGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:

    explicit IcarusGeometryHelper(fhicl::ParameterSet const& pset);

  private:

    ChannelMapAlgPtr_t
    doConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                             std::string const& detectorName) const override;

    fhicl::ParameterSet fPset; ///< copy of configuration parameter set
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(Icarus::IcarusGeometryHelper,
                                   geo::ExptGeoHelperInterface,
                                   SHARED)

#endif // Icarus_ExptGeoHelperInterface_h
