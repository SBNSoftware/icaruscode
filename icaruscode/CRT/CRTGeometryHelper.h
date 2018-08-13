///////////////////////////////////////////////////////////////////////////////
/// \file CRTGeometryHelper.h
/// \brief Auxiliary detector geometry helper service for SBND geometries.
///
/// Handles SBND-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class.
///
/// Originally ported from LArIATAuxDetGeometryHelper.h (Author: brebel@fnal.gov)
/// and modified for SBND (Author: mastbaum@uchicago.edu)
///
/// \version $Id: 0.0 $
/// \author chilge@rams.colostate.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef ICARUS_CRTExptGeoHelperInterface_h
#define ICARUS_CRTExptGeoHelperInterface_h

#include "larcore/Geometry/AuxDetExptGeoHelperInterface.h"
#include "larcorealg/Geometry/AuxDetChannelMapAlg.h"
#include "icaruscode/CRT/CRTChannelMapAlg.h"
#include <memory>

namespace geo
{
  class CRTChannelMapAlg;
}

namespace icarus {
namespace crt{

  class CRTGeometryHelper : public geo::AuxDetExptGeoHelperInterface {
  public:

    CRTGeometryHelper(fhicl::ParameterSet const & pset,
                      art::ActivityRegistry &);

  private:

    void doConfigureAuxDetChannelMapAlg(
        fhicl::ParameterSet const& sortingParameters,
        geo::AuxDetGeometryCore* geom) override;

    AuxDetChannelMapAlgPtr_t doGetAuxDetChannelMapAlg() const override;

    fhicl::ParameterSet fPset; ///< Copy of configuration parameter set
    std::shared_ptr<geo::CRTChannelMapAlg> fCRTChannelMap; ///< Channel map

  };

}  //namespace crt
}  // namespace icarus

DECLARE_ART_SERVICE_INTERFACE_IMPL(icarus::crt::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface, LEGACY)

#endif  // ICARUS_CRTExptGeoHelperInterface_h

