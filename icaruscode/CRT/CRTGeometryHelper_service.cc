///////////////////////////////////////////////////////////////////////////////
/// \file CRTGeometryHelper_service.cc
///
/// Originally ported from LArIATAuxDetGeometryHelper.h (Author: brebel@fnal.gov)
/// and modified for SBND (Author: mastbaum@uchicago.edu)
///
/// \version $Id: 0.0 $
/// \author chilge@rams.colostate.edu
///////////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/AuxDetChannelMapAlg.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "icaruscode/CRT/CRTGeometryHelper.h"
#include "icaruscode/CRT/CRTChannelMapAlg.h"
#include <memory>

namespace icarus {
namespace crt{

  //------------------------------------------------------------------------
  CRTGeometryHelper::CRTGeometryHelper(
      fhicl::ParameterSet const& pset)
  : fPset(pset) {}

  //------------------------------------------------------------------------
  CRTGeometryHelper::AuxDetChannelMapAlgPtr_t
  CRTGeometryHelper::doConfigureAuxDetChannelMapAlg(
      fhicl::ParameterSet const& sortingParameters) const {

    return std::make_unique<geo::CRTChannelMapAlg>(sortingParameters);

  }

}  // namespace crt
}  // namespace icarus

DEFINE_ART_SERVICE_INTERFACE_IMPL(icarus::crt::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface)
