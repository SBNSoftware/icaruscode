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

  //------------------------------------------------------------------------
  CRTGeometryHelper::CRTGeometryHelper(
      fhicl::ParameterSet const& pset, art::ActivityRegistry &)
  : fPset(pset), fChannelMap() {}

  //------------------------------------------------------------------------
  void CRTGeometryHelper::doConfigureAuxDetChannelMapAlg(
      fhicl::ParameterSet const& sortingParameters,
      geo::AuxDetGeometryCore* geom) {
    fChannelMap = \
      std::make_shared<geo::CRTChannelMapAlg>(sortingParameters);

    if (fChannelMap) {
      geom->ApplyChannelMap(fChannelMap);
    }
  }

  //------------------------------------------------------------------------
  CRTGeometryHelper::AuxDetChannelMapAlgPtr_t
  CRTGeometryHelper::doGetAuxDetChannelMapAlg() const {
    return fChannelMap;
  }

}  // namespace icarus

DEFINE_ART_SERVICE_INTERFACE_IMPL(icarus::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface)

