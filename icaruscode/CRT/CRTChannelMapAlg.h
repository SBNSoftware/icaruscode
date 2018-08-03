///////////////////////////////////////////////////////////////////////////////
/// \file CRTChannelMapAlg.h
/// \brief Algorithm class for SBND auxiliary detector channel mapping
///
/// Originally ported from AuxDetChannelMapLArIATAlg.cxx (Author: brebel@fnal.gov)
/// and modified for SBND (Author: mastbaum@uchicago.edu)
/// \version $Id: 0.0 $
/// \author chilge@rams.colostate.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef ICARUS_CRTChannelMapAlg_h
#define ICARUS_CRTChannelMapAlg_h

#include "larcorealg/Geometry/AuxDetChannelMapAlg.h"
#include "icaruscode/CRT/CRTGeoObjectSorter.h"
#include "fhiclcpp/ParameterSet.h"
#include "TVector3.h"
#include <vector>

//namespace icarus{
//namespace crt {
namespace geo{

  class CRTChannelMapAlg : public AuxDetChannelMapAlg {
  public:
    CRTChannelMapAlg(fhicl::ParameterSet const& p);

    void Initialize(AuxDetGeometryData_t& geodata) override;

    void Uninitialize() override;

    uint32_t PositionToAuxDetChannel(
        double const worldLoc[3],
        std::vector<geo::AuxDetGeo*> const& auxDets,
        size_t& ad,
        size_t& sv) const override;

    const TVector3 AuxDetChannelToPosition(
        uint32_t const& channel,
        std::string const& auxDetName,
        std::vector<geo::AuxDetGeo*> const& auxDets) const override;

  private:
    geo::CRTGeoObjectSorter fSorter; ///< Class to sort geo objects
  };

}  // namespace crt
//} //namespace icarus

#endif  // ICARUS_CRTChannelMapAlg_h

