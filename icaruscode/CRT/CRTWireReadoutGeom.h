///////////////////////////////////////////////////////////////////////////////
/// \file CRTWireReadoutGeom.h
/// \brief Algorithm class for SBND auxiliary detector channel mapping
///
/// Originally ported from AuxDetChannelMapLArIATAlg.cxx (Author: brebel@fnal.gov)
/// and modified for SBND (Author: mastbaum@uchicago.edu)
/// \version $Id: 0.0 $
/// \author chilge@rams.colostate.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef ICARUS_CRTWireReadoutGeom_h
#define ICARUS_CRTWireReadoutGeom_h

#include "larcorealg/Geometry/AuxDetWireReadoutGeom.h"
#include "fhiclcpp/ParameterSet.h"
#include "TVector3.h"
#include <vector>

//namespace icarus{
//namespace crt {
namespace geo{

  class CRTWireReadoutGeom : {
  public:
    CRTWireReadoutGeom(fhicl::ParameterSet const& p);

    void Initialize(AuxDetGeometryData_t& geodata) override;

    void Uninitialize() override;

    uint32_t PositionToAuxDetChannel(
        geo::Point_t const& worldLoc,
        std::vector<geo::AuxDetGeo> const& auxDets,
        size_t& ad,
        size_t& sv) const override;

    geo::Point_t AuxDetChannelToPosition(
        uint32_t channel,
        std::string const& auxDetName,
        std::vector<geo::AuxDetGeo> const& auxDets) const override;


  private:
  };

}  // namespace crt
//} //namespace icarus

#endif  // ICARUS_CRTWireReadoutGeom_h
