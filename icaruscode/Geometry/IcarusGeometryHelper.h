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
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/AuxDetGeo.h"

#include <memory>
#include <vector>

// Forward declarations
//
class TString;

namespace geo
{
  class ChannelMapAlg;
}

// Declaration
//
namespace Icarus
{
  class IcarusGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
  
    IcarusGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry &reg );
    ~IcarusGeometryHelper() throw();

    // Public interface for ExptGeoHelperInterface (for reference purposes)
    //
    // Configure and initialize the channel map.
    //
    // void  ConfigureChannelMapAlg( const TString & detectorName, 
    //                               fhicl::ParameterSet const & sortingParam,
    //                               std::vector<geo::CryostatGeo*> & c,
    //				     std::vector<geo::AuxDetGeo*>   & ad );
    //
    // Returns null pointer if the initialization failed
    // NOTE:  the sub-class owns the ChannelMapAlg object
    //
    // std::shared_ptr<const geo::ChannelMapAlg> & GetChannelMapAlg() const;
  
  private:

    void doConfigureChannelMapAlg( fhicl::ParameterSet const & sortingParameters, geo::GeometryCore* geom ) override;
    ChannelMapAlgPtr_t doGetChannelMapAlg() const override;
    
    fhicl::ParameterSet fPset; ///< copy of configuration parameter set
    //art::ActivityRegistry fReg; ///< copy of activity registry
    std::shared_ptr<geo::ChannelMapAlg> fChannelMap;
  
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(Icarus::IcarusGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // Icarus_ExptGeoHelperInterface_h
