////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterICARUS.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "icaruscode/Geometry/GeoObjectSorterICARUS.h"

#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

namespace geo{

  //----------------------------------------------------------------------------
  // Define sort order for cryostats in standard configuration
  static bool sortAuxDetStandard(const AuxDetGeo* ad1, const AuxDetGeo* ad2)
  {
    std::string type1 = "", type2 = "";
    switch (ad1->NSensitiveVolume()) {
        case 20 : type1 = "MINOS"; break;
        case 16 : type1 = "CERN"; break;
        case 64 : type1 = "DC"; break;
    }
    switch (ad2->NSensitiveVolume()) {
        case 20 : type2 = "MINOS"; break;
        case 16 : type2 = "CERN"; break;
        case 64 : type2 = "DC"; break;
    }

    // sort based off of GDML name, module number
    std::string ad1name = (ad1->TotalVolume())->GetName();
    std::string ad2name = (ad2->TotalVolume())->GetName();
    // assume volume name is "volAuxDet_<type>_module_###_<region>"
    std::string base1 = "volAuxDet_"+type1+"_module_";
    std::string base2 = "volAuxDet_"+type2+"_module_";

    int ad1Num = atoi( ad1name.substr( base1.size(), 3).c_str() );
    int ad2Num = atoi( ad2name.substr( base2.size(), 3).c_str() );  
  
    return ad1Num < ad2Num;
   
  }

  //----------------------------------------------------------------------------
  // Define sort order for cryostats in standard configuration
  static bool sortAuxDetSensitiveStandard(const AuxDetSensitiveGeo* ad1, const AuxDetSensitiveGeo* ad2)
  {
    std::string type1 = "", type2 = "";

    // sort based off of GDML name, assuming ordering is encoded
    std::string ad1name = (ad1->TotalVolume())->GetName();
    std::string ad2name = (ad2->TotalVolume())->GetName();

    if ( ad1name.find("MINOS") != std::string::npos ) type1 = "MINOS";
    if ( ad1name.find("CERN") != std::string::npos ) type1 = "CERN";
    if ( ad1name.find("DC") != std::string::npos ) type1 = "DC";
    if ( ad2name.find("MINOS") != std::string::npos ) type2 = "MINOS";
    if ( ad2name.find("CERN") != std::string::npos ) type2 = "CERN";
    if ( ad2name.find("DC") != std::string::npos ) type2 = "DC";

    // assume volume name is "volAuxDetSensitive_<type>_module_###_strip_##"
    std::string baseMod1 = "volAuxDetSensitive_"+type1+"module_";
    std::string baseStr1 = "volAuxDetSensitive_"+type1+"module_###_strip_";
    std::string baseMod2 = "volAuxDetSensitive_"+type2+"module_";
    std::string baseStr2 = "volAuxDetSensitive_"+type2+"module_###_strip_";

    int ad1Num = atoi( ad1name.substr( baseMod1.size(), 3).c_str() );
    int ad2Num = atoi( ad2name.substr( baseMod2.size(), 3).c_str() );

    if(ad1Num!=ad2Num) return ad1Num < ad2Num;

    ad1Num = atoi( ad1name.substr( baseStr1.size(), 2).c_str() );
    ad2Num = atoi( ad2name.substr( baseStr2.size(), 2).c_str() );

    
    return ad1Num < ad2Num;
   
  }

  //----------------------------------------------------------------------------
  // Define sort order for cryostats in standard configuration
  static bool sortCryoStandard(const CryostatGeo* c1, const CryostatGeo* c2)
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.}; 
    c1->LocalToWorld(local, xyz1);
    c2->LocalToWorld(local, xyz2);

    return xyz1[0] < xyz2[0];   
  }


  //----------------------------------------------------------------------------
  // Define sort order for tpcs in standard configuration.
  static bool sortTPCStandard(const TPCGeo* t1, const TPCGeo* t2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    // sort TPCs according to x
    if(xyz1[0] < xyz2[0]) return true;

    return false;
  }

  const double EPSILON = 0.000001;
  
  //----------------------------------------------------------------------------
  // Define sort order for planes in standard configuration
  static bool sortPlaneStandard(const PlaneGeo* p1, const PlaneGeo* p2) 
  {
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double local[3] = {0.};
    p1->LocalToWorld(local, xyz1);
    p2->LocalToWorld(local, xyz2);

    //if the planes are in the same drift coordinate, lower Z is first plane
    if( std::abs(xyz1[0] - xyz2[0]) < EPSILON)
      return xyz1[2] < xyz2[2];

    //else
    // drift direction is negative, plane number increases in drift direction
    return xyz1[0] > xyz2[0];
  }


  //----------------------------------------------------------------------------
  bool sortWireStandard(WireGeo* w1, WireGeo* w2){
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};

    w1->GetCenter(xyz1); w2->GetCenter(xyz2);

    //we have horizontal wires...
    if( std::abs(xyz1[2]-xyz2[2]) < EPSILON)
      return xyz1[1] < xyz2[1];

    //in the other cases...
    return xyz1[2] < xyz2[2];
  }

  //----------------------------------------------------------------------------
  GeoObjectSorterICARUS::GeoObjectSorterICARUS(fhicl::ParameterSet const& p)
  {
  }

  //----------------------------------------------------------------------------
  GeoObjectSorterICARUS::~GeoObjectSorterICARUS()
  {
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterICARUS::SortAuxDets(std::vector<geo::AuxDetGeo*> & adgeo) const
  {
    std::sort(adgeo.begin(), adgeo.end(), sortAuxDetStandard);
    
    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterICARUS::SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo) const
  {
    std::sort(adsgeo.begin(), adsgeo.end(), sortAuxDetSensitiveStandard);
    
    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterICARUS::SortCryostats(std::vector<geo::CryostatGeo*> & cgeo) const
  {
    std::sort(cgeo.begin(), cgeo.end(), sortCryoStandard);
    
    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterICARUS::SortTPCs(std::vector<geo::TPCGeo*>  & tgeo) const
  {
    
    std::sort(tgeo.begin(), tgeo.end(), sortTPCStandard);

    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterICARUS::SortPlanes(std::vector<geo::PlaneGeo*> & pgeo,
					   geo::DriftDirection_t  const& driftDir) const
  {
    // sort the planes to increase in drift direction
    // The drift direction has to be set before this method is called.  It is set when
    // the CryostatGeo objects are sorted by the CryostatGeo::SortSubVolumes method
    if     (driftDir == geo::kPosX) std::sort(pgeo.rbegin(), pgeo.rend(), sortPlaneStandard);
    else if(driftDir == geo::kNegX) std::sort(pgeo.begin(),  pgeo.end(),  sortPlaneStandard);
    else if(driftDir == geo::kUnknownDrift)
      throw cet::exception("TPCGeo") << "Drift direction is unknown, can't sort the planes\n";

    return;
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterICARUS::SortWires(std::vector<geo::WireGeo*> & wgeo) const
  {
    std::sort(wgeo.begin(), wgeo.end(), sortWireStandard);

    return;
  }

}
