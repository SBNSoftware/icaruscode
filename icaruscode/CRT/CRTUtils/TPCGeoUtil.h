#ifndef IC_CRT_TPCGEOUTIL_H
#define IC_CRT_TPCGEOUTIL_H
// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"

// c++
#include <vector>

namespace icarus {
  namespace TPCGeoUtil {
    int DetectedInTPC(std::vector<art::Ptr<recob::Hit>> const& hits, geo::GeometryCore const& geometry);
    // Work out the drift limits for a collection of hits
    std::pair<double, double> XLimitsFromHits(geo::GeometryCore const& geometry,
                                              geo::WireReadoutGeom const& wireReadout,
                                              std::vector<art::Ptr<recob::Hit>> const& hits);
    // Is point inside given TPC
    bool InsideTPC(geo::Point_t point, const geo::TPCGeo& tpc, double buffer);
    int DriftDirectionFromHits(const geo::GeometryCore& geometry, std::vector<art::Ptr<recob::Hit>> const& hits);
  } // namespace TPCGeoUtil
} // namespace icarus
#endif
