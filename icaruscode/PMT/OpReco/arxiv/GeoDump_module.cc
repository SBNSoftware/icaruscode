////////////////////////////////////////////////////////////////////////
// Class:       GeoDump
// Plugin Type: analyzer (art v3_01_01)
// File:        GeoDump_module.cc
//
// Generated at Tue Feb 12 03:06:11 2019 by Kazuhiro Terao using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"

class GeoDump;


class GeoDump : public art::EDAnalyzer {
public:
  explicit GeoDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GeoDump(GeoDump const&) = delete;
  GeoDump(GeoDump&&) = delete;
  GeoDump& operator=(GeoDump const&) = delete;
  GeoDump& operator=(GeoDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


GeoDump::GeoDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void GeoDump::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  auto geop = lar::providerFrom<geo::Geometry>(); 

  std::cout << geop->Ncryostats() << " cryostats" << std::endl;
  //for(auto const& cryostat : geop->Cryostats()) {
  for(size_t c=0; c<geop->Ncryostats(); ++c) {
    auto const& cryostat = geop->Cryostat(c);
    std::cout<<"Cryostat " << cryostat.ID() << " ... " << cryostat.NTPC() << " TPCs and " << cryostat.NOpDet() << " opdets" << std::endl;
    auto const& cryobox  = cryostat.Boundaries();
    std::cout<<"  (" << cryobox.MinX() << "," << cryobox.MinY() << "," << cryobox.MinZ() << ")" 
	     << " => ("
	     << cryobox.MaxX() << "," << cryobox.MaxY() << "," << cryobox.MaxZ() << ")" << std::endl;

    for(size_t t=0; t<geop->NTPC(); ++t) {
      // Why there's no CryostatGeo::TPCs()?
      if(!cryostat.HasTPC(t)) continue;
      auto const& tpc = cryostat.TPC(t);
      std::cout<<"    TPC ID=" << t << " ... " << tpc.Nplanes() << " planes" << std::endl;
      auto const& tpcbox = tpc.BoundingBox();
      auto const& tpcabox = tpc.ActiveBoundingBox();
      std::cout<<"    BB (" << tpcbox.MinX() << "," << tpcbox.MinY() << "," << tpcbox.MinZ() << ")" 
	       << " => ("
	       << tpcbox.MaxX() << "," << tpcbox.MaxY() << "," << tpcbox.MaxZ() << ")" << std::endl;
      std::cout<<"    Active BB (" << tpcabox.MinX() << "," << tpcabox.MinY() << "," << tpcabox.MinZ() << ")" 
	       << " => ("
	       << tpcabox.MaxX() << "," << tpcabox.MaxY() << "," << tpcabox.MaxZ() << ")" << std::endl;
      
      for(size_t p=0; p<tpc.Nplanes(); ++p) {
	auto const& plane = tpc.Plane(p);
	std::cout<<"      Plane ID=" << p << " ... " << plane.Nwires() << " wires, thetaZ=" << plane.ThetaZ() << std::endl;
      }
    }
    for(size_t o=0; o<cryostat.NOpDet(); ++o) {
      auto const& opdet = cryostat.OpDet(o);
      std::cout << "OpDet ID="<<o<< " ... (" << opdet.GetCenter().x() << "," << opdet.GetCenter().y() << "," << opdet.GetCenter().z() << ")" << std::endl;
    }
    for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {
      std::cout << "OpChannel " << opch << " => OpDet " << geop->OpDetFromOpChannel(opch) << std::endl;
    }
  }
  std::cout<<std::endl;

}

DEFINE_ART_MODULE(GeoDump)
