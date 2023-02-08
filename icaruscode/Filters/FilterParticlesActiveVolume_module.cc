////////////////////////////////////////////////////////////////////////
/// \file  FilterNoDirtNeutrinos_module.cc
/// \brief Simple EDFilter to require neutrino interaction in TPC
///
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef FILTER_FILTERPARTICLESACTIVEVOLUME_H
#define FILTER_FILTERPARTICLESACTIVEVOLUME_H

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// LArSoft Includes
#include "nug4/ParticleNavigation/ParticleList.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/sim.h"
#include "larcore/Geometry/Geometry.h"

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

namespace simb{
  class MCTruth;
}

namespace sim{
  class ParticleList;
}

///Geant4 interface 
namespace simfilter {  
 
  class FilterParticlesActiveVolume : public art::EDFilter
  {  
  // explicit EDFilter(ParameterSet const&)  
  public:

    explicit FilterParticlesActiveVolume(fhicl::ParameterSet const &pset);
    virtual ~FilterParticlesActiveVolume();
    
    bool filter(art::Event&) ;
    virtual void reconfigure(fhicl::ParameterSet const&)  ;
      
    virtual void beginJob()  ;
    /*
    virtual void endJob()  ;
    virtual bool beginRun(art::Run &)  ;
    virtual bool endRun(art::Run &)  ;
    virtual bool beginSubRun(art::SubRun &)  ;
    virtual bool endSubRun(art::SubRun &)  ;
    */
    private:

      double fXmax;
      double fXmin;
      double fZmax;
      double fZmin;
      double fYmax;
      double fYmin;
      int finActive;
      int filtpart;
  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterParticlesActiveVolume::FilterParticlesActiveVolume(fhicl::ParameterSet const& pset) : EDFilter{pset}
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  FilterParticlesActiveVolume::~FilterParticlesActiveVolume()
  {
  }

  //-----------------------------------------------------------------------
  void FilterParticlesActiveVolume::beginJob()
  {
    art::ServiceHandle<geo::Geometry> geo;
  
  }

  //-----------------------------------------------------------------------
  void FilterParticlesActiveVolume::reconfigure(fhicl::ParameterSet const& p)
  {
      finActive        = p.get< int >("inActive");
      if (finActive==0) {
      fXmax            = p.get< double >("Xmax");
      fYmax            = p.get< double >("Ymax");
      fZmax            = p.get< double >("Zmax");
      fXmin            = p.get< double >("Xmin");
      fYmin            = p.get< double >("Ymin");
      fZmin            = p.get< double >("Zmin");
      }
      if (finActive==1) {
          fXmax=-999;
          fYmax=-999;
          fZmax=-999;
          fXmin=-999;
          fYmin=-999;
          fZmin=-999;
      }
      filtpart       =  p.get< int >("filterpart");
      
    return;
  }

  //-----------------------------------------------------------------------
  bool FilterParticlesActiveVolume::filter(art::Event& evt)
  {
    bool interactionDesired(false);
    //get the list of particles from this event
    art::ServiceHandle<geo::Geometry> geom;

      
    // * MC truth information

      //std::vector< art::Handle< std::vector<simb::MCTruth> > > allmclists;
      //evt.getManyByType(allmclists);
      auto allmclists = evt.getMany< std::vector<simb::MCTruth> >();


      std::cout << fXmin << " " << fXmax << " " << fYmin << " " <<fYmax << " " << fZmin << " " << fZmax << std::endl;
 
      geo::CryostatGeo const& cryo0 = geom->Cryostat(geo::CryostatID{0});
      geo::CryostatGeo const& cryo1 = geom->Cryostat(geo::CryostatID{1});

      geo::TPCGeo const& tpc00 = cryo0.TPC(0);
      geo::Point_t xyzcenter00 = tpc00.GetActiveVolumeCenter();
      std::cout << xyzcenter00.X() << " " << xyzcenter00.Y() << " " << xyzcenter00.Z() << std::endl;

      geo::TPCGeo const& tpc01 = cryo0.TPC(1);
      geo::Point_t xyzcenter01 = tpc01.GetActiveVolumeCenter();
      std::cout << xyzcenter01.X() << " " << xyzcenter01.Y() << " " << xyzcenter01.Z() << std::endl;

      geo::TPCGeo const& tpc10 = cryo1.TPC(0);
      geo::Point_t xyzcenter10 = tpc10.GetActiveVolumeCenter();
      std::cout << xyzcenter10.X() << " " << xyzcenter10.Y() << " " << xyzcenter10.Z() << std::endl;

      geo::TPCGeo const& tpc11 = cryo1.TPC(1);
      geo::Point_t xyzcenter11 = tpc11.GetActiveVolumeCenter();
      std::cout << xyzcenter11.X() << " " << xyzcenter11.Y() << " " << xyzcenter11.Z() << std::endl;

      double h00=tpc00.ActiveHalfHeight();
      double h01=tpc01.ActiveHalfHeight();
      double h10=tpc10.ActiveHalfHeight();
      double h11=tpc11.ActiveHalfHeight();
      double w00=tpc00.ActiveHalfWidth();
      double w01=tpc01.ActiveHalfWidth();
      double w10=tpc10.ActiveHalfWidth();
      double w11=tpc11.ActiveHalfWidth();
      double l00=tpc00.ActiveLength();
      double l01=tpc01.ActiveLength();
      double l10=tpc10.ActiveLength();
      double l11=tpc11.ActiveLength();
      
      std::cout << h00 << " " << w00 << " " << l00 << std::endl;
      std::cout << h01 << " " << w01 << " " << l01 << std::endl;
      std::cout << h10 << " " << w10 << " " << l10 << std::endl;
      std::cout << h11 << " " << w11 << " " << l11 << std::endl;
      
      for(size_t mcl = 0; mcl < allmclists.size(); ++mcl){
          art::Handle< std::vector<simb::MCTruth> > mclistHandle = allmclists[mcl];
          for(size_t m = 0; m < mclistHandle->size(); ++m){
              art::Ptr<simb::MCTruth> mct(mclistHandle, m);
              for(int ipart=0;ipart<mct->NParticles();ipart++){
                  int pdg=mct->GetParticle(ipart).PdgCode();
                  double xx=mct->GetParticle(ipart).Vx();
                  double yy=mct->GetParticle(ipart).Vy();
                  double zz=mct->GetParticle(ipart).Vz();
      
		  if (finActive==1 && pdg==filtpart)
		    {
                      if (xx>(xyzcenter00.X()-w00) && xx<(xyzcenter00.X()+w00) && yy>(xyzcenter00.Y()-h00) && yy<(xyzcenter00.Y()+h00) && zz>(xyzcenter00.Z()-l00/2) && zz<(xyzcenter00.Z()+l00/2))
			{
			  interactionDesired = true;
			}
                      if (xx>(xyzcenter01.X()-w01) && xx<(xyzcenter01.X()+w01) && yy>(xyzcenter01.Y()-h01) && yy<(xyzcenter01.Y()+h01) && zz>(xyzcenter01.Z()-l01/2) && zz<(xyzcenter01.Z()+l01/2))
			{
			  interactionDesired = true;
			}
                      if (xx>(xyzcenter10.X()-w10) && xx<(xyzcenter10.X()+w10) && yy>(xyzcenter10.Y()-h10) && yy<(xyzcenter10.Y()+h10) && zz>(xyzcenter10.Z()-l10/2) && zz<(xyzcenter10.Z()+l10/2))
			{
			  interactionDesired = true;
			}
                      if (xx>(xyzcenter11.X()-w11) && xx<(xyzcenter11.X()+w11) && yy>(xyzcenter11.Y()-h11) && yy<(xyzcenter11.Y()+h11) && zz>(xyzcenter11.Z()-l11/2) && zz<(xyzcenter11.Z()+l11/2))
			{
			  interactionDesired = true;
			}
                        
                      if (finActive==0 && pdg==filtpart)
                        {
			  if (xx>fXmin && xx<fXmax && yy>fYmin && yy<fYmax && zz>fZmin && zz<fZmax)
			    {
			      interactionDesired = true;
			    }
                        }
        
		    }
	  //	std::cout << "FilterNoDirtParticles: i is " << i << std::endl ;
	  // Now walk through trajectory and see if it enters the TPC
	   // trajectory loop
	 // end Genie particle
	      }
          }
      }// loop on MCPHandle

  return interactionDesired;
    
  } // end FilterNoDirtParticles()function

} // namespace simfilter

namespace simfilter {

  DEFINE_ART_MODULE(FilterParticlesActiveVolume)

} // namespace simfilter

#endif // FILTER_FILTERNEUTRINOSACTIVEVOLUME_H
