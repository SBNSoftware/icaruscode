////////////////////////////////////////////////////////////////////////
/// \file  FilterNoDirtNeutrinos_module.cc
/// \brief Simple EDFilter to require neutrino interaction in TPC
///
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef FILTER_FILTERNEUTRINOSACTIVEVOLUME_H
#define FILTER_FILTERNEUTRINOSACTIVEVOLUME_H

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// LArSoft Includes
#include "nutools/ParticleNavigation/ParticleList.h"
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
 
  class FilterNeutrinosActiveVolume : public art::EDFilter
  {  
  // explicit EDFilter(ParameterSet const&)  
  public:

    explicit FilterNeutrinosActiveVolume(fhicl::ParameterSet const &pset);
    virtual ~FilterNeutrinosActiveVolume();
    
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
      int filtype;
      int typetofilt;
  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterNeutrinosActiveVolume::FilterNeutrinosActiveVolume(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  FilterNeutrinosActiveVolume::~FilterNeutrinosActiveVolume()
  {
  }

  //-----------------------------------------------------------------------
  void FilterNeutrinosActiveVolume::beginJob()
  {
    //    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<geo::Geometry> geo;
  
  }

  //-----------------------------------------------------------------------
  void FilterNeutrinosActiveVolume::reconfigure(fhicl::ParameterSet const& p)
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
      filtype       =  p.get< int >("filterType");
      if (filtype==1) {
          typetofilt            = p.get< int >("Interaction");
      }
      if (finActive==0) {
          typetofilt=-999;
      }
      
    return;
  }

  //-----------------------------------------------------------------------
  bool FilterNeutrinosActiveVolume::filter(art::Event& evt)
  {
    bool interactionDesired(false);
    //get the list of particles from this event
    art::ServiceHandle<geo::Geometry> geom;

      
    // * MC truth information

      std::vector< art::Handle< std::vector<simb::MCTruth> > > allmclists;
      evt.getManyByType(allmclists);


      std::cout << fXmin << " " << fXmax << " " << fYmin << " " <<fYmax << " " << fZmin << " " << fZmax << std::endl;
/*
      double xmin = 67.1;
     double xmax = 368.2;
     double ymin = -173;
     double ymax = 143;
     double zmin = -913;
     double zmax = 883;
*/
 
      geo::CryostatGeo const& cryo0 = geom->Cryostat(0);
      geo::CryostatGeo const& cryo1 = geom->Cryostat(1);

      geo::TPCGeo const& tpc00 = cryo0.TPC(0);
      TVector3 xyzcenter00 = tpc00.GetActiveVolumeCenter();
      std::cout << xyzcenter00[0] << " " << xyzcenter00[1] << " " << xyzcenter00[2] << std::endl;

      geo::TPCGeo const& tpc01 = cryo0.TPC(1);
      TVector3 xyzcenter01 = tpc01.GetActiveVolumeCenter();
      std::cout << xyzcenter01[0] << " " << xyzcenter01[1] << " " << xyzcenter01[2] << std::endl;

      geo::TPCGeo const& tpc10 = cryo1.TPC(0);
      TVector3 xyzcenter10 = tpc10.GetActiveVolumeCenter();
      std::cout << xyzcenter10[0] << " " << xyzcenter10[1] << " " << xyzcenter10[2] << std::endl;

      geo::TPCGeo const& tpc11 = cryo1.TPC(1);
      TVector3 xyzcenter11 = tpc11.GetActiveVolumeCenter();
      std::cout << xyzcenter11[0] << " " << xyzcenter11[1] << " " << xyzcenter11[2] << std::endl;

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
      
      /*

       
       std::cout << ((geom->Cryostat(0)).TPC(0)).GetActiveVolumeCenter() << std::endl;
      std::cout << ((geom->Cryostat(1)).TPC(0)).GetActiveVolumeCenter() << std::endl;
      std::cout << ((geom->Cryostat(0)).TPC(1)).GetActiveVolumeCenter() << std::endl;
      std::cout << ((geom->Cryostat(1)).TPC(1)).GetActiveVolumeCenter() << std::endl;
*/
      
      for(size_t mcl = 0; mcl < allmclists.size(); ++mcl){
          art::Handle< std::vector<simb::MCTruth> > mclistHandle = allmclists[mcl];
          for(size_t m = 0; m < mclistHandle->size(); ++m){
              art::Ptr<simb::MCTruth> mct(mclistHandle, m);
              std::cout << mct->GetNeutrino().InteractionType() << " INTERACTION TYPE " << std::endl;
              int value_of_interaction=mct->GetNeutrino().InteractionType();
              for(int ipart=0;ipart<mct->NParticles();ipart++){
                  int pdg=mct->GetParticle(ipart).PdgCode();
                  double xx=mct->GetParticle(ipart).Vx();
                  double yy=mct->GetParticle(ipart).Vy();
                  double zz=mct->GetParticle(ipart).Vz();
      
                  if ((filtype==1 && value_of_interaction==typetofilt) || filtype==0) {
                      
                  if (pdg==12 || pdg==14) {
                      if (finActive==1)
                        {
                          if (xx>(xyzcenter00[0]-w00) && xx<(xyzcenter00[0]+w00) && yy>(xyzcenter00[1]-h00) && yy<(xyzcenter00[1]+h00) && zz>(xyzcenter00[2]-l00/2) && zz<(xyzcenter00[2]+l00/2))
                          {
                              interactionDesired = true;
                          }
                            if (xx>(xyzcenter01[0]-w01) && xx<(xyzcenter01[0]+w01) && yy>(xyzcenter01[1]-h01) && yy<(xyzcenter01[1]+h01) && zz>(xyzcenter01[2]-l01/2) && zz<(xyzcenter01[2]+l01/2))
                            {
                                interactionDesired = true;
                            }
                            if (xx>(xyzcenter10[0]-w10) && xx<(xyzcenter10[0]+w10) && yy>(xyzcenter10[1]-h10) && yy<(xyzcenter10[1]+h10) && zz>(xyzcenter10[2]-l10/2) && zz<(xyzcenter10[2]+l10/2))
                            {
                                interactionDesired = true;
                            }
                            if (xx>(xyzcenter11[0]-w11) && xx<(xyzcenter11[0]+w11) && yy>(xyzcenter11[1]-h11) && yy<(xyzcenter11[1]+h11) && zz>(xyzcenter11[2]-l11/2) && zz<(xyzcenter11[2]+l11/2))
                            {
                                interactionDesired = true;
                            }
                        }
                      if (finActive==0)
                        {
                            if (xx>fXmin && xx<fXmax && yy>fYmin && yy<fYmax && zz>fZmin && zz<fZmax)
                                {
                                    interactionDesired = true;
                                }
                        }
        }
                  }
	  //	std::cout << "FilterNoDirtNeutrinos: i is " << i << std::endl ;
	  // Now walk through trajectory and see if it enters the TPC
	   // trajectory loop
	 // end Genie particle
    }
          }
      }// loop on MCPHandle

  return interactionDesired;
    
  } // end FilterNoDirtNeutrinos()function

} // namespace simfilter

namespace simfilter {

  DEFINE_ART_MODULE(FilterNeutrinosActiveVolume)

} // namespace simfilter

#endif // FILTER_FILTERNEUTRINOSACTIVEVOLUME_H

