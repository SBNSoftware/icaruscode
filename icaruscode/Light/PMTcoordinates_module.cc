////////////////////////////////////////////////////////////////////////
// Class:       PMTcoordinates
// Plugin Type: analyzer (art v2_07_03)
// File:        PMTcoordinates_module.cc
//
// Generated at Mon Sep 25 15:10:31 2017 by Andrea Falcone using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"


// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

const int nPMTs = 360;
const int PMTs_per_TPC = 90;

namespace icarus {
  class PMTcoordinates;
}


class icarus::PMTcoordinates : public art::EDAnalyzer 
{
public:
  explicit PMTcoordinates(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  PMTcoordinates(PMTcoordinates const &) = delete;
  PMTcoordinates(PMTcoordinates &&) = delete;
  PMTcoordinates & operator = (PMTcoordinates const &) = delete;
  PMTcoordinates & operator = (PMTcoordinates &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  //void beginJob();
  //void reconfigure(fhicl::ParameterSet const & p);

private:

double coordinate [3][nPMTs];

int PMT [nPMTs];

float photons [nPMTs];
};


icarus::PMTcoordinates::PMTcoordinates(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

}

void icarus::PMTcoordinates::analyze(art::Event const & e)
{
 art::ServiceHandle<geo::Geometry> geom;
 art::ServiceHandle<opdet::SimPhotonCounter> optical;

 optical->TotalPhotonVector(total);	


  //auto event = evt.id().event();

for (int channel=0; channel < nPMTs; channel++)
{
 double xyz[3];
 float total;
 geom->OpDetGeoFromOpChannel(channel).GetCenter(xyz);

 coordinate[0][channel] = xyz[0];
 coordinate[1][channel] = xyz[1];
 coordinate[2][channel] = xyz[2];

 std::cout<< total << '\t'<< coordinate[0][channel] << '\t' << coordinate[1][channel] << '\t' << coordinate[2][channel] << std::endl;
}

}

DEFINE_ART_MODULE(icarus::PMTcoordinates)
