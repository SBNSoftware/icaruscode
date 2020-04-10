////////////////////////////////////////////////////////////////////////
// Class:       TPCDemo
// Plugin Type: analyzer (art v2_07_03)
// File:        TPCDemo_module.cc
//
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
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

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


#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "nusimdata/SimulationBase/MCTrajectory.h"

#include "nug4/ParticleNavigation/ParticleList.h"

#include "lardataobj/Simulation/sim.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"


namespace icarus {
  class TPCDemo;
}


class icarus::TPCDemo : public art::EDAnalyzer
{
public:
explicit TPCDemo(fhicl::ParameterSet const & p);

TPCDemo(TPCDemo const &) = delete;
TPCDemo(TPCDemo &&) = delete;
TPCDemo & operator = (TPCDemo const &) = delete;
TPCDemo & operator = (TPCDemo &&) = delete;

void analyze(art::Event const & e) override;

void beginJob() override;

private:

 
TTree* fTree;
TTree* recotrack_tree;

int event;

int event_type;

int gt_0;
int gt_1;

double vertex_x;
double vertex_y;
double vertex_z;
double tr_vertex[3];
double tr_end[3];
double tr_length;
double tr_theta;
double tr_phi;

art::InputTag typoLabel;

};


icarus::TPCDemo::TPCDemo(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  typoLabel  (p.get<art::InputTag>("tiponi", "generator"))
 // More initializers here.
{
    std::cout << " Let's try to call some TPC information " << std::endl;
}

void icarus::TPCDemo::analyze(art::Event const & evt)
{
////////////////////////////////// Create the LArsoft services and service handle//////////////////////////////

art::ServiceHandle<geo::Geometry> geom;
 

////////////////////////////////// Event number//////////////////////////////

      event = evt.id().event();

////////////////////////////////// Accessing MCTruth //////////////////////////////
      std::vector< art::Handle< std::vector<simb::MCTruth> > > type;
      evt.getManyByType(type);
      art::Handle< std::vector<simb::MCParticle> > particleVecHandle;
      evt.getByLabel("largeant", particleVecHandle);

      std::vector<const simb::MCParticle*> particleVec;
      if (particleVecHandle.isValid())
      {
          for(size_t idx = 0; idx < particleVecHandle->size(); idx++) particleVec.push_back(&particleVecHandle->at(idx)); //
      }

////////////////////////////////// Accesing Reconstructed tracks//////////////////////////////
    art::Handle< std::vector<recob::Track> > trackVecHandle;
    evt.getByLabel("pandoraTrackICARUSCryo0", trackVecHandle);

    
    std::vector<const recob::Track*> trackVec;
    if (trackVecHandle.isValid())
    {
        for(size_t idx = 0; idx < trackVecHandle->size(); idx++) trackVec.push_back(&trackVecHandle->at(idx)); //
    }


////////////////////////////////// Looping over the True MC information//////////////////////////////


for(size_t mcl = 0; mcl < type.size(); ++mcl)
{	
	art::Handle< std::vector<simb::MCTruth> > mclistHandle = type[mcl];
	
	for(size_t m = 0; m < mclistHandle->size(); ++m)
	{
	  art::Ptr<simb::MCTruth> mct(mclistHandle, m);	
	  for(int ipart=0;ipart<mct->NParticles();ipart++)
	  {	

			event_type=mct->GetParticle(0).PdgCode();	
			vertex_x=mct->GetParticle(0).Vx();	
			vertex_y=mct->GetParticle(0).Vy();
                        vertex_z=mct->GetParticle(0).Vz();

	   }
	}

}

fTree->Fill();
std::cout << " after filling " << fTree << std::endl;

////////////////////////////////// Looping over the TPC Reconstructed Tracks  MC information//////////////////////////////
//
std::cout<<" track "<<trackVec.size()<<std::endl;

for(size_t iTrk = 0; iTrk < trackVec.size(); iTrk++){
        
          tr_vertex[0] =trackVec[iTrk]->Start().X();
          tr_vertex[1] =trackVec[iTrk]->Start().Y();
          tr_vertex[2] =trackVec[iTrk]->Start().Z();
          tr_end[0] =trackVec[iTrk]->End().X();
          tr_end[1] =trackVec[iTrk]->End().Y();
          tr_end[2] =trackVec[iTrk]->End().Z();
          tr_length    =trackVec[iTrk]->Length();
          tr_theta     =trackVec[iTrk]->Theta();
          tr_phi      =trackVec[iTrk]->Phi();           
}

recotrack_tree->Fill();
std::cout << " after filling " << fTree << std::endl;
}

void icarus::TPCDemo::beginJob()
{
    std::cout << " TPCDemo beginjob " << std::endl;
    
art::ServiceHandle<art::TFileService> tfs;
fTree = tfs->make<TTree>("TrueTree","tree for the True GENIE information");
recotrack_tree= tfs->make<TTree>("RecoTrack","tree Recon information");
fTree->Branch("event",&event,"event/I");
fTree->Branch("event_type",&event_type,"event_type/I");
fTree->Branch("vertex_x",&vertex_x,"vertex_x/D");
fTree->Branch("vertex_y",&vertex_y,"vertex_y/D");
fTree->Branch("vertex_z",&vertex_z,"vertex_z/D");

recotrack_tree->Branch("tr_vertex",              &tr_vertex,              "tr_vertex[3]/D");
recotrack_tree->Branch("tr_end",                 &tr_end,                 "tr_end[3]/D");
recotrack_tree->Branch("tr_theta",              &tr_theta,              "tr_theta/D");
recotrack_tree->Branch("tr_phi",              &tr_phi,              "tr_phi/D");
recotrack_tree->Branch("tr_length",              &tr_length,              "tr_length/D");
}

DEFINE_ART_MODULE(icarus::TPCDemo)
