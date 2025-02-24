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
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/WireReadout.h"
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

#include "larana/OpticalDetector/SimPhotonCounter.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "lardataobj/Simulation/SimChannel.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

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

////////////////////////////////// Define some constant variable //////////////////////////////
const int nPMTs = 360;
//const int PMTs_per_TPC = 90;
const int MaxPhotons = 10000;

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
void beginJob() override;
//void reconfigure(fhicl::ParameterSet const & p);

private:

//TRandom* Ran;
 
TTree* fTree;

////////////////////////////////// Variable in th tree//////////////////////////////

int event;

int event_type;

int is_Neutrino;
int Neutrino_Interaction;

int noPMT[nPMTs];

int turned_PMT;

double PMTx[nPMTs];
double PMTy[nPMTs];
double PMTz[nPMTs];

int Cryostat[nPMTs];	
int TPC[nPMTs];

float photons_collected[nPMTs];
float QE_photons_collected[nPMTs];
float photon_time[nPMTs][MaxPhotons];

float firstphoton_time[nPMTs];

float true_barycentre_x;
float true_barycentre_y;
float true_barycentre_z;

double vertex_x;
double vertex_y;
double vertex_z;

float reco_barycentre_y;
float reco_barycentre_z;

float total_quenched_energy;
float total_coll_photons;

float PMT_error_y;
float PMT_error_z;
float PMT_total_error;
  
art::InputTag photonLabel;
art::InputTag chargeLabel;
art::InputTag typoLabel;

};


icarus::PMTcoordinates::PMTcoordinates(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  photonLabel(p.get<art::InputTag>("fottoni", "largeant")),
  chargeLabel(p.get<art::InputTag>("carconi", "largeant")),
  typoLabel  (p.get<art::InputTag>("tiponi", "generator"))
 // More initializers here.
{
    std::cout << " PMT coordinates constructor " << std::endl;
}

void icarus::PMTcoordinates::analyze(art::Event const & evt)
{
////////////////////////////////// Create the LArsoft services and service handle//////////////////////////////

geo::WireReadoutGeom const& wireReadoutAlg = art::ServiceHandle<geo::WireReadout const>()->Get();
 
std::vector<sim::SimPhotons> const& optical  = *(evt.getValidHandle<std::vector<sim::SimPhotons>>(photonLabel));
std::vector<sim::SimChannel> const& charge   = *(evt.getValidHandle<std::vector<sim::SimChannel>>(chargeLabel));
//std::vector<simb::MCTruth> const& type    = *(evt.getValidHandle<std::vector<simb::MCTruth>>(typoLabel));

////////////////////////////////// Event number//////////////////////////////

event = evt.id().event();

//std::vector< art::Handle< std::vector<simb::MCTruth> > > type;
//evt.getManyByType(type);
auto type = evt.getMany< std::vector<simb::MCTruth> >();

for(size_t mcl = 0; mcl < type.size(); ++mcl)
{	
	art::Handle< std::vector<simb::MCTruth> > mclistHandle = type[mcl];
	
	for(size_t m = 0; m < mclistHandle->size(); ++m)
	{
		art::Ptr<simb::MCTruth> mct(mclistHandle, m);	
//		for(int ipart=0;ipart<mct->NParticles();ipart++)
//		{	
//			int pdg=mct->GetParticle(ipart).PdgCode();	
//			double xx=mct->GetParticle(ipart).Vx();	
//			double yy=mct->GetParticle(ipart).Vy();
//                	double zz=mct->GetParticle(ipart).Vz();

			event_type=mct->GetParticle(0).PdgCode();	
			vertex_x=mct->GetParticle(0).Vx();	
			vertex_y=mct->GetParticle(0).Vy();
                  	vertex_z=mct->GetParticle(0).Vz();

			if (event_type==12||event_type==-12||event_type==14||event_type==-14||event_type==16||event_type==-16)
			{
				is_Neutrino=1;
				Neutrino_Interaction=mct->GetNeutrino().InteractionType();
			}	
			else
			{
				is_Neutrino=0;
				Neutrino_Interaction=-9999;			
			}		

//		}
	}

}


////////////////////////////////// Putting at 0 all the variables//////////////////////////////

for (int g=0; g<MaxPhotons; g++)
{
	for (int u=0; u<360; u++)
	{
		photon_time[u][g]=0;
	}
} 

true_barycentre_x =0;
true_barycentre_y =0;
true_barycentre_z =0;

total_quenched_energy =0;

////////////////////////////////// Charge part: identify the baricentre of the event //////////////////////////////

for (std::size_t chargechannel = 0;  chargechannel<charge.size(); ++chargechannel) //loop on SimChannel
{ 	
	auto const& channeltdcide = charge.at(chargechannel).TDCIDEMap();
	
	for (std::size_t TDCnu = 0;  TDCnu<channeltdcide.size(); ++TDCnu) 	//loop on TDC
	{

		sim::TDCIDE const& tdcide = channeltdcide.at(TDCnu);

		for (std::size_t IDEnu = 0;  IDEnu<tdcide.second.size(); ++IDEnu) 	//loop on IDE
		{
			sim::IDE const& ida = tdcide.second.at(IDEnu);

//			std::cout << "IDA     " << ida.x << '\t' << ida.y << '\t' << ida.z << std::endl;

			true_barycentre_x = true_barycentre_x + ida.x*ida.energy;
			true_barycentre_y = true_barycentre_y + ida.y*ida.energy;
			true_barycentre_z = true_barycentre_z + ida.z*ida.energy;
			total_quenched_energy      = total_quenched_energy + ida.energy;

		}	//loop on IDE
		
	} 	//loop on TDC

}//loop on SimChannel

true_barycentre_x = true_barycentre_x/total_quenched_energy;
true_barycentre_y = true_barycentre_y/total_quenched_energy;
true_barycentre_z = true_barycentre_z/total_quenched_energy;

total_quenched_energy = total_quenched_energy/3; 

////////////////////////////////// Light part //////////////////////////////////////////////////
turned_PMT=0;

reco_barycentre_y=0;
reco_barycentre_z=0;

total_coll_photons=0;

for (std::size_t channel = 0; channel < optical.size(); ++channel) {

	sim::SimPhotons const& photon_vec = optical[channel];

	noPMT[channel] = channel;	

	photons_collected[channel]= photon_vec.size();
    std::cout << " channel " << channel << " photons collected " << photons_collected[channel] << std::endl;
//	double media = photons_collected[channel]*QE;

//	QE_photons_collected[channel]= Ran.Poisson(media);

	QE_photons_collected[channel]= 0.06*photons_collected[channel];

	if (photons_collected[channel]>0){

	turned_PMT++;
	}

        auto const xyz = wireReadoutAlg.OpDetGeoFromOpChannel(channel).GetCenter();

        PMTx[channel] = xyz.X();
        PMTy[channel] = xyz.Y();
        PMTz[channel] = xyz.Z();

	reco_barycentre_y = reco_barycentre_y + PMTy[channel]*photons_collected[channel];
	reco_barycentre_z = reco_barycentre_z + PMTz[channel]*photons_collected[channel];
	total_coll_photons= total_coll_photons + photons_collected[channel];
        std::cout << " channel " << channel << " total photons  " << total_coll_photons << std::endl;
	firstphoton_time[channel] = 100000000;

	if (photons_collected[channel]>0)
	{
		for (size_t i = 0; i<photon_vec.size() && int(i)< MaxPhotons; ++i)
		{
			photon_time[channel][i]= photon_vec.at(i).Time;

			if (photon_time[channel][i]<firstphoton_time[channel])
			{				
				firstphoton_time[channel]=photon_time[channel][i];
			}
		}

	}


//	std::cout << PMTx[channel] << '\t' << PMTy[channel] << '\t' << PMTz[channel] << std::endl;

	if (PMTx[channel]<0){Cryostat[channel]=0;}
	if (PMTx[channel]>0){Cryostat[channel]=1;}

	if (PMTx[channel]<-200){TPC[channel]=0;}
	if (PMTx[channel]>-200 && PMTx[channel]<0){TPC[channel]=1;}
	if (PMTx[channel]<200 && PMTx[channel]>0){TPC[channel]=2;}
	if (PMTx[channel]>200){TPC[channel]=3;}
    	
}

//total_coll_photons = total_coll_photons;

std::cout << " fotoni finale = " <<total_coll_photons <<std::endl; 

reco_barycentre_y = reco_barycentre_y/total_coll_photons;
reco_barycentre_z = reco_barycentre_z/total_coll_photons;


PMT_error_y = reco_barycentre_y-true_barycentre_y;
PMT_error_z = reco_barycentre_z-true_barycentre_z;
PMT_total_error = sqrt((PMT_error_y*PMT_error_y)+(PMT_error_z*PMT_error_z));

fTree->Fill();
    std::cout << " after filling " << fTree << std::endl;
}

void icarus::PMTcoordinates::beginJob()
{
    std::cout << " PMTcoordinates beginjob " << std::endl;
    
art::ServiceHandle<art::TFileService> tfs;
fTree = tfs->make<TTree>("lighttree","tree for the light response");

fTree->Branch("event",&event,"event/I");
fTree->Branch("event_type",&event_type,"event_type/I");
fTree->Branch("is_Neutrino",&is_Neutrino,"is_Neutrino/I");
fTree->Branch("Neutrino_Interaction",&Neutrino_Interaction,"Neutrino_Interaction/I");
fTree->Branch("total_quenched_energy",&total_quenched_energy,"total_quenched_energy");
fTree->Branch("Cryostat",&Cryostat,("Cryostat[" + std::to_string(nPMTs) + "]/I").c_str());
fTree->Branch("TPC",&TPC,("TPC[" + std::to_string(nPMTs) + "]/I").c_str());
fTree->Branch("noPMT",&noPMT,("noPMT[" + std::to_string(nPMTs) + "]/I").c_str());
fTree->Branch("PMTx",&PMTx,("PMTx[" + std::to_string(nPMTs) + "]/D").c_str());
fTree->Branch("PMTy",&PMTy,("PMTy[" + std::to_string(nPMTs) + "]/D").c_str());
fTree->Branch("PMTz",&PMTz,("PMTz[" + std::to_string(nPMTs) + "]/D").c_str());
fTree->Branch("turned_PMT",&turned_PMT,"turned_PMT/I");
fTree->Branch("total_coll_photons",&total_coll_photons,"total_coll_photons/F");
fTree->Branch("photons_colleted",&photons_collected,("photons_collected[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("QE_photons_colleted",&QE_photons_collected,("QE_photons_collected[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_time",&firstphoton_time,("firstphoton_time[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("photon_time",&photon_time,"photon_time[360][10000]/F");
fTree->Branch("vertex_x",&vertex_x,"vertex_x/D");
fTree->Branch("vertex_y",&vertex_y,"vertex_y/D");
fTree->Branch("vertex_z",&vertex_z,"vertex_z/D");
fTree->Branch("true_barycentre_x",&true_barycentre_x,"true_barycentre_x/F");
fTree->Branch("true_barycentre_y",&true_barycentre_y,"true_barycentre_y/F");
fTree->Branch("true_barycentre_z",&true_barycentre_z,"true_barycentre_z/F");
fTree->Branch("reco_barycentre_y",&reco_barycentre_y,"reco_barycentre_y/F");
fTree->Branch("reco_barycentre_z",&reco_barycentre_z,"reco_barycentre_z/F");
fTree->Branch("PMT_error_y",&PMT_error_y,"PMT_error_y/F");
fTree->Branch("PMT_error_z",&PMT_error_z,"PMT_error_z/F");
fTree->Branch("PMT_total_error",&PMT_total_error,"PMT_total_error/F");
}

DEFINE_ART_MODULE(icarus::PMTcoordinates)
