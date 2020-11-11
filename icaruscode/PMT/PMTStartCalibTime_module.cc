////////////////////////////////////////////////////////////////////////
// Class:       PMTStartCalibTime
// Plugin Type: analyzer (art v2_07_03)
// File:        PMTStartCalibTime_module.cc
//
// Prepared by Christian Farnese writing down few ideas for Calibration in time for PMt
// starting from PMTcoordinates module.
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

#include "larana/OpticalDetector/SimPhotonCounter.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "lardataobj/Simulation/SimChannel.h"

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

////////////////////////////////// Define some constant variable //////////////////////////////
const int nPMTs = 360;
//const int PMTs_per_TPC = 90;
const int MaxPhotons = 10000;

namespace icarus {
  class PMTStartCalibTime;
}


class icarus::PMTStartCalibTime : public art::EDAnalyzer
{
public:
explicit PMTStartCalibTime(fhicl::ParameterSet const & p);

// Plugins should not be copied or assigned.
PMTStartCalibTime(PMTStartCalibTime const &) = delete;
PMTStartCalibTime(PMTStartCalibTime &&) = delete;
PMTStartCalibTime & operator = (PMTStartCalibTime const &) = delete;
PMTStartCalibTime & operator = (PMTStartCalibTime &&) = delete;

// Required functions.
void analyze(art::Event const & e) override;

// Selected optional functions.
void beginJob() override;
//void reconfigure(fhicl::ParameterSet const & p);

private:

//TRandom* Ran;
 
TTree* fTree;

////////////////////////////////// Variable in th tree//////////////////////////////

int run;
int subrun;
int event;
int nparticles = 0;
int event_type;

int gt_0;
int gt_1;

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
float firstphoton_time_expected[nPMTs];
float firstphoton_tp_expected[nPMTs];
float firstphoton_xp_expected[nPMTs];
float firstphoton_yp_expected[nPMTs];
float firstphoton_zp_expected[nPMTs];


    float firstphoton_time_expected_r[nPMTs];

float firstphoton_tp_expected_r[nPMTs];
float firstphoton_xp_expected_r[nPMTs];
float firstphoton_yp_expected_r[nPMTs];
float firstphoton_zp_expected_r[nPMTs];

float firstphoton_tp_expected_v2[nPMTs];
float firstphoton_xp_expected_v2[nPMTs];
float firstphoton_yp_expected_v2[nPMTs];
float firstphoton_zp_expected_v2[nPMTs];
float firstphoton_tp_expected_v3[nPMTs];
float firstphoton_xp_expected_v3[nPMTs];
float firstphoton_yp_expected_v3[nPMTs];
float firstphoton_zp_expected_v3[nPMTs];


float firstphoton_time_expected_v2[nPMTs];
float firstphoton_time_expected_v3[nPMTs];
float firstphoton_time_expected_d2[nPMTs];
float firstphoton_time_expected_d3[nPMTs];

float true_barycentre_x;
float true_barycentre_y;
float true_barycentre_z;

double vertex_x;
double vertex_y;
double vertex_z;

float total_quenched_energy;
float total_coll_photons;
  
art::InputTag photonLabel;
art::InputTag chargeLabel;
art::InputTag typoLabel;

};


icarus::PMTStartCalibTime::PMTStartCalibTime(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  photonLabel(p.get<art::InputTag>("fottoni", "largeant")),
  chargeLabel(p.get<art::InputTag>("carconi", "largeant")),
  typoLabel  (p.get<art::InputTag>("tiponi", "generator"))
 // More initializers here.
{
    std::cout << " Let's try to calib in time the PMTs... " << std::endl;
}

void icarus::PMTStartCalibTime::analyze(art::Event const & evt)
{
////////////////////////////////// Create the LArsoft services and service handle//////////////////////////////

art::ServiceHandle<geo::Geometry> geom;
 
std::vector<sim::SimPhotons> const& optical  = *(evt.getValidHandle<std::vector<sim::SimPhotons>>(photonLabel));
std::vector<sim::SimChannel> const& charge   = *(evt.getValidHandle<std::vector<sim::SimChannel>>(chargeLabel));
//std::vector<simb::MCTruth> const& type    = *(evt.getValidHandle<std::vector<simb::MCTruth>>(typoLabel));

////////////////////////////////// Event number//////////////////////////////

run = evt.id().run();
subrun = evt.id().subRun();
event = evt.id().event();

std::vector< art::Handle< std::vector<simb::MCTruth> > > type;
evt.getManyByType(type);
      art::Handle< std::vector<simb::MCParticle> > particleVecHandle;
      evt.getByLabel("largeant", particleVecHandle);

      std::vector<const simb::MCParticle*> particleVec;
      if (particleVecHandle.isValid())
      {
        for(size_t idx = 0; idx < particleVecHandle->size(); idx++)
        {
            const simb::MCParticle* particle = &particleVecHandle->at(idx);
            particleVec.push_back(particle); //

            // Count cosmics in the event
            if( particle->Mother() == 0 ){
                nparticles++;
            }

        }
      }
    art::Handle< std::vector<recob::Track> > trackVecHandle;
    evt.getByLabel("pandoraTrackICARUSCryo0", trackVecHandle);//at the moment I am calling here directly a particular kind of track: this should be generalized, putting the label of the track type in the fhicl...
    
    std::vector<const recob::Track*> trackVec;
    if (trackVecHandle.isValid())
    {
        for(size_t idx = 0; idx < trackVecHandle->size(); idx++) trackVec.push_back(&trackVecHandle->at(idx)); //
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

total_coll_photons=0;

for (std::size_t channel = 0; channel < optical.size(); ++channel) {

	sim::SimPhotons const& photon_vec = optical[channel];

	noPMT[channel] = channel;	

	photons_collected[channel]= photon_vec.size();
    //std::cout << " channel " << channel << " photons collected " << photons_collected[channel] << std::endl;
//	double media = photons_collected[channel]*QE;

//	QE_photons_collected[channel]= Ran.Poisson(media);

	QE_photons_collected[channel]= 0.06*photons_collected[channel];

	if (photons_collected[channel]>0){

	turned_PMT++;
	}

	double xyz[3];

	geom->OpDetGeoFromOpChannel(channel).GetCenter(xyz);

	PMTx[channel] = xyz[0];
	PMTy[channel] = xyz[1];
	PMTz[channel] = xyz[2];

	total_coll_photons= total_coll_photons + photons_collected[channel];
        //std::cout << " channel " << channel << " total photons  " << total_coll_photons << std::endl;
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



for(int ijk=0;ijk<360;ijk++)
{
    firstphoton_time_expected[ijk]=100000000;
    firstphoton_time_expected_r[ijk]=100000000;
    firstphoton_tp_expected[ijk]=100000000;
    firstphoton_tp_expected_r[ijk]=100000000;
    firstphoton_xp_expected[ijk]=100000000;
    firstphoton_xp_expected_r[ijk]=100000000;
    firstphoton_yp_expected[ijk]=100000000;
    firstphoton_yp_expected_r[ijk]=100000000;
    firstphoton_zp_expected[ijk]=100000000;
    firstphoton_zp_expected_r[ijk]=100000000;

    firstphoton_tp_expected_v2[ijk]=100000000;
    firstphoton_tp_expected_v3[ijk]=100000000;
    firstphoton_xp_expected_v2[ijk]=100000000;
    firstphoton_xp_expected_v3[ijk]=100000000;
    firstphoton_yp_expected_v2[ijk]=100000000;
    firstphoton_yp_expected_v3[ijk]=100000000;
    firstphoton_zp_expected_v2[ijk]=100000000;
    firstphoton_zp_expected_v3[ijk]=100000000;

    firstphoton_time_expected_v3[ijk]=100000000;
    firstphoton_time_expected_v2[ijk]=100000000;
    firstphoton_time_expected_d2[ijk]=100000000;
    firstphoton_time_expected_d3[ijk]=100000000;
}

for(size_t mcl = 0; mcl < type.size(); ++mcl)
{	
	art::Handle< std::vector<simb::MCTruth> > mclistHandle = type[mcl];
	
	for(size_t m = 0; m < mclistHandle->size(); ++m)
	{
	  art::Ptr<simb::MCTruth> mct(mclistHandle, m);	
	  for(int ipart=0;ipart<mct->NParticles();ipart++)
	  {	
		  //int pdg=mct->GetParticle(ipart).PdgCode();	
		  //double xx=mct->GetParticle(ipart).Vx();	
		  //double yy=mct->GetParticle(ipart).Vy();
		  //double zz=mct->GetParticle(ipart).Vz();

			event_type=mct->GetParticle(0).PdgCode();	
			vertex_x=mct->GetParticle(0).Vx();	
			vertex_y=mct->GetParticle(0).Vy();
            vertex_z=mct->GetParticle(0).Vz();



	   }
	}

}

    gt_0=1;
    gt_1=1;

    
    for(const auto& tracknow3 : trackVec)
    {
        std::cout << tracknow3->NumberTrajectoryPoints() << std::endl;
                int ntraj = tracknow3->NumberTrajectoryPoints();
        int last_good_point=ntraj;
        for ( size_t ijk=0;ijk<tracknow3->NumberTrajectoryPoints();ijk++)
        {
            auto trajPoint = tracknow3->TrajectoryPoint(ijk);
            if(trajPoint.position.X()<-600)
            {
                last_good_point=ijk;
                ijk=ntraj+1;
            }
            //std::cout << "Start of trajectory at " << trajPoint.position.X() << " " <<  trajPoint.position.Y() << " " << trajPoint.position.Z()  << " cm " << std::endl;
        }



	float y_max_value=-150;
        int point_y_max=-1;

        for ( size_t ijk=0;ijk<tracknow3->NumberTrajectoryPoints();ijk++)
        {
            auto trajPoint = tracknow3->TrajectoryPoint(ijk);
            if(trajPoint.position.Y()>y_max_value)
            {
                point_y_max=ijk;
                y_max_value=trajPoint.position.Y();
            }
            //std::cout << "Start of trajectory at " << trajPoint.position.X() << " " <<  trajPoint.position.Y() << " " << trajPoint.position.Z()  << " cm " << std::endl;
         }


        float c_f_x;
        float c_f_y;
        float c_f_z;

          
if(point_y_max>-1)
{
        auto first_point=tracknow3->TrajectoryPoint(point_y_max);
            c_f_x=first_point.position.X();
            c_f_y=first_point.position.Y();
            c_f_z=first_point.position.Z();
/*
        auto first_point=tracknow3->TrajectoryPoint(0);
        auto last_point=tracknow3->TrajectoryPoint(last_good_point);
        float coordinates_y_0=first_point.position.Y();
        float coordinates_y_l=last_point.position.Y();
        float c_f_x;
        float c_f_y;
        float c_f_z;
        if (coordinates_y_0>coordinates_y_l) {
            c_f_x=first_point.position.X();
            c_f_y=first_point.position.Y();
            c_f_z=first_point.position.Z();
        }
        if (coordinates_y_0<=coordinates_y_l) {
            c_f_x=last_point.position.X();
            c_f_y=last_point.position.Y();
            c_f_z=last_point.position.Z();
        }
*/
std::cout << "Start of trajectory at " << first_point.position.X() << " " <<  first_point.position.Y() << " " << first_point.position.Z()  << " cm " << std::endl;
        for ( int ijk2=0;ijk2<last_good_point;ijk2++)
        {
        auto trajPoint = tracknow3->TrajectoryPoint(ijk2);
        //std::cout << "Start of trajectory at " << trajPoint.position.X() << " " <<  trajPoint.position.Y() << " " << trajPoint.position.Z()  << " cm " << std::endl;
        float xpr=trajPoint.position.X();
        float ypr=trajPoint.position.Y();
        float zpr=trajPoint.position.Z();
        float distance_from_first_point=sqrt((xpr-c_f_x)*(xpr-c_f_x)+(ypr-c_f_y)*(ypr-c_f_y)+(zpr-c_f_z)*(zpr-c_f_z))/100;    
        float tpr=distance_from_first_point/0.3;
        int factor=0;
        if(xpr>0)factor=180;
            if(xpr>-400)
            {
        for(int ijk=0;ijk<180;ijk++)
        {
            float distance_x=xpr-PMTx[ijk+factor];
            float distance_y=ypr-PMTy[ijk+factor];
            float distance_z=zpr-PMTz[ijk+factor];
            float total_distance=sqrt(distance_x*distance_x+distance_y*distance_y+distance_z*distance_z)/100;
            float time_travel_light=total_distance/(0.1);//velocita' e' 3*10^8 m/s cioe' 0.3 m/ns ma sto usando c/n e ho trovato n=1.5 da Petrillo al collaboration meeting
            //std::cout << time_travel_light+tp << " " << tp << " " << time_travel_light << std::endl;
            if(firstphoton_time_expected_r[ijk+factor]>(tpr+time_travel_light))
{
firstphoton_time_expected_r[ijk+factor]=(tpr+time_travel_light);
           firstphoton_tp_expected_r[ijk+factor]=(tpr);
           firstphoton_xp_expected_r[ijk+factor]=(xpr);
           firstphoton_yp_expected_r[ijk+factor]=(ypr);
           firstphoton_zp_expected_r[ijk+factor]=(zpr);
}
        }
        }
        }
}
        
        /*
        for ( recob::Trajectory::const_iterator trajectory = tracknow3->Trajectory().begin(); trajectory != tracknow3->Trajectory().end(); ++trajectory)
        {
            float xp=(*trajectory).first.X();
            float yp=(*trajectory).first.Y();
            float zp=(*trajectory).first.Z();
            float tp=(*trajectory).first.T();
            std::cout << xp << " " << yp << " " << zp << " " << tp << std::endl;
        }
         */
    }
    
      for(const auto& particlenow3 : particleVec)
      {

for ( simb::MCTrajectory::const_iterator trajectory = particlenow3->Trajectory().begin(); trajectory != particlenow3->Trajectory().end(); ++trajectory)
{
	  float xp=(*trajectory).first.X(); 
	  float yp=(*trajectory).first.Y();   
	  float zp=(*trajectory).first.Z(); 
	  float tp=(*trajectory).first.T(); 
	  //std::cout << xp << " " << yp << " " << zp << " " << tp << std::endl;
	  int factor=0;
    if (yp>-185 && yp<140 && fabs(zp)<910)
    {
        if (xp<-343 || xp>-247) {
            gt_0=0;
        }
        if (xp<-193 || xp>-98) {
            gt_1=0;
        }
    }
	  if(xp>0)factor=180;
          if(yp>-185 && yp<140 && fabs(zp)<910 && fabs(xp)<370 && fabs(xp)>70){
	  for(int ijk=0;ijk<180;ijk++)
	    {
	      float distance_x=xp-PMTx[ijk+factor];
	      float distance_y=yp-PMTy[ijk+factor];
	      float distance_z=zp-PMTz[ijk+factor];
	      float total_distance=sqrt(distance_x*distance_x+distance_y*distance_y+distance_z*distance_z)/100;
	      float time_travel_light_v2=total_distance/(0.2);//velocity used is 20 cm/ns
	      float time_travel_light=total_distance/(0.1);//velocity used is 10 cm/ns
         float time_travel_light_v3=total_distance/(0.05);//velocity used is 5 cm/ns
	      if(firstphoton_time_expected[ijk+factor]>(tp+time_travel_light))
{
firstphoton_time_expected[ijk+factor]=(tp+time_travel_light);
firstphoton_tp_expected[ijk+factor]=(tp);
firstphoton_xp_expected[ijk+factor]=(xp);
firstphoton_yp_expected[ijk+factor]=(yp);
firstphoton_zp_expected[ijk+factor]=(zp);
}
              if(firstphoton_time_expected_v2[ijk+factor]>(tp+time_travel_light_v2))
{
firstphoton_time_expected_v2[ijk+factor]=(tp+time_travel_light_v2);
firstphoton_tp_expected_v2[ijk+factor]=(tp);
firstphoton_xp_expected_v2[ijk+factor]=(xp);
firstphoton_yp_expected_v2[ijk+factor]=(yp);
firstphoton_zp_expected_v2[ijk+factor]=(zp);
}
              if(firstphoton_time_expected_v3[ijk+factor]>(tp+time_travel_light_v3))
{
firstphoton_time_expected_v3[ijk+factor]=(tp+time_travel_light_v3);
firstphoton_tp_expected_v3[ijk+factor]=(tp);
firstphoton_xp_expected_v3[ijk+factor]=(xp);
firstphoton_yp_expected_v3[ijk+factor]=(yp);
firstphoton_zp_expected_v3[ijk+factor]=(zp);
}

	    }
    for(int ijk=0;ijk<180;ijk++)
    {
        float distance_x=(xp-25)-PMTx[ijk+factor];
        float distance_y=yp-PMTy[ijk+factor];
        float distance_z=zp-PMTz[ijk+factor];
        float total_distance=sqrt(distance_x*distance_x+distance_y*distance_y+distance_z*distance_z)/100;
        float time_travel_light_d2=total_distance/(0.1);//velocita' e' 3*10^8 m/s cioe' 0.3 m/ns ma sto usando c/n e ho trovato n=1.5 da Petrillo al collaboration meeting
        if(firstphoton_time_expected_d2[ijk+factor]>(tp+time_travel_light_d2))firstphoton_time_expected_d2[ijk+factor]=(tp+time_travel_light_d2);
    }
    for(int ijk=0;ijk<180;ijk++)
    {
        float distance_x=(xp+25)-PMTx[ijk+factor];
        float distance_y=yp-PMTy[ijk+factor];
        float distance_z=zp-PMTz[ijk+factor];
        float total_distance=sqrt(distance_x*distance_x+distance_y*distance_y+distance_z*distance_z)/100;
        float time_travel_light_d3=total_distance/(0.1);//velocita' e' 3*10^8 m/s cioe' 0.3 m/ns ma sto usando c/n e ho trovato n=1.5 da Petrillo al collaboration meeting
        if(firstphoton_time_expected_d3[ijk+factor]>(tp+time_travel_light_d3))firstphoton_time_expected_d3[ijk+factor]=(tp+time_travel_light_d3);
    }
}
    
    
 }
      }



fTree->Fill();

nparticles=0;
    std::cout << " after filling " << fTree << std::endl;
}

void icarus::PMTStartCalibTime::beginJob()
{
    std::cout << " PMTStartCalibTime beginjob " << std::endl;
    
art::ServiceHandle<art::TFileService> tfs;
fTree = tfs->make<TTree>("lighttree","tree for the light response");

fTree->Branch("run",&run,"run/I");
fTree->Branch("subrun",&subrun,"subrun/I");
fTree->Branch("event",&event,"event/I");
fTree->Branch("nparticles",&nparticles,"nparticles/I");
fTree->Branch("event_type",&event_type,"event_type/I");
fTree->Branch("gt_0",&gt_0,"gt_0/I");
fTree->Branch("gt_1",&gt_1,"gt_1/I");
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

fTree->Branch("firstphoton_time_expected",&firstphoton_time_expected,("firstphoton_time_expected[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_time_expected_r",&firstphoton_time_expected_r,("firstphoton_time_expected_r[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_tp_expected",&firstphoton_tp_expected,("firstphoton_tp_expected[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_tp_expected_r",&firstphoton_tp_expected_r,("firstphoton_tp_expected_r[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_xp_expected",&firstphoton_xp_expected,("firstphoton_xp_expected[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_xp_expected_r",&firstphoton_xp_expected_r,("firstphoton_xp_expected_r[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_yp_expected",&firstphoton_yp_expected,("firstphoton_yp_expected[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_yp_expected_r",&firstphoton_yp_expected_r,("firstphoton_yp_expected_r[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_zp_expected",&firstphoton_zp_expected,("firstphoton_zp_expected[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_zp_expected_r",&firstphoton_zp_expected_r,("firstphoton_zp_expected_r[" + std::to_string(nPMTs) + "]/F").c_str());

fTree->Branch("firstphoton_tp_expected_v2",&firstphoton_tp_expected_v2,("firstphoton_tp_expected_v2[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_tp_expected_v3",&firstphoton_tp_expected_v3,("firstphoton_tp_expected_v3[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_xp_expected_v2",&firstphoton_xp_expected_v2,("firstphoton_xp_expected_v2[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_xp_expected_v3",&firstphoton_xp_expected_v3,("firstphoton_xp_expected_v3[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_yp_expected_v2",&firstphoton_yp_expected_v2,("firstphoton_yp_expected_v3[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_yp_expected_v3",&firstphoton_yp_expected_v3,("firstphoton_yp_expected_v2[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_zp_expected_v2",&firstphoton_zp_expected_v2,("firstphoton_zp_expected_v2[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_zp_expected_v3",&firstphoton_zp_expected_v3,("firstphoton_zp_expected_v3[" + std::to_string(nPMTs) + "]/F").c_str());


fTree->Branch("firstphoton_time_expected_v2",&firstphoton_time_expected_v2,("firstphoton_time_expected_v2[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_time_expected_v3",&firstphoton_time_expected_v3,("firstphoton_time_expected_v3[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_time_expected_d2",&firstphoton_time_expected_d2,("firstphoton_time_expected_d2[" + std::to_string(nPMTs) + "]/F").c_str());
fTree->Branch("firstphoton_time_expected_d3",&firstphoton_time_expected_d3,("firstphoton_time_expected_d3[" + std::to_string(nPMTs) + "]/F").c_str());
//fTree->Branch("photon_time",&photon_time,"photon_time[360][10000]/F");


fTree->Branch("vertex_x",&vertex_x,"vertex_x/D");
fTree->Branch("vertex_y",&vertex_y,"vertex_y/D");
fTree->Branch("vertex_z",&vertex_z,"vertex_z/D");
fTree->Branch("true_barycentre_x",&true_barycentre_x,"true_barycentre_x/F");
fTree->Branch("true_barycentre_y",&true_barycentre_y,"true_barycentre_y/F");
fTree->Branch("true_barycentre_z",&true_barycentre_z,"true_barycentre_z/F");
}

DEFINE_ART_MODULE(icarus::PMTStartCalibTime)
